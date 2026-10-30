"""
Landslide Candidate Patch Detection Pipeline
---------------------------------------------
Publication-grade geospatial analysis script for detecting potential landslide 
source regions (Downslope Motion) and toe uplift zones (Toe Uplift) from InSAR displacement GeoTIFF data.

Outputs QGIS-compatible vector files (.gpkg and .geojson) containing labeled patch polygons
with spatial statistics.
"""

import os
import numpy as np
import rasterio
from rasterio.features import shapes
import xarray as xr
import rioxarray
import geopandas as gpd
from shapely.geometry import shape
from scipy.ndimage import median_filter, label, generate_binary_structure


def load_and_align_data(disp_path: str = "displacement.tif", dem_path: str = "dem.tif"):
    """
    Loads displacement raster and DEM dataset, aligning DEM spatial extent 
    and CRS to match displacement.tif pixel grid. Calculates DEM slope (degrees).
    """
    if not os.path.exists(disp_path):
        raise FileNotFoundError(f"Displacement raster not found: {disp_path}")

    # Load displacement with rasterio & rioxarray
    disp_da = rioxarray.open_rasterio(disp_path).squeeze()
    disp_arr = disp_da.values.astype(np.float32)

    # Spatial parameters from displacement raster
    crs = disp_da.rio.crs
    transform = disp_da.rio.transform()
    pixel_width = abs(transform[0])
    pixel_height = abs(transform[4])
    pixel_area_m2 = pixel_width * pixel_height

    # Load DEM if available and align
    slope_deg_arr = None
    dem_arr = None

    if os.path.exists(dem_path):
        try:
            dem_ds = xr.open_dataset(dem_path)
            dem_z = dem_ds['z'].rio.write_crs("EPSG:4326")
            dem_reprojected = dem_z.rio.reproject_match(disp_da)
            dem_arr = dem_reprojected.values.astype(np.float32)

            # Compute slope in degrees using spatial gradient
            dy, dx = np.gradient(dem_arr, pixel_height, pixel_width)
            slope_rad = np.arctan(np.sqrt(dx**2 + dy**2))
            slope_deg_arr = np.degrees(slope_rad)
        except Exception as e:
            print(f"[WARNING] Could not process dem.tif: {e}. Slope calculation skipped.")

    return {
        "disp_da": disp_da,
        "disp_arr": disp_arr,
        "dem_arr": dem_arr,
        "slope_arr": slope_deg_arr,
        "crs": crs,
        "transform": transform,
        "pixel_width": pixel_width,
        "pixel_height": pixel_height,
        "pixel_area_m2": pixel_area_m2,
    }


def extract_landslide_patches(
    disp_path: str = "displacement.tif",
    dem_path: str = "dem.tif",
    x_thresh_mm: float = 5.0,
    min_area_m2: float = 50.0,
    min_slope_deg: float = 0.0,
    k_sigma: float = 0.0,
    output_gpkg: str = "candidate_patches.gpkg",
    output_geojson: str = "candidate_patches.geojson"
):
    """
    Executes candidate patch extraction for both Downslope Motion and Toe Uplift.

    Parameters:
    ----------
    disp_path : str
        Path to displacement.tif (values in mm)
    dem_path : str
        Path to dem.tif (elevation dataset)
    x_thresh_mm : float
        Displacement threshold in mm (default = 5.0 mm)
    min_area_m2 : float
        Minimum patch area threshold in sq. meters (default = 50.0 m^2, user-controlled)
    min_slope_deg : float
        Minimum topographic slope threshold in degrees (default = 0.0 deg)
    k_sigma : float
        Multiplier for statistical background noise thresholding (0.0 to use exact x_thresh_mm)
    output_gpkg : str
        Output QGIS GeoPackage file path
    output_geojson : str
        Output QGIS GeoJSON file path
    """
    # 1. Load data
    data = load_and_align_data(disp_path, dem_path)
    disp_arr = data["disp_arr"]
    dem_arr = data["dem_arr"]
    slope_arr = data["slope_arr"]
    crs = data["crs"]
    transform = data["transform"]

    # 2. Mask invalid / NaN values
    valid_mask = ~np.isnan(disp_arr)

    # 3. Compute adaptive or fixed threshold
    if k_sigma > 0.0:
        background_pixels = disp_arr[valid_mask & (np.abs(disp_arr) < 15.0)]
        mu_bg = float(np.mean(background_pixels)) if len(background_pixels) > 0 else 0.0
        sigma_bg = float(np.std(background_pixels)) if len(background_pixels) > 0 else 1.0
        effective_thresh = max(x_thresh_mm, abs(mu_bg) + k_sigma * sigma_bg)
    else:
        effective_thresh = float(x_thresh_mm)

    print(f"[INFO] Effective Threshold Applied: ±{effective_thresh:.2f} mm")

    # 4. Spatial 3x3 Median Filter to reduce speckle noise
    filtered_disp = median_filter(np.nan_to_num(disp_arr, nan=0.0), size=3)

    # Topographic slope mask
    if slope_arr is not None and min_slope_deg > 0.0:
        slope_mask = valid_mask & (slope_arr >= min_slope_deg)
    else:
        slope_mask = valid_mask

    # 5. Define Binary Masks for Downslope vs. Toe Uplift
    downslope_mask = (filtered_disp < -effective_thresh) & slope_mask
    uplift_mask = (filtered_disp > +effective_thresh) & slope_mask

    # 6. Process each motion type
    motion_types = [
        ("Downslope", downslope_mask, "#e41a1c"),   # Red marking for Downslope (Range lengthening)
        ("Toe_Uplift", uplift_mask, "#377eb8")      # Blue marking for Toe Uplift (Range shortening)
    ]

    structure_8 = generate_binary_structure(2, 2)  # 8-neighbor connectivity
    patch_records = []
    global_patch_counter = 1

    with rasterio.open(disp_path) as src:
        for motion_name, mask, color in motion_types:
            labeled_arr, num_features = label(mask, structure=structure_8)
            print(f"[INFO] Identified {num_features} raw candidate clusters for {motion_name}")

            for geom_dict, val in shapes(labeled_arr.astype(np.int32), mask=mask, transform=transform):
                if val == 0:
                    continue

                poly = shape(geom_dict)
                area_m2 = poly.area

                # User-controlled minimum area filter (y m^2)
                if area_m2 < min_area_m2:
                    continue

                perimeter_m = poly.length
                compactness = (4 * np.pi * area_m2) / (perimeter_m**2) if perimeter_m > 0 else 0.0

                # Extract internal raster statistics with orientation-robust pixel indexing
                minx, miny, maxx, maxy = poly.bounds
                r1, c1 = src.index(minx, maxy)
                r2, c2 = src.index(maxx, miny)

                r_s, r_e = min(r1, r2), max(r1, r2) + 1
                c_s, c_e = min(c1, c2), max(c1, c2) + 1

                r_s, r_e = max(0, r_s), min(disp_arr.shape[0], r_e)
                c_s, c_e = max(0, c_s), min(disp_arr.shape[1], c_e)

                patch_pixels = disp_arr[r_s:r_e, c_s:c_e]
                patch_valid = patch_pixels[~np.isnan(patch_pixels)]

                if len(patch_valid) == 0:
                    continue

                mean_disp = float(np.mean(patch_valid))
                max_disp = float(np.min(patch_valid) if motion_name == "Downslope" else np.max(patch_valid))
                std_disp = float(np.std(patch_valid))

                # Extract DEM elevation & slope stats if present
                mean_elev = float(np.nanmean(dem_arr[r_s:r_e, c_s:c_e])) if dem_arr is not None else 0.0
                mean_slope = float(np.nanmean(slope_arr[r_s:r_e, c_s:c_e])) if slope_arr is not None else 0.0

                patch_records.append({
                    "Patch_ID": f"PATCH_{global_patch_counter:03d}",
                    "Motion_Type": motion_name,
                    "Area_m2": round(area_m2, 2),
                    "Perimeter_m": round(perimeter_m, 2),
                    "Compactness": round(compactness, 3),
                    "Mean_Displ_mm": round(mean_disp, 2),
                    "Max_Displ_mm": round(max_disp, 2),
                    "Std_Displ_mm": round(std_disp, 2),
                    "Mean_Elev_m": round(mean_elev, 1),
                    "Mean_Slope_deg": round(mean_slope, 1),
                    "Pixel_Count": len(patch_valid),
                    "QGIS_Color": color,
                    "geometry": poly
                })

                global_patch_counter += 1

    # 7. Create GeoDataFrame
    if len(patch_records) == 0:
        print("[WARNING] No patches met the specified area and displacement thresholds.")
        gdf = gpd.GeoDataFrame(columns=[
            "Patch_ID", "Motion_Type", "Area_m2", "Perimeter_m", "Compactness",
            "Mean_Displ_mm", "Max_Displ_mm", "Std_Displ_mm", "Mean_Elev_m", 
            "Mean_Slope_deg", "Pixel_Count", "QGIS_Color", "geometry"
        ], crs=crs)
    else:
        gdf = gpd.GeoDataFrame(patch_records, crs=crs)

    # 8. Export directly importable QGIS formats
    gdf.to_file(output_gpkg, driver="GPKG", layer="candidate_patches")
    gdf.to_file(output_geojson, driver="GeoJSON")

    print(f"[SUCCESS] Exported {len(gdf)} Candidate Patches:")
    print(f"  -> QGIS GeoPackage: {os.path.abspath(output_gpkg)}")
    print(f"  -> QGIS GeoJSON:    {os.path.abspath(output_geojson)}")

    return gdf


if __name__ == "__main__":
    # Standard script entrypoint with default parameters (x = 5 mm, y = 50 m^2)
    extract_landslide_patches(
        disp_path="displacement.tif",
        dem_path="dem.tif",
        x_thresh_mm=5.0,
        min_area_m2=50.0,
        output_gpkg="candidate_patches.gpkg",
        output_geojson="candidate_patches.geojson"
    )
