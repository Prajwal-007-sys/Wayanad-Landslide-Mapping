import os
import sys
import asf_search as asf
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds

# =====================================================================
# NASA EARTHDATA LOGIN CREDENTIALS
# Please replace 'your_username_here' and 'your_password_here' with
# your actual Earthdata credentials.
# =====================================================================
EARTHDATA_USERNAME = "prajwal_007"
EARTHDATA_PASSWORD = "Prominous@123"
# =====================================================================

# Query Parameters (Choorlamala, Wayanad)
BBOX = [76.12683, 11.45700, 76.14517, 11.47500]  # [West, South, East, North]
DATE_RANGE = "2024-01-01T00:00:00Z/2024-08-01T23:59:59Z"
OUT_DIR = "output_asf"
POLARIZATIONS = ["VV"]
TARGET_BURST_ID = "165_353270_IW3"

# Helper to convert W,S,E,N BBOX to WKT Polygon for ASF Search API
def bbox_to_wkt(bbox):
    w, s, e, n = bbox
    return f"POLYGON(({w} {s}, {e} {s}, {e} {n}, {w} {n}, {w} {s}))"

# Helper to crop all downloaded TIFFs to target BBOX
def crop_downloaded_bursts(directory, bbox):
    import glob
    tiff_files = glob.glob(os.path.join(directory, "*.tiff")) + glob.glob(os.path.join(directory, "*.tif"))
    for filepath in tiff_files:
        print(f"  Cropping {os.path.basename(filepath)}...")
        temp_filepath = filepath + ".tmp"
        try:
            with rasterio.open(filepath) as src:
                bounds = transform_bounds("EPSG:4326", src.crs, *bbox)
                window = from_bounds(*bounds, transform=src.transform)
                
                # Read all bands to preserve real/imaginary values
                data = src.read(window=window)
                transform = src.window_transform(window)
                
                profile = src.profile.copy()
                bands = data.shape[0] if len(data.shape) == 3 else 1
                height = data.shape[1] if len(data.shape) == 3 else data.shape[0]
                width = data.shape[2] if len(data.shape) == 3 else data.shape[1]
                
                profile.update(
                    height=height,
                    width=width,
                    transform=transform,
                    count=bands
                )
            
            with rasterio.open(temp_filepath, "w", **profile) as dst:
                dst.write(data)
            
            # Replace original file with cropped file
            os.replace(temp_filepath, filepath)
        except Exception as e:
            print(f"    [WARNING] Failed to crop {os.path.basename(filepath)}: {e}")
            if os.path.exists(temp_filepath):
                os.remove(temp_filepath)

def main():
    print("--------------------------------------------------")
    print("      ASF Sentinel-1 SLC Burst Extractor          ")
    print("--------------------------------------------------")
    
    os.makedirs(OUT_DIR, exist_ok=True)
    
    wkt_polygon = bbox_to_wkt(BBOX)
    start_date, end_date = DATE_RANGE.split("/")
    
    print(f"Area of Interest (WKT): {wkt_polygon}")
    print(f"Temporal Window: {start_date} to {end_date}")
    print("Querying ASF catalog...")
    
    # 1. Perform Search (Public search, does not require authentication)
    try:
        results = asf.search(
            dataset=asf.DATASET.SLC_BURST,
            intersectsWith=wkt_polygon,
            start=start_date,
            end=end_date,
            polarization=POLARIZATIONS,
            fullBurstID=TARGET_BURST_ID
        )
    except Exception as e:
        print(f"[ERROR] Failed to query ASF API: {e}")
        sys.exit(1)
        
    total_found = len(results)
    print(f"\n[SUCCESS] Found {total_found} SLC burst files (filtered to {TARGET_BURST_ID} and VV).")
    
    if total_found == 0:
        print("No matching bursts found for the given criteria. Exiting.")
        sys.exit(0)
        
    # Group bursts by fullBurstID to show geographic coverage frames
    unique_burst_frames = set()
    for r in results:
        burst_meta = r.properties.get("burst")
        if burst_meta:
            unique_burst_frames.add(burst_meta.get("fullBurstID"))
            
    print(f"These products cover {len(unique_burst_frames)} unique burst frames:")
    for frame_id in sorted(unique_burst_frames):
        print(f"  - {frame_id}")
        
    # 2. Authentication check for downloading
    if EARTHDATA_USERNAME == "your_username_here" or EARTHDATA_PASSWORD == "your_password_here":
        print("\n" + "="*80)
        print("  CREDENTIALS REQUIRED TO DOWNLOAD IMAGERY")
        print("="*80)
        print(f"  Search succeeded! We found {total_found} items.")
        print("  To download these bursts (including both TIFF imagery and XML metadata):")
        print(f"  1. Open this file: {os.path.abspath(__file__)}")
        print("  2. Fill in 'EARTHDATA_USERNAME' and 'EARTHDATA_PASSWORD' at the top.")
        print("  3. Run the script again.")
        print("="*80 + "\n")
        sys.exit(0)
        
    print("\nAuthenticating with NASA Earthdata Login...")
    try:
        session = asf.ASFSession().auth_with_creds(EARTHDATA_USERNAME, EARTHDATA_PASSWORD)
    except Exception as e:
        print(f"[ERROR] Authentication failed. Please verify credentials. Details: {e}")
        sys.exit(1)
        
    print(f"Authentication successful! Starting download of {total_found} assets to '{OUT_DIR}'...")
    print("Note: Downloading both TIFF imagery and XML metadata for each burst.")
    
    # 3. Download both image and metadata files
    try:
        results.download(
            path=OUT_DIR,
            session=session,
            fileType=asf.FileDownloadType.ALL_FILES
        )
        print(f"\n[SUCCESS] Finished downloading. Files saved to '{OUT_DIR}'.")
        
        # Crop and organize downloaded bursts
        print("\nCropping and organizing downloaded burst files...")
        crop_downloaded_bursts(OUT_DIR, BBOX)
        
        # Organize into standard PyGMTSAR/insardev directory layout
        annotation_dir = os.path.join(OUT_DIR, TARGET_BURST_ID, "annotation")
        measurement_dir = os.path.join(OUT_DIR, TARGET_BURST_ID, "measurement")
        os.makedirs(annotation_dir, exist_ok=True)
        os.makedirs(measurement_dir, exist_ok=True)
        
        import glob
        import shutil
        for f in glob.glob(os.path.join(OUT_DIR, "*.xml")):
            shutil.move(f, os.path.join(annotation_dir, os.path.basename(f)))
        for f in glob.glob(os.path.join(OUT_DIR, "*.tiff")) + glob.glob(os.path.join(OUT_DIR, "*.tif")):
            shutil.move(f, os.path.join(measurement_dir, os.path.basename(f)))
            
        print("[SUCCESS] All files cropped and organized into subfolders.")
    except Exception as e:
        print(f"\n[ERROR] An error occurred during download: {e}")
        print("Please ensure your Earthdata account has authorized the 'Alaska Satellite Facility' application.")

if __name__ == "__main__":
    main()
