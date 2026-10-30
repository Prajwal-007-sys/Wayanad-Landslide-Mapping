# Bundle for InSAR Processing

This bundle contains the necessary notebook, scripts, data, and libraries to run the Mundakkai Hill (Wayanad, Kerala - 2024)Subsidence analysis.

## Structure
- `Imperial_Valley_Subsidence.ipynb` - Main notebook
- `stac_asf.py` - Helper script to query/download data from ASF
- `output_asf_pre/` - Pre-disaster data directory (before July 30, 2024 — contains DEM, orbits, and Sentinel-1 SLC bursts)
- `output_asf_post/` - Post-disaster data directory (after July 30, 2024)
- `InSARdev/` - Local modified python packages

## Setup Instructions

1. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Install Local InSARdev packages**:
   These packages contain custom local edits. Install them from the local folder:
   ```bash
   pip install ./InSARdev/insardev_toolkit
   pip install ./InSARdev/insardev_pygmtsar
   pip install ./InSARdev/insardev
   ```

3. **Start Jupyter**:
   ```bash
   jupyter notebook Imperial_Valley_Subsidence.ipynb
   ```
