# Analysis scripts

## Google Earth Engine

1. `gee/01_endpoint_landcover_retained_cores.js` uses the archived endpoint images for retained-core construction and endpoint-threshold sensitivity. It also contains an optional classifier run against the current Sentinel-2 catalogue; the original confusion matrices are supplied separately.
2. `gee/02_endpoint_landsat_dlst.js` documents the matched-pixel endpoint LST workflow and S0–S3 support scenarios.
3. `gee/03_scene_meteorology_annual_scl.js` creates the 75-acquisition audit, ERA5-Land variables, forest diagnostics, annual vegetation-support series, and SCL-based p90/p95/p97 screens calculated from all fully supported pixels.
4. `gee/04_core_acquisition_exports.js` exports the year-split Landsat/ERA5-Land core–acquisition tables used to construct the public 60-acquisition panel.

The GEE scripts retain private asset identifiers for provenance. They cannot be rerun without authorised access to the inputs listed in `docs/ASSET_MANIFEST.md`.

## R

`r/01_endpoint_statistics.R` is an independent base-R check of the cell-level Welch test, Wilcoxon check, and Hedges’ g for the 47 UU and 17 FF cores.

## Python

Run from the repository root:

```bash
python -m pip install -r requirements.txt
python scripts/python/01_validate_release.py
python scripts/python/02_weather_adjusted_models.py
python scripts/python/03_morphology_analysis.py
python scripts/python/04_footprint_audit_analysis.py
```

The Python scripts read the released derived tables and write new outputs under `results/reproduced_*`. They rerun the endpoint, weather-adjusted, morphology, and footprint-audit analyses, but not the private upstream geospatial classification.

See `../docs/RUNNING_THE_CODE.md` for Windows commands, GEE flags, and expected checkpoints.
