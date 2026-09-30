# Running the released code

Run commands from the repository root. Generated reproduction outputs are written under `results/reproduced_*` and are excluded from version control.

## Python on Windows PowerShell

```powershell
python --version
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts\python\01_validate_release.py
.\.venv\Scripts\python.exe scripts\python\02_weather_adjusted_models.py
.\.venv\Scripts\python.exe scripts\python\03_morphology_analysis.py
.\.venv\Scripts\python.exe scripts\python\04_footprint_audit_analysis.py
```

The validation script should report 1,000 validation observations per endpoint, 83 retained cores, 40 morphology cells, 60 acquisitions, and 3,823 core–acquisition observations.

The primary weather-adjusted model should reproduce an estimate of approximately 0.791474 °C, with a 95% confidence interval of 0.036583–1.546365 °C and p = 0.040443.

## Base R endpoint check on Windows PowerShell

No contributed R packages are required.

```powershell
Rscript --version
Rscript .\scripts\r\01_endpoint_statistics.R
```

If `Rscript` is not on the system path:

```powershell
$rscript = Get-ChildItem "C:\Program Files\R" -Filter Rscript.exe -Recurse | Sort-Object FullName -Descending | Select-Object -First 1 -ExpandProperty FullName
& $rscript --version
& $rscript .\scripts\r\01_endpoint_statistics.R
```

Expected rounded results are 2.144 °C for UU, 1.050 °C for FF, a difference of 1.094 °C, a Welch 95% confidence interval of 0.714–1.474 °C, p = 8.50 × 10⁻⁷, and Hedges’ g = 1.386.

## Google Earth Engine

The GEE scripts require authorised access to the private assets listed in `ASSET_MANIFEST.md`.

Run the scripts in numerical order:

1. `01_endpoint_landcover_retained_cores.js`
2. `02_endpoint_landsat_dlst.js`
3. `03_scene_meteorology_annual_scl.js`
4. `04_core_acquisition_exports.js`

Script 01 uses the archived endpoint images for retained-core screening. `RUN_CURRENT_CATALOG_RECONSTRUCTION_DIAGNOSTIC` is `false` by default because a new classification against the current catalogue does not recreate the original validation run. The original confusion matrices are supplied as a derived CSV.

Script 02 suppresses the large threshold-by-support Console preview by default. Set `PRINT_THRESHOLD_THERMAL_ROWS` to `true` only if the large server-side collection must be inspected; the export definition remains available.

Script 03 uses all 74,431 fully supported anchor pixels and sets `maxRaw = 100000`, forcing direct percentile calculation rather than a histogram approximation. Expected retained UU counts are 13 at p90, 40 at p95, and 44 at p97.

Script 04 defaults to `EXPORT_YEARS = [2017]` for a two-task test. A successful test produces six scene-meteorology rows and 384 unfiltered core–scene rows. Change the year list only when the remaining year-split exports are required.

Do not start GEE export tasks unless new Drive outputs are intentionally required.
