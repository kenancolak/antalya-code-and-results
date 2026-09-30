# Antalya retained-core LST study: code and derived results

This repository accompanies the study **“Urban–Forest Land Surface Temperature Change under Stability and Purity Constraints: Antalya, 2017–2025”** by Kenan Çolak and Ziya Gençel.

The repository contains the Google Earth Engine (GEE) workflow, an independent R check of the endpoint statistics, Python scripts for the released-table analyses, and the derived tables used for the endpoint, meteorological, annual-stability, morphology, and footprint-audit results.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── CHANGELOG.md
├── VERSION
├── LICENSE_CODE
├── LICENSE_DATA
├── requirements.txt
├── scripts/
│   ├── gee/
│   ├── r/
│   └── python/
├── data/
│   └── derived/
├── results/
└── docs/
    ├── REPOSITORY_SCOPE.md
    ├── ASSET_MANIFEST.md
    ├── DATA_DICTIONARY.md
    ├── RUNNING_THE_CODE.md
    └── VALIDATION_REPORT.md
```

## Analysis sample

- Initial rectangular grid: 3,060 cells
- Eligible terrestrial population after Urban Atlas, coverage, sea, and inland-water exclusions: 2,190 cells
- Endpoint-retained cores: 47 urban (UU), 17 forest (FF), 4 agriculture (AA), 15 bare surface (BB), and 0 greenhouse (GG); 83 total
- Final SCL-p95 morphology population: 40 UU cells
- Landsat audit: 75 June–August acquisitions during 2017–2025
- Meteorological model: 60 acquisitions, 64 cores, and 3,823 core–acquisition observations

`endpoint_grid_stats_pre_inland_water_2301.csv` contains the 2,301 cells remaining after the Urban Atlas, coverage, and sea exclusions. The ArcGIS inland-water overlay subsequently removed 111 cells. None of the 83 retained cores intersected those additional exclusions. The file name records this intermediate stage; the final eligible terrestrial population contains 2,190 cells.

## Main analytical settings

- Study period: 2017–2025
- Grid unit: 300 × 300 m
- Projected analysis CRS: EPSG:32636
- Sentinel-2 collection: `COPERNICUS/S2_SR_HARMONIZED`
- Endpoint classification predictors: B2, B3, B4, B5, B6, B7, B8, B11, B12, NDVI, NDBI, NDWI, BSI, and TCB
- Random Forest trees: 500
- Random Forest bag fraction: 0.632
- Random seed: 42
- Endpoint valid-data threshold: 0.99
- Common endpoint support threshold: 0.98
- Endpoint change threshold: 0.10
- Endpoint class-purity threshold: 0.85
- Final annual cloud mask: SCL classes 2, 4, 5, and 6
- Annual spectral calibration: empirical p95 thresholds from 74,431 fully supported pixels
- Annual unstable-pixel limit: 0.10
- Landsat WRS path/row: 178/034
- Landsat collections: `LANDSAT/LC08/C02/T1_L2` and `LANDSAT/LC09/C02/T1_L2`

## Results reported in the article

- Overall classification accuracy: 85.2% in 2017 and 86.6% in 2025
- Cohen’s kappa: 0.815 in 2017 and 0.833 in 2025
- Mean endpoint ΔLST: 2.144 °C for UU and 1.050 °C for FF
- Raw UU–FF mean ΔLST contrast: 1.094 °C
- Cell-level Welch 95% CI: 0.714–1.474 °C; p = 8.50 × 10⁻⁷
- Primary weather-adjusted divergence: 0.791 °C
- Weather-adjusted 95% CI: 0.037–1.546 °C; p = 0.040
- Current BCR–ΔLST Spearman correlation in the 40-cell SCL-p95 population: ρ = 0.569
- Baseline-proxy BCR–ΔLST correlation: ρ = 0.547
- Net BCR-change correlation: ρ = −0.042
- Current BCR–ΔLST correlation among 17 cells with no detected footprint change: ρ = 0.738

The archived endpoint confusion matrices underlying the accuracy and kappa checkpoints are provided in `data/derived/classification_validation_confusion_matrices.csv`.

## What can be run from this repository

- The Python scripts rerun the endpoint, meteorological, morphology, and footprint-audit calculations from the released derived tables.
- The base-R script provides an independent check of the cell-level UU–FF endpoint statistics.
- The GEE scripts document the upstream geospatial processing and retain the private asset identifiers used in the study. Those stages require authorised access to the non-redistributed inputs.
- Script 01 also includes an optional classifier run against the current Sentinel-2 catalogue. It is provided as a diagnostic because catalogue changes mean that it may not reproduce the original map/validation pairing. The original confusion matrices are included as a derived table.

## Execution

The Python reproductions can be run directly from the repository root:

```bash
python -m pip install -r requirements.txt
python scripts/python/01_validate_release.py
python scripts/python/02_weather_adjusted_models.py
python scripts/python/03_morphology_analysis.py
python scripts/python/04_footprint_audit_analysis.py
```

The independent base-R endpoint check can be run with:

```bash
Rscript scripts/r/01_endpoint_statistics.R
```

The GEE scripts require access to non-redistributed study-created assets. Their identifiers are retained as provenance information, but public access to the repository does not confer access to those assets.

Detailed Windows PowerShell commands and expected checkpoints are provided in `docs/RUNNING_THE_CODE.md`.

## Checks performed before release

The release was checked on 29 September 2026. The Python scripts returned 83 retained cores, the 40-cell morphology population, 60 model acquisitions, and 3,823 observations. Runs made with access to the study GEE assets returned:

- 47 UU, 17 FF, 4 AA, 15 BB, and 0 GG endpoint cores;
- 83 default retained-core rows and 332 rows across four endpoint-support scenarios;
- 75 JJA Landsat 8/9 candidate acquisitions;
- 74,431 fully supported annual-calibration pixels;
- 13, 40, and 44 UU cores under the SCL-p90, p95, and p97 screens; and
- six 2017 scene-meteorology rows and 384 unfiltered 2017 core–scene rows in the year-split export test.

## Archive scope

The repository does not redistribute Sentinel-2, Landsat, ERA5-Land, SRTM, Urban Atlas, Microsoft Global Building Footprints, or Esri World Imagery Wayback source data. It also does not redistribute training and validation geometries, high-resolution reference imagery, restricted municipal building data, source building-footprint geometries, or the original boundary and water-mask layers.

The released tables are derived analysis outputs. They allow the reported statistical analyses to be rerun, but the repository is not a self-contained copy of every upstream geospatial input.

## Citation

Citation metadata are provided in `CITATION.cff`. The all-versions Zenodo DOI is:

`https://doi.org/10.5281/zenodo.21836778`

## Licences

- Code in `scripts/`: MIT License (`LICENSE_CODE`)
- Derived tables and documentation: CC BY 4.0 (`LICENSE_DATA`)

Third-party source datasets remain governed by their respective providers’ terms.

## Contact

**Kenan Çolak**  
Department of City and Regional Planning, Faculty of Architecture  
Akdeniz University, Antalya, Türkiye  
Email: kenancolaktr@gmail.com
