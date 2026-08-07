# Antalya retained-core LST study: code and derived results

This repository accompanies the study **“Reducing Land Cover Ambiguity in Long Term Surface Urban Heat Assessment: A Retained Core Grid Approach in Antalya, Türkiye”** by Kenan Çolak and Ziya Gençel.

It archives the Google Earth Engine (GEE) scripts used to document the analytical workflow and a selected set of derived tables supporting the reported results. The repository is intended for **methodological transparency, inspection of analytical settings, and access to the archived outputs**.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── scripts/
│   ├── 01_landcover_stability.js
│   └── 02_thermal_dlst.js
├── data/
│   └── derived/
│       ├── endpoint_core_grid_stats.csv
│       ├── endpoint_threshold_sensitivity.csv
│       ├── annual_grid_stats.csv
│       ├── annual_threshold_calibration.csv
│       ├── annual_threshold_sensitivity.csv
│       ├── thermal_primary_retained_cores.csv
│       ├── thermal_scenarios_retained_cores.csv
│       ├── endpoint_threshold_dlst_summary.csv
│       ├── landsat_scene_metadata.csv
│       └── reference_class_lulc_lst_summary.csv
└── docs/
    ├── REPOSITORY_SCOPE.md
    ├── ASSET_MANIFEST.md
    └── DATA_DICTIONARY.md
```

## Analysis scripts

### `scripts/01_landcover_stability.js`

Documents the Sentinel-2 land-cover classification workflow, accuracy assessment, endpoint retained-core screening, endpoint-threshold sensitivity, and annual spectral diagnostics.

### `scripts/02_thermal_dlst.js`

Documents the Landsat Collection 2 Level 2 surface-temperature workflow, common-pixel endpoint comparison, thermal-support scenarios, and thermal output tables.

The scripts retain the study asset identifiers, run identifiers, thresholds, scene dates, class codes, and export logic used in the archived workflow.

## Derived result tables

| File | Rows | Contents |
|---|---:|---|
| `endpoint_core_grid_stats.csv` | 2,301 | Endpoint class fractions, data-support metrics, purity/change metrics, and retained-core status |
| `endpoint_threshold_sensitivity.csv` | 24 | Retained-class counts across endpoint change/purity threshold combinations |
| `annual_grid_stats.csv` | 2,301 | Annual spectral diagnostics and annual urban-stability flags |
| `annual_threshold_calibration.csv` | 3 | Percentile-based annual spectral threshold calibration |
| `annual_threshold_sensitivity.csv` | 12 | Annual-stability sensitivity combinations and retained counts |
| `thermal_primary_retained_cores.csv` | 83 | Primary LST results for the default retained-core set |
| `thermal_scenarios_retained_cores.csv` | 332 | Four thermal-support scenarios for the 83 retained cores |
| `endpoint_threshold_dlst_summary.csv` | 96 | Thermal summaries and urban–forest contrasts across endpoint thresholds and thermal scenarios |
| `landsat_scene_metadata.csv` | 6 | Metadata for the six Landsat scenes in the primary three-scene-per-endpoint comparison |
| `reference_class_lulc_lst_summary.csv` | 3 | Reference-class LST summaries used for contextual comparison |

The CSV files are archived analysis outputs. Public file names were shortened for readability; their numeric contents were not recalculated during repository preparation.

## Main analytical settings

- Study period: 2017–2025
- Grid unit: 300 × 300 m
- Projected analysis CRS: EPSG:32636
- Sentinel-2 collection: `COPERNICUS/S2_SR_HARMONIZED`
- Random Forest trees: 500
- Random Forest bag fraction: 0.632
- Random seed: 42
- Endpoint valid-data threshold: 0.99
- Common endpoint support threshold: 0.98
- Endpoint change threshold: 0.10
- Endpoint class-purity threshold: 0.85
- Default annual spectral calibration percentile: 95th percentile
- Landsat WRS path/row: 178/34
- Landsat collections: `LANDSAT/LC08/C02/T1_L2` and `LANDSAT/LC09/C02/T1_L2`

## Archived result checkpoints

The archived tables contain the following principal checkpoints:

- Overall classification accuracy: 0.852 (2017) and 0.866 (2025)
- Cohen's kappa: 0.815 (2017) and 0.8325 (2025)
- Default retained cores: 47 UU, 17 FF, 4 AA, 0 GG, and 15 BB; 83 total
- Urban endpoint cores retained after the annual spectral diagnostic: 43
- Primary mean ΔLST: 2.144 °C for UU and 1.050 °C for FF
- Primary UU − FF mean ΔLST contrast: 1.094 °C
- Welch 95% confidence interval for the primary UU − FF contrast: 0.714 to 1.474 °C
- Welch p-value: 8.50 × 10⁻⁷

## Scope of the archive

The repository does not redistribute public Sentinel-2 or Landsat imagery. Dataset identifiers and processing logic are recorded in the scripts, and the primary Landsat scene metadata are included in `data/derived/landsat_scene_metadata.csv`.

Several study-created inputs referenced by the scripts—including the study boundary, fixed grid, training and validation data, and canonical endpoint land-cover assets—are not included in this archive. The original GEE paths are retained as provenance information. Accordingly, the repository documents the workflow and preserves selected outputs; it is not presented as a self-contained executable replication package.

Additional documentation is provided in:

- `docs/REPOSITORY_SCOPE.md` — archive purpose and boundaries;
- `docs/ASSET_MANIFEST.md` — principal input dependencies; and
- `docs/DATA_DICTIONARY.md` — definitions of the main released fields.

## Citation

Citation metadata are provided in `CITATION.cff`. When a versioned archival DOI is available for this repository, that DOI identifies the corresponding release.

## Contact

**Kenan Çolak**  
Department of City and Regional Planning, Faculty of Architecture  
Akdeniz University, Antalya, Türkiye  
Email: kenancolaktr@gmail.com
