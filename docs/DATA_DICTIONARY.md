# Data dictionary

## Classification validation

### `classification_validation_confusion_matrices.csv`

Each row is one reference class from an independent, class-balanced endpoint validation sample. The five `predicted_*` columns form the corresponding confusion-matrix row. Row totals are 200; each year contains 1,000 observations.

| Field | Meaning |
|---|---|
| `year` | Endpoint validation year |
| `reference_class_code`, `reference_class_label` | Reference class |
| `predicted_urban` … `predicted_bare` | Counts assigned to each mapped class |
| `reference_total`, `correct` | Row total and diagonal count |
| `producer_accuracy_pct` | Correct divided by the reference-class total |
| `overall_accuracy_pct`, `kappa` | Year-level balanced-sample checkpoints repeated across the five class rows |

## Population and endpoint tables

### `eligible_population_construction.csv`

| Field | Meaning |
|---|---|
| `stage` | Sequential eligibility step |
| `cells_excluded` | Cells removed at the step |
| `cells_remaining` | Cells remaining after the step |

### `endpoint_grid_stats_pre_inland_water_2301.csv`

This table predates the 111-cell inland-water exclusion.

| Field | Meaning | Unit/coding |
|---|---|---|
| `GridID` | Stable analysis-cell identifier | text |
| `u17`, `f17`, `a17`, `g17`, `b17` | 2017 urban, forest, agriculture, greenhouse, and bare fractions | 0–1 |
| `u25`, `f25`, `a25`, `g25`, `b25` | Corresponding 2025 class fractions | 0–1 |
| `dom17`, `dom25` | Dominant class code | 1 U; 2 F; 3 A; 4 G; 5 B |
| `changeRatio` | Changed common-valid pixels divided by all common-valid pixels | 0–1 |
| `purity17`, `purity25` | Persistent dominant-class shares at the endpoints | 0–1 |
| `validFrac17`, `validFrac25` | Valid classification support | 0–1 |
| `commonFrac_full`, `commonFrac_union` | Common classification-support metrics | 0–1 |
| `keep_default` | Passes the default endpoint screen | 0/1 |
| `core_default` | Retained-core class | UU, FF, AA, GG, BB or missing |

### `endpoint_thermal_retained_cores_83.csv`

| Field | Meaning | Unit/coding |
|---|---|---|
| `core_default` | Retained-core group | UU, FF, AA or BB |
| `lst2017`, `lst2025` | Mean matched-support endpoint LST | °C |
| `dlst` | `lst2025 − lst2017` | °C |
| `commonFrac_landsat` | Common valid thermal support | 0–1 |
| `x_utm`, `y_utm` | Cell-centroid coordinates | metres; EPSG:32636 |

## Annual SCL tables

### `annual_scl_core_stats_64.csv`

| Field | Meaning | Unit/coding |
|---|---|---|
| `maxdNDVI`, `maxdNDBI`, `maxdTCB` | Maximum absolute consecutive-year change | index units |
| `annual_pair_count` | Number of fully supported consecutive-year comparisons | count; maximum 8 |
| `annual_full_support` | Share of 10 m pixels with all eight comparisons | 0–1 |
| `unstableRatio_p90`, `unstableRatio_p95`, `unstableRatio_p97` | Cell share exceeding at least one percentile-specific spectral threshold | 0–1 |
| `annual_stable_urban_p90`, `annual_stable_urban_p95`, `annual_stable_urban_p97` | Urban-core retention indicators | 0/1 |

### `annual_scl_thresholds.csv`

| Field | Meaning |
|---|---|
| `anchor_percentile` | Empirical percentile used for threshold calibration |
| `anchor_valid_pixel_count` | Fully supported pixels used in calibration |
| `threshold_maxdNDVI`, `threshold_maxdNDBI`, `threshold_maxdTCB` | Index-specific instability thresholds |
| `threshold_method` | Direct all-valid-pixel percentile method; no random pixel subsample |

### `annual_core_vegetation_support_2017_2025.csv`

`vegetation_support_frac` is the fraction of valid pixels with NDVI ≥0.20. `IVS_year` equals one minus that fraction.

## Acquisition-level tables

### `landsat_scene_meteorology_audit_75.csv`

| Field | Meaning | Unit/coding |
|---|---|---|
| `scene_id` | Landsat product identifier | text |
| `date_utc`, `acquisition_utc` | Acquisition date/time | UTC |
| `sensor` | Landsat 8 or Landsat 9 | text |
| `model_scene_keep_core95` | Meets model scene-eligibility rules | 0/1 |
| `t2m_C_aoi`, `td2m_C_aoi` | Interpolated 2 m air/dew-point temperature | °C |
| `vpd_kPa_aoi` | Vapour-pressure deficit | kPa |
| `wind10_mps_aoi` | 10 m wind speed | m s⁻¹ |
| `ssrd_prev3h_Wm2_aoi` | Antecedent three-hour shortwave radiation | W m⁻² |
| `precip_prev7d_mm_aoi`, `precip_prev30d_mm_aoi` | Antecedent precipitation | mm |
| `sm_root_0_100_ff_background` | Forest-background 0–100 cm soil water | m³ m⁻³ |

### `core_acquisition_panel_60acq_3823obs.csv`

| Field | Meaning |
|---|---|
| `GridID`, `scene_id` | Cell and acquisition identifiers |
| `LST_C` | Cell-level daytime LST for the acquisition |
| `NDVI`, `NDMI` | Acquisition-specific optical diagnostics |
| `core_scene_valid_frac` | Valid thermal coverage within the cell |
| `core_scene_optical_valid_frac` | Valid optical coverage within the cell |
| `block_600`, `block_1000`, `block_1500`, `block_2000` | UTM-aligned spatial-block identifiers |

## Morphology and footprint tables

### `morphology_scl_p95_40.csv`

| Field | Meaning | Unit |
|---|---|---|
| `BCR_2025_pct` | 2025-aligned building-footprint area divided by cell area | % |
| `p75_storeys` | Within-cell 75th percentile of building-level storey counts | storeys |
| `IVS_median_2017_2025` | Median annual inverse vegetation support | 0–1 |
| `delta_LST_C` | Matched-support endpoint LST change | °C |
| `functional_label_en` | Descriptive local functional context | text |

### `manual_footprint_audit_47.csv`

| Field | Meaning | Unit |
|---|---|---|
| `added_m2`, `removed_m2` | Manually audited footprint additions/removals | m² |
| `net_change_m2` | Added minus removed footprint | m² |

### `morphology_footprint_audit_scl_p95_40.csv`

| Field | Meaning | Unit |
|---|---|---|
| `BCR_2025_pct` | Current building coverage ratio | % |
| `delta_BCR_pp` | Net footprint change divided by 90,000 m² | percentage points |
| `BCR_baseline_proxy_pct` | Current BCR minus net BCR change | % |
| `turnover_pp` | Added plus removed area divided by 90,000 m² | percentage points |
| `any_footprint_change` | At least one addition or removal detected | Boolean |

## Geometry and system fields

Some GEE-exported tables retain `.geo`, `system:index`, run identifiers, and duplicate source fields for provenance. They are not independent analytical variables.
