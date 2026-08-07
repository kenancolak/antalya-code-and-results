# Data dictionary

This dictionary defines the principal fields needed to interpret the released result tables. Some CSV files retain additional Earth Engine or GIS export fields for provenance.

## Identifier and endpoint-classification fields

| Field | Meaning | Unit / coding |
|---|---|---|
| `GridID` | Grid-cell identifier used in the thermal result tables | text |
| `run_uid`, `runid` | Analysis-run identifier retained in the original output | text |
| `u17`, `f17`, `a17`, `g17`, `b17` | 2017 fractions assigned to urban, forest, agriculture, greenhouse, or bare surface | 0–1 |
| `u25`, `f25`, `a25`, `g25`, `b25` | Equivalent 2025 class fractions | 0–1 |
| `dom17`, `dom25` | Dominant endpoint class code | 1 urban; 2 forest; 3 agriculture; 4 greenhouse; 5 bare surface |
| `dom17_lbl`, `dom25_lbl` | Short dominant-class label | U, F, A, G, B |
| `purity17`, `purity25` | Dominant-class fraction used by the archived endpoint screen | 0–1 |
| `changeRatio` | Within-cell endpoint class-change metric used by the archived workflow | 0–1 |
| `validFrac17`, `validFrac25` | Valid classified fraction at each endpoint | 0–1 |
| `commonFrac_full` | Fraction valid at both endpoints relative to the full cell | 0–1 |
| `commonFrac_union` | Common-valid fraction relative to the valid union | 0–1 |
| `eligible` | Passes the endpoint data-support conditions in the archived workflow | 0/1 |
| `domSame` | Dominant endpoint class agrees between 2017 and 2025 | 0/1 |
| `keep_default` | Passes the default endpoint screen | 0/1 |
| `core_default` | Default retained-core class | UU, FF, AA, GG, BB; blank for non-retained cells |

## Annual spectral fields

| Field | Meaning | Unit / coding |
|---|---|---|
| `maxdNDVI`, `maxdNDBI`, `maxdTCB` | Maximum absolute consecutive-year change in the named spectral metric | index units |
| `anchor_pct`, `default_anchor_pct` | Percentile used for annual threshold calibration | percentile |
| `thr_ndvi`, `thr_ndbi`, `thr_tcb` | Calibrated threshold values in the calibration table | index units |
| `default_thr_ndvi`, `default_thr_ndbi`, `default_thr_tcb` | Default threshold values attached to grid-level outputs | index units |
| `unstableRatio_default` | Fraction flagged by at least one default annual spectral diagnostic | 0–1 |
| `annualKeep_default` | Passes the default annual urban diagnostic | 0/1 |
| `confirmedStableUrban_default` | Default UU cell retained by the annual diagnostic | 0/1 |
| `nAnnualStableUrban` | Number of cells passing an annual-stability sensitivity combination | count |
| `nIntersection_UU` | Number also belonging to the endpoint UU core set | count |

## Endpoint-threshold sensitivity fields

| Field | Meaning | Unit / coding |
|---|---|---|
| `changeThr` | Endpoint within-cell change threshold represented by the row | proportion |
| `purityThr` | Endpoint dominant-class purity threshold represented by the row | proportion |
| `nRetained`, `nTotal` | Number of retained cells for the represented threshold/scenario | count |
| `nUU`, `nFF`, `nAA`, `nGG`, `nBB` | Retained class counts | count |
| `feasible` | Indicates whether the archived minimum UU/FF sample-size criterion is met | 0/1 |

## Thermal fields

| Field | Meaning | Unit / coding |
|---|---|---|
| `scenario` | Landsat scene-support scenario | text |
| `lst2017`, `lst2025` | Mean endpoint land surface temperature on common valid pixels | °C |
| `dlst` | `lst2025 - lst2017` | °C |
| `commonFrac_landsat` | Fraction with valid thermal support at both endpoints | 0–1 |
| `x_utm`, `y_utm` | Grid-centroid coordinates | metres, EPSG:32636 |
| `lon`, `lat` | Grid-centroid coordinates | decimal degrees |

## Thermal summary and inference fields

| Field | Meaning | Unit / coding |
|---|---|---|
| `mean_dlst_UU`, `mean_dlst_FF`, `mean_dlst_AA`, `mean_dlst_BB` | Mean ΔLST for the retained class under the represented threshold/scenario | °C |
| `sd_dlst_UU`, `sd_dlst_FF`, `sd_dlst_AA`, `sd_dlst_BB` | Standard deviation of class-specific ΔLST | °C |
| `UU_minus_FF` | Difference between the UU and FF mean ΔLST values | °C |
| `welch_se` | Standard error for the Welch UU–FF comparison | °C |
| `welch_df` | Welch–Satterthwaite degrees of freedom | numeric |
| `welch_ci_low`, `welch_ci_high` | Lower and upper limits of the reported Welch confidence interval | °C |
| `welch_p` | Welch comparison p-value | probability |

## Landsat scene metadata

`landsat_scene_metadata.csv` records the metadata retained for the six primary Landsat scenes. Important fields include `date_acquired`, `landsat_product_id`, `landsat_scene_id`, `processing_level`, `spacecraft_id`, `sensor_id`, `scene_center_time_utc`, `cloud_cover`, `wrs_path`, and `wrs_row`.

## Reference-class summary

`reference_class_lulc_lst_summary.csv` contains the archived reference-class LST summaries. Important fields include `reference_name`, `reference_year`, `class_code`, `class_name`, `lst2017_mean_c`, `lst2025_mean_c`, `dlst_mean_c`, `common_valid_fraction`, and `n_landsat_pixels_common`.

## Export-only fields

Fields such as `system:index`, `.geo`, and legacy `Shape_*` attributes are retained from the original exports. They are auxiliary provenance/GIS fields rather than primary analytical variables.
