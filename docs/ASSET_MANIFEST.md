# Asset manifest

## Public source collections

| Source | Identifier or product | Role |
|---|---|---|
| Sentinel-2 Level-2A | `COPERNICUS/S2_SR_HARMONIZED` | Endpoint classification and annual SCL-based spectral screening |
| Landsat 8 Level-2 | `LANDSAT/LC08/C02/T1_L2` | Daytime surface temperature and optical diagnostics |
| Landsat 9 Level-2 | `LANDSAT/LC09/C02/T1_L2` | Daytime surface temperature and optical diagnostics |
| ERA5-Land hourly | `ECMWF/ERA5_LAND/HOURLY` | Acquisition-specific meteorology and forest-background hydroclimate |
| Urban Atlas 2018 | DOI `10.2909/fb4dffa1-6ceb-4cc0-8372-1ed354c285e6` | Analysis frame, coverage/water controls, and functional context |
| SRTMGL1 v003 | DOI `10.5067/MEaSUREs/SRTM/SRTMGL1.003` | Grid-mean elevation |
| Microsoft Global Building Footprints | Provider repository | Candidate/current building footprints |
| Esri World Imagery Wayback | Provider service | Baseline high-resolution image audit |

## Study-created or restricted inputs not redistributed

- Study boundary and fixed 300 × 300 m grid geometries
- Training polygons and endpoint-specific validation samples
- Classified endpoint image assets used in the analysis
- Inland-water exclusion layer used in ArcGIS
- Antalya Metropolitan Municipality institutional building/base-map data
- Source Microsoft and manually reconciled building-footprint geometries
- High-resolution reference imagery and street-level imagery
- ArcGIS deep-learning candidate-footprint layer

## Public derived replacements

The `data/derived/` directory provides tabular outputs needed to inspect the reported results. These tables do not grant or imply redistribution rights for the underlying source imagery or geometries.

`classification_validation_confusion_matrices.csv` preserves the independently reported historical endpoint validation counts without redistributing the validation-point geometries or reference imagery.

## GEE identifiers

Private GEE identifiers remain visible in the scripts to document the processing chain. Upstream processing requires authorised access to those assets or replacement inputs. The endpoint image identifiers are used for retained-core screening; the optional classifier branch in Script 01 runs against the current catalogue and is intended as a diagnostic.
