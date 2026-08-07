# Input asset manifest

This manifest records the principal inputs referenced by the archived scripts. Study-created GEE assets are listed to document the input chain; they are not distributed in this repository.

## Study-created GEE assets

| Asset path in archived scripts | Type | Role | Included here |
|---|---|---|---|
| `projects/tikkenan/assets/sinir` | FeatureCollection | Study boundary | No |
| `projects/tikkenan/assets/denemefish5` | FeatureCollection | Fixed 300 × 300 m grid | No |
| `projects/tikkenan/assets/train3a_balanced_canonical_v01` | FeatureCollection | Canonical balanced training data | No |
| `projects/tikkenan/assets/valpoints200_17` | FeatureCollection | 2017 validation data | No |
| `projects/tikkenan/assets/valpoints200_25` | FeatureCollection | 2025 validation data | No |
| `projects/tikkenan/assets/final_lulc2017raw` | Image | Canonical 2017 endpoint land-cover map referenced by the thermal script | No |
| `projects/tikkenan/assets/final_lulc2025raw` | Image | Canonical 2025 endpoint land-cover map referenced by the thermal script | No |

## Public Earth Engine datasets

| Dataset | GEE collection | Role |
|---|---|---|
| Sentinel-2 Level-2A harmonised surface reflectance | `COPERNICUS/S2_SR_HARMONIZED` | Endpoint land-cover classification and annual spectral diagnostics |
| Landsat 8 Collection 2 Level 2 | `LANDSAT/LC08/C02/T1_L2` | Endpoint LST |
| Landsat 9 Collection 2 Level 2 | `LANDSAT/LC09/C02/T1_L2` | 2025 thermal support and sensitivity scenarios |

Copies of the public satellite products are not redistributed. The scripts record the collection identifiers and processing operations; `data/derived/landsat_scene_metadata.csv` contains the primary Landsat scene metadata archived with the study outputs.
