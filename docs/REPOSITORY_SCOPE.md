# Repository scope

## Purpose

This release contains the analytical code and derived outputs used for the revised study. The statistical analyses can be rerun from the released tables, while the upstream geospatial work can be inspected in the GEE scripts.

## What can be reproduced directly

Using only the released files, a reader can:

- inspect the archived endpoint confusion matrices and verify the reported overall accuracy and kappa values;
- verify the 47 UU and 17 FF endpoint populations and the raw 1.094 °C contrast;
- reproduce the cell-level endpoint statistics;
- inspect all 75 audited Landsat acquisitions and the 60-acquisition model sample;
- rerun the acquisition-level meteorological models from 3,823 observations;
- verify the final 40-cell SCL-p95 morphology population;
- rerun the primary BCR, p75-storey, and IVS analyses;
- rerun the manual footprint-audit correlations and HC3 model; and
- inspect the reported threshold, acquisition-support, clustering, and leave-one-out outputs.

## What cannot be reproduced from this archive alone

The complete upstream geospatial workflow cannot be rerun without non-redistributed study-created and third-party inputs. These include the training and validation geometries, original study boundary and water masks, restricted municipal building data, source building-footprint geometries, and high-resolution reference imagery.

The GEE scripts preserve the original asset identifiers and processing settings as provenance. Access to those identifiers is not equivalent to public access to the underlying assets.

Script 01 separates retained-core screening from an optional Random Forest diagnostic. Screening uses the endpoint image assets that produced the released tables. The optional branch runs against the current Sentinel-2 catalogue and may differ from the original map/validation pairing. The original confusion matrices are included as a derived table.

## Code coverage

The base-R script independently checks the endpoint UU–FF statistics. The Python scripts rerun the released-table statistical chain, including the weather-adjusted models, morphology analysis, and footprint audit. They do not replace the non-redistributed upstream geometries and imagery.

## The 2,301- and 2,190-cell distinction

`endpoint_grid_stats_pre_inland_water_2301.csv` represents the grid after Urban Atlas/coverage and sea exclusions. A separate conservative ArcGIS inland-water overlay then removed 111 cells, producing the final 2,190-cell eligible terrestrial population. None of the 83 retained cores intersected the additional inland-water exclusions. Therefore, the step changed the eligible-population denominator but not retained-core membership or the reported thermal analyses.

## Outcome and inference limits

The outcome is daytime Landsat radiometric surface temperature at satellite overpass. The released results concern criterion-defined retained populations, not Antalya as a whole. They do not establish annual temperature trends, night-time behaviour, thermal comfort, health effects, or causal effects of building coverage.
