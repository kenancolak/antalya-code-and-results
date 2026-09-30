# Changelog

## 2.0.1 — 2026-09-29

- Separated endpoint screening from the optional current-catalogue Random Forest diagnostic in GEE Script 01.
- Added the archived 2017 and 2025 endpoint confusion matrices and validation checks.
- Suppressed the memory-intensive threshold-by-support Console preview in GEE Script 02 while retaining its export definition.
- Replaced the sampled annual-percentile calculation with a direct 74,431-pixel calculation in GEE Script 03 and forced raw percentile evaluation with `maxRaw = 100000`.
- Added Windows execution instructions and clarified the roles of the R, Python, and GEE scripts.
- Added the results of the GEE checks for the endpoint populations, support scenarios, annual SCL screens, and 2017 year-split exports.

## 2.0.0 — 2026-09-29

- Updated the repository title and documentation to match the revised article.
- Replaced the historical 43-cell annual morphology checkpoint with the final, temporally consistent SCL-p95 population of 40 urban cores.
- Added the 75-acquisition Landsat/ERA5-Land audit and the final 60-acquisition, 3,823-observation core–acquisition panel.
- Added the weather-adjusted fixed-effects model and its principal sensitivity outputs.
- Added the manual footprint audit for all 47 retained urban cores and the baseline-proxy BCR analyses for the 40-cell morphology population.
- Added public R and Python analysis scripts with repository-relative paths.
- Renamed the 2,301-row endpoint grid table to state explicitly that it predates the 111-cell inland-water exclusion.
- Updated the data dictionary, asset manifest, archive scope, licences, and validation checks.

## 1.0.0 — 2026-08-05

- Initial public archive of the endpoint classification, retained-core, and Landsat ΔLST workflow.
