# Release checks

Validation was completed on 29 September 2026 before packaging version 2.0.1.

## Local released-table checks

- Classification validation: 1,000 observations per endpoint; 852 correct in 2017 and 866 correct in 2025.
- Endpoint population: 47 UU, 17 FF, 4 AA, 15 BB, and 0 GG; 83 total.
- Annual screening: 13, 40, and 44 UU cores at p90, p95, and p97.
- Annual threshold calibration: 74,431 fully supported pixels.
- Acquisition panel: 60 acquisitions, 64 cores, and 3,823 observations.
- Morphology population: 40 SCL-p95 UU cores.
- No-footprint-change subset: 17 UU cores.

The four Python scripts completed successfully in a fresh Windows Python 3.12 virtual environment. The primary weather-adjusted model reproduced 0.791474 °C (95% CI 0.036583–1.546365; p = 0.040443).

## GEE checks

| Script | Result checked |
|---|---|
| 01 | The analysis endpoint assets are used for retained-core screening; current-catalogue classification is labelled diagnostic. |
| 02 | 47 UU, 17 FF, 4 AA, 15 BB, 0 GG; 83 default rows; 332 retained-core rows across four supports; 9,204 full-grid support rows. |
| 03 | 47 UU, 17 FF; 75 candidate acquisitions; nine annual S2 composites; 576 annual UU/FF vegetation-support rows; 74,431 calibration pixels; 13/40/44 p90/p95/p97 UU counts. |
| 04 | 47 UU, 17 FF, and 75 acquisitions; 2017 test exports completed successfully. |

The 2017 GEE export test produced six scene-meteorology rows and 384 unfiltered core–scene rows. Every scene contained 47 UU and 17 FF rows. The six scene identifiers and all overlapping analytical values matched the archived audit/panel data; two core–scene rows are absent from the final 60-acquisition panel because the later row-validity filter was applied.

## Scope of these checks

The private training and validation geometries and the endpoint rasters are not redistributed. The released confusion matrices preserve the reported validation counts. Because the optional classifier run uses the current catalogue, its results can differ from the original map/validation pairing.
