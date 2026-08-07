# Repository scope

## Purpose

This repository archives the principal Google Earth Engine code and selected derived result tables associated with the Antalya retained-core land-surface-temperature study.

Its purpose is to make the analytical logic, parameter choices, dataset identifiers, and reported output tables inspectable in a stable public repository.

## Included

- two GEE JavaScript workflows;
- analysis parameters and dataset identifiers recorded in those scripts;
- selected grid-level and summary result tables;
- primary Landsat scene metadata;
- definitions of the principal released fields; and
- citation metadata for the repository.

## Outside the archive

The repository does not include:

- private study-created boundary and grid assets;
- training and validation geometries;
- study-created endpoint land-cover image assets;
- redistributed Sentinel-2 or Landsat imagery;
- every intermediate export produced during the analysis; or
- external GIS/statistical project files not required to interpret the archived tables.

## Interpretation

The JavaScript files are archived methodological code documenting the workflow used in the study. The CSV files are archived outputs associated with the reported analyses.

Because several study-created GEE inputs are outside the archive, the repository is not presented as a self-contained executable replication package. Private GEE asset paths are retained in the scripts only as provenance information.
