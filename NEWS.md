# scop

## scop 0.1.3

*   **New Functionality**: Added `TACSPlot()` for creating FACS-like plots.

## scop 0.0.9

*   **Bug Fixes**: Corrected an issue in `PrepareEnv()`. The default Python version is now set to `3.10`, which ensures that Python-dependent functions like `RunPAGA()` and `RunSCVELO()` run correctly.

## scop 0.0.6

*   **New Functionality**: Introduced `GetAssayData5()`, a reimplementation of `GetAssayData()`, for compatibility with Seurat v5 `Assay` objects.
*   **Enhancements**: Updated `GetFeaturesData()` and `AddFeaturesData()` to support retrieving and adding feature metadata.

## scop 0.0.5

*   **Data**: Updated the `pancreas_sub` and `panc8_sub` test datasets to the Seurat v5 object format.
