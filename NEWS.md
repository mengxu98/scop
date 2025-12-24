# scop

# scop 0.7.6

* **func**:
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.5

* **bugs**:
  * `CellScoring()`: Fixed bug where the function failed to build results. This issue reported in [#98](https://github.com/mengxu98/scop/issues/98).

* **func**:
  * `RunDEtest()`: Fixed compatibility issue with [SeuratObject](https://satijalab.github.io/seurat-object/) 5.0.0+ by replacing deprecated `Assays()` `slot` argument with `LayerData()`. This issue reported in [#100](https://github.com/mengxu98/scop/issues/100).
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.3

* **func**:
  * `RunCellTypist()`: New function for cell type annotation using the CellTypist method.
  * `CellTypistModels()`: New function for downloading and managing CellTypist pre-trained models.

# scop 0.7.2

* **func**:
  * `RunCellRank()`: Performance optimizations and code improvements.

# scop 0.7.1

* **bugs**:
  * `CellDimPlot()`: Fixed issue where NA values appeared in labels. This issue reported in [#93](https://github.com/mengxu98/scop/issues/93).

# scop 0.7.0

* **func**:
  * `PrepareEnv()`: Integrated `uv` as the primary *Python* package installer for improved installation speed.
  * `check_python()`: Now uses `uv` as the primary installation tool with `pip` as fallback, significantly improving package installation speed.
  * Added `find_uv()` and `install_uv()` internal functions for managing `uv` package manager installation and detection.

# scop 0.6.6

* **docs**:
  * Unified documentation format across all R functions:
    * Standardized return value tags: Changed all `@returns` to `@return` for consistency.
    * Unified parameter documentation: Replaced all `\code{value}` with Markdown backticks `` `value` `` format.
    * Standardized default value descriptions.
    * Added `@md` tags: Added `@md` tags to all functions using Markdown syntax in documentation.
    * Enhanced cross-references: Added `@seealso` links to related functions where appropriate.

# scop 0.6.5

* **func**:
  * `PrepareEnv()`:
    * Added comprehensive environment variable configuration to prevent crashes when calling *Python* functions, including setting thread limits for OMP, OPENBLAS, MKL, NUMBA, and other libraries. This improves stability on all platforms, especially Apple silicon Macs.
    * Added `accept_conda_tos()` function to automatically accept conda Terms of Service for required channels, improving the conda environment setup process.
    * Fixed conda Terms of Service acceptance issue in `PrepareEnv()`. The function now automatically accepts conda Terms of Service for required channels, eliminating the need for manual acceptance. This addresses the issue reported in [#85](https://github.com/mengxu98/scop/issues/85).
  * Multiple Python-based functions (`RunPAGA`, `RunSCVELO`, `RunPalantir`, `RunCellRank`, `RunWOT`, `RunPHATE`, `RunPaCMAP`, `RunTriMap`): Enhanced message formatting and code improvements.
  * `PrepareSCExplorer()`: Fixed package version dependency issues with `shiny` and `bslib` compatibility. The function now properly handles `bslib` theme configuration to work with both `shiny` 1.6.0 and 1.7.0+, addressing compatibility errors reported in [#87](https://github.com/mengxu98/scop/issues/87).

* **refactor**:
  * Improved code formatting and consistency across multiple functions.
  * Enhanced Python functions in `inst/python/functions.py` with better error handling and message formatting.

* **docs**:
  * Updated documentation for multiple functions to reflect code improvements.

# scop 0.6.2

* **func**:
  * `CellChatPlot()`: Adjusted the size of saved figures for better file size optimization.

* **docs**:
  * Updated README.md to remove references to Monocle2 and Monocle3 (deprecated functions).

# scop 0.6.1

* **func**:
  * `PrepareEnv()`: Improved message formatting and simplified log output for better user experience.
  * Added `get_conda_envs_dir()` helper function to centralize conda environment directory retrieval.
  * `integration_scop()`: Enhanced `integration_method` parameter definition with explicit method list for better code clarity.

* **refactor**:
  * Moved `exist_python_pkgs()` function to `check_package.R` for better code organization.
  * Replaced direct `conda_info()$envs_dirs[1]` calls with `get_conda_envs_dir()` helper function for consistency.
  * `RunSCExplorer()`: Updated to use `thisplot::palette_list` and `thisplot::slim_data()` instead of `scop::palette_list` and `scop::slim_data()`.
  * Added `thisplot` to dependency checks in `RunSCExplorer()`.

* **docs**:
  * Updated documentation across multiple functions.

# scop 0.6.0

* **func**:
  * `PrepareEnv()`: Enhanced with environment caching mechanism to avoid redundant environment preparation. Improved message formatting and error handling.
  * Python-based functions (`RunPAGA()`, `RunSCVELO()`, `RunPalantir()`, `RunCellRank()`, `RunWOT()`) now automatically call `PrepareEnv()` internally, eliminating the need for users to manually prepare the Python environment before using these functions.
  * `cluster_within_group2()`: New function for clustering within groups.
  * Multiple plotting functions: Replaced `geom_sankey()` with `ggsankey::geom_sankey()` for better Sankey diagram support.
  * Multiple functions: Replaced `:::` operator with `get_namespace_fun()` for safer namespace access.

* **refactor**:
  * Removed `RunMonocle()` function and related documentation (`RunMonocle2.Rd`, `RunMonocle3.Rd`).
  * Removed `projection_functions.R` file (functions moved to other locations).
  * Replaced custom theme functions with `thisplot::theme_this()` (exported as `theme_scop()`).
  * Replaced direct `log_message()` calls with `thisutils::log_message()` for consistency.
  * Removed `palette_list` data object.

* **deps**:
  * Moved `cli` from `Suggests` to `Imports` for better message formatting support.
  * Added `thisplot` to `Imports` for theme and utility functions.
  * Added `ggsankey` to `Suggests` for Sankey diagram support.
  * Added remote dependencies: `theislab/destiny`, `mengxu98/thisplot`.
  * Removed unused dependencies: `Biobase`, `BiocGenerics`, `concaveman`, `DDRTree`, `glmGamPoi`, `hexbin`, `monocle`, `png`, `ragg`, `tidyr`.

* **docs**:
  * Updated documentation across multiple functions to reflect code refactoring.
  * Improved code organization and maintainability.

# scop 0.5.5

* **bugs**:
  * Fixed `VelocityPlot()` function error in `plot_type = "grid"` mode: replaced vectorized arrow length with fixed-length arrows (using mean length) to resolve `vapply()` error that occurred when `grid::arrow()` received a vector instead of a single value, see [#72](https://github.com/mengxu98/scop/issues/72), [#74](https://github.com/mengxu98/scop/issues/74).

# scop 0.5.4

* **bugs**:
  * Fixed parameter name error in `CheckDataType()` function calls: changed `data` parameter to `object` in `RunKNNMap()`, `RunSymphonyMap()`, `RunScmap()`, `RunPCAMap()`, and `RunSingleR()` functions, see [#68](https://github.com/mengxu98/scop/issues/68).
  * Fixed `SingleCellExperiment` object creation in `RunScmap()` and `RunSingleR()` functions: changed from coercing `SummarizedExperiment` to directly constructing `SingleCellExperiment` objects.

# scop 0.5.3

* **func**:
  * `PrepareDB()`: Changed default `Ensembl_version` parameter from `103` to `NULL` for more flexible version handling.
  * Added *Python* version `log_message()` for Python-based functions (`RunSCVELO()`, `RunPAGA()`, `RunPalantir()`, `RunCellRank()`, `RunWOT()`) and added `verbose` parameter inheritance and improved message formatting using cli-style formatting.

* **refactor**:
  * Delete `harmonizomeapi.py` file.
  * Move `scop_analysis.py` into a single `functions.py` file in `inst/python/` for better code organization and maintainability.

* **docs**:
  * Improved parameter documentation consistency.

# scop 0.5.1

* **docs**:
  * Improved reference formatting and consistency across multiple functions.
  * Enhanced documentation clarity and readability.

# scop 0.5.0

* **func**:
  * `RunCellChat()`: New function to perform CellChat analysis for investigating cell-to-cell communication with support for human, mouse, and zebrafish species.
  * `CellChatPlot()`: New function to visualize CellChat analysis results with various plot types and customization options.
  * Multiple integration functions: Improved error messages and message formatting for better user experience.

* **deps**:
  * Added [CellChat](https://github.com/jinworks/CellChat) package dependency with remote repository `jinworks/CellChat`.

* **docs**:
  * Updated README.md with improved code formatting and examples.
  * Enhanced documentation for cell communication analysis functions.
  * Improved error messages and user guidance across integration functions.

* **refactor**:
  * Removed some example figures to optimize package installation size.

# scop 0.4.0

* **func**:
  * `RunProportionTest()`: New function to perform Monte-carlo permutation test for quantifying cell proportion differences between conditions.
  * `ProportionTestPlot()`: New function to generate proportion test plots with customizable significance thresholds and visualization options.
  * Multiple *Python*-based functions: add `\dontrun{}` blocks for Github workfolw checking.

* **docs**:
  * Added comprehensive documentation for new proportion testing functions.
  * Enhanced example documentation across multiple functions.
  * Updated package documentation and examples.

# scop 0.3.4

* **docs**:
  * Updated workflow examples and function documentation.

# scop 0.3.3

* **func**:
  * Multiple functions: Improved parameter documentation formatting and consistency across the package.

# scop 0.3.2

* **func**:
  * `GetFeaturesData()` and `AddFeaturesData()`: Enhanced argument clarity, added input validation, and standardized return values for `Seurat`, `Assay`, and `Assay5` objects.
  * `CellCorHeatmap()`: 
    - Renamed parameters: `query_cell_annotation` → `query_annotation`, `ref_cell_annotation` → `ref_annotation`.
    - Improved error message formatting using cli-style formatting.
    - Simplified variable assignments and improved readability.

* **docs**:
  * Comprehensive documentation updates across multiple functions including `AnnotateFeatures`, `CellDimPlot`, `CellStatPlot`, `FeatureStatPlot`, `GroupHeatmap`, `RunCellQC`, and others.
  * Improved parameter descriptions and function clarity.

# scop 0.3.1

* **func**:
  * `EnrichmentPlot()` and `GSEAPlot()`: Removed conditional font face styling (`face = ifelse()` logic) for better text rendering consistency. Set the default value of `lineheight` from `0.5` to `0.7`.
  * Updated `check_r()` function for improved package checking functionality.
  * Updated reexports functionality.

* **docs**:
  * Updated documentation formatting and consistency.

# scop 0.3.0

* **func**:
  * Fixed `segmentation faults` and `R crashes` on *M-series* MacBook when running *Python* functions.
  * `RunPAGA()`: Enhanced with *M-series* MacBook detection and automatic environment configuration.
  * `RunSCVELO()`: Added ARM64-specific optimizations to prevent crashes and ensure stable execution.
  * `RunCellRank()`: Implemented *M-series* compatibility with proper NUMBA configuration.
  * `RunPalantir()`: Added ARM64 support with single-threaded execution mode.
  * `RunWOT()`: Enhanced with *M-series* MacBook environment variable settings.
  * `RunTriMap()`: Added *M-series* MacBook compatibility for dimensionality reduction.
  * `RunPaCMAP()`: Implemented ARM64-specific environment configuration.
  * `RunPHATE()`: Added *M-series* MacBook support for non-linear dimensionality reduction.
  * `RunCellQC()`: Enhanced both scrublet and doubletdetection functions with ARM64 compatibility.

* **docs**:
  * Updated function documentation to reflect *M-series* MacBook compatibility.
  * Added technical notes about ARM64 architecture considerations.

# scop 0.2.9

* **bugs**:
  * Fix bug for `RunSCVELO()`.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.7

* **func**:
  * Added an internal function `.check_pkg_status()` to check if an *R* package is installed.
  * Update function `CheckDataType()` to *S4* class function.
  * Update function `standard_scop()`, make it more efficient.

* **data**:
  * Delete `lifemap` data, including: `lifemap_cell`, `lifemap_compartment` and `lifemap_organ`.
  * Reconstructed the sample data for both `panc8_sub` and `pancreas_sub`, retaining only the basic `Seurat` object.

# scop 0.2.6

* **func**:
  * Added `remove_r()` function for easy remove *R* packages.
  * Rename function: `RemovePackages()` to `remove_python()`.
  * Removed other methods of installing *R* packages from the `check_r()` function, only retaining [pak::pak](https://pak.r-lib.org/reference/pak.html). 
  * Delete useless import packages: `BBmisc`, `BiocManager`, `covr`, `devtools`, `promises` and `withr`.
  * Optimize the structure of `_pkgdown.yml` file.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.5

* **func**:
  * Rename function: `palette_scop()` to `palette_colors()`.

# scop 0.2.4

* **func**:
  * Rename functions: `check_srt_merge()` to `CheckDataMerge()`, `check_srt_list()` to `CheckDataList` and `check_data_type()` to `CheckDataType`.

# scop 0.2.2

* **func**:
  * Replace all `BiocParallel::bplapply()` with `thisutils::parallelize_fun()`.

* **bugs**:
  * Fix bugs in `RunSingleR()`.

# scop 0.2.0

* **func**:
  * Added `remove_python()` function for easy remove *Python* packages.

* **bugs**:
  * Corrected an issue in `py_to_r2()` function (intrinsic function), which ensures that Python-dependent functions like `RunPAGA()` and `RunSCVELO()` function run correctly.

# scop 0.1.9

* **func**:
  * Update `CellScoring()` and `AddModuleScore2()` functions. Now, new parameters `cores` and `verbose` have been added.
  * `AddModuleScore2()` function no longer uses the `BiocParallel::bpparam()` function to enable parallelization, but [thisutils::parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html), and the `cores` parameter is used to control the number of cores in [thisutils::parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).

# scop 0.1.5

* **bugs**:
  * Fix error for `RunPalantir()` function, see [#23](https://github.com/mengxu98/scop/issues/23).

# scop 0.1.4

* **func**:
  * Update `.onAttach()`, now `.onAttach()` will print more information about *conda* and *Python*.
  * Update `PrepareEnv()` function for easy add or update a conda environments and install *Python* packages.
  * Added `ListEnv()` and `RemoveEnv()` functions for easy management of *conda* environment and *Python* packages.

# scop 0.1.3

* **func**:
  * Added `TACSPlot()` function for creating FACS-like plots. Please refer to [Kernfeld et al. paper](https://doi.org/10.1016/j.immuni.2018.04.015) and [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271) for specific information.

# scop 0.0.9

* **bugs**:
  * Fix a bug for `PrepareEnv()` function.
  * The default *Python* version is now set to `3.10-1`.

# scop 0.0.6

* **func**:
  * Add `GetAssayData5()` function, a reimplementation of `GetAssayData()`, for compatibility with Seurat v5 `Assay` objects.
  * Updated `GetFeaturesData()` and `AddFeaturesData()` function to support retrieving and adding feature metadata.

# scop 0.0.5

* **data**:
  * Updated the `pancreas_sub` and `panc8_sub` test datasets to the Seurat v5 object format.
