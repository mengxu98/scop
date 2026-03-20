# scop

# scop 0.8.6

* **feat**:
  * Added `RunGSVA()` for gene set variation analysis and `GSVAPlot()` for visualization of [GSVA](https://github.com/rcastelo/GSVA) results. Related issue #146 (@mengxu98).
  * `RunMetabolism()`
    * Added `RunMetabolism()` for single-cell metabolism pathway scoring with support for `AUCell`, `GSVA`, `ssGSEA`, and optional `VISION`, and `MetabolismPlot()` for visualization
    * `RunMetabolism()` now uses [scMetabolism](https://github.com/wu-yc/scMetabolism) KEGG / Reactome metabolism pathway definitions to identify metabolism-related terms, then rebuilds pathway gene sets from updated `PrepareDB()` annotations, enabling cached database updates and cross-species metabolism scoring.
    Related issue #146 (@mengxu98).
  * Added `RunDecontX()` and integrated optional `decontX` ambient RNA decontamination into `RunCellQC()`, including contamination metadata, optional decontaminated assay output, and threshold-based QC filtering. Related issue #147 (@mengxu98).

* **fix**:
  * `FeatureStatPlot()`: Fixed duplicated X/Y axis titles and main title when `stack = TRUE` with `theme_args` or `title`. Stack assembly now strips axis and plot titles from non-edge panels and draws a single shared ylab/title via gtable; theme styling from `theme_args` (e.g. `axis.title.y`, `title`/`plot.title`) is applied to these shared grobs. Related issue #145 (@PanSX-Dr).

* **docs**:
  * Updated the README badge/logo.

* **data**:
  * Removed the example datasets `ifnb_sub`, `ref_scHCL`, and `ref_scZCL`.

# scop 0.8.5

* **feat**:
  * Default palette changed from `"Paired"` to `"Chinese"`. See `thisplot::ChineseColors()` for details.
  * Unified the `cores` parameter across multiple functions (`DynamicHeatmap()`, `FeatureHeatmap()`, `GroupHeatmap()`, and `heatmap_enrichment()`), ensuring `cores` is correctly threaded through to `RunEnrichment()`.
  * `RunMonocle2()` and `RunMonocle3()`: Added `xlab` and `ylab` parameters, passed through to internal `CellDimPlot()` and `FeatureDimPlot()` calls.
  * `ListDB()`: Now supports multiple species in the `species` parameter simultaneously and adds `Species` and `DB` columns to the output data frame for clearer identification.
  * Moved `is_outlier` to `thisutils::is_outlier()`.
  * `configure_apple_silicon_env()` (*Python* function): Added OpenMP compatibility handling on macOS arm64 by prepending environment `lib` paths to `DYLD_FALLBACK_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` and preloading `libomp.dylib`/`libiomp5.dylib` via `ctypes.CDLL(..., RTLD_GLOBAL)` to reduce `scanpy`/`python-igraph` import conflicts.
  * `scVI_integrate()`: Unified parameter naming by renaming `num_threads` to `cores` for consistency across integration functions.
  * Added `find_neighbors_and_clusters()` and `run_nonlinear_reduction()` helper functions to unify operations across integration functions and reduce redundant code.

* **fix**:
  * `PrepareDB()`: 
    * TF database now uses [AnimalTFDB4](https://github.com/mengxu98/AnimalTFDB4) as the data source.
    * MP database - the Web Archive URL year is now dynamically determined from the current date, and the archive snapshot date for all MP-related file downloads (`VOC_MammalianPhenotype.rpt`, `MGI_Gene_Model_Coord.rpt`, `MGI_GenePheno.rpt`) is extracted from the server index page to ensure consistent versioning.
    * hTFtarget data download URL has been changed from `"http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"` to `"https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"`.
    * CSPA data download URL has been changed from `"https://wlab.ethz.ch/cspa/data/S1_File.xlsx"` to `"https://raw.githubusercontent.com/mengxu98/CSPA/main/S1_File.xlsx"`.

    Related issus #76 (@hwa2Hu), #139 (@pengding774-dot), #140 (@mengxu98).
  * `LIGER_integrate()`
    * Migrated to the `rliger` 2.x workflow (`rliger::runIntegration()` + `rliger::quantileNorm()` on `Seurat` object) and now prepares/uses the `ligerScaleData` layer via `rliger::scaleNotCenter()` before integration.
    * Updated argument naming/style from `LIGER_dims_use` to `liger_dims_use`, and removed legacy quantile-normalization parameter compatibility mapping (`ref_dataset`), keeping `reference` as the supported interface.
  * Optimized the installation of some *Python* packages on Apple Silicon devices.

* **docs**:
  * Updated README example code and visualizations. After regenerating figures, the package size was reduced by ~3 MB.

# scop 0.8.4

* **fix**:
  * `FeatureStatPlot()` / `ExpressionStatPlot()`: Fixed box and violin x-axis misalignment when `add_box = TRUE` with `split.by`. Groups with fewer than 2 observations are now filtered before violin density estimation (with a warning), and the violin layer uses a consistent `position_dodge(width = 0.9)` to match the boxplot. Related issue #123 (@oranges7). 

# scop 0.8.3

* **feat**:
  * `RunDynamicFeatures()`: Added [PreTSA](https://github.com/haotian-zhuang/PreTSA/) method for dynamic feature fitting. The [PreTSA](https://github.com/haotian-zhuang/PreTSA/) algorithm in `scop` has been re implemented to support parallelization for higher performance. Original research: [PreTSA: computationally efficient modeling of temporal and spatial gene expression patterns](https://doi.org/10.1186/s13059-026-03994-3). Use `fit_method = "pretsa"` for B-spline-based piecewise truncated spline analysis; `fit_method = "gam"` (default) keeps generalized additive models. PreTSA supports `knot` (`0` or `"auto"`) and `max_knot_allowed` when `knot = "auto"`. Relate issue #133.
  * `CellDimPlot()` and `FeatureDimPlot()`: Added `legend.title` parameter (default `NULL`) to control the legend title. When `NULL`, default titles are used (e.g. group name for `CellDimPlot`, feature or empty for `FeatureDimPlot`).
  * `ExpressionStatPlot()` and `FeatureStatPlot()`: Added `legend.title` parameter (default `NULL`) for single-legend plots. When `NULL`, the default title (e.g. `keynm` or feature/group name) is used. `FeatureStatPlot()` forwards `legend.title` to `ExpressionStatPlot()`.
  * Moved `StatPlot` function to `thisplot::StatPlot`.

* **fix**:
  * `DynamicHeatmap()` / `heatmap_enrichment()`: Fixed incorrect `db` handling when using custom `TERM2GENE`/`TERM2NAME`. Enrichment results with `Database = "custom"` could be incorrectly filtered by default `db` values (e.g. `"GO_BP"`), causing false "No term enriched using the threshold" warnings even when enrichment succeeded. Relate issue #133 (@1228849000).

# scop 0.8.2

* **feat**:
  * Add the `DEtestPlot()` function, which calls the original `VolcanoPlot()` and adds two plot types, Manhattan and Ring, controlled by the `plot_type` parameter (`c("volcano", "manhattan", "ring")`, default `"volcano"`). Add standalone functions `DEtestManhattanPlot()` and `DEtestRingPlot()` for direct use. Relate issue #121 (@ericavalentini).
  * Differential expression visualization (`DEtestPlot()`, `VolcanoPlot()`, `DEtestManhattanPlot()`, `DEtestRingPlot()`): added `res` parameter to accept existing DE results (data.frame). When `res` is provided, `srt` is ignored. Data processing supports: (1) `group1` or `cluster` column for grouped plots; (2) no grouping column for a single panel; (3) gene names from row names when `gene` column is missing. Relate issue #129 (@mengxu98).
  * `PrepareEnv()`: Update `version` parameter to specify the *Python* version of the conda environment. Default is `"3.11-1"` on *Windows* and `"3.10-1"` on *macOS* and *Unix*. Relate issue #103 (@PanSX-Dr).
  * `RunNMF()`: Add the `cores` parameter for `RunNMF()` and optimize the printed message.
  * Removed the setting in *Python* functions that prevents drawing functions from causing *R* crashes.

* **fix**:
  * `RunPalantir()`: Fixed unused `plot_format` parameter error. The parameter is now properly excluded from arguments passed to Python functions. Relate issue #126 (@christinejay990202-dev).

# scop 0.8.1

* **fix**:
  * `FeatureHeatmap()`: Fixed `group_palcolor` when a named vector is passed: the function previously used only the first color for all groups because `group_palcolor[[1]]` on a vector returns a single element. Now when `group.by` has length 1, a vector is automatically wrapped as a list; when `within_groups = TRUE`, `group_palcolor` is expanded in line with `group_palette`.
  * `GroupHeatmap()`: Same `group_palcolor` fix as `FeatureHeatmap()`: support for named-vector input and correct expansion when `within_groups = TRUE`.

* **docs**:
  * `group_by` modified to `group.by` to unify documentation, involving functions: `RunCellRank()`, `RunDEtest()`, `VolcanoPlot()`, `RunEnrichment()`, `EnrichmentPlot()`, `RunGSEA()`, `GSEAPlot()`, `RunPalantir()`, `RunPAGA()`, `RunWOT()`, `RunSCVELO()`, releated to [#120](https://github.com/mengxu98/scop/issues/120).

# scop 0.8.0

* **feat**:
  * `RunMonocle2()`: New function for performing [Monocle2](https://github.com/mengxu98/monocle) trajectory analysis with support for various dimensionality reduction methods (DDRTree, ICA, tSNE, SimplePPT, L1-graph, SGL-tree). Uses the fixed version of monocle2 from [mengxu98/monocle](https://github.com/mengxu98/monocle).
  * `RunMonocle3()`: New function for performing [Monocle3](https://github.com/cole-trapnell-lab/monocle3) trajectory analysis with support for cell ordering, trajectory learning, and pseudotime computation.
  * `RunCytoTRACE()`: New function for running [CytoTRACE 2](https://github.com/digitalcytometry/cytotrace2) analysis to predict cellular potency scores and categories (Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent) with support for human and mouse species.
  * `CytoTRACEPlot()`: New function for visualizing [CytoTRACE 2](https://github.com/digitalcytometry/cytotrace2) analysis results.

# scop 0.7.9

* **docs**:
  * Optimized `@inheritParams` usage to reduce redundant parameter definitions.
  * Unified parameter naming: renamed `group_by` to `group.by` for consistency with Seurat conventions in `RunMonocle3()`, `RunPalantir()`, `RunSCVELO()`, `RunWOT()`, `CellCorHeatmap()`, and `RunDynamicFeatures()`.

# scop 0.7.8

* **fix**:
  * `RunPalantir()`: Fixed unused `plot_format` parameter error. The parameter is now properly excluded from arguments passed to Python functions. Relate issue #114 (@Moonerss).
  * `RunSCVELO()`: Fixed PAGA computation error by replacing `scv.tl.paga` with `sc.tl.paga` (scanpy implementation) for better stability. The function now uses the same PAGA implementation as `RunPAGA()` function.
  * `RunCellRank()`: Fixed GPCCA Schur decomposition error by adding fallback mechanism. When `brandts` method fails with "subspace_angles" error, the function automatically tries `krylov` method. If both methods fail, it automatically switches to CFLARE estimator for more robust computation.
  * `RunCellRank()`: Fixed `recover_dynamics` error by ensuring `velocity_graph` and `velocity_graph_neg` are properly set before calling `scv.tl.recover_dynamics()` for latent time computation.

* **feat**:
  * `adata_to_srt()`: Removed automatic removal of "X_" prefix from dimensionality reduction names in `obsm` keys. The function now preserves original reduction names as they are stored in AnnData objects.

* **data**:
  * Reducing the size of `pancreas_sub` example dataset.

# scop 0.7.7

* **feat**:
  * `adata_to_srt()`: Enhanced to support multiple AnnData object types including *Python* AnnData objects (from scanpy/reticulate), R6 AnnData objects from the `anndata` package (AnnDataR6), and R6 AnnData objects from the `anndataR` package (InMemoryAnnData). Added internal helper functions `get_adata_element()` and `get_adata_names()` for better compatibility. Relate issues #67 (@lisch7), #91 (@mengxu98) and [commit91#issuecomment](https://github.com/mengxu98/scop/issues/91#issuecomment-3659404993).

* **fix**:
  * `RunDEtest()`: Fixed error when comparing one cluster against multiple clusters using `group1` and `group2` parameters. Relate issue #111 (@zhaoxiaoyan9225).
  * `AnnotateFeatures()`: Fixed bug where the function would fail when processing GTF file annotations due to column name matching issues during data naming. The function now correctly handles column name intersections when merging annotation data.

# scop 0.7.6

* **feat**:
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.5

* **fix**:
  * `CellScoring()`: Fixed bug where the function failed to build results. Relate issue #98 (@SuperrNaruto).

* **feat**:
  * `RunDEtest()`: Fixed compatibility issue with [SeuratObject](https://satijalab.github.io/seurat-object/) 5.0.0+ by replacing deprecated `Assays()` `slot` argument with `LayerData()`. Relate issue #100 (@mattizecos).
  * `RunDM()`: Added automatic PCA-based dimensionality reduction when using many features (>1000) to speed up diffusion map computation. The `npcs` parameter can be used to control the number of principal components used for pre-processing.

# scop 0.7.3

* **feat**:
  * `RunCellTypist()`: New function for cell type annotation using the CellTypist method.
  * `CellTypistModels()`: New function for downloading and managing CellTypist pre-trained models.

# scop 0.7.2

* **feat**:
  * `RunCellRank()`: Performance optimizations and code improvements.

# scop 0.7.1

* **fix**:
  * `CellDimPlot()`: Fixed issue where NA values appeared in labels. Relate issue #93 (@12345nkjil).

# scop 0.7.0

* **feat**:
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

* **feat**:
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

* **feat**:
  * `CellChatPlot()`: Adjusted the size of saved figures for better file size optimization.

* **docs**:
  * Updated README.md to remove references to Monocle2 and Monocle3 (deprecated functions).

# scop 0.6.1

* **feat**:
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

* **feat**:
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

* **fix**:
  * Fixed `VelocityPlot()` function error in `plot_type = "grid"` mode: replaced vectorized arrow length with fixed-length arrows (using mean length) to resolve `vapply()` error that occurred when `grid::arrow()` received a vector instead of a single value, see [#72](https://github.com/mengxu98/scop/issues/72), [#74](https://github.com/mengxu98/scop/issues/74).

# scop 0.5.4

* **fix**:
  * Fixed parameter name error in `CheckDataType()` function calls: changed `data` parameter to `object` in `RunKNNMap()`, `RunSymphonyMap()`, `RunScmap()`, `RunPCAMap()`, and `RunSingleR()` functions, see [#68](https://github.com/mengxu98/scop/issues/68).
  * Fixed `SingleCellExperiment` object creation in `RunScmap()` and `RunSingleR()` functions: changed from coercing `SummarizedExperiment` to directly constructing `SingleCellExperiment` objects.

# scop 0.5.3

* **feat**:
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

* **feat**:
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

* **feat**:
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

* **feat**:
  * Multiple functions: Improved parameter documentation formatting and consistency across the package.

# scop 0.3.2

* **feat**:
  * `GetFeaturesData()` and `AddFeaturesData()`: Enhanced argument clarity, added input validation, and standardized return values for `Seurat`, `Assay`, and `Assay5` objects.
  * `CellCorHeatmap()`: 
    - Renamed parameters: `query_cell_annotation` → `query_annotation`, `ref_cell_annotation` → `ref_annotation`.
    - Improved error message formatting using cli-style formatting.
    - Simplified variable assignments and improved readability.

* **docs**:
  * Comprehensive documentation updates across multiple functions including `AnnotateFeatures`, `CellDimPlot`, `CellStatPlot`, `FeatureStatPlot`, `GroupHeatmap`, `RunCellQC`, and others.
  * Improved parameter descriptions and function clarity.

# scop 0.3.1

* **feat**:
  * `EnrichmentPlot()` and `GSEAPlot()`: Removed conditional font face styling (`face = ifelse()` logic) for better text rendering consistency. Set the default value of `lineheight` from `0.5` to `0.7`.
  * Updated `check_r()` function for improved package checking functionality.
  * Updated reexports functionality.

* **docs**:
  * Updated documentation formatting and consistency.

# scop 0.3.0

* **feat**:
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

* **fix**:
  * Fix bug for `RunSCVELO()`.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.7

* **feat**:
  * Added an internal function `.check_pkg_status()` to check if an *R* package is installed.
  * Update function `CheckDataType()` to *S4* class function.
  * Update function `standard_scop()`, make it more efficient.

* **data**:
  * Delete `lifemap` data, including: `lifemap_cell`, `lifemap_compartment` and `lifemap_organ`.
  * Reconstructed the sample data for both `panc8_sub` and `pancreas_sub`, retaining only the basic `Seurat` object.

# scop 0.2.6

* **feat**:
  * Added `remove_r()` function for easy remove *R* packages.
  * Rename function: `RemovePackages()` to `remove_python()`.
  * Removed other methods of installing *R* packages from the `check_r()` function, only retaining [pak::pak](https://pak.r-lib.org/reference/pak.html). 
  * Delete useless import packages: `BBmisc`, `BiocManager`, `covr`, `devtools`, `promises` and `withr`.
  * Optimize the structure of `_pkgdown.yml` file.

* **docs**:
  * Updated documentation for some functions.

# scop 0.2.5

* **feat**:
  * Rename function: `palette_scop()` to `palette_colors()`.

# scop 0.2.4

* **feat**:
  * Rename functions: `check_srt_merge()` to `CheckDataMerge()`, `check_srt_list()` to `CheckDataList` and `check_data_type()` to `CheckDataType`.

# scop 0.2.2

* **feat**:
  * Replace all `BiocParallel::bplapply()` with `thisutils::parallelize_fun()`.

* **fix**:
  * Fix bugs in `RunSingleR()`.

# scop 0.2.0

* **feat**:
  * Added `remove_python()` function for easy remove *Python* packages.

* **fix**:
  * Corrected an issue in `py_to_r2()` function (intrinsic function), which ensures that Python-dependent functions like `RunPAGA()` and `RunSCVELO()` function run correctly.

# scop 0.1.9

* **feat**:
  * Update `CellScoring()` and `AddModuleScore2()` functions. Now, new parameters `cores` and `verbose` have been added.
  * `AddModuleScore2()` function no longer uses the `BiocParallel::bpparam()` function to enable parallelization, but `thisutils::parallelize_fun`, and the `cores` parameter is used to control the number of cores in `thisutils::parallelize_fun`.

# scop 0.1.5

* **fix**:
  * Fix error for `RunPalantir()` function, see [#23](https://github.com/mengxu98/scop/issues/23).

# scop 0.1.4

* **feat**:
  * Update `.onAttach()`, now `.onAttach()` will print more information about *conda* and *Python*.
  * Update `PrepareEnv()` function for easy add or update a conda environments and install *Python* packages.
  * Added `ListEnv()` and `RemoveEnv()` functions for easy management of *conda* environment and *Python* packages.

# scop 0.1.3

* **feat**:
  * Added `TACSPlot()` function for creating FACS-like plots. Please refer to [Kernfeld et al. paper](https://doi.org/10.1016/j.immuni.2018.04.015) and [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271) for specific information.

# scop 0.0.9

* **fix**:
  * Fix a bug for `PrepareEnv()` function.
  * The default *Python* version is now set to `3.10-1`.

# scop 0.0.6

* **feat**:
  * Add `GetAssayData5()` function, a reimplementation of `GetAssayData()`, for compatibility with Seurat v5 `Assay` objects.
  * Updated `GetFeaturesData()` and `AddFeaturesData()` function to support retrieving and adding feature metadata.

# scop 0.0.5

* **data**:
  * Updated the `pancreas_sub` and `panc8_sub` test datasets to the Seurat v5 object format.

# scop 0.0.1

* **Initial version**
