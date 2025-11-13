# Changelog

## scop 0.6.0

- **func**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Enhanced with environment caching mechanism to avoid redundant
    environment preparation. Improved message formatting and error
    handling.
  - Python-based functions
    ([`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    `RunCellRank()`,
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md))
    now automatically call
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    internally, eliminating the need for users to manually prepare the
    Python environment before using these functions.
  - [`cluster_within_group2()`](https://mengxu98.github.io/scop/reference/cluster_within_group2.md):
    New function for clustering within groups.
  - Multiple plotting functions: Replaced `geom_sankey()` with
    [`ggsankey::geom_sankey()`](https://rdrr.io/pkg/ggsankey/man/geom_sankey.html)
    for better Sankey diagram support.
  - Multiple functions: Replaced `:::` operator with
    `get_namespace_fun()` for safer namespace access.
- **refactor**:
  - Removed `RunMonocle()` function and related documentation
    (`RunMonocle2.Rd`, `RunMonocle3.Rd`).
  - Removed `projection_functions.R` file (functions moved to other
    locations).
  - Replaced custom theme functions with
    [`thisplot::theme_this()`](https://mengxu98.github.io/thisplot/reference/theme_this.html)
    (exported as `theme_scop()`).
  - Replaced direct `log_message()` calls with
    [`thisutils::log_message()`](https://mengxu98.github.io/thisutils/reference/log_message.html)
    for consistency.
  - Removed `palette_list` data object.
- **deps**:
  - Moved `cli` from `Suggests` to `Imports` for better message
    formatting support.
  - Added `thisplot` to `Imports` for theme and utility functions.
  - Added `ggsankey` to `Suggests` for Sankey diagram support.
  - Added remote dependencies: `theislab/destiny`, `mengxu98/thisplot`.
  - Removed unused dependencies: `Biobase`, `BiocGenerics`,
    `concaveman`, `DDRTree`, `glmGamPoi`, `hexbin`, `monocle`, `png`,
    `ragg`, `tidyr`.
- **docs**:
  - Updated documentation across multiple functions to reflect code
    refactoring.
  - Improved code organization and maintainability.

## scop 0.5.5

- **bugs**:
  - Fixed
    [`VelocityPlot()`](https://mengxu98.github.io/scop/reference/VelocityPlot.md)
    function error in `plot_type = "grid"` mode: replaced vectorized
    arrow length with fixed-length arrows (using mean length) to resolve
    [`vapply()`](https://rdrr.io/r/base/lapply.html) error that occurred
    when [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html) received a
    vector instead of a single value, see
    [\#72](https://github.com/mengxu98/scop/issues/72),
    [\#74](https://github.com/mengxu98/scop/issues/74).

## scop 0.5.4

- **bugs**:
  - Fixed parameter name error in
    [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
    function calls: changed `data` parameter to `object` in
    [`RunKNNMap()`](https://mengxu98.github.io/scop/reference/RunKNNMap.md),
    [`RunSymphonyMap()`](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md),
    [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md),
    [`RunPCAMap()`](https://mengxu98.github.io/scop/reference/RunPCAMap.md),
    and
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
    functions, see [\#68](https://github.com/mengxu98/scop/issues/68).
  - Fixed `SingleCellExperiment` object creation in
    [`RunScmap()`](https://mengxu98.github.io/scop/reference/RunScmap.md)
    and
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md)
    functions: changed from coercing `SummarizedExperiment` to directly
    constructing `SingleCellExperiment` objects.

## scop 0.5.3

- **func**:
  - [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md):
    Changed default `Ensembl_version` parameter from `103` to `NULL` for
    more flexible version handling.
  - Added *Python* version `log_message()` for Python-based functions
    ([`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    `RunCellRank()`,
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md))
    and added `verbose` parameter inheritance and improved message
    formatting using cli-style formatting.
- **refactor**:
  - Delete `harmonizomeapi.py` file.
  - Move `scop_analysis.py` into a single `functions.py` file in
    `inst/python/` for better code organization and maintainability.
- **docs**:
  - Improved parameter documentation consistency.

## scop 0.5.1

- **docs**:
  - Improved reference formatting and consistency across multiple
    functions.
  - Enhanced documentation clarity and readability.

## scop 0.5.0

- **func**:
  - [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md):
    New function to perform CellChat analysis for investigating
    cell-to-cell communication with support for human, mouse, and
    zebrafish species.
  - [`CellChatPlot()`](https://mengxu98.github.io/scop/reference/CellChatPlot.md):
    New function to visualize CellChat analysis results with various
    plot types and customization options.
  - Multiple integration functions: Improved error messages and message
    formatting for better user experience.
- **deps**:
  - Added [CellChat](https://github.com/jinworks/CellChat) package
    dependency with remote repository `jinworks/CellChat`.
- **docs**:
  - Updated README.md with improved code formatting and examples.
  - Enhanced documentation for cell communication analysis functions.
  - Improved error messages and user guidance across integration
    functions.
- **refactor**:
  - Removed some example figures to optimize package installation size.

## scop 0.4.0

- **func**:
  - [`RunProportionTest()`](https://mengxu98.github.io/scop/reference/RunProportionTest.md):
    New function to perform Monte-carlo permutation test for quantifying
    cell proportion differences between conditions.
  - [`ProportionTestPlot()`](https://mengxu98.github.io/scop/reference/ProportionTestPlot.md):
    New function to generate proportion test plots with customizable
    significance thresholds and visualization options.
  - Multiple *Python*-based functions: add `\dontrun{}` blocks for
    Github workfolw checking.
- **docs**:
  - Added comprehensive documentation for new proportion testing
    functions.
  - Enhanced example documentation across multiple functions.
  - Updated package documentation and examples.

## scop 0.3.4

- **docs**:
  - Updated workflow examples and function documentation.

## scop 0.3.3

- **func**:
  - Multiple functions: Improved parameter documentation formatting and
    consistency across the package.

## scop 0.3.2

- **func**:
  - [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
    and
    [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md):
    Enhanced argument clarity, added input validation, and standardized
    return values for `Seurat`, `Assay`, and `Assay5` objects.
  - [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md):
    - Renamed parameters: `query_cell_annotation` → `query_annotation`,
      `ref_cell_annotation` → `ref_annotation`.
    - Improved error message formatting using cli-style formatting.
    - Simplified variable assignments and improved readability.
- **docs**:
  - Comprehensive documentation updates across multiple functions
    including `AnnotateFeatures`, `CellDimPlot`, `CellStatPlot`,
    `FeatureStatPlot`, `GroupHeatmap`, `RunCellQC`, and others.
  - Improved parameter descriptions and function clarity.

## scop 0.3.1

- **func**:
  - [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
    and
    [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md):
    Removed conditional font face styling (`face = ifelse()` logic) for
    better text rendering consistency. Set the default value of
    `lineheight` from `0.5` to `0.7`.
  - Updated
    [`check_r()`](https://mengxu98.github.io/scop/reference/check_r.md)
    function for improved package checking functionality.
  - Updated reexports functionality.
- **docs**:
  - Updated documentation formatting and consistency.

## scop 0.3.0

- **func**:
  - Fixed `segmentation faults` and `R crashes` on *M-series* MacBook
    when running *Python* functions.
  - [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md):
    Enhanced with *M-series* MacBook detection and automatic environment
    configuration.
  - [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md):
    Added ARM64-specific optimizations to prevent crashes and ensure
    stable execution.
  - `RunCellRank()`: Implemented *M-series* compatibility with proper
    NUMBA configuration.
  - [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md):
    Added ARM64 support with single-threaded execution mode.
  - [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md):
    Enhanced with *M-series* MacBook environment variable settings.
  - [`RunTriMap()`](https://mengxu98.github.io/scop/reference/RunTriMap.md):
    Added *M-series* MacBook compatibility for dimensionality reduction.
  - [`RunPaCMAP()`](https://mengxu98.github.io/scop/reference/RunPaCMAP.md):
    Implemented ARM64-specific environment configuration.
  - [`RunPHATE()`](https://mengxu98.github.io/scop/reference/RunPHATE.md):
    Added *M-series* MacBook support for non-linear dimensionality
    reduction.
  - [`RunCellQC()`](https://mengxu98.github.io/scop/reference/RunCellQC.md):
    Enhanced both scrublet and doubletdetection functions with ARM64
    compatibility.
- **docs**:
  - Updated function documentation to reflect *M-series* MacBook
    compatibility.
  - Added technical notes about ARM64 architecture considerations.

## scop 0.2.9

- **bugs**:
  - Fix bug for
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md).
- **docs**:
  - Updated documentation for some functions.

## scop 0.2.7

- **func**:
  - Added an internal function `.check_pkg_status()` to check if an *R*
    package is installed.
  - Update function
    [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
    to *S4* class function.
  - Update function
    [`standard_scop()`](https://mengxu98.github.io/scop/reference/standard_scop.md),
    make it more efficient.
- **data**:
  - Delete `lifemap` data, including: `lifemap_cell`,
    `lifemap_compartment` and `lifemap_organ`.
  - Reconstructed the sample data for both `panc8_sub` and
    `pancreas_sub`, retaining only the basic `Seurat` object.

## scop 0.2.6

- **func**:
  - Added
    [`remove_r()`](https://mengxu98.github.io/scop/reference/remove_r.md)
    function for easy remove *R* packages.
  - Rename function: `RemovePackages()` to
    [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md).
  - Removed other methods of installing *R* packages from the
    [`check_r()`](https://mengxu98.github.io/scop/reference/check_r.md)
    function, only retaining
    [pak::pak](https://pak.r-lib.org/reference/pak.html).
  - Delete useless import packages: `BBmisc`, `BiocManager`, `covr`,
    `devtools`, `promises` and `withr`.
  - Optimize the structure of `_pkgdown.yml` file.
- **docs**:
  - Updated documentation for some functions.

## scop 0.2.5

- **func**:
  - Rename function: `palette_scop()` to `palette_colors()`.

## scop 0.2.4

- **func**:
  - Rename functions: `check_srt_merge()` to
    [`CheckDataMerge()`](https://mengxu98.github.io/scop/reference/CheckDataMerge.md),
    `check_srt_list()` to `CheckDataList` and `check_data_type()` to
    `CheckDataType`.

## scop 0.2.2

- **func**:
  - Replace all
    [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
    with
    [`thisutils::parallelize_fun()`](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).
- **bugs**:
  - Fix bugs in
    [`RunSingleR()`](https://mengxu98.github.io/scop/reference/RunSingleR.md).

## scop 0.2.0

- **func**:
  - Added
    [`remove_python()`](https://mengxu98.github.io/scop/reference/remove_python.md)
    function for easy remove *Python* packages.
- **bugs**:
  - Corrected an issue in `py_to_r2()` function (intrinsic function),
    which ensures that Python-dependent functions like
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md)
    and
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md)
    function run correctly.

## scop 0.1.9

- **func**:
  - Update
    [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md)
    and `AddModuleScore2()` functions. Now, new parameters `cores` and
    `verbose` have been added.
  - `AddModuleScore2()` function no longer uses the
    [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
    function to enable parallelization, but
    [thisutils::parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html),
    and the `cores` parameter is used to control the number of cores in
    [thisutils::parallelize_fun](https://mengxu98.github.io/thisutils/reference/parallelize_fun.html).

## scop 0.1.5

- **bugs**:
  - Fix error for
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md)
    function, see [\#23](https://github.com/mengxu98/scop/issues/23).

## scop 0.1.4

- **func**:
  - Update `.onAttach()`, now `.onAttach()` will print more information
    about *conda* and *Python*.
  - Update
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    function for easy add or update a conda environments and install
    *Python* packages.
  - Added
    [`ListEnv()`](https://mengxu98.github.io/scop/reference/ListEnv.md)
    and
    [`RemoveEnv()`](https://mengxu98.github.io/scop/reference/RemoveEnv.md)
    functions for easy management of *conda* environment and *Python*
    packages.

## scop 0.1.3

- **func**:
  - Added
    [`TACSPlot()`](https://mengxu98.github.io/scop/reference/TACSPlot.md)
    function for creating FACS-like plots. Please refer to [Kernfeld et
    al. paper](https://doi.org/10.1016/j.immuni.2018.04.015) and
    [Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271)
    for specific information.

## scop 0.0.9

- **bugs**:
  - Fix a bug for
    [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)
    function.
  - The default *Python* version is now set to `3.10-1`.

## scop 0.0.6

- **func**:
  - Add
    [`GetAssayData5()`](https://mengxu98.github.io/scop/reference/GetAssayData5.md)
    function, a reimplementation of
    [`GetAssayData()`](https://satijalab.github.io/seurat-object/reference/AssayData.html),
    for compatibility with Seurat v5 `Assay` objects.
  - Updated
    [`GetFeaturesData()`](https://mengxu98.github.io/scop/reference/GetFeaturesData.md)
    and
    [`AddFeaturesData()`](https://mengxu98.github.io/scop/reference/AddFeaturesData.md)
    function to support retrieving and adding feature metadata.

## scop 0.0.5

- **data**:
  - Updated the `pancreas_sub` and `panc8_sub` test datasets to the
    Seurat v5 object format.
