# Changelog

## scop 0.8.1

- **bugs**:
  - [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md):
    Fixed `group_palcolor` when a named vector is passed: the function
    previously used only the first color for all groups because
    `group_palcolor[[1]]` on a vector returns a single element. Now when
    `group.by` has length 1, a vector is automatically wrapped as a
    list; when `within_groups = TRUE`, `group_palcolor` is expanded in
    line with `group_palette`.
  - [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md):
    Same `group_palcolor` fix as
    [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md):
    support for named-vector input and correct expansion when
    `within_groups = TRUE`.
- **docs**:
  - `group_by` modified to `group.by` to unify documentation, involving
    functions:
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
    [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md),
    [`VolcanoPlot()`](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
    [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
    [`EnrichmentPlot()`](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
    [`RunGSEA()`](https://mengxu98.github.io/scop/reference/RunGSEA.md),
    [`GSEAPlot()`](https://mengxu98.github.io/scop/reference/GSEAPlot.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    releated to [\#120](https://github.com/mengxu98/scop/issues/120).

## scop 0.8.0

- **func**:
  - [`RunMonocle2()`](https://mengxu98.github.io/scop/reference/RunMonocle2.md):
    New function for performing
    [Monocle2](https://github.com/mengxu98/monocle) trajectory analysis
    with support for various dimensionality reduction methods (DDRTree,
    ICA, tSNE, SimplePPT, L1-graph, SGL-tree). Uses the fixed version of
    monocle2 from
    [mengxu98/monocle](https://github.com/mengxu98/monocle).
  - [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md):
    New function for performing
    [Monocle3](https://github.com/cole-trapnell-lab/monocle3) trajectory
    analysis with support for cell ordering, trajectory learning, and
    pseudotime computation.
  - [`RunCytoTRACE()`](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md):
    New function for running [CytoTRACE
    2](https://github.com/digitalcytometry/cytotrace2) analysis to
    predict cellular potency scores and categories (Differentiated,
    Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent) with
    support for human and mouse species.
  - [`CytoTRACEPlot()`](https://mengxu98.github.io/scop/reference/CytoTRACEPlot.md):
    New function for visualizing [CytoTRACE
    2](https://github.com/digitalcytometry/cytotrace2) analysis results.

## scop 0.7.9

- **docs**:
  - Optimized `@inheritParams` usage to reduce redundant parameter
    definitions.
  - Unified parameter naming: renamed `group_by` to `group.by` for
    consistency with Seurat conventions in
    [`RunMonocle3()`](https://mengxu98.github.io/scop/reference/RunMonocle3.md),
    [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md),
    [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
    [`RunWOT()`](https://mengxu98.github.io/scop/reference/RunWOT.md),
    [`CellCorHeatmap()`](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md),
    and
    [`RunDynamicFeatures()`](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md).

## scop 0.7.8

- **bugs**:
  - [`RunPalantir()`](https://mengxu98.github.io/scop/reference/RunPalantir.md):
    Fixed unused `plot_format` parameter error. The parameter is now
    properly excluded from arguments passed to Python functions. This
    issue reported in
    [\#114](https://github.com/mengxu98/scop/issues/114).
  - [`RunSCVELO()`](https://mengxu98.github.io/scop/reference/RunSCVELO.md):
    Fixed PAGA computation error by replacing `scv.tl.paga` with
    `sc.tl.paga` (scanpy implementation) for better stability. The
    function now uses the same PAGA implementation as
    [`RunPAGA()`](https://mengxu98.github.io/scop/reference/RunPAGA.md)
    function.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Fixed GPCCA Schur decomposition error by adding fallback mechanism.
    When `brandts` method fails with “subspace_angles” error, the
    function automatically tries `krylov` method. If both methods fail,
    it automatically switches to CFLARE estimator for more robust
    computation.
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Fixed `recover_dynamics` error by ensuring `velocity_graph` and
    `velocity_graph_neg` are properly set before calling
    `scv.tl.recover_dynamics()` for latent time computation.
- **func**:
  - [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md):
    Removed automatic removal of “X\_” prefix from dimensionality
    reduction names in `obsm` keys. The function now preserves original
    reduction names as they are stored in AnnData objects.
- **data**:
  - Reducing the size of `pancreas_sub` example dataset.

## scop 0.7.7

- **func**:
  - [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md):
    Enhanced to support multiple AnnData object types including *Python*
    AnnData objects (from scanpy/reticulate), R6 AnnData objects from
    the `anndata` package (AnnDataR6), and R6 AnnData objects from the
    `anndataR` package (InMemoryAnnData). Added internal helper
    functions `get_adata_element()` and `get_adata_names()` for better
    compatibility. This enhancement addresses the issue reported in
    [\#67](https://github.com/mengxu98/scop/issues/67),
    [\#91](https://github.com/mengxu98/scop/issues/91) and
    [91#issuecomment](https://github.com/mengxu98/scop/issues/91#issuecomment-3659404993).
- **bugs**:
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md):
    Fixed error when comparing one cluster against multiple clusters
    using `group1` and `group2` parameters. This issue reported in
    [\#111](https://github.com/mengxu98/scop/issues/111).
  - [`AnnotateFeatures()`](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md):
    Fixed bug where the function would fail when processing GTF file
    annotations due to column name matching issues during data naming.
    The function now correctly handles column name intersections when
    merging annotation data.

## scop 0.7.6

- **func**:
  - [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md):
    Added automatic PCA-based dimensionality reduction when using many
    features (\>1000) to speed up diffusion map computation. The `npcs`
    parameter can be used to control the number of principal components
    used for pre-processing.

## scop 0.7.5

- **bugs**:
  - [`CellScoring()`](https://mengxu98.github.io/scop/reference/CellScoring.md):
    Fixed bug where the function failed to build results. This issue
    reported in [\#98](https://github.com/mengxu98/scop/issues/98).
- **func**:
  - [`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md):
    Fixed compatibility issue with
    [SeuratObject](https://satijalab.github.io/seurat-object/) 5.0.0+ by
    replacing deprecated
    [`Assays()`](https://satijalab.github.io/seurat-object/reference/ObjectAccess.html)
    `slot` argument with
    [`LayerData()`](https://satijalab.github.io/seurat-object/reference/Layers.html).
    This issue reported in
    [\#100](https://github.com/mengxu98/scop/issues/100).
  - [`RunDM()`](https://mengxu98.github.io/scop/reference/RunDM.md):
    Added automatic PCA-based dimensionality reduction when using many
    features (\>1000) to speed up diffusion map computation. The `npcs`
    parameter can be used to control the number of principal components
    used for pre-processing.

## scop 0.7.3

- **func**:
  - [`RunCellTypist()`](https://mengxu98.github.io/scop/reference/RunCellTypist.md):
    New function for cell type annotation using the CellTypist method.
  - [`CellTypistModels()`](https://mengxu98.github.io/scop/reference/CellTypistModels.md):
    New function for downloading and managing CellTypist pre-trained
    models.

## scop 0.7.2

- **func**:
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Performance optimizations and code improvements.

## scop 0.7.1

- **bugs**:
  - [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md):
    Fixed issue where NA values appeared in labels. This issue reported
    in [\#93](https://github.com/mengxu98/scop/issues/93).

## scop 0.7.0

- **func**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Integrated `uv` as the primary *Python* package installer for
    improved installation speed.
  - [`check_python()`](https://mengxu98.github.io/scop/reference/check_python.md):
    Now uses `uv` as the primary installation tool with `pip` as
    fallback, significantly improving package installation speed.
  - Added `find_uv()` and `install_uv()` internal functions for managing
    `uv` package manager installation and detection.

## scop 0.6.6

- **docs**:
  - Unified documentation format across all R functions:
    - Standardized return value tags: Changed all `@returns` to
      `@return` for consistency.
    - Unified parameter documentation: Replaced all `\code{value}` with
      Markdown backticks `` `value` `` format.
    - Standardized default value descriptions.
    - Added `@md` tags: Added `@md` tags to all functions using Markdown
      syntax in documentation.
    - Enhanced cross-references: Added `@seealso` links to related
      functions where appropriate.

## scop 0.6.5

- **func**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    - Added comprehensive environment variable configuration to prevent
      crashes when calling *Python* functions, including setting thread
      limits for OMP, OPENBLAS, MKL, NUMBA, and other libraries. This
      improves stability on all platforms, especially Apple silicon
      Macs.
    - Added `accept_conda_tos()` function to automatically accept conda
      Terms of Service for required channels, improving the conda
      environment setup process.
    - Fixed conda Terms of Service acceptance issue in
      [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md).
      The function now automatically accepts conda Terms of Service for
      required channels, eliminating the need for manual acceptance.
      This addresses the issue reported in
      [\#85](https://github.com/mengxu98/scop/issues/85).
  - Multiple Python-based functions (`RunPAGA`, `RunSCVELO`,
    `RunPalantir`, `RunCellRank`, `RunWOT`, `RunPHATE`, `RunPaCMAP`,
    `RunTriMap`): Enhanced message formatting and code improvements.
  - [`PrepareSCExplorer()`](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md):
    Fixed package version dependency issues with `shiny` and `bslib`
    compatibility. The function now properly handles `bslib` theme
    configuration to work with both `shiny` 1.6.0 and 1.7.0+, addressing
    compatibility errors reported in
    [\#87](https://github.com/mengxu98/scop/issues/87).
- **refactor**:
  - Improved code formatting and consistency across multiple functions.
  - Enhanced Python functions in `inst/python/functions.py` with better
    error handling and message formatting.
- **docs**:
  - Updated documentation for multiple functions to reflect code
    improvements.

## scop 0.6.2

- **func**:
  - [`CellChatPlot()`](https://mengxu98.github.io/scop/reference/CellChatPlot.md):
    Adjusted the size of saved figures for better file size
    optimization.
- **docs**:
  - Updated README.md to remove references to Monocle2 and Monocle3
    (deprecated functions).

## scop 0.6.1

- **func**:
  - [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md):
    Improved message formatting and simplified log output for better
    user experience.
  - Added `get_conda_envs_dir()` helper function to centralize conda
    environment directory retrieval.
  - [`integration_scop()`](https://mengxu98.github.io/scop/reference/integration_scop.md):
    Enhanced `integration_method` parameter definition with explicit
    method list for better code clarity.
- **refactor**:
  - Moved `exist_python_pkgs()` function to `check_package.R` for better
    code organization.
  - Replaced direct `conda_info()$envs_dirs[1]` calls with
    `get_conda_envs_dir()` helper function for consistency.
  - [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md):
    Updated to use
    [`thisplot::palette_list`](https://mengxu98.github.io/thisplot/reference/palette_list.html)
    and
    [`thisplot::slim_data()`](https://mengxu98.github.io/thisplot/reference/slim_data.html)
    instead of `scop::palette_list` and `scop::slim_data()`.
  - Added `thisplot` to dependency checks in
    [`RunSCExplorer()`](https://mengxu98.github.io/scop/reference/RunSCExplorer.md).
- **docs**:
  - Updated documentation across multiple functions.

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
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
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
    [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md),
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
  - [`RunCellRank()`](https://mengxu98.github.io/scop/reference/RunCellRank.md):
    Implemented *M-series* compatibility with proper NUMBA
    configuration.
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
