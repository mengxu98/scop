# Run Palantir analysis

Run Palantir analysis

## Usage

``` r
RunPalantir(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group.by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  n_pcs = 30,
  n_neighbors = 30,
  dm_n_components = 10,
  dm_alpha = 0,
  dm_n_eigs = NULL,
  early_group = NULL,
  early_cell = NULL,
  terminal_cells = NULL,
  terminal_groups = NULL,
  num_waypoints = 1200,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  adjust_early_cell = FALSE,
  adjust_terminal_cells = FALSE,
  max_iterations = 25,
  magic_impute = FALSE,
  magic_layer = "MAGIC_imputed_data",
  cores = 1,
  point_size = 20,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "on data",
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "palantir",
  dirpath = "./",
  envname = NULL,
  conda = "auto",
  backend = c("python", "cpp"),
  allow_approximate = FALSE,
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. If provided, `adata` will be ignored.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.

- layer_x:

  Layer name for assay_x in the Seurat object.

- assay_y:

  Assays to convert as layers in the anndata object.

- layer_y:

  Layer names for the assay_y in the Seurat object.

- adata:

  An anndata object.

- group.by:

  Metadata column(s) used to color cells.

- linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

- basis:

  The basis to use for reduction, e.g., `"UMAP"`.

- n_pcs:

  Number of principal components to use for linear reduction.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph.

- dm_n_components:

  The number of diffusion components to calculate.

- dm_alpha:

  Normalization parameter for the diffusion operator.

- dm_n_eigs:

  Number of eigen vectors to use.

- early_group:

  Name of the group to start Palantir analysis from.

- early_cell:

  Name of the cell to start Palantir analysis from.

- terminal_cells:

  Character vector specifying terminal cells for Palantir analysis.

- terminal_groups:

  Character vector specifying terminal groups for Palantir analysis.

- num_waypoints:

  Number of waypoints to be included.

- scale_components:

  Should the cell fate probabilities be scaled for each component
  independently?

- use_early_cell_as_start:

  Should the starting cell for each terminal group be set as early_cell?

- adjust_early_cell:

  Whether to adjust the early cell to the cell with the minimum
  pseudotime value.

- adjust_terminal_cells:

  Whether to adjust the terminal cells to the cells with the maximum
  pseudotime value for each terminal group.

- max_iterations:

  Maximum number of iterations for pseudotime convergence.

- magic_impute:

  Whether to calculate a Palantir MAGIC expression layer. This layer is
  for visualization and trend fitting, not count-based testing.

- magic_layer:

  Name of the MAGIC layer to store.

- cores:

  The number of cores to use for `cellrank`.

- point_size:

  The point size for plotting.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc.

- show_plot:

  Whether to show the plot.

- save_plot:

  Whether to save plots to files.

- plot_format:

  Format for saved plots: `"pdf"`, `"png"`, or `"svg"`.

- plot_dpi:

  Resolution (DPI) for saved plots.

- plot_prefix:

  Prefix for saved plot filenames.

- dirpath:

  The directory to save the plots.

- envname:

  Optional Python environment name. `NULL` uses the current SCOP
  environment selection.

- conda:

  Conda-compatible executable used by
  [PrepareEnv](https://mengxu98.github.io/scop/reference/PrepareEnv.md).

- backend:

  Backend used to compute Palantir. `"python"` keeps the reference
  Palantir workflow and remains the default. `"cpp"` uses an approximate
  C++ workflow and stores results in `srt@tools[["Palantir"]]`. The C++
  backend never switches to Python implicitly.

- allow_approximate:

  Whether to allow the approximate C++ workflow. This must be `TRUE`
  when `backend = "cpp"`.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[PalantirTrajectoryPlot](https://mengxu98.github.io/scop/reference/PalantirTrajectoryPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:30:11] Start standard processing workflow...
#> ℹ [2026-09-06 22:30:12] Checking a list of <Seurat>...
#> ! [2026-09-06 22:30:12] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:30:12] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:30:12] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:30:12] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:30:12] Number of available HVF: 2000
#> ℹ [2026-09-06 22:30:12] Finished check
#> ℹ [2026-09-06 22:30:12] Perform `ScaleData()`
#> ℹ [2026-09-06 22:30:12] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:30:13] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:30:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:30:13] Reorder clusters...
#> ℹ [2026-09-06 22:30:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:30:13] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:30:21] Standard processing workflow completed
pancreas_sub <- RunPalantir(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  early_group = "Ductal",
  terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon"),
  backend = "cpp",
  allow_approximate = TRUE
)
#> ℹ [2026-09-06 22:30:21] Computing Palantir KNN graph with BiocNeighbors...
#> ✔ [2026-09-06 22:30:21] Palantir cpp backend completed

FeatureDimPlot(
  pancreas_sub,
  c("palantir_pseudotime", "palantir_diff_potential")
)


FeatureDimPlot(
  pancreas_sub,
  grep(
    "TerminalState_.*_diff_potential$",
    colnames(pancreas_sub@meta.data),
    value = TRUE
  )
)


PalantirTrajectoryPlot(
  pancreas_sub,
  reduction = "UMAP",
  pseudotime_interval = c(0, 0.9)
)


PalantirTrajectoryPlot(
  pancreas_sub,
  reduction = "UMAP",
  cell_color = "branch_selection",
  pseudotime_interval = c(0, 0.9)
)
```
