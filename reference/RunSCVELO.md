# Run scVelo workflow

scVelo is a scalable toolkit for RNA velocity analysis in single cells.
This function runs an enhanced scVelo workflow on a Seurat object with
improved error handling, version compatibility, and modular design.

## Usage

``` r
RunSCVELO(
  srt = NULL,
  adata = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  group_by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  mode = "stochastic",
  fitting_by = "stochastic",
  magic_impute = FALSE,
  knn = 5,
  t = 2,
  min_shared_counts = 30,
  n_pcs = 30,
  n_neighbors = 30,
  filter_genes = TRUE,
  min_counts = 3,
  min_counts_u = 3,
  normalize_per_cell = TRUE,
  log_transform = TRUE,
  use_raw = FALSE,
  diff_kinetics = FALSE,
  stream_smooth = NULL,
  stream_density = 2,
  arrow_length = 5,
  arrow_size = 5,
  arrow_density = 0.5,
  denoise = FALSE,
  denoise_topn = 3,
  kinetics = FALSE,
  kinetics_topn = 100,
  calculate_velocity_genes = FALSE,
  compute_velocity_confidence = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_paga = TRUE,
  top_n = 6,
  cores = 1,
  palette = "Paired",
  palcolor = NULL,
  legend.position = "on data",
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "scvelo",
  dirpath = "./scvelo",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. Default is `NULL`. If provided, `adata` will be
  ignored.

- adata:

  An anndata object. Default is `NULL`.

- assay_x:

  Assay to convert in the anndata object.

- layer_x:

  Layer name for `assay_x` in the Seurat object.

- assay_y:

  Assay to convert in the anndata object.

- layer_y:

  Layer names for the `assay_y` in the Seurat object.

- group_by:

  Variable to use for grouping cells in the Seurat object.

- linear_reduction:

  Linear reduction method to use, e.g., `"PCA"`.

- nonlinear_reduction:

  Non-linear reduction method to use, e.g., `"UMAP"`.

- basis:

  The basis to use for reduction, e.g., `"UMAP"`.

- mode:

  Velocity estimation models to use. Can be a vector containing
  `"deterministic"`, `"stochastic"`, and/or `"dynamical"`.

- fitting_by:

  Method used to fit gene velocities for dynamical modeling, e.g.,
  "stochastic".

- magic_impute:

  Flag indicating whether to perform magic imputation.

- knn:

  The number of nearest neighbors for `magic.MAGIC`.

- t:

  power to which the diffusion operator is powered for `magic.MAGIC`.

- min_shared_counts:

  Minimum number of counts (both unspliced and spliced) required for a
  gene.

- n_pcs:

  Number of principal components to use for linear reduction. Default is
  `30`.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph. Default is
  `30`.

- filter_genes:

  Whether to filter genes based on minimum counts.

- min_counts:

  Minimum counts for gene filtering.

- min_counts_u:

  Minimum unspliced counts for gene filtering.

- normalize_per_cell:

  Whether to normalize counts per cell.

- log_transform:

  Whether to apply log transformation.

- use_raw:

  Whether to use raw data for dynamical modeling.

- diff_kinetics:

  Whether to use differential kinetics.

- stream_smooth:

  Multiplication factor for scale in Gaussian kernel around grid point.

- stream_density:

  Controls the closeness of streamlines. When density = 2 (default), the
  domain is divided into a 60x60 grid, whereas density linearly scales
  this grid. Each cell in the grid can have, at most, one traversing
  streamline.

- arrow_length:

  Length of arrows.

- arrow_size:

  Size of arrows.

- arrow_density:

  Amount of velocities to show.

- denoise:

  Boolean flag indicating whether to denoise.

- denoise_topn:

  Number of genes with highest likelihood selected to infer velocity
  directions.

- kinetics:

  Boolean flag indicating whether to estimate RNA kinetics.

- kinetics_topn:

  Number of genes with highest likelihood selected to infer velocity
  directions.

- calculate_velocity_genes:

  Boolean flag indicating whether to calculate velocity genes.

- compute_velocity_confidence:

  Whether to compute velocity confidence metrics.

- compute_terminal_states:

  Whether to compute terminal states (root and end points).

- compute_pseudotime:

  Whether to compute velocity pseudotime.

- compute_paga:

  Whether to compute PAGA (Partition-based graph abstraction).

- top_n:

  The number of top features to plot.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Paired"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc. Default is `"on data"`.

- show_plot:

  Whether to show the plot. Default is `FALSE`.

- save_plot:

  Whether to save plots to files. Default is `FALSE`.

- plot_format:

  Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.

- plot_dpi:

  Resolution (DPI) for saved plots. Default is `300`.

- plot_prefix:

  Prefix for saved plot filenames. Default is "cellrank".

- dirpath:

  The directory to save the plots. Default is `"./cellrank"`.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.
  Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md),
[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[RunPAGA](https://mengxu98.github.io/scop/reference/RunPAGA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSCVELO(
  pancreas_sub,
  assay_x = "RNA",
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)

FeatureDimPlot(
  pancreas_sub,
  c(
    "stochastic_length",
    "stochastic_confidence",
    "stochastic_pseudotime"
  )
)

VelocityPlot(
  pancreas_sub,
  reduction = "UMAP",
  plot_type = "stream"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = NA,
  velocity = "stochastic"
)

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSCVELO(
  pancreas_sub,
  assay_x = "RNA",
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  mode = c("deterministic", "stochastic"),
  filter_genes = TRUE,
  min_counts = 5,
  compute_velocity_confidence = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_paga = TRUE
)
} # }
```
