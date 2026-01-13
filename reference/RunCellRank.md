# Run CellRank analysis with kernel-estimator architecture

CellRank is a powerful toolkit for studying cellular dynamics using
Markov state modeling. This function implements the modern
kernel-estimator architecture recommended by CellRank, which provides
more flexibility and advanced features compared to the legacy API.

## Usage

``` r
RunCellRank(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group_by = NULL,
  cores = 1,
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
  stream_smooth = NULL,
  stream_density = 2,
  arrow_size = 5,
  arrow_length = 5,
  arrow_density = 0.5,
  calculate_velocity_genes = FALSE,
  denoise = FALSE,
  kinetics = FALSE,
  kernel_type = c("velocity", "pseudotime", "cytotrace"),
  time_key = "dpt_pseudotime",
  estimator_type = c("GPCCA", "CFLARE"),
  use_connectivity_kernel = TRUE,
  velocity_weight = 0.8,
  connectivity_weight = 0.2,
  softmax_scale = 4,
  n_macrostates = NULL,
  schur_method = c("krylov", "brandts"),
  n_cells_terminal = 10,
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "cellrank",
  legend.position = "on data",
  palette = "Paired",
  palcolor = NULL,
  dirpath = "./cellrank",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. Default is `NULL`. If provided, `adata` will be
  ignored.

- assay_x:

  Assay to convert in the anndata object.

- layer_x:

  Layer name for `assay_x` in the Seurat object.

- assay_y:

  Assay to convert in the anndata object.

- layer_y:

  Layer names for the `assay_y` in the Seurat object.

- adata:

  An anndata object. Default is `NULL`.

- group_by:

  Variable to use for grouping cells in the Seurat object.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- linear_reduction:

  Linear reduction method to use, e.g., `"PCA"`.

- nonlinear_reduction:

  Non-linear reduction method to use, e.g., `"UMAP"`.

- basis:

  The basis to use for reduction, e.g., `"UMAP"`.

- mode:

  Velocity estimation models to use. Can be `"deterministic"`,
  `"stochastic"`, or `"dynamical"`.

- fitting_by:

  Method used to fit gene velocities for dynamical modeling. Default is
  `"stochastic"`.

- magic_impute:

  Flag indicating whether to perform magic imputation. Default is
  `FALSE`.

- knn:

  The number of nearest neighbors for `magic.MAGIC`. Default is `5`.

- t:

  power to which the diffusion operator is powered for `magic.MAGIC`.
  Default is `2`.

- min_shared_counts:

  Minimum number of counts (both unspliced and spliced) required for a
  gene. Default is `30`.

- n_pcs:

  Number of principal components to use for linear reduction. Default is
  `30`.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph. Default is
  `30`.

- stream_smooth:

  Multiplication factor for scale in Gaussian kernel around grid point.

- stream_density:

  Controls the closeness of streamlines. Default is `2`.

- arrow_size:

  Size of arrows. Default is `5`.

- arrow_length:

  Length of arrows.

- arrow_density:

  Amount of velocities to show.

- calculate_velocity_genes:

  Boolean flag indicating whether to calculate velocity genes.

- denoise:

  Boolean flag indicating whether to denoise.

- kinetics:

  Boolean flag indicating whether to estimate RNA kinetics.

- kernel_type:

  Type of kernel to use: `"velocity"` (default, requires
  spliced/unspliced), `"pseudotime"` (requires pre-computed pseudotime
  or auto-computes DPT), or `"cytotrace"` (auto-computes CytoTRACE
  score, suitable for RNA-only data).

- time_key:

  Key in metadata for pseudotime. Used when
  `kernel_type = "pseudotime"`. If the key doesn't exist, DPT pseudotime
  will be computed automatically. Default is `"dpt_pseudotime"`.

- estimator_type:

  Type of estimator to use: `"GPCCA"` (default) or `"CFLARE"`. GPCCA
  provides coarse-grained analysis and Schur decomposition.

- use_connectivity_kernel:

  Whether to combine the main kernel with ConnectivityKernel. Default is
  `TRUE`.

- velocity_weight:

  Weight for the VelocityKernel when combining with ConnectivityKernel.
  Default is `0.8`.

- connectivity_weight:

  Weight for the ConnectivityKernel when combining with VelocityKernel.
  Default is `0.2`. Weights are automatically normalized to sum to
  `1.0`.

- softmax_scale:

  Scaling parameter for softmax transformation of velocity kernel.
  Default is `4`.

- n_macrostates:

  Number of macrostates to compute. If `NULL` (default), automatically
  determined based on eigenvalue spectrum.

- schur_method:

  Method for Schur decomposition: `"krylov"` or `"brandts"`. Only used
  for GPCCA estimator.

- n_cells_terminal:

  Minimum number of cells required for a state to be considered
  terminal. Default is `10`.

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

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc. Default is `"on data"`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Paired"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- dirpath:

  The directory to save the plots. Default is `"./cellrank"`.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.
  Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Returns a Seurat object if `return_seurat = TRUE` or an anndata object
with CellRank results stored in `obsm`, `obs`, and `varm` slots. The
estimator and kernel objects are stored in `srt@misc$cellrank`.

Lineage pseudotime columns are automatically added to metadata with
prefix `"Lineage_"` (e.g., `Lineage_Alpha`, `Lineage_Beta`), compatible
with
[LineagePlot](https://mengxu98.github.io/scop/reference/LineagePlot.md),
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
and
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md).
These are computed by combining base pseudotime (`cytotrace_pseudotime`,
`dpt_pseudotime`, or `latent_time`) with fate probabilities.

## See also

[RunSCVELO](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
[RunPAGA](https://mengxu98.github.io/scop/reference/RunPAGA.md),
[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunCellRank(
  srt = pancreas_sub,
  group_by = "SubCellType",
  cores = 6
)

CellDimPlot(
  pancreas_sub,
  group.by = "term_states_fwd",
  reduction = "umap",
  label = TRUE
)

FeatureDimPlot(
  pancreas_sub,
  features = "latent_time",
  reduction = "umap"
)

FeatureDimPlot(
  pancreas_sub,
  features = c("stochastic_confidence", "stochastic_length"),
  reduction = "umap"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = "cellrank_pseudotime",
  lineages_span = 0.1,
  lineages_trim = c(0.05, 0.95)
)

DynamicPlot(
  pancreas_sub,
  lineages = "cellrank_pseudotime",
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType"
)
} # }
```
