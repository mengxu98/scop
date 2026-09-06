# Run CellRank analysis

CellRank is a toolkit for studying cellular dynamics using Markov state
modeling.

## Usage

``` r
RunCellRank(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group.by = NULL,
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
  kernel_type = c("velocity", "pseudotime", "cytotrace", "wot"),
  time_key = "dpt_pseudotime",
  time_field = "Time",
  growth_iters = 3L,
  tmap_out = "tmaps/tmap_out",
  recalculate = FALSE,
  estimator_type = c("GPCCA", "CFLARE"),
  use_connectivity_kernel = TRUE,
  velocity_weight = 0.8,
  connectivity_weight = 0.2,
  softmax_scale = 4,
  n_macrostates = NULL,
  schur_method = c("brandts", "krylov"),
  schur_n_components = NULL,
  n_cells_terminal = 10,
  terminal_states = NULL,
  terminal_state_agg = c("top_n", "union"),
  driver_lineages = NULL,
  compute_lineage_drivers = TRUE,
  backward = FALSE,
  backend = c("python", "cpp"),
  allow_approximate = FALSE,
  max_dense_gib = 8,
  show_plot = TRUE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "cellrank",
  legend.position = "on data",
  palette = "Chinese",
  palcolor = NULL,
  dirpath = "./cellrank",
  envname = NULL,
  conda = "auto",
  recompute_neighbors = TRUE,
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

- cores:

  The number of cores to use for `cellrank`.

- linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

- basis:

  The basis to use for reduction, e.g., `"UMAP"`.

- mode:

  Velocity estimation models to use. Can be `"deterministic"`,
  `"stochastic"`, or `"dynamical"`. If the corresponding velocity
  reduction (e.g. `stochastic_umap`) is not found, falls back to the
  scVelo convention (e.g. `velocity_umap`), which is what an object
  converted from an AnnData (via
  [adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md))
  contains.

- fitting_by:

  Method used to fit gene velocities for dynamical modeling.

- magic_impute:

  Flag indicating whether to perform magic imputation.

- knn:

  The number of nearest neighbors for `magic.MAGIC`.

- t:

  Power to which the diffusion operator is powered for `magic.MAGIC`.

- min_shared_counts:

  Minimum number of counts (both unspliced and spliced) required for a
  gene.

- n_pcs:

  Number of principal components to use for linear reduction.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph.

- stream_smooth:

  Multiplication factor for scale in Gaussian kernel around grid point.

- stream_density:

  Controls the closeness of streamlines. When density = 2 (default), the
  domain is divided into a 60x60 grid, whereas density linearly scales
  this grid. Each cell in the grid can have, at most, one traversing
  streamline.

- arrow_size:

  Size of arrows.

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
  or auto-computes DPT), `"cytotrace"` (auto-computes CytoTRACE score,
  suitable for RNA-only data), or `"wot"` (uses Waddington-OT transport
  maps through CellRank's RealTimeKernel).

- time_key:

  Key in metadata for pseudotime. Used when
  `kernel_type = "pseudotime"`. If the key doesn't exist, DPT pseudotime
  will be computed automatically.

- time_field:

  Key in metadata for experimental time. Used when
  `kernel_type = "wot"`.

- growth_iters:

  Number of growth iterations passed to `wot.ot.OTModel`.

- tmap_out:

  Directory used to store or read Waddington-OT transport maps.

- recalculate:

  Whether to recompute Waddington-OT transport maps even when `tmap_out`
  already exists.

- estimator_type:

  Type of estimator to use: `"GPCCA"` (default) or `"CFLARE"`. GPCCA
  provides coarse-grained analysis and Schur decomposition.

- use_connectivity_kernel:

  Whether to combine the main kernel with ConnectivityKernel.

- velocity_weight:

  Weight for the VelocityKernel when combining with ConnectivityKernel.

- connectivity_weight:

  Weight for the ConnectivityKernel when combining with VelocityKernel.
  Weights are automatically normalized to sum to `1.0`.

- softmax_scale:

  Scaling parameter for softmax transformation of velocity kernel.

- n_macrostates:

  Number of macrostates to compute. If `NULL` (default), automatically
  determined based on eigenvalue spectrum.

- schur_method:

  Method for Schur decomposition: `"brandts"` (default) or `"krylov"`.
  Only used for the GPCCA estimator. The `"krylov"` method requires the
  optional Python packages `petsc4py` and `slepc4py`.

- schur_n_components:

  Number of Schur components. If `NULL`, retain the existing size
  heuristic.

- n_cells_terminal:

  Minimum number of cells required for a state to be considered
  terminal.

- terminal_states:

  Optional macrostate names to set as terminal states. Names are
  validated against the computed macrostates before fate estimation.

- terminal_state_agg:

  Aggregation policy for combined terminal states.

- driver_lineages:

  Optional lineage names for driver-gene computation.

- compute_lineage_drivers:

  Whether to compute and store lineage drivers.

- backward:

  Whether to compute backward transitions.

- backend:

  Backend for computation: `"python"` (default) or `"cpp"`. The C++ path
  is an explicitly opted-in approximation and does not reproduce the
  complete CellRank estimator. The selected backend is stable: the C++
  path never switches to Python because of an argument value.

- allow_approximate:

  Whether to allow the approximate C++ path. This must be `TRUE` when
  `backend = "cpp"`.

- max_dense_gib:

  Maximum estimated GiB allowed for the dense cell-by-cell working
  matrices used by the C++ path.

- show_plot:

  Whether to show the plot.

- save_plot:

  Whether to save plots to files.

- plot_format:

  Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.

- plot_dpi:

  Resolution (DPI) for saved plots.

- plot_prefix:

  Prefix for saved plot filenames.

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- dirpath:

  The directory to save the plots.

- envname:

  Optional Python environment name. `NULL` uses the current SCOP
  environment selection.

- conda:

  Conda-compatible executable used by
  [PrepareEnv](https://mengxu98.github.io/scop/reference/PrepareEnv.md).

- recompute_neighbors:

  Whether to rebuild Scanpy neighbors before a pseudotime kernel. If
  `FALSE`, the existing graph must be complete and match the current
  AnnData cell order.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Returns a Seurat object if `return_seurat = TRUE` or an anndata object
with CellRank results stored in `obsm`, `obs`, and `varm` slots. The
estimator and kernel objects are stored in `srt@misc$cellrank`.

## See also

[RunSCVELO](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
[RunPAGA](https://mengxu98.github.io/scop/reference/RunPAGA.md),
[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md),
[RunCellRankTrends](https://mengxu98.github.io/scop/reference/RunCellRankTrends.md),
[RunCellRankEnrichment](https://mengxu98.github.io/scop/reference/RunCellRankEnrichment.md),
[CellRankPlot](https://mengxu98.github.io/scop/reference/CellRankPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunCellRank(
  srt = pancreas_sub,
  group.by = "SubCellType"
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
