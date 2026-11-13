# Run CellRank analysis

CellRank is a toolkit for studying cellular dynamics using Markov state
modeling.

## Usage

``` r
RunCellRank(
  object = NULL,
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
  kernel_type = c("velocity", "pseudotime", "cytotrace", "wot", "moscot"),
  time_key = "dpt_pseudotime",
  time_field = "Time",
  time_values = NULL,
  moscot_args = list(),
  realtime_args = list(),
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
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A Seurat object. If provided, `adata` will be ignored.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.

- layer_x:

  Layer name for assay_x in the Seurat object.

- assay_y:

  Assays, or extra layers of `assay_x`, to convert as layers in the
  anndata object. Each name is matched against the assays of the object
  first, then against the layers of `assay_x`, which covers velocity
  matrices such as `spliced` and `unspliced` stored as layers of a
  Seurat v5 `RNA` assay.

- layer_y:

  Layer names for the assays in `assay_y`. It is ignored for names that
  are matched as layers of `assay_x`, where the layer name itself is
  used.

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
  maps through CellRank's RealTimeKernel), or `"moscot"` (uses moscot's
  TemporalProblem through RealTimeKernel).

- time_key:

  Key in metadata for pseudotime. Used when
  `kernel_type = "pseudotime"`. If the key doesn't exist, DPT pseudotime
  will be computed automatically.

- time_field:

  Key in metadata for experimental time. Used when
  `kernel_type %in% c("wot", "moscot")`.

- time_values:

  Optional named mapping from non-numeric experimental-time labels to
  numeric values used to order moscot time points. Numeric columns do
  not need a mapping.

- moscot_args:

  Named lists for moscot's `growth`, `prepare`, and `solve` stages. Only
  used when `kernel_type = "moscot"`.

- realtime_args:

  Named lists for `RealTimeKernel`'s `from_moscot` and `transition`
  stages. Only used when `kernel_type = "moscot"`.

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
  is an approximate implementation and issues a warning that results may
  not be completely identical to the reference workflow. The selected
  computation path is stable: the C++ path never switches to Python
  because of an argument value.

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

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

Returns a Seurat object if `return_seurat = TRUE` or an anndata object
with CellRank results stored in `obsm`, `obs`, and `varm` slots. The
estimator and kernel objects are stored in `srt@misc$cellrank`.

## Details

Use
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
for stored fate probabilities and
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
for terminal states.
[CellRankPlot](https://mengxu98.github.io/scop/reference/CellRankPlot.md)
adds circular fate coordinates, projected transition vectors, random
walks, and experimental-time probability flow. A flow plot requires an
experimental-time result from the moscot workflow.

A stored pseudotime column can also be used with
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
and
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md),
as for other trajectory methods. These functions fit expression along
pseudotime using the common dynamic analysis workflow.
[RunCellRankTrends](https://mengxu98.github.io/scop/reference/RunCellRankTrends.md)
instead fits CellRank's GAM with lineage-probability weights; its stored
curves are a separate analysis.
[RunCellRankEnrichment](https://mengxu98.github.io/scop/reference/RunCellRankEnrichment.md)
tests over-representation in those gene modules, whereas
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)
fits changing gene-set scores along pseudotime.

## See also

[RunSCVELO](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
[RunPAGA](https://mengxu98.github.io/scop/reference/RunPAGA.md),
[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md),
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md),
[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
[RunCellRankTrends](https://mengxu98.github.io/scop/reference/RunCellRankTrends.md),
[RunCellRankEnrichment](https://mengxu98.github.io/scop/reference/RunCellRankEnrichment.md),
[CellRankPlot](https://mengxu98.github.io/scop/reference/CellRankPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-11-13 10:23:16] Start standard processing workflow...
#> ℹ [2026-11-13 10:23:16] Checking a list of <Seurat>...
#> ! [2026-11-13 10:23:16] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-13 10:23:16] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-13 10:23:16] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-13 10:23:16] Use the separate HVF from `srt_list`
#> ℹ [2026-11-13 10:23:16] Number of available HVF: 2000
#> ℹ [2026-11-13 10:23:16] Finished check
#> ℹ [2026-11-13 10:23:16] Perform `ScaleData()`
#> ℹ [2026-11-13 10:23:16] Perform pca linear dimension reduction
#> ℹ [2026-11-13 10:23:17] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-13 10:23:17] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-13 10:23:17] Reorder clusters...
#> ℹ [2026-11-13 10:23:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-13 10:23:17] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-13 10:23:25] Standard processing workflow completed
reduction <- DefaultReduction(pancreas_sub, pattern = "UMAP")
pancreas_sub <- RunSCVELO(
  object = pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = reduction,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  backend = "cpp",
  show_plot = FALSE
)
#> ℹ [2026-11-13 10:23:25] Running scanpy-compatible preprocessing (15998 features -> filter + normalize)...
#> ℹ [2026-11-13 10:23:30] Running scVelo "stochastic" mode with `backend = 'cpp'` (9699 features)
#> ✔ [2026-11-13 10:23:34] scVelo "stochastic" mode completed
#> ✔ [2026-11-13 10:23:34] scVelo cpp backend completed
pancreas_sub <- RunCellRank(
  object = pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = reduction,
  kernel_type = "pseudotime",
  time_key = "stochastic_pseudotime",
  backend = "cpp",
  show_plot = FALSE
)
#> ! [2026-11-13 10:23:34] RunCellRank(backend = "cpp") uses an approximate C++ implementation; results may not be completely identical to the reference python workflow.
#> ✔ [2026-11-13 10:23:37] CellRank cpp backend ("gpcca" + "pseudotime_connectivity_combined" kernel) completed

CellDimPlot(
  pancreas_sub,
  group.by = "cellrank_terminal_states",
  reduction = reduction
)


FeatureDimPlot(
  pancreas_sub,
  features = "cellrank_fate_Lineage1",
  reduction = reduction
)


CellRankPlot(
  pancreas_sub,
  plot_type = "circular"
)
#> Scale for colour is already present.
#> Adding another scale for colour, which will replace the existing scale.


CellRankPlot(
  pancreas_sub,
  plot_type = "projection",
  reduction = reduction,
  group.by = "SubCellType"
)


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = "stochastic_pseudotime",
  features = c("Rbp4", "Sst", "Pyy"),
  assay = "RNA"
)
#> ℹ [2026-11-13 10:23:38] Start find dynamic features
#> ℹ [2026-11-13 10:23:38] Data type is raw counts
#> ℹ [2026-11-13 10:23:38] Number of candidate features (union): 3
#> ℹ [2026-11-13 10:23:38] Data type is raw counts
#> ℹ [2026-11-13 10:23:38] Calculating dynamic features for "stochastic_pseudotime"...
#> ℹ [2026-11-13 10:23:38] Using 1 core
#> ⠙ [2026-11-13 10:23:38] Running for 1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-11-13 10:23:38] Completed 3 tasks in 111ms
#> 
#> ℹ [2026-11-13 10:23:38] Building results
#> ✔ [2026-11-13 10:23:38] Find dynamic features done

DynamicHeatmap(
  pancreas_sub,
  lineages = "stochastic_pseudotime",
  features = c("Rbp4", "Sst", "Pyy"),
  use_fitted = TRUE,
  assay = "RNA",
  cell_annotation = "SubCellType"
)
#> $plot

#> 
#> $matrix
#>      TCAATCTGTCAAAGATstochastic_pseudotime
#> Pyy                             -0.7947626
#> Rbp4                            -0.7885953
#> Sst                             -0.3342073
#>      CTTGGCTAGACAGACCstochastic_pseudotime
#> Pyy                             -0.7946591
#> Rbp4                            -0.7878362
#> Sst                             -0.3341013
#>      GTATTCTCACTTGGATstochastic_pseudotime
#> Pyy                             -0.7944789
#> Rbp4                            -0.7864836
#> Sst                             -0.3339125
#>      AGATTGCAGCTGTCTAstochastic_pseudotime
#> Pyy                             -0.7944413
#> Rbp4                            -0.7861966
#> Sst                             -0.3338725
#>      CCATTCGCATGCCTAAstochastic_pseudotime
#> Pyy                             -0.7944270
#> Rbp4                            -0.7860874
#> Sst                             -0.3338573
#>      GTATTCTGTCCGTCAGstochastic_pseudotime
#> Pyy                             -0.7943137
#> Rbp4                            -0.7852113
#> Sst                             -0.3337353
#>      TGGACGCTCAGGTAAAstochastic_pseudotime
#> Pyy                             -0.7942202
#> Rbp4                            -0.7844782
#> Sst                             -0.3336334
#>      ATCATCTTCTGTCAAGstochastic_pseudotime
#> Pyy                             -0.7942075
#> Rbp4                            -0.7843778
#> Sst                             -0.3336194
#>      TTCTACATCTGCAGTAstochastic_pseudotime
#> Pyy                             -0.7941766
#> Rbp4                            -0.7841332
#> Sst                             -0.3335854
#>      TATCTCAGTTCGCTAAstochastic_pseudotime
#> Pyy                             -0.7941388
#> Rbp4                            -0.7838327
#> Sst                             -0.3335437
#>      CAGCCGAAGTAGGTGCstochastic_pseudotime
#> Pyy                             -0.7940973
#> Rbp4                            -0.7835008
#> Sst                             -0.3334976
#>      AGTCTTTCACCAGGTCstochastic_pseudotime
#> Pyy                             -0.7940717
#> Rbp4                            -0.7832951
#> Sst                             -0.3334690
#>      ATTACTCAGTTGCAGGstochastic_pseudotime
#> Pyy                             -0.7940451
#> Rbp4                            -0.7830814
#> Sst                             -0.3334394
#>      CAGCAGCCACACTGCGstochastic_pseudotime
#> Pyy                             -0.7940392
#> Rbp4                            -0.7830337
#> Sst                             -0.3334328
#>      AGTGAGGGTTACGGAGstochastic_pseudotime
#> Pyy                             -0.7939977
#> Rbp4                            -0.7826979
#> Sst                             -0.3333862
#>      TGAGCATTCGGGAGTAstochastic_pseudotime
#> Pyy                             -0.7939542
#> Rbp4                            -0.7823444
#> Sst                             -0.3333372
#>      AACTCAGGTTGGAGGTstochastic_pseudotime
#> Pyy                             -0.7939306
#> Rbp4                            -0.7821514
#> Sst                             -0.3333104
#>      GGACGTCGTTTAGCTGstochastic_pseudotime
#> Pyy                             -0.7939279
#> Rbp4                            -0.7821298
#> Sst                             -0.3333074
#>      CTAGCCTGTTCGGGCTstochastic_pseudotime
#> Pyy                             -0.7939158
#> Rbp4                            -0.7820309
#> Sst                             -0.3332937
#>      TAGTGGTCAACTGCGCstochastic_pseudotime
#> Pyy                             -0.7939141
#> Rbp4                            -0.7820165
#> Sst                             -0.3332917
#>      GAATGAACAGCCTTGGstochastic_pseudotime
#> Pyy                             -0.7938604
#> Rbp4                            -0.7815756
#> Sst                             -0.3332306
#>      TGGCCAGAGCTACCTAstochastic_pseudotime
#> Pyy                             -0.7938385
#> Rbp4                            -0.7813951
#> Sst                             -0.3332057
#>      ATTACTCTCAGGCGAAstochastic_pseudotime
#> Pyy                             -0.7937882
#> Rbp4                            -0.7809795
#> Sst                             -0.3331481
#>      GATGAAAAGTTGTAGAstochastic_pseudotime
#> Pyy                             -0.7937525
#> Rbp4                            -0.7806827
#> Sst                             -0.3331071
#>      CCGTACTTCAATAAGGstochastic_pseudotime
#> Pyy                             -0.7937501
#> Rbp4                            -0.7806634
#> Sst                             -0.3331044
#>      TTAACTCCAAAGAATCstochastic_pseudotime
#> Pyy                             -0.7937429
#> Rbp4                            -0.7806033
#> Sst                             -0.3330961
#>      TATGCCCTCTTACCTAstochastic_pseudotime
#> Pyy                             -0.7937354
#> Rbp4                            -0.7805404
#> Sst                             -0.3330874
#>      CCATTCGGTAGGGACTstochastic_pseudotime
#> Pyy                             -0.7937172
#> Rbp4                            -0.7803883
#> Sst                             -0.3330664
#>      TAGGCATTCCAAAGTCstochastic_pseudotime
#> Pyy                             -0.7934874
#> Rbp4                            -0.7784514
#> Sst                             -0.3327990
#>      CTACACCGTCGCCATGstochastic_pseudotime
#> Pyy                             -0.7934787
#> Rbp4                            -0.7783773
#> Sst                             -0.3327888
#>      AACTGGTGTGCCTGTGstochastic_pseudotime
#> Pyy                             -0.7934122
#> Rbp4                            -0.7778095
#> Sst                             -0.3327106
#>      TGCGCAGGTACGCACCstochastic_pseudotime
#> Pyy                             -0.7933471
#> Rbp4                            -0.7772512
#> Sst                             -0.3326338
#>      CTCATTACAAGCGCTCstochastic_pseudotime
#> Pyy                             -0.7933020
#> Rbp4                            -0.7768636
#> Sst                             -0.3325804
#>      TGACTTTAGCGTTGCCstochastic_pseudotime
#> Pyy                             -0.7932959
#> Rbp4                            -0.7768112
#> Sst                             -0.3325732
#>      GATCGTAGTCCTAGCGstochastic_pseudotime
#> Pyy                             -0.7932740
#> Rbp4                            -0.7766219
#> Sst                             -0.3325472
#>      GGACAGATCCTGCAGGstochastic_pseudotime
#> Pyy                             -0.7932658
#> Rbp4                            -0.7765512
#> Sst                             -0.3325375
#>      TTTATGCGTCGAACAGstochastic_pseudotime
#> Pyy                             -0.7931538
#> Rbp4                            -0.7755825
#> Sst                             -0.3324045
#>      CGAACATAGCCAGTAGstochastic_pseudotime
#> Pyy                             -0.7931500
#> Rbp4                            -0.7755492
#> Sst                             -0.3323999
#>      AAGCCGCAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.7931280
#> Rbp4                            -0.7753581
#> Sst                             -0.3323737
#>      TTTATGCTCCATTCTAstochastic_pseudotime
#> Pyy                             -0.7931065
#> Rbp4                            -0.7751712
#> Sst                             -0.3323481
#>      TAGAGCTGTTGTCGCGstochastic_pseudotime
#> Pyy                             -0.7930835
#> Rbp4                            -0.7749717
#> Sst                             -0.3323207
#>      TAGAGCTGTATGAAACstochastic_pseudotime
#> Pyy                             -0.7930165
#> Rbp4                            -0.7743888
#> Sst                             -0.3322409
#>      CGTGTAATCCATGCTCstochastic_pseudotime
#> Pyy                             -0.7929849
#> Rbp4                            -0.7741134
#> Sst                             -0.3322032
#>      TGACTAGTCCACGTGGstochastic_pseudotime
#> Pyy                             -0.7929741
#> Rbp4                            -0.7740192
#> Sst                             -0.3321903
#>      TGCGCAGTCTGCCAGGstochastic_pseudotime
#> Pyy                             -0.7929658
#> Rbp4                            -0.7739468
#> Sst                             -0.3321804
#>      TGAGAGGGTAAGGGCTstochastic_pseudotime
#> Pyy                             -0.7928873
#> Rbp4                            -0.7732616
#> Sst                             -0.3320867
#>      AACTCAGTCTGGTTCCstochastic_pseudotime
#> Pyy                             -0.7928576
#> Rbp4                            -0.7730019
#> Sst                             -0.3320513
#>      CTCTGGTAGTCCCACGstochastic_pseudotime
#> Pyy                             -0.7928099
#> Rbp4                            -0.7725859
#> Sst                             -0.3319945
#>      ACAGCTATCTGGGCCAstochastic_pseudotime
#> Pyy                             -0.7928008
#> Rbp4                            -0.7725061
#> Sst                             -0.3319836
#>      CGGTTAAGTAAACACAstochastic_pseudotime
#> Pyy                             -0.7927946
#> Rbp4                            -0.7724520
#> Sst                             -0.3319762
#>      AGGCCGTGTACCAGTTstochastic_pseudotime
#> Pyy                             -0.7927938
#> Rbp4                            -0.7724451
#> Sst                             -0.3319753
#>      CGCCAAGCACCGGAAAstochastic_pseudotime
#> Pyy                             -0.7927434
#> Rbp4                            -0.7720055
#> Sst                             -0.3319154
#>      CGGTTAAAGAGCTGCAstochastic_pseudotime
#> Pyy                             -0.7927217
#> Rbp4                            -0.7718164
#> Sst                             -0.3318896
#>      ATAAGAGCATCGATTGstochastic_pseudotime
#> Pyy                             -0.7926670
#> Rbp4                            -0.7713396
#> Sst                             -0.3318247
#>      GTAACTGTCAGTGCATstochastic_pseudotime
#> Pyy                             -0.7926310
#> Rbp4                            -0.7710273
#> Sst                             -0.3317822
#>      CCTAGCTTCGATGAGGstochastic_pseudotime
#> Pyy                             -0.7926006
#> Rbp4                            -0.7707626
#> Sst                             -0.3317463
#>      GATGCTAGTAGAAGGAstochastic_pseudotime
#> Pyy                             -0.7925867
#> Rbp4                            -0.7706418
#> Sst                             -0.3317299
#>      TGTGGTACAAGACGTGstochastic_pseudotime
#> Pyy                             -0.7925766
#> Rbp4                            -0.7705544
#> Sst                             -0.3317180
#>      TCGTACCCACGCATCGstochastic_pseudotime
#> Pyy                             -0.7925713
#> Rbp4                            -0.7705087
#> Sst                             -0.3317118
#>      CTGCGGATCAGCATGTstochastic_pseudotime
#> Pyy                             -0.7925529
#> Rbp4                            -0.7703495
#> Sst                             -0.3316902
#>      CAGCCGAAGCGATATAstochastic_pseudotime
#> Pyy                             -0.7925506
#> Rbp4                            -0.7703290
#> Sst                             -0.3316874
#>      CAGCAGCAGACAAAGGstochastic_pseudotime
#> Pyy                             -0.7925450
#> Rbp4                            -0.7702808
#> Sst                             -0.3316809
#>      TCTGAGAGTCAAACTCstochastic_pseudotime
#> Pyy                             -0.7924519
#> Rbp4                            -0.7694780
#> Sst                             -0.3315721
#>      CTAGAGTCACTGTTAGstochastic_pseudotime
#> Pyy                             -0.7924476
#> Rbp4                            -0.7694406
#> Sst                             -0.3315670
#>      ACGGGCTGTGCGATAGstochastic_pseudotime
#> Pyy                             -0.7923887
#> Rbp4                            -0.7689357
#> Sst                             -0.3314987
#>      ACCGTAAGTTATCACGstochastic_pseudotime
#> Pyy                             -0.7923716
#> Rbp4                            -0.7687896
#> Sst                             -0.3314790
#>      AAGACCTCAGGTTTCAstochastic_pseudotime
#> Pyy                             -0.7923520
#> Rbp4                            -0.7686228
#> Sst                             -0.3314565
#>      TGTATTCCAATGACCTstochastic_pseudotime
#> Pyy                             -0.7923258
#> Rbp4                            -0.7684005
#> Sst                             -0.3314265
#>      CCTAGCTCATACTACGstochastic_pseudotime
#> Pyy                             -0.7923167
#> Rbp4                            -0.7683235
#> Sst                             -0.3314161
#>      CACACTCAGATCCCATstochastic_pseudotime
#> Pyy                             -0.7923034
#> Rbp4                            -0.7682104
#> Sst                             -0.3314009
#>      ACTGATGGTGAGGGTTstochastic_pseudotime
#> Pyy                             -0.7922760
#> Rbp4                            -0.7679795
#> Sst                             -0.3313697
#>      CAGCTAAAGTCCATACstochastic_pseudotime
#> Pyy                             -0.7922398
#> Rbp4                            -0.7676755
#> Sst                             -0.3313288
#>      TGCTGCTGTCCGTCAGstochastic_pseudotime
#> Pyy                             -0.7922387
#> Rbp4                            -0.7676663
#> Sst                             -0.3313276
#>      CACATAGAGCGTGTCCstochastic_pseudotime
#> Pyy                             -0.7922358
#> Rbp4                            -0.7676419
#> Sst                             -0.3313243
#>      AGGGAGTAGGGTATCGstochastic_pseudotime
#> Pyy                             -0.7922257
#> Rbp4                            -0.7675578
#> Sst                             -0.3313130
#>      CAGCGACCACAGGAGTstochastic_pseudotime
#> Pyy                             -0.7922159
#> Rbp4                            -0.7674755
#> Sst                             -0.3313019
#>      AGAGCTTGTTAAAGACstochastic_pseudotime
#> Pyy                             -0.7922028
#> Rbp4                            -0.7673669
#> Sst                             -0.3312874
#>      GCTGCTTGTCTTCAAGstochastic_pseudotime
#> Pyy                             -0.7921841
#> Rbp4                            -0.7672110
#> Sst                             -0.3312664
#>      TGGCGCATCACAACGTstochastic_pseudotime
#> Pyy                             -0.7921464
#> Rbp4                            -0.7669002
#> Sst                             -0.3312247
#>      GTCAAGTGTGGTCCGTstochastic_pseudotime
#> Pyy                             -0.7921446
#> Rbp4                            -0.7668848
#> Sst                             -0.3312227
#>      ACGTCAACAGACAAATstochastic_pseudotime
#> Pyy                             -0.7921061
#> Rbp4                            -0.7665694
#> Sst                             -0.3311804
#>      GACTAACTCAGTCCCTstochastic_pseudotime
#> Pyy                             -0.7921026
#> Rbp4                            -0.7665408
#> Sst                             -0.3311766
#>      TTGACTTAGACAAGCCstochastic_pseudotime
#> Pyy                             -0.7919467
#> Rbp4                            -0.7652994
#> Sst                             -0.3310111
#>      CGATGGCTCTGTGCAAstochastic_pseudotime
#> Pyy                             -0.7919284
#> Rbp4                            -0.7651577
#> Sst                             -0.3309923
#>      GAGCAGAAGCGTGTCCstochastic_pseudotime
#> Pyy                             -0.7919091
#> Rbp4                            -0.7650094
#> Sst                             -0.3309727
#>      ACAGCCGAGTGGTCCCstochastic_pseudotime
#> Pyy                             -0.7918623
#> Rbp4                            -0.7646549
#> Sst                             -0.3309258
#>      AGAATAGCAGTACACTstochastic_pseudotime
#> Pyy                             -0.7918144
#> Rbp4                            -0.7643015
#> Sst                             -0.3308793
#>      AGAGCGATCTCTTGATstochastic_pseudotime
#> Pyy                             -0.7917567
#> Rbp4                            -0.7638882
#> Sst                             -0.3308251
#>      ATTATCCAGTGAACATstochastic_pseudotime
#> Pyy                             -0.7916505
#> Rbp4                            -0.7631756
#> Sst                             -0.3307327
#>      CACACTCTCATGTGGTstochastic_pseudotime
#> Pyy                             -0.7916072
#> Rbp4                            -0.7629072
#> Sst                             -0.3306983
#>      CGTTAGAAGTATCTCGstochastic_pseudotime
#> Pyy                             -0.7915897
#> Rbp4                            -0.7628034
#> Sst                             -0.3306851
#>      TCTATTGCACAGGAGTstochastic_pseudotime
#> Pyy                             -0.7914409
#> Rbp4                            -0.7620581
#> Sst                             -0.3305933
#>      CTTCTCTAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.7914011
#> Rbp4                            -0.7619191
#> Sst                             -0.3305777
#>      CCTACACTCCTTTACAstochastic_pseudotime
#> Pyy                             -0.7913237
#> Rbp4                            -0.7617790
#> Sst                             -0.3305665
#>      CTGCCTAGTCTCCATCstochastic_pseudotime
#> Pyy                             -0.7912836
#> Rbp4                            -0.7618188
#> Sst                             -0.3305772
#>      GGAAAGCGTCTTGATGstochastic_pseudotime
#> Pyy                             -0.7912790
#> Rbp4                            -0.7618316
#> Sst                             -0.3305797
#>      CCGTGGATCCAAGCCGstochastic_pseudotime
#> Pyy                             -0.7912727
#> Rbp4                            -0.7618523
#> Sst                             -0.3305835
#>      ATCCACCGTACCATCAstochastic_pseudotime
#> Pyy                             -0.7912558
#> Rbp4                            -0.7619328
#> Sst                             -0.3305974
#>      AGCTTGACAGTAAGATstochastic_pseudotime
#> Pyy                             -0.7912317
#> Rbp4                            -0.7621424
#> Sst                             -0.3306311
#>      GTGCATAAGCTCCTCTstochastic_pseudotime
#> Pyy                             -0.7912120
#> Rbp4                            -0.7625115
#> Sst                             -0.3306878
#>      CAACCAATCGGGAGTAstochastic_pseudotime
#> Pyy                             -0.7912062
#> Rbp4                            -0.7627300
#> Sst                             -0.3307206
#>      CCTCTGAGTATAAACGstochastic_pseudotime
#> Pyy                             -0.7912027
#> Rbp4                            -0.7629853
#> Sst                             -0.3307585
#>      ACAGCTATCAGCATGTstochastic_pseudotime
#> Pyy                             -0.7912017
#> Rbp4                            -0.7631754
#> Sst                             -0.3307865
#>      GTGTTAGCACCATCCTstochastic_pseudotime
#> Pyy                             -0.7912020
#> Rbp4                            -0.7634158
#> Sst                             -0.3308218
#>      CATCAGACATGCATGTstochastic_pseudotime
#> Pyy                             -0.7912020
#> Rbp4                            -0.7634344
#> Sst                             -0.3308245
#>      CACACAACAGCCAGAAstochastic_pseudotime
#> Pyy                             -0.7912043
#> Rbp4                            -0.7637187
#> Sst                             -0.3308660
#>      GGGCATCCATGACATCstochastic_pseudotime
#> Pyy                             -0.7912054
#> Rbp4                            -0.7638135
#> Sst                             -0.3308798
#>      GGAATAAGTTCGGGCTstochastic_pseudotime
#> Pyy                             -0.7912078
#> Rbp4                            -0.7639866
#> Sst                             -0.3309049
#>      CATATGGCATGTTGACstochastic_pseudotime
#> Pyy                             -0.7912202
#> Rbp4                            -0.7646132
#> Sst                             -0.3309954
#>      CTACACCGTCGAATCTstochastic_pseudotime
#> Pyy                             -0.7912544
#> Rbp4                            -0.7657749
#> Sst                             -0.3311618
#>      CCCATACCACATCCAAstochastic_pseudotime
#> Pyy                             -0.7912670
#> Rbp4                            -0.7661295
#> Sst                             -0.3312123
#>      GGAGCAAAGCTTCGCGstochastic_pseudotime
#> Pyy                             -0.7913015
#> Rbp4                            -0.7670027
#> Sst                             -0.3313364
#>      CTAGCCTGTCTCTTATstochastic_pseudotime
#> Pyy                             -0.7913230
#> Rbp4                            -0.7674993
#> Sst                             -0.3314067
#>      GAGTCCGCAATTGCTGstochastic_pseudotime
#> Pyy                             -0.7913445
#> Rbp4                            -0.7679696
#> Sst                             -0.3314731
#>      CGGAGCTCATTGGGCCstochastic_pseudotime
#> Pyy                             -0.7913526
#> Rbp4                            -0.7681428
#> Sst                             -0.3314975
#>      AGGGAGTAGTGTCTCAstochastic_pseudotime
#> Pyy                             -0.7913922
#> Rbp4                            -0.7689522
#> Sst                             -0.3316108
#>      AGCGGTCGTGTATGGGstochastic_pseudotime
#> Pyy                             -0.7914226
#> Rbp4                            -0.7695455
#> Sst                             -0.3316931
#>      GCGCAACTCACAACGTstochastic_pseudotime
#> Pyy                             -0.7914495
#> Rbp4                            -0.7700507
#> Sst                             -0.3317626
#>      TTCTCAAAGCCGCCTAstochastic_pseudotime
#> Pyy                             -0.7915053
#> Rbp4                            -0.7710614
#> Sst                             -0.3318999
#>      GAACATCGTAATCGTCstochastic_pseudotime
#> Pyy                             -0.7915167
#> Rbp4                            -0.7712606
#> Sst                             -0.3319267
#>      CGAGCACAGCTGAAATstochastic_pseudotime
#> Pyy                             -0.7915251
#> Rbp4                            -0.7714076
#> Sst                             -0.3319464
#>      TGGTTCCGTCACCTAAstochastic_pseudotime
#> Pyy                             -0.7915308
#> Rbp4                            -0.7715068
#> Sst                             -0.3319596
#>      CGTCTACAGATCTGCTstochastic_pseudotime
#> Pyy                             -0.7915428
#> Rbp4                            -0.7717122
#> Sst                             -0.3319870
#>      ACGTCAACACTGCCAGstochastic_pseudotime
#> Pyy                             -0.7915919
#> Rbp4                            -0.7725403
#> Sst                             -0.3320962
#>      CTGAAACTCGCCGTGAstochastic_pseudotime
#> Pyy                             -0.7916197
#> Rbp4                            -0.7729966
#> Sst                             -0.3321555
#>      CTACCCAGTGCGCTTGstochastic_pseudotime
#> Pyy                             -0.7916273
#> Rbp4                            -0.7731216
#> Sst                             -0.3321716
#>      CTTAGGACACAGGTTTstochastic_pseudotime
#> Pyy                             -0.7916309
#> Rbp4                            -0.7731787
#> Sst                             -0.3321789
#>      CGATGGCGTTACGGAGstochastic_pseudotime
#> Pyy                             -0.7917315
#> Rbp4                            -0.7747631
#> Sst                             -0.3323787
#>      TACGGATAGAGCTGGTstochastic_pseudotime
#> Pyy                             -0.7917324
#> Rbp4                            -0.7747782
#> Sst                             -0.3323806
#>      TGGGAAGCATGACATCstochastic_pseudotime
#> Pyy                             -0.7917684
#> Rbp4                            -0.7753256
#> Sst                             -0.3324475
#>      TAGCCGGTCAACACACstochastic_pseudotime
#> Pyy                             -0.7917885
#> Rbp4                            -0.7756265
#> Sst                             -0.3324838
#>      TGACGGCAGGAATGGAstochastic_pseudotime
#> Pyy                             -0.7917933
#> Rbp4                            -0.7756991
#> Sst                             -0.3324925
#>      CGGACACAGTGTCTCAstochastic_pseudotime
#> Pyy                             -0.7918140
#> Rbp4                            -0.7760062
#> Sst                             -0.3325292
#>      TCTCTAAGTTACCGATstochastic_pseudotime
#> Pyy                             -0.7918243
#> Rbp4                            -0.7761572
#> Sst                             -0.3325470
#>      TGGCCAGGTCGGCATCstochastic_pseudotime
#> Pyy                             -0.7918928
#> Rbp4                            -0.7771522
#> Sst                             -0.3326622
#>      ACTGATGTCCCTAACCstochastic_pseudotime
#> Pyy                             -0.7918958
#> Rbp4                            -0.7771944
#> Sst                             -0.3326670
#>      TTCCCAGGTATTCTCTstochastic_pseudotime
#> Pyy                             -0.7919333
#> Rbp4                            -0.7777277
#> Sst                             -0.3327268
#>      CCTAGCTCAGACTCGCstochastic_pseudotime
#> Pyy                             -0.7919344
#> Rbp4                            -0.7777429
#> Sst                             -0.3327285
#>      AAGGCAGAGCTGAAATstochastic_pseudotime
#> Pyy                             -0.7919382
#> Rbp4                            -0.7777969
#> Sst                             -0.3327344
#>      CGGGTCAAGGCTCAGAstochastic_pseudotime
#> Pyy                             -0.7919502
#> Rbp4                            -0.7779656
#> Sst                             -0.3327530
#>      ATCATGGTCCAGAGGAstochastic_pseudotime
#> Pyy                             -0.7920186
#> Rbp4                            -0.7789120
#> Sst                             -0.3328542
#>      CTGATCCTCAGTTGACstochastic_pseudotime
#> Pyy                             -0.7920275
#> Rbp4                            -0.7790345
#> Sst                             -0.3328670
#>      AGATTGCCACCAGGCTstochastic_pseudotime
#> Pyy                             -0.7920306
#> Rbp4                            -0.7790771
#> Sst                             -0.3328714
#>      TCTGGAACACACGCTGstochastic_pseudotime
#> Pyy                             -0.7920694
#> Rbp4                            -0.7796031
#> Sst                             -0.3329249
#>      CAGCAGCGTTCACGGCstochastic_pseudotime
#> Pyy                             -0.7920700
#> Rbp4                            -0.7796110
#> Sst                             -0.3329257
#>      TTCTCCTCAATGGTCTstochastic_pseudotime
#> Pyy                             -0.7920772
#> Rbp4                            -0.7797074
#> Sst                             -0.3329353
#>      TTGACTTTCCAGTAGTstochastic_pseudotime
#> Pyy                             -0.7921173
#> Rbp4                            -0.7802443
#> Sst                             -0.3329876
#>      GGGCACTGTGGCAAACstochastic_pseudotime
#> Pyy                             -0.7921193
#> Rbp4                            -0.7802706
#> Sst                             -0.3329901
#>      CGAGCACTCTACCTGCstochastic_pseudotime
#> Pyy                             -0.7921447
#> Rbp4                            -0.7806066
#> Sst                             -0.3330217
#>      TTGACTTTCAGCAACTstochastic_pseudotime
#> Pyy                             -0.7921524
#> Rbp4                            -0.7807092
#> Sst                             -0.3330312
#>      GGTATTGGTCCCGACAstochastic_pseudotime
#> Pyy                             -0.7921700
#> Rbp4                            -0.7809393
#> Sst                             -0.3330522
#>      ACACCGGCATGTCCTCstochastic_pseudotime
#> Pyy                             -0.7921707
#> Rbp4                            -0.7809485
#> Sst                             -0.3330530
#>      CTGAAACAGACTAGGCstochastic_pseudotime
#> Pyy                             -0.7921726
#> Rbp4                            -0.7809743
#> Sst                             -0.3330553
#>      GGACAAGGTGATAAGTstochastic_pseudotime
#> Pyy                             -0.7921860
#> Rbp4                            -0.7811489
#> Sst                             -0.3330709
#>      TACCTTAAGTACATGAstochastic_pseudotime
#> Pyy                             -0.7921951
#> Rbp4                            -0.7812686
#> Sst                             -0.3330814
#>      GGGAATGAGTGTACGGstochastic_pseudotime
#> Pyy                             -0.7922091
#> Rbp4                            -0.7814496
#> Sst                             -0.3330971
#>      GAAATGACAACTGGCCstochastic_pseudotime
#> Pyy                             -0.7922216
#> Rbp4                            -0.7816116
#> Sst                             -0.3331108
#>      TTTGCGCAGACAGACCstochastic_pseudotime
#> Pyy                             -0.7922276
#> Rbp4                            -0.7816896
#> Sst                             -0.3331174
#>      CACATTTGTTTCGCTCstochastic_pseudotime
#> Pyy                             -0.7922332
#> Rbp4                            -0.7817629
#> Sst                             -0.3331235
#>      TTTGGTTCACCAGATTstochastic_pseudotime
#> Pyy                             -0.7922367
#> Rbp4                            -0.7818069
#> Sst                             -0.3331271
#>      TAAGTGCAGGGTTTCTstochastic_pseudotime
#> Pyy                             -0.7922482
#> Rbp4                            -0.7819555
#> Sst                             -0.3331392
#>      TGCTACCCAGGGCATAstochastic_pseudotime
#> Pyy                             -0.7922620
#> Rbp4                            -0.7821334
#> Sst                             -0.3331534
#>      ATTGGTGTCCAGTATGstochastic_pseudotime
#> Pyy                             -0.7922653
#> Rbp4                            -0.7821751
#> Sst                             -0.3331567
#>      ATTTCTGGTAGTGAATstochastic_pseudotime
#> Pyy                             -0.7922724
#> Rbp4                            -0.7822672
#> Sst                             -0.3331639
#>      CTAGCCTTCAGCGATTstochastic_pseudotime
#> Pyy                             -0.7922725
#> Rbp4                            -0.7822674
#> Sst                             -0.3331639
#>      CACCTTGAGAGTAATCstochastic_pseudotime
#> Pyy                             -0.7922951
#> Rbp4                            -0.7825567
#> Sst                             -0.3331858
#>      GAAATGACACATGACTstochastic_pseudotime
#> Pyy                             -0.7923259
#> Rbp4                            -0.7829487
#> Sst                             -0.3332139
#>      CGGACTGAGGTAGCTGstochastic_pseudotime
#> Pyy                             -0.7923260
#> Rbp4                            -0.7829501
#> Sst                             -0.3332140
#>      CTCAGAAGTCCTCTTGstochastic_pseudotime
#> Pyy                             -0.7923280
#> Rbp4                            -0.7829757
#> Sst                             -0.3332158
#>      AACTGGTTCTTACCGCstochastic_pseudotime
#> Pyy                             -0.7923338
#> Rbp4                            -0.7830498
#> Sst                             -0.3332209
#>      CTCGTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -0.7923401
#> Rbp4                            -0.7831290
#> Sst                             -0.3332262
#>      CACAGTAGTTCCTCCAstochastic_pseudotime
#> Pyy                             -0.7923418
#> Rbp4                            -0.7831503
#> Sst                             -0.3332276
#>      GGGTTGCCAAGCTGAGstochastic_pseudotime
#> Pyy                             -0.7923471
#> Rbp4                            -0.7832180
#> Sst                             -0.3332321
#>      TACGGTACATGTTCCCstochastic_pseudotime
#> Pyy                             -0.7923620
#> Rbp4                            -0.7834066
#> Sst                             -0.3332442
#>      TGAAAGAGTTATCACGstochastic_pseudotime
#> Pyy                             -0.7923664
#> Rbp4                            -0.7834622
#> Sst                             -0.3332477
#>      TACAGTGCAATCTACGstochastic_pseudotime
#> Pyy                             -0.7923752
#> Rbp4                            -0.7835730
#> Sst                             -0.3332545
#>      GGACGTCCACGCATCGstochastic_pseudotime
#> Pyy                             -0.7923893
#> Rbp4                            -0.7837511
#> Sst                             -0.3332650
#>      CTAACTTCATTGGGCCstochastic_pseudotime
#> Pyy                             -0.7923902
#> Rbp4                            -0.7837622
#> Sst                             -0.3332656
#>      AAACCTGAGCCTTGATstochastic_pseudotime
#> Pyy                             -0.7923927
#> Rbp4                            -0.7837935
#> Sst                             -0.3332674
#>      CCTTCCCAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.7923949
#> Rbp4                            -0.7838213
#> Sst                             -0.3332690
#>      TCACGAATCCTATGTTstochastic_pseudotime
#> Pyy                             -0.7923965
#> Rbp4                            -0.7838423
#> Sst                             -0.3332701
#>      AAAGCAACACATGGGAstochastic_pseudotime
#> Pyy                             -0.7923984
#> Rbp4                            -0.7838656
#> Sst                             -0.3332715
#>      GAAGCAGTCAAGGCTTstochastic_pseudotime
#> Pyy                             -0.7924022
#> Rbp4                            -0.7839137
#> Sst                             -0.3332741
#>      AGACGTTTCTTACCTAstochastic_pseudotime
#> Pyy                             -0.7924048
#> Rbp4                            -0.7839458
#> Sst                             -0.3332759
#>      TAGACCAAGGAATCGCstochastic_pseudotime
#> Pyy                             -0.7924111
#> Rbp4                            -0.7840253
#> Sst                             -0.3332801
#>      TGAGCATTCACCCGAGstochastic_pseudotime
#> Pyy                             -0.7924112
#> Rbp4                            -0.7840270
#> Sst                             -0.3332802
#>      GGGACCTTCGAATGCTstochastic_pseudotime
#> Pyy                             -0.7924130
#> Rbp4                            -0.7840492
#> Sst                             -0.3332814
#>      GTACTCCCACGTCTCTstochastic_pseudotime
#> Pyy                             -0.7924226
#> Rbp4                            -0.7841704
#> Sst                             -0.3332876
#>      ACTGAACAGCCTCGTGstochastic_pseudotime
#> Pyy                             -0.7924265
#> Rbp4                            -0.7842197
#> Sst                             -0.3332900
#>      CTTAACTCAAGCGTAGstochastic_pseudotime
#> Pyy                             -0.7924270
#> Rbp4                            -0.7842260
#> Sst                             -0.3332903
#>      ACACCGGCATCACGTAstochastic_pseudotime
#> Pyy                             -0.7924346
#> Rbp4                            -0.7843212
#> Sst                             -0.3332949
#>      CAAGAAATCGTACCGGstochastic_pseudotime
#> Pyy                             -0.7924372
#> Rbp4                            -0.7843538
#> Sst                             -0.3332965
#>      CGAACATCAGACGTAGstochastic_pseudotime
#> Pyy                             -0.7924386
#> Rbp4                            -0.7843715
#> Sst                             -0.3332973
#>      ACGGCCACACGTGAGAstochastic_pseudotime
#> Pyy                             -0.7924506
#> Rbp4                            -0.7845229
#> Sst                             -0.3333041
#>      GAATGAACAATCGAAAstochastic_pseudotime
#> Pyy                             -0.7924543
#> Rbp4                            -0.7845694
#> Sst                             -0.3333061
#>      CTAAGACAGCCACGCTstochastic_pseudotime
#> Pyy                             -0.7924547
#> Rbp4                            -0.7845743
#> Sst                             -0.3333063
#>      TAGTTGGGTACCGTATstochastic_pseudotime
#> Pyy                             -0.7924586
#> Rbp4                            -0.7846234
#> Sst                             -0.3333083
#>      CCTTCGACAGCTGTTAstochastic_pseudotime
#> Pyy                             -0.7924606
#> Rbp4                            -0.7846486
#> Sst                             -0.3333093
#>      AGCAGCCTCTTTAGTCstochastic_pseudotime
#> Pyy                             -0.7924656
#> Rbp4                            -0.7847114
#> Sst                             -0.3333118
#>      CGTCACTCATAGTAAGstochastic_pseudotime
#> Pyy                             -0.7924663
#> Rbp4                            -0.7847195
#> Sst                             -0.3333122
#>      TACTTACCAACTGCGCstochastic_pseudotime
#> Pyy                             -0.7924710
#> Rbp4                            -0.7847794
#> Sst                             -0.3333145
#>      CATGACAGTAGCTTGTstochastic_pseudotime
#> Pyy                             -0.7924777
#> Rbp4                            -0.7848640
#> Sst                             -0.3333176
#>      CCTTCCCCAAAGTGCGstochastic_pseudotime
#> Pyy                             -0.7924808
#> Rbp4                            -0.7849030
#> Sst                             -0.3333189
#>      GGGTTGCTCGAACGGAstochastic_pseudotime
#> Pyy                             -0.7924824
#> Rbp4                            -0.7849221
#> Sst                             -0.3333196
#>      CTGGTCTCAGGAATGCstochastic_pseudotime
#> Pyy                             -0.7924827
#> Rbp4                            -0.7849265
#> Sst                             -0.3333197
#>      TGAGGGACAGACTCGCstochastic_pseudotime
#> Pyy                             -0.7924854
#> Rbp4                            -0.7849605
#> Sst                             -0.3333209
#>      GGTGTTAGTTCCTCCAstochastic_pseudotime
#> Pyy                             -0.7924855
#> Rbp4                            -0.7849611
#> Sst                             -0.3333209
#>      TGGTTCCCACCGCTAGstochastic_pseudotime
#> Pyy                             -0.7924920
#> Rbp4                            -0.7850436
#> Sst                             -0.3333236
#>      GACACGCCAAAGGAAGstochastic_pseudotime
#> Pyy                             -0.7924920
#> Rbp4                            -0.7850438
#> Sst                             -0.3333236
#>      CACAGTATCCTCCTAGstochastic_pseudotime
#> Pyy                             -0.7924941
#> Rbp4                            -0.7850693
#> Sst                             -0.3333244
#>      CGTAGGCAGCGTAGTGstochastic_pseudotime
#> Pyy                             -0.7924946
#> Rbp4                            -0.7850757
#> Sst                             -0.3333246
#>      AACGTTGGTGCCTTGGstochastic_pseudotime
#> Pyy                             -0.7924947
#> Rbp4                            -0.7850770
#> Sst                             -0.3333246
#>      CAAGGCCAGTGTGGCAstochastic_pseudotime
#> Pyy                             -0.7925067
#> Rbp4                            -0.7852289
#> Sst                             -0.3333289
#>      CTCTGGTCAGTTCATGstochastic_pseudotime
#> Pyy                             -0.7925132
#> Rbp4                            -0.7853113
#> Sst                             -0.3333310
#>      TATCTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -0.7925136
#> Rbp4                            -0.7853160
#> Sst                             -0.3333311
#>      GCTGCAGGTTACCGATstochastic_pseudotime
#> Pyy                             -0.7925159
#> Rbp4                            -0.7853447
#> Sst                             -0.3333318
#>      TACGGATAGTTTAGGAstochastic_pseudotime
#> Pyy                             -0.7925215
#> Rbp4                            -0.7854155
#> Sst                             -0.3333333
#>      ACTTTCAAGGATATACstochastic_pseudotime
#> Pyy                             -0.7925219
#> Rbp4                            -0.7854207
#> Sst                             -0.3333334
#>      ATAGACCTCTAACTGGstochastic_pseudotime
#> Pyy                             -0.7925238
#> Rbp4                            -0.7854454
#> Sst                             -0.3333340
#>      GACTACAGTTCGTGATstochastic_pseudotime
#> Pyy                             -0.7925270
#> Rbp4                            -0.7854849
#> Sst                             -0.3333347
#>      AAACGGGAGATATGGTstochastic_pseudotime
#> Pyy                             -0.7925292
#> Rbp4                            -0.7855132
#> Sst                             -0.3333353
#>      CTCGTCATCGGTCTAAstochastic_pseudotime
#> Pyy                             -0.7925294
#> Rbp4                            -0.7855161
#> Sst                             -0.3333353
#>      TCTCTAAGTCGCATCGstochastic_pseudotime
#> Pyy                             -0.7925304
#> Rbp4                            -0.7855289
#> Sst                             -0.3333355
#>      GAAATGAGTCGGATCCstochastic_pseudotime
#> Pyy                             -0.7925340
#> Rbp4                            -0.7855735
#> Sst                             -0.3333363
#>      GACGGCTAGGTGATATstochastic_pseudotime
#> Pyy                             -0.7925370
#> Rbp4                            -0.7856122
#> Sst                             -0.3333369
#>      CGTCACTCACGAAATAstochastic_pseudotime
#> Pyy                             -0.7925431
#> Rbp4                            -0.7856890
#> Sst                             -0.3333380
#>      TATCAGGCACCAGATTstochastic_pseudotime
#> Pyy                             -0.7925464
#> Rbp4                            -0.7857321
#> Sst                             -0.3333385
#>      TGACTTTAGCGTTTACstochastic_pseudotime
#> Pyy                             -0.7925492
#> Rbp4                            -0.7857674
#> Sst                             -0.3333389
#>      GGGCACTCATTTCAGGstochastic_pseudotime
#> Pyy                             -0.7925502
#> Rbp4                            -0.7857804
#> Sst                             -0.3333390
#>      GAGTCCGAGTTAAGTGstochastic_pseudotime
#> Pyy                             -0.7925510
#> Rbp4                            -0.7857902
#> Sst                             -0.3333391
#>      GTTCTCGCAGACACTTstochastic_pseudotime
#> Pyy                             -0.7925527
#> Rbp4                            -0.7858122
#> Sst                             -0.3333393
#>      GTATTCTCACATTAGCstochastic_pseudotime
#> Pyy                             -0.7925554
#> Rbp4                            -0.7858460
#> Sst                             -0.3333396
#>      AACTCCCTCGAGAGCAstochastic_pseudotime
#> Pyy                             -0.7925579
#> Rbp4                            -0.7858783
#> Sst                             -0.3333398
#>      CGTCAGGAGCCTATGTstochastic_pseudotime
#> Pyy                             -0.7925607
#> Rbp4                            -0.7859141
#> Sst                             -0.3333400
#>      GTTCATTAGAAGATTCstochastic_pseudotime
#> Pyy                             -0.7925686
#> Rbp4                            -0.7860149
#> Sst                             -0.3333404
#>      TCAGCAACACCAGGCTstochastic_pseudotime
#> Pyy                             -0.7925708
#> Rbp4                            -0.7860425
#> Sst                             -0.3333404
#>      CGTCAGGCATGGGAACstochastic_pseudotime
#> Pyy                             -0.7925767
#> Rbp4                            -0.7861177
#> Sst                             -0.3333404
#>      CCTACACAGATCCGAGstochastic_pseudotime
#> Pyy                             -0.7925768
#> Rbp4                            -0.7861195
#> Sst                             -0.3333404
#>      CAGAATCGTGTAAGTAstochastic_pseudotime
#> Pyy                             -0.7925770
#> Rbp4                            -0.7861225
#> Sst                             -0.3333404
#>      TGCGTGGAGGTAAACTstochastic_pseudotime
#> Pyy                             -0.7925840
#> Rbp4                            -0.7862116
#> Sst                             -0.3333401
#>      TGATTTCAGTGACTCTstochastic_pseudotime
#> Pyy                             -0.7925855
#> Rbp4                            -0.7862313
#> Sst                             -0.3333400
#>      TCAATCTCATCTCCCAstochastic_pseudotime
#> Pyy                             -0.7925864
#> Rbp4                            -0.7862428
#> Sst                             -0.3333399
#>      CCTTACGAGATCACGGstochastic_pseudotime
#> Pyy                             -0.7925909
#> Rbp4                            -0.7863007
#> Sst                             -0.3333395
#>      TGAGCCGTCATGTCCCstochastic_pseudotime
#> Pyy                             -0.7925947
#> Rbp4                            -0.7863487
#> Sst                             -0.3333391
#>      CCACCTATCACCAGGCstochastic_pseudotime
#> Pyy                             -0.7926028
#> Rbp4                            -0.7864538
#> Sst                             -0.3333378
#>      GCAATCATCCAACCAAstochastic_pseudotime
#> Pyy                             -0.7926028
#> Rbp4                            -0.7864538
#> Sst                             -0.3333378
#>      AGAGTGGGTAGCACGAstochastic_pseudotime
#> Pyy                             -0.7926031
#> Rbp4                            -0.7864578
#> Sst                             -0.3333377
#>      AGTAGTCGTTCTGTTTstochastic_pseudotime
#> Pyy                             -0.7926032
#> Rbp4                            -0.7864592
#> Sst                             -0.3333377
#>      GCGGGTTGTTACTGACstochastic_pseudotime
#> Pyy                             -0.7926038
#> Rbp4                            -0.7864671
#> Sst                             -0.3333376
#>      GCATGATTCAAACGGGstochastic_pseudotime
#> Pyy                             -0.7926058
#> Rbp4                            -0.7864926
#> Sst                             -0.3333372
#>      TCATTACGTCAGAAGCstochastic_pseudotime
#> Pyy                             -0.7926150
#> Rbp4                            -0.7866111
#> Sst                             -0.3333351
#>      GATCTAGGTGCCTTGGstochastic_pseudotime
#> Pyy                             -0.7926173
#> Rbp4                            -0.7866414
#> Sst                             -0.3333345
#>      AGTAGTCCACATGACTstochastic_pseudotime
#> Pyy                             -0.7926190
#> Rbp4                            -0.7866641
#> Sst                             -0.3333340
#>      CGTTAGAGTGCAACGAstochastic_pseudotime
#> Pyy                             -0.7926298
#> Rbp4                            -0.7868038
#> Sst                             -0.3333306
#>      CTCGTACCAATGGATAstochastic_pseudotime
#> Pyy                             -0.7926306
#> Rbp4                            -0.7868147
#> Sst                             -0.3333303
#>      TTGCCGTTCTCTGTCGstochastic_pseudotime
#> Pyy                             -0.7926406
#> Rbp4                            -0.7869452
#> Sst                             -0.3333262
#>      ACATACGAGTAGGCCAstochastic_pseudotime
#> Pyy                             -0.7926464
#> Rbp4                            -0.7870210
#> Sst                             -0.3333236
#>      TTTACTGGTGCTAGCCstochastic_pseudotime
#> Pyy                             -0.7926519
#> Rbp4                            -0.7870934
#> Sst                             -0.3333208
#>      AGCTTGACAATAGAGTstochastic_pseudotime
#> Pyy                             -0.7926547
#> Rbp4                            -0.7871301
#> Sst                             -0.3333193
#>      CTCAGAATCAGGTTCAstochastic_pseudotime
#> Pyy                             -0.7926642
#> Rbp4                            -0.7872559
#> Sst                             -0.3333137
#>      CGTCACTAGGGTGTGTstochastic_pseudotime
#> Pyy                             -0.7926713
#> Rbp4                            -0.7873499
#> Sst                             -0.3333090
#>      CTTGGCTAGGACACCAstochastic_pseudotime
#> Pyy                             -0.7926719
#> Rbp4                            -0.7873584
#> Sst                             -0.3333086
#>      TGTCCCAAGTTTAGGAstochastic_pseudotime
#> Pyy                             -0.7926739
#> Rbp4                            -0.7873844
#> Sst                             -0.3333072
#>      TTCTACAAGATCGATAstochastic_pseudotime
#> Pyy                             -0.7926821
#> Rbp4                            -0.7874948
#> Sst                             -0.3333010
#>      AGCTCCTGTGTGCCTGstochastic_pseudotime
#> Pyy                             -0.7926906
#> Rbp4                            -0.7876076
#> Sst                             -0.3332941
#>      GCCTCTATCGCTTAGAstochastic_pseudotime
#> Pyy                             -0.7926940
#> Rbp4                            -0.7876536
#> Sst                             -0.3332910
#>      GACGTTAAGACGCAACstochastic_pseudotime
#> Pyy                             -0.7926963
#> Rbp4                            -0.7876845
#> Sst                             -0.3332889
#>      GGGTCTGCATGAAGTAstochastic_pseudotime
#> Pyy                             -0.7926980
#> Rbp4                            -0.7877077
#> Sst                             -0.3332873
#>      CTTTGCGAGGAATGGAstochastic_pseudotime
#> Pyy                             -0.7927038
#> Rbp4                            -0.7877859
#> Sst                             -0.3332817
#>      TACGGGCAGCGTCAAGstochastic_pseudotime
#> Pyy                             -0.7927045
#> Rbp4                            -0.7877959
#> Sst                             -0.3332810
#>      CCGTTCACACATGACTstochastic_pseudotime
#> Pyy                             -0.7927053
#> Rbp4                            -0.7878071
#> Sst                             -0.3332802
#>      GACGCGTCAGTCGTGCstochastic_pseudotime
#> Pyy                             -0.7927057
#> Rbp4                            -0.7878122
#> Sst                             -0.3332798
#>      ACTTACTAGCCGGTAAstochastic_pseudotime
#> Pyy                             -0.7927097
#> Rbp4                            -0.7878660
#> Sst                             -0.3332757
#>      GATGAAACATGGTTGTstochastic_pseudotime
#> Pyy                             -0.7927158
#> Rbp4                            -0.7879489
#> Sst                             -0.3332690
#>      AGGCCACAGCTGTCTAstochastic_pseudotime
#> Pyy                             -0.7927173
#> Rbp4                            -0.7879693
#> Sst                             -0.3332673
#>      CACATTTTCGCAGGCTstochastic_pseudotime
#> Pyy                             -0.7927181
#> Rbp4                            -0.7879808
#> Sst                             -0.3332664
#>      CGAGCACTCTTCGAGAstochastic_pseudotime
#> Pyy                             -0.7927189
#> Rbp4                            -0.7879920
#> Sst                             -0.3332654
#>      CAAGGCCAGCTGTTCAstochastic_pseudotime
#> Pyy                             -0.7927219
#> Rbp4                            -0.7880324
#> Sst                             -0.3332619
#>      TTGAACGGTGATGCCCstochastic_pseudotime
#> Pyy                             -0.7927232
#> Rbp4                            -0.7880506
#> Sst                             -0.3332603
#>      ACTGAGTAGTGACTCTstochastic_pseudotime
#> Pyy                             -0.7927235
#> Rbp4                            -0.7880544
#> Sst                             -0.3332600
#>      GTATCTTCAGTGGGATstochastic_pseudotime
#> Pyy                             -0.7927266
#> Rbp4                            -0.7880975
#> Sst                             -0.3332561
#>      TCTGGAATCGGTGTCGstochastic_pseudotime
#> Pyy                             -0.7927274
#> Rbp4                            -0.7881078
#> Sst                             -0.3332552
#>      CATGACACAACAACCTstochastic_pseudotime
#> Pyy                             -0.7927275
#> Rbp4                            -0.7881099
#> Sst                             -0.3332550
#>      GAGGTGAGTGACAAATstochastic_pseudotime
#> Pyy                             -0.7927282
#> Rbp4                            -0.7881192
#> Sst                             -0.3332541
#>      AACTTTCGTGACCAAGstochastic_pseudotime
#> Pyy                             -0.7927299
#> Rbp4                            -0.7881425
#> Sst                             -0.3332520
#>      TACACGACAACGATCTstochastic_pseudotime
#> Pyy                             -0.7927360
#> Rbp4                            -0.7882268
#> Sst                             -0.3332439
#>      CCACTACAGAGTACCGstochastic_pseudotime
#> Pyy                             -0.7927372
#> Rbp4                            -0.7882444
#> Sst                             -0.3332421
#>      TCTTTCCTCGTCTGCTstochastic_pseudotime
#> Pyy                             -0.7927407
#> Rbp4                            -0.7882926
#> Sst                             -0.3332372
#>      GGCCGATTCACTCTTAstochastic_pseudotime
#> Pyy                             -0.7927408
#> Rbp4                            -0.7882943
#> Sst                             -0.3332371
#>      GTAACGTGTCACTGGCstochastic_pseudotime
#> Pyy                             -0.7927427
#> Rbp4                            -0.7883200
#> Sst                             -0.3332344
#>      CGGTTAATCGTCTGAAstochastic_pseudotime
#> Pyy                             -0.7927448
#> Rbp4                            -0.7883501
#> Sst                             -0.3332312
#>      TCAGGTAAGTCAATAGstochastic_pseudotime
#> Pyy                             -0.7927473
#> Rbp4                            -0.7883845
#> Sst                             -0.3332275
#>      TGGCCAGCACGAAAGCstochastic_pseudotime
#> Pyy                             -0.7927502
#> Rbp4                            -0.7884259
#> Sst                             -0.3332230
#>      GATGAGGTCTGGGCCAstochastic_pseudotime
#> Pyy                             -0.7927540
#> Rbp4                            -0.7884787
#> Sst                             -0.3332170
#>      GCAGTTATCTTGAGACstochastic_pseudotime
#> Pyy                             -0.7927543
#> Rbp4                            -0.7884831
#> Sst                             -0.3332165
#>      TTCTACAGTACTCTCCstochastic_pseudotime
#> Pyy                             -0.7927574
#> Rbp4                            -0.7885272
#> Sst                             -0.3332114
#>      CGCTGGAGTCCATCCTstochastic_pseudotime
#> Pyy                             -0.7927586
#> Rbp4                            -0.7885436
#> Sst                             -0.3332095
#>      CGTGTAAGTCTTGATGstochastic_pseudotime
#> Pyy                             -0.7927622
#> Rbp4                            -0.7885955
#> Sst                             -0.3332032
#>      TGCGGGTTCTCAAACGstochastic_pseudotime
#> Pyy                             -0.7927690
#> Rbp4                            -0.7886919
#> Sst                             -0.3331912
#>      TACGGATGTGTTCGATstochastic_pseudotime
#> Pyy                             -0.7927717
#> Rbp4                            -0.7887303
#> Sst                             -0.3331862
#>      ATCCGAAGTGACTCATstochastic_pseudotime
#> Pyy                             -0.7927832
#> Rbp4                            -0.7888960
#> Sst                             -0.3331634
#>      CTCTGGTCAACACCTAstochastic_pseudotime
#> Pyy                             -0.7927877
#> Rbp4                            -0.7889625
#> Sst                             -0.3331538
#>      AAACGGGTCAGCTCTCstochastic_pseudotime
#> Pyy                             -0.7927882
#> Rbp4                            -0.7889687
#> Sst                             -0.3331528
#>      TGAGCCGTCCGAACGCstochastic_pseudotime
#> Pyy                             -0.7927982
#> Rbp4                            -0.7891167
#> Sst                             -0.3331300
#>      AGAGCTTAGCCGATTTstochastic_pseudotime
#> Pyy                             -0.7928078
#> Rbp4                            -0.7892587
#> Sst                             -0.3331065
#>      CGAATGTAGTTCCACAstochastic_pseudotime
#> Pyy                             -0.7928084
#> Rbp4                            -0.7892678
#> Sst                             -0.3331050
#>      CTTAGGAGTCCGAATTstochastic_pseudotime
#> Pyy                             -0.7928102
#> Rbp4                            -0.7892957
#> Sst                             -0.3331002
#>      ATTGGACAGTAATCCCstochastic_pseudotime
#> Pyy                             -0.7928204
#> Rbp4                            -0.7894494
#> Sst                             -0.3330724
#>      GGTGTTAAGGATCGCAstochastic_pseudotime
#> Pyy                             -0.7928209
#> Rbp4                            -0.7894572
#> Sst                             -0.3330709
#>      TGGACGCAGTGTTGAAstochastic_pseudotime
#> Pyy                             -0.7928223
#> Rbp4                            -0.7894786
#> Sst                             -0.3330669
#>      TCGTACCCATATACGCstochastic_pseudotime
#> Pyy                             -0.7928233
#> Rbp4                            -0.7894940
#> Sst                             -0.3330640
#>      CGACTTCCAGCCTATAstochastic_pseudotime
#> Pyy                             -0.7928240
#> Rbp4                            -0.7895051
#> Sst                             -0.3330619
#>      AAGGCAGTCCGAGCCAstochastic_pseudotime
#> Pyy                             -0.7928285
#> Rbp4                            -0.7895733
#> Sst                             -0.3330486
#>      ACGAGGAAGTTACCCAstochastic_pseudotime
#> Pyy                             -0.7928351
#> Rbp4                            -0.7896771
#> Sst                             -0.3330275
#>      TCAGCTCAGTCGTACTstochastic_pseudotime
#> Pyy                             -0.7928422
#> Rbp4                            -0.7897890
#> Sst                             -0.3330038
#>      CGTCAGGTCCTGCCATstochastic_pseudotime
#> Pyy                             -0.7928857
#> Rbp4                            -0.7905123
#> Sst                             -0.3328195
#>      TGATTTCAGTACGCCCstochastic_pseudotime
#> Pyy                             -0.7928892
#> Rbp4                            -0.7905735
#> Sst                             -0.3328012
#>      AGAGTGGGTCACCTAAstochastic_pseudotime
#> Pyy                             -0.7928914
#> Rbp4                            -0.7906124
#> Sst                             -0.3327894
#>      CGATCGGCACGCCAGTstochastic_pseudotime
#> Pyy                             -0.7929033
#> Rbp4                            -0.7908303
#> Sst                             -0.3327198
#>      TGCCCATCACTGTGTAstochastic_pseudotime
#> Pyy                             -0.7929122
#> Rbp4                            -0.7909978
#> Sst                             -0.3326622
#>      AACTTTCCACTGCCAGstochastic_pseudotime
#> Pyy                             -0.7929238
#> Rbp4                            -0.7912280
#> Sst                             -0.3325770
#>      GTCCTCATCATTTGGGstochastic_pseudotime
#> Pyy                             -0.7929405
#> Rbp4                            -0.7915827
#> Sst                             -0.3324314
#>      GTGCGGTAGTAGGCCAstochastic_pseudotime
#> Pyy                             -0.7929441
#> Rbp4                            -0.7916633
#> Sst                             -0.3323957
#>      CGACCTTAGCGCTTATstochastic_pseudotime
#> Pyy                             -0.7929675
#> Rbp4                            -0.7922475
#> Sst                             -0.3321085
#>      GTCTCGTAGGCTAGACstochastic_pseudotime
#> Pyy                             -0.7929705
#> Rbp4                            -0.7923332
#> Sst                             -0.3320621
#>      GCATGATAGCCTTGATstochastic_pseudotime
#> Pyy                             -0.7929807
#> Rbp4                            -0.7926454
#> Sst                             -0.3318840
#>      ACTTTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -0.7929850
#> Rbp4                            -0.7927905
#> Sst                             -0.3317968
#>      CCAATCCAGAGTAAGGstochastic_pseudotime
#> Pyy                             -0.7929855
#> Rbp4                            -0.7928078
#> Sst                             -0.3317862
#>      ATTATCCAGAAACCTAstochastic_pseudotime
#> Pyy                             -0.7929917
#> Rbp4                            -0.7930414
#> Sst                             -0.3316400
#>      TCGGGACAGAGAACAGstochastic_pseudotime
#> Pyy                             -0.7929979
#> Rbp4                            -0.7933097
#> Sst                             -0.3314654
#>      AGCGGTCGTACCGAGAstochastic_pseudotime
#> Pyy                             -0.7930010
#> Rbp4                            -0.7934584
#> Sst                             -0.3313666
#>      TTCCCAGAGTTACCCAstochastic_pseudotime
#> Pyy                             -0.7930018
#> Rbp4                            -0.7934971
#> Sst                             -0.3313407
#>      CTCGAGGTCTGTTGAGstochastic_pseudotime
#> Pyy                             -0.7930067
#> Rbp4                            -0.7937692
#> Sst                             -0.3311582
#>      ATTATCCAGACCCACCstochastic_pseudotime
#> Pyy                             -0.7930075
#> Rbp4                            -0.7938169
#> Sst                             -0.3311263
#>      CCTCAGTGTTACCAGTstochastic_pseudotime
#> Pyy                             -0.7930117
#> Rbp4                            -0.7940855
#> Sst                             -0.3309486
#>      GCCTCTAAGTCCGGTCstochastic_pseudotime
#> Pyy                             -0.7930195
#> Rbp4                            -0.7947977
#> Sst                             -0.3305011
#>      TTCGAAGAGTACGCCCstochastic_pseudotime
#> Pyy                             -0.7930216
#> Rbp4                            -0.7951540
#> Sst                             -0.3302960
#>      GCTGCGAAGCTAGCCCstochastic_pseudotime
#> Pyy                             -0.7930216
#> Rbp4                            -0.7951638
#> Sst                             -0.3302906
#>      GCGGGTTAGACTTGAAstochastic_pseudotime
#> Pyy                             -0.7930212
#> Rbp4                            -0.7958583
#> Sst                             -0.3299457
#>      GTATCTTCATCGATGTstochastic_pseudotime
#> Pyy                             -0.7930191
#> Rbp4                            -0.7961440
#> Sst                             -0.3298304
#>      ACTTACTAGAGGGCTTstochastic_pseudotime
#> Pyy                             -0.7930141
#> Rbp4                            -0.7965348
#> Sst                             -0.3297045
#>      GGGTCTGCACAGGCCTstochastic_pseudotime
#> Pyy                             -0.7930134
#> Rbp4                            -0.7965804
#> Sst                             -0.3296924
#>      GGTATTGGTACTCGCGstochastic_pseudotime
#> Pyy                             -0.7930082
#> Rbp4                            -0.7968492
#> Sst                             -0.3296340
#>      TGTGGTATCCAGTAGTstochastic_pseudotime
#> Pyy                             -0.7930075
#> Rbp4                            -0.7968796
#> Sst                             -0.3296288
#>      CTCAGAACACTAGTACstochastic_pseudotime
#> Pyy                             -0.7929985
#> Rbp4                            -0.7972191
#> Sst                             -0.3295923
#>      CAGAATCGTCATTAGCstochastic_pseudotime
#> Pyy                             -0.7929971
#> Rbp4                            -0.7972661
#> Sst                             -0.3295905
#>      ACCAGTATCGTTTAGGstochastic_pseudotime
#> Pyy                             -0.7929964
#> Rbp4                            -0.7972879
#> Sst                             -0.3295899
#>      TTTACTGCAGCATACTstochastic_pseudotime
#> Pyy                             -0.7929923
#> Rbp4                            -0.7974099
#> Sst                             -0.3295903
#>      GCAGCCACAGCGTTCGstochastic_pseudotime
#> Pyy                             -0.7929864
#> Rbp4                            -0.7975678
#> Sst                             -0.3295997
#>      CTGTGCTCATTACCTTstochastic_pseudotime
#> Pyy                             -0.7929858
#> Rbp4                            -0.7975828
#> Sst                             -0.3296011
#>      CTCGAAAGTAAGAGGAstochastic_pseudotime
#> Pyy                             -0.7929726
#> Rbp4                            -0.7978807
#> Sst                             -0.3296503
#>      TATGCCCGTCTAGCCGstochastic_pseudotime
#> Pyy                             -0.7929547
#> Rbp4                            -0.7982068
#> Sst                             -0.3297528
#>      GGGCATCGTGGTCTCGstochastic_pseudotime
#> Pyy                             -0.7929055
#> Rbp4                            -0.7988538
#> Sst                             -0.3301311
#>      CATGACACATCGATTGstochastic_pseudotime
#> Pyy                             -0.7929044
#> Rbp4                            -0.7988659
#> Sst                             -0.3301407
#>      AACTTTCCACAGTCGCstochastic_pseudotime
#> Pyy                             -0.7928807
#> Rbp4                            -0.7990970
#> Sst                             -0.3303426
#>      GAAATGAGTCACCCAGstochastic_pseudotime
#> Pyy                             -0.7928797
#> Rbp4                            -0.7991061
#> Sst                             -0.3303513
#>      ACTTACTCAATCAGAAstochastic_pseudotime
#> Pyy                             -0.7928657
#> Rbp4                            -0.7992262
#> Sst                             -0.3304724
#>      TGCTGCTAGGATGCGTstochastic_pseudotime
#> Pyy                             -0.7928619
#> Rbp4                            -0.7992577
#> Sst                             -0.3305059
#>      AACTCAGGTGCCTGTGstochastic_pseudotime
#> Pyy                             -0.7928505
#> Rbp4                            -0.7993459
#> Sst                             -0.3306042
#>      GTATCTTGTGGCCCTAstochastic_pseudotime
#> Pyy                             -0.7928407
#> Rbp4                            -0.7994180
#> Sst                             -0.3306891
#>      GTCCTCATCTGGAGCCstochastic_pseudotime
#> Pyy                             -0.7928092
#> Rbp4                            -0.7996245
#> Sst                             -0.3309575
#>      TGACTAGTCCTATTCAstochastic_pseudotime
#> Pyy                             -0.7927771
#> Rbp4                            -0.7998031
#> Sst                             -0.3312226
#>      TGAGCATGTTAAAGACstochastic_pseudotime
#> Pyy                             -0.7927385
#> Rbp4                            -0.7999842
#> Sst                             -0.3315273
#>      TGCCAAACATGACATCstochastic_pseudotime
#> Pyy                             -0.7927117
#> Rbp4                            -0.8000924
#> Sst                             -0.3317295
#>      AAGTCTGTCCGAGCCAstochastic_pseudotime
#> Pyy                             -0.7926934
#> Rbp4                            -0.8001585
#> Sst                             -0.3318618
#>      ACATCAGTCCTTGGTCstochastic_pseudotime
#> Pyy                             -0.7926728
#> Rbp4                            -0.8002269
#> Sst                             -0.3320070
#>      GTGCGGTCAGCGTAAGstochastic_pseudotime
#> Pyy                             -0.7926554
#> Rbp4                            -0.8002798
#> Sst                             -0.3321255
#>      GTCTCGTCATAGTAAGstochastic_pseudotime
#> Pyy                             -0.7926389
#> Rbp4                            -0.8003260
#> Sst                             -0.3322344
#>      GATGAGGAGTGGACGTstochastic_pseudotime
#> Pyy                             -0.7926368
#> Rbp4                            -0.8003317
#> Sst                             -0.3322480
#>      AGAGCTTAGGTTACCTstochastic_pseudotime
#> Pyy                             -0.7926274
#> Rbp4                            -0.8003561
#> Sst                             -0.3323082
#>      GCTGGGTGTCTAGAGGstochastic_pseudotime
#> Pyy                             -0.7926245
#> Rbp4                            -0.8003636
#> Sst                             -0.3323269
#>      GGGATGAGTTGTGGAGstochastic_pseudotime
#> Pyy                             -0.7926233
#> Rbp4                            -0.8003666
#> Sst                             -0.3323346
#>      AGCGGTCGTCGAATCTstochastic_pseudotime
#> Pyy                             -0.7926108
#> Rbp4                            -0.8003967
#> Sst                             -0.3324122
#>      AGGCCGTCACGACGAAstochastic_pseudotime
#> Pyy                             -0.7925575
#> Rbp4                            -0.8005066
#> Sst                             -0.3327262
#>      CGAGAAGCACGCATCGstochastic_pseudotime
#> Pyy                             -0.7925571
#> Rbp4                            -0.8005072
#> Sst                             -0.3327283
#>      GTACTCCTCCAACCAAstochastic_pseudotime
#> Pyy                             -0.7925515
#> Rbp4                            -0.8005170
#> Sst                             -0.3327594
#>      ATGTGTGCACAAGTAAstochastic_pseudotime
#> Pyy                             -0.7925497
#> Rbp4                            -0.8005201
#> Sst                             -0.3327692
#>      TGACTAGGTGTATGGGstochastic_pseudotime
#> Pyy                             -0.7925039
#> Rbp4                            -0.8005885
#> Sst                             -0.3330109
#>      AGGGATGGTTAGGGTGstochastic_pseudotime
#> Pyy                             -0.7925004
#> Rbp4                            -0.8005930
#> Sst                             -0.3330287
#>      CGTGTAATCATTATCCstochastic_pseudotime
#> Pyy                             -0.7924579
#> Rbp4                            -0.8006391
#> Sst                             -0.3332328
#>      CTGAAACTCGATAGAAstochastic_pseudotime
#> Pyy                             -0.7924357
#> Rbp4                            -0.8006575
#> Sst                             -0.3333331
#>      TTGAACGAGGAGTTTAstochastic_pseudotime
#> Pyy                             -0.7924324
#> Rbp4                            -0.8006599
#> Sst                             -0.3333474
#>      AGCGTCGCACAGGCCTstochastic_pseudotime
#> Pyy                             -0.7924286
#> Rbp4                            -0.8006626
#> Sst                             -0.3333640
#>      CCACGGAGTTACGACTstochastic_pseudotime
#> Pyy                             -0.7924113
#> Rbp4                            -0.8006735
#> Sst                             -0.3334379
#>      CTGAAACTCTGTTTGTstochastic_pseudotime
#> Pyy                             -0.7924107
#> Rbp4                            -0.8006738
#> Sst                             -0.3334403
#>      AGATTGCAGTTCCACAstochastic_pseudotime
#> Pyy                             -0.7924015
#> Rbp4                            -0.8006787
#> Sst                             -0.3334787
#>      GTATCTTGTCGCGGTTstochastic_pseudotime
#> Pyy                             -0.7923919
#> Rbp4                            -0.8006832
#> Sst                             -0.3335179
#>      GAGTCCGAGCGATGACstochastic_pseudotime
#> Pyy                             -0.7923916
#> Rbp4                            -0.8006833
#> Sst                             -0.3335191
#>      AACCATGCACAGGAGTstochastic_pseudotime
#> Pyy                             -0.7923812
#> Rbp4                            -0.8006874
#> Sst                             -0.3335606
#>      ATGTGTGTCAAAGACAstochastic_pseudotime
#> Pyy                             -0.7923771
#> Rbp4                            -0.8006888
#> Sst                             -0.3335767
#>      CCGTACTAGGAGTACCstochastic_pseudotime
#> Pyy                             -0.7923653
#> Rbp4                            -0.8006922
#> Sst                             -0.3336225
#>      CCGTACTAGTCCGTATstochastic_pseudotime
#> Pyy                             -0.7923627
#> Rbp4                            -0.8006929
#> Sst                             -0.3336323
#>      GACGTTAAGACGCTTTstochastic_pseudotime
#> Pyy                             -0.7923260
#> Rbp4                            -0.8006971
#> Sst                             -0.3337672
#>      CCACGGAGTCAACATCstochastic_pseudotime
#> Pyy                             -0.7922870
#> Rbp4                            -0.8006922
#> Sst                             -0.3338996
#>      CATGGCGAGAGTCTGGstochastic_pseudotime
#> Pyy                             -0.7922816
#> Rbp4                            -0.8006907
#> Sst                             -0.3339173
#>      GAAATGAAGCTCCTCTstochastic_pseudotime
#> Pyy                             -0.7922627
#> Rbp4                            -0.8006844
#> Sst                             -0.3339771
#>      ACGCCAGTCCCAAGTAstochastic_pseudotime
#> Pyy                             -0.7922524
#> Rbp4                            -0.8006801
#> Sst                             -0.3340086
#>      CGGTTAAGTCTAGCCGstochastic_pseudotime
#> Pyy                             -0.7922414
#> Rbp4                            -0.8006747
#> Sst                             -0.3340417
#>      GTCTCGTTCACTTCATstochastic_pseudotime
#> Pyy                             -0.7922191
#> Rbp4                            -0.8006617
#> Sst                             -0.3341065
#>      GACCAATTCAACGGCCstochastic_pseudotime
#> Pyy                             -0.7921873
#> Rbp4                            -0.8006382
#> Sst                             -0.3341936
#>      CTCTGGTAGCAATCTCstochastic_pseudotime
#> Pyy                             -0.7920903
#> Rbp4                            -0.8005311
#> Sst                             -0.3344260
#>      TCTCTAAGTGGCAAACstochastic_pseudotime
#> Pyy                             -0.7920781
#> Rbp4                            -0.8005138
#> Sst                             -0.3344519
#>      CATCGGGTCTTCGGTCstochastic_pseudotime
#> Pyy                             -0.7920603
#> Rbp4                            -0.8004872
#> Sst                             -0.3344885
#>      CAACCAACACCTCGTTstochastic_pseudotime
#> Pyy                             -0.7920104
#> Rbp4                            -0.8004031
#> Sst                             -0.3345834
#>      CTGGTCTTCGCATGGCstochastic_pseudotime
#> Pyy                             -0.7919847
#> Rbp4                            -0.8003543
#> Sst                             -0.3346282
#>      ACTGTCCTCGCTAGCGstochastic_pseudotime
#> Pyy                             -0.7919768
#> Rbp4                            -0.8003385
#> Sst                             -0.3346414
#>      CTAAGACAGCACGCCTstochastic_pseudotime
#> Pyy                             -0.7919435
#> Rbp4                            -0.8002681
#> Sst                             -0.3346945
#>      CGGCTAGGTACTCGCGstochastic_pseudotime
#> Pyy                             -0.7919141
#> Rbp4                            -0.8002004
#> Sst                             -0.3347379
#>      TCTCATAGTGGACGATstochastic_pseudotime
#> Pyy                             -0.7918853
#> Rbp4                            -0.8001293
#> Sst                             -0.3347773
#>      GGATGTTTCCAGTATGstochastic_pseudotime
#> Pyy                             -0.7918574
#> Rbp4                            -0.8000555
#> Sst                             -0.3348128
#>      CTTAGGATCACTCTTAstochastic_pseudotime
#> Pyy                             -0.7918442
#> Rbp4                            -0.8000191
#> Sst                             -0.3348285
#>      TACGGATTCCAAGTACstochastic_pseudotime
#> Pyy                             -0.7918203
#> Rbp4                            -0.7999499
#> Sst                             -0.3348558
#>      TGAGAGGAGCCGTCGTstochastic_pseudotime
#> Pyy                             -0.7918080
#> Rbp4                            -0.7999131
#> Sst                             -0.3348690
#>      CTTGGCTCAGCCACCAstochastic_pseudotime
#> Pyy                             -0.7918000
#> Rbp4                            -0.7998886
#> Sst                             -0.3348774
#>      AATCGGTCATATGGTCstochastic_pseudotime
#> Pyy                             -0.7917809
#> Rbp4                            -0.7998284
#> Sst                             -0.3348965
#>      TACTTACCATGGGACAstochastic_pseudotime
#> Pyy                             -0.7917708
#> Rbp4                            -0.7997956
#> Sst                             -0.3349062
#>      TGAGGGAAGAGGTACCstochastic_pseudotime
#> Pyy                             -0.7917417
#> Rbp4                            -0.7996972
#> Sst                             -0.3349323
#>      TCCCGATCAACGATGGstochastic_pseudotime
#> Pyy                             -0.7916820
#> Rbp4                            -0.7994766
#> Sst                             -0.3349780
#>      GTCACAACACGTCTCTstochastic_pseudotime
#> Pyy                             -0.7916139
#> Rbp4                            -0.7991927
#> Sst                             -0.3350175
#>      TTTGGTTCACGAAGCAstochastic_pseudotime
#> Pyy                             -0.7916133
#> Rbp4                            -0.7991902
#> Sst                             -0.3350178
#>      GGGAATGTCACAAACCstochastic_pseudotime
#> Pyy                             -0.7915943
#> Rbp4                            -0.7991046
#> Sst                             -0.3350264
#>      ACAGCCGCATCGTCGGstochastic_pseudotime
#> Pyy                             -0.7915841
#> Rbp4                            -0.7990574
#> Sst                             -0.3350306
#>      CGTGAGCTCAGCGATTstochastic_pseudotime
#> Pyy                             -0.7915453
#> Rbp4                            -0.7988696
#> Sst                             -0.3350441
#>      ATCATGGGTCCGCTGAstochastic_pseudotime
#> Pyy                             -0.7915382
#> Rbp4                            -0.7988338
#> Sst                             -0.3350461
#>      GGCAATTCATCTCCCAstochastic_pseudotime
#> Pyy                             -0.7915086
#> Rbp4                            -0.7986796
#> Sst                             -0.3350530
#>      CCTATTATCATCGCTCstochastic_pseudotime
#> Pyy                             -0.7914951
#> Rbp4                            -0.7986065
#> Sst                             -0.3350554
#>      CAAGGCCCACATTTCTstochastic_pseudotime
#> Pyy                             -0.7914523
#> Rbp4                            -0.7983635
#> Sst                             -0.3350594
#>      CTCACACGTCACACGCstochastic_pseudotime
#> Pyy                             -0.7914491
#> Rbp4                            -0.7983448
#> Sst                             -0.3350595
#>      AACCATGGTTCGGCACstochastic_pseudotime
#> Pyy                             -0.7914448
#> Rbp4                            -0.7983192
#> Sst                             -0.3350596
#>      CTGCGGACAGACTCGCstochastic_pseudotime
#> Pyy                             -0.7914420
#> Rbp4                            -0.7983026
#> Sst                             -0.3350597
#>      ACGCAGCGTAATCACCstochastic_pseudotime
#> Pyy                             -0.7914356
#> Rbp4                            -0.7982636
#> Sst                             -0.3350596
#>      CCATGTCGTCGGCTCAstochastic_pseudotime
#> Pyy                             -0.7914324
#> Rbp4                            -0.7982442
#> Sst                             -0.3350596
#>      TTGCGTCCACCGATATstochastic_pseudotime
#> Pyy                             -0.7914319
#> Rbp4                            -0.7982413
#> Sst                             -0.3350596
#>      CACAAACAGTAGCGGTstochastic_pseudotime
#> Pyy                             -0.7914314
#> Rbp4                            -0.7982381
#> Sst                             -0.3350595
#>      ACTGATGCATCAGTACstochastic_pseudotime
#> Pyy                             -0.7913943
#> Rbp4                            -0.7980041
#> Sst                             -0.3350567
#>      TAGTTGGCACGTAAGGstochastic_pseudotime
#> Pyy                             -0.7913878
#> Rbp4                            -0.7979621
#> Sst                             -0.3350557
#>      CAGCTAAGTACCGAGAstochastic_pseudotime
#> Pyy                             -0.7913807
#> Rbp4                            -0.7979148
#> Sst                             -0.3350546
#>      CGTCCATTCGCGCCAAstochastic_pseudotime
#> Pyy                             -0.7913749
#> Rbp4                            -0.7978761
#> Sst                             -0.3350536
#>      TTCGGTCAGCTAGCCCstochastic_pseudotime
#> Pyy                             -0.7913516
#> Rbp4                            -0.7977162
#> Sst                             -0.3350484
#>      GTTAAGCAGGATTCGGstochastic_pseudotime
#> Pyy                             -0.7913502
#> Rbp4                            -0.7977060
#> Sst                             -0.3350480
#>      GGACAGATCCTTTCGGstochastic_pseudotime
#> Pyy                             -0.7913480
#> Rbp4                            -0.7976906
#> Sst                             -0.3350474
#>      AGTTGGTTCACATAGCstochastic_pseudotime
#> Pyy                             -0.7913129
#> Rbp4                            -0.7974359
#> Sst                             -0.3350361
#>      TTCTCAACAATAGAGTstochastic_pseudotime
#> Pyy                             -0.7913011
#> Rbp4                            -0.7973469
#> Sst                             -0.3350314
#>      GTAGTCACATCGATGTstochastic_pseudotime
#> Pyy                             -0.7912935
#> Rbp4                            -0.7972886
#> Sst                             -0.3350281
#>      GCTGCGAAGTATTGGAstochastic_pseudotime
#> Pyy                             -0.7912873
#> Rbp4                            -0.7972398
#> Sst                             -0.3350253
#>      CGTCAGGGTAGGCATGstochastic_pseudotime
#> Pyy                             -0.7912824
#> Rbp4                            -0.7972015
#> Sst                             -0.3350230
#>      AACTCAGTCCTCGCATstochastic_pseudotime
#> Pyy                             -0.7912781
#> Rbp4                            -0.7971673
#> Sst                             -0.3350210
#>      AAATGCCCACCAACCGstochastic_pseudotime
#> Pyy                             -0.7912701
#> Rbp4                            -0.7971038
#> Sst                             -0.3350170
#>      CCGTTCAAGGACCACAstochastic_pseudotime
#> Pyy                             -0.7912472
#> Rbp4                            -0.7969175
#> Sst                             -0.3350044
#>      AATCGGTTCCAAACTGstochastic_pseudotime
#> Pyy                             -0.7912291
#> Rbp4                            -0.7967680
#> Sst                             -0.3349936
#>      TCACGAAGTCCGAGTCstochastic_pseudotime
#> Pyy                             -0.7911376
#> Rbp4                            -0.7960150
#> Sst                             -0.3349313
#>      GGCTGGTAGTCACGCCstochastic_pseudotime
#> Pyy                             -0.7911285
#> Rbp4                            -0.7959419
#> Sst                             -0.3349247
#>      GGACGTCCATGGATGGstochastic_pseudotime
#> Pyy                             -0.7910436
#> Rbp4                            -0.7952923
#> Sst                             -0.3348635
#>      CAGCCGATCGTTTATCstochastic_pseudotime
#> Pyy                             -0.7910355
#> Rbp4                            -0.7952336
#> Sst                             -0.3348578
#>      TTGACTTGTTCCATGAstochastic_pseudotime
#> Pyy                             -0.7910054
#> Rbp4                            -0.7950185
#> Sst                             -0.3348366
#>      GGCTGGTTCATGCATGstochastic_pseudotime
#> Pyy                              -0.790955
#> Rbp4                             -0.794672
#> Sst                              -0.334802
#>      ACACCCTAGGACAGCTstochastic_pseudotime
#> Pyy                             -0.7909422
#> Rbp4                            -0.7945868
#> Sst                             -0.3347934
#>      CGATGTAGTACTCTCCstochastic_pseudotime
#> Pyy                             -0.7909380
#> Rbp4                            -0.7945585
#> Sst                             -0.3347905
#>      TCAGGTAAGGCGTACAstochastic_pseudotime
#> Pyy                             -0.7907460
#> Rbp4                            -0.7933846
#> Sst                             -0.3346702
#>      AGCAGCCGTGTAACGGstochastic_pseudotime
#> Pyy                             -0.7905451
#> Rbp4                            -0.7923024
#> Sst                             -0.3345589
#>      CCTACCACATCGATTGstochastic_pseudotime
#> Pyy                             -0.7905436
#> Rbp4                            -0.7922951
#> Sst                             -0.3345582
#>      AAGACCTAGTAGGCCAstochastic_pseudotime
#> Pyy                             -0.7904733
#> Rbp4                            -0.7919408
#> Sst                             -0.3345221
#>      GGACAAGTCAAAGACAstochastic_pseudotime
#> Pyy                             -0.7904082
#> Rbp4                            -0.7916219
#> Sst                             -0.3344899
#>      GACTAACGTCATCCCTstochastic_pseudotime
#> Pyy                             -0.7903488
#> Rbp4                            -0.7913379
#> Sst                             -0.3344613
#>      CGCCAAGCAAAGTCAAstochastic_pseudotime
#> Pyy                             -0.7903259
#> Rbp4                            -0.7912298
#> Sst                             -0.3344505
#>      CCGGGATAGATCTGAAstochastic_pseudotime
#> Pyy                             -0.7901479
#> Rbp4                            -0.7904161
#> Sst                             -0.3343701
#>      TACACGAGTGCAACTTstochastic_pseudotime
#> Pyy                             -0.7900801
#> Rbp4                            -0.7901167
#> Sst                             -0.3343410
#>      CGTTCTGTCACAATGCstochastic_pseudotime
#> Pyy                             -0.7898968
#> Rbp4                            -0.7893317
#> Sst                             -0.3342658
#>      CAGCCGATCCTATTCAstochastic_pseudotime
#> Pyy                             -0.7897578
#> Rbp4                            -0.7887560
#> Sst                             -0.3342119
#>      AGGGATGTCGCATGGCstochastic_pseudotime
#> Pyy                             -0.7894817
#> Rbp4                            -0.7876538
#> Sst                             -0.3341115
#>      GTAGTCAGTCCGACGTstochastic_pseudotime
#> Pyy                             -0.7892904
#> Rbp4                            -0.7869168
#> Sst                             -0.3340465
#>      AGCAGCCGTCGCATCGstochastic_pseudotime
#> Pyy                             -0.7891969
#> Rbp4                            -0.7865632
#> Sst                             -0.3340159
#>      TGCTACCCATAAGACAstochastic_pseudotime
#> Pyy                             -0.7881570
#> Rbp4                            -0.7828520
#> Sst                             -0.3337173
#>      AGCTCCTTCAATACCGstochastic_pseudotime
#> Pyy                             -0.7881250
#> Rbp4                            -0.7827428
#> Sst                             -0.3337092
#>      GAATGAAAGAATGTGTstochastic_pseudotime
#> Pyy                             -0.7873312
#> Rbp4                            -0.7801074
#> Sst                             -0.3335208
#>      TTGTAGGTCAAAGTAGstochastic_pseudotime
#> Pyy                             -0.7870326
#> Rbp4                            -0.7791466
#> Sst                             -0.3334564
#>      AAACCTGGTAAGTGGCstochastic_pseudotime
#> Pyy                             -0.7869933
#> Rbp4                            -0.7790212
#> Sst                             -0.3334481
#>      AAATGCCCATGCTAGTstochastic_pseudotime
#> Pyy                             -0.7869354
#> Rbp4                            -0.7788368
#> Sst                             -0.3334361
#>      GACTGCGAGAATCTCCstochastic_pseudotime
#> Pyy                             -0.7868838
#> Rbp4                            -0.7786731
#> Sst                             -0.3334254
#>      CGATTGAAGTCAAGGCstochastic_pseudotime
#> Pyy                             -0.7867654
#> Rbp4                            -0.7782988
#> Sst                             -0.3334013
#>      CCATGTCTCCGCATAAstochastic_pseudotime
#> Pyy                             -0.7867443
#> Rbp4                            -0.7782323
#> Sst                             -0.3333970
#>      CAAGAAAAGAGTACCGstochastic_pseudotime
#> Pyy                             -0.7866669
#> Rbp4                            -0.7779890
#> Sst                             -0.3333816
#>      GGCGTGTCACCAGGTCstochastic_pseudotime
#> Pyy                             -0.7864514
#> Rbp4                            -0.7773160
#> Sst                             -0.3333394
#>      CCGTGGATCTTAGAGCstochastic_pseudotime
#> Pyy                             -0.7861919
#> Rbp4                            -0.7765138
#> Sst                             -0.3332904
#>      TGCGGGTAGACGCACAstochastic_pseudotime
#> Pyy                             -0.7859384
#> Rbp4                            -0.7757382
#> Sst                             -0.3332443
#>      AAGCCGCCAATGGACGstochastic_pseudotime
#> Pyy                             -0.7855317
#> Rbp4                            -0.7745095
#> Sst                             -0.3331735
#>      ACACTGAGTCCATGATstochastic_pseudotime
#> Pyy                             -0.7847699
#> Rbp4                            -0.7722547
#> Sst                             -0.3330509
#>      AGCTTGAAGACACGACstochastic_pseudotime
#> Pyy                             -0.7838104
#> Rbp4                            -0.7694901
#> Sst                             -0.3329120
#>      TCTATTGTCGTATCAGstochastic_pseudotime
#> Pyy                             -0.7837198
#> Rbp4                            -0.7692329
#> Sst                             -0.3328997
#>      TCGCGAGAGTACGCGAstochastic_pseudotime
#> Pyy                             -0.7833030
#> Rbp4                            -0.7680582
#> Sst                             -0.3328445
#>      TCAACGATCAACCATGstochastic_pseudotime
#> Pyy                             -0.7812479
#> Rbp4                            -0.7624382
#> Sst                             -0.3326056
#>      CACCAGGCATAAAGGTstochastic_pseudotime
#> Pyy                             -0.7812348
#> Rbp4                            -0.7624031
#> Sst                             -0.3326042
#>      TGAGCCGTCCGCTGTTstochastic_pseudotime
#> Pyy                             -0.7802239
#> Rbp4                            -0.7597313
#> Sst                             -0.3325033
#>      ACTGTCCCAGTTTACGstochastic_pseudotime
#> Pyy                             -0.7800411
#> Rbp4                            -0.7592540
#> Sst                             -0.3324861
#>      GATCGATAGAAAGTGGstochastic_pseudotime
#> Pyy                             -0.7794114
#> Rbp4                            -0.7576222
#> Sst                             -0.3324287
#>      GATCGTAAGTAGGTGCstochastic_pseudotime
#> Pyy                             -0.7789320
#> Rbp4                            -0.7563926
#> Sst                             -0.3323871
#>      CACAGGCCAAGGCTCCstochastic_pseudotime
#> Pyy                             -0.7778795
#> Rbp4                            -0.7537298
#> Sst                             -0.3323013
#>      TGTGTTTGTATGAATGstochastic_pseudotime
#> Pyy                             -0.7776866
#> Rbp4                            -0.7532471
#> Sst                             -0.3322863
#>      GGCTCGAAGATCACGGstochastic_pseudotime
#> Pyy                             -0.7772472
#> Rbp4                            -0.7521531
#> Sst                             -0.3322532
#>      GGGTCTGTCATGTGGTstochastic_pseudotime
#> Pyy                             -0.7763982
#> Rbp4                            -0.7500610
#> Sst                             -0.3321921
#>      TACAGTGGTAAACGCGstochastic_pseudotime
#> Pyy                             -0.7762973
#> Rbp4                            -0.7498141
#> Sst                             -0.3321851
#>      TCTCTAATCCCAACGGstochastic_pseudotime
#> Pyy                             -0.7760029
#> Rbp4                            -0.7490963
#> Sst                             -0.3321650
#>      AGTTGGTGTTGTCGCGstochastic_pseudotime
#> Pyy                             -0.7756266
#> Rbp4                            -0.7481833
#> Sst                             -0.3321399
#>      CCAGCGACAGACGTAGstochastic_pseudotime
#> Pyy                             -0.7752984
#> Rbp4                            -0.7473913
#> Sst                             -0.3321186
#>      TCAATCTCATCCGCGAstochastic_pseudotime
#> Pyy                             -0.7747885
#> Rbp4                            -0.7461685
#> Sst                             -0.3320864
#>      TCAGATGCACTCGACGstochastic_pseudotime
#> Pyy                             -0.7740744
#> Rbp4                            -0.7444709
#> Sst                             -0.3320433
#>      GTTACAGTCCGATATGstochastic_pseudotime
#> Pyy                             -0.7732052
#> Rbp4                            -0.7424272
#> Sst                             -0.3319935
#>      TGGACGCTCTTGTTTGstochastic_pseudotime
#> Pyy                             -0.7729897
#> Rbp4                            -0.7419241
#> Sst                             -0.3319816
#>      CTGCGGACATTAGCCAstochastic_pseudotime
#> Pyy                             -0.7705612
#> Rbp4                            -0.7363518
#> Sst                             -0.3318583
#>      CACCACTGTGTGTGCCstochastic_pseudotime
#> Pyy                             -0.7674059
#> Rbp4                            -0.7293513
#> Sst                             -0.3317229
#>      GTAACTGGTCACTGGCstochastic_pseudotime
#> Pyy                             -0.7656879
#> Rbp4                            -0.7256415
#> Sst                             -0.3316586
#>      GGCGTGTCACCGATATstochastic_pseudotime
#> Pyy                             -0.7648361
#> Rbp4                            -0.7238267
#> Sst                             -0.3316289
#>      AGAGTGGCACCTCGGAstochastic_pseudotime
#> Pyy                             -0.7634883
#> Rbp4                            -0.7209870
#> Sst                             -0.3315844
#>      CAGCTAACATGCCTAAstochastic_pseudotime
#> Pyy                             -0.7633028
#> Rbp4                            -0.7205990
#> Sst                             -0.3315785
#>      ATCCGAAGTACACCGCstochastic_pseudotime
#> Pyy                             -0.7632614
#> Rbp4                            -0.7205125
#> Sst                             -0.3315772
#>      CAAGGCCTCTCAACTTstochastic_pseudotime
#> Pyy                             -0.7632464
#> Rbp4                            -0.7204813
#> Sst                             -0.3315767
#>      ACTGATGCATGTCTCCstochastic_pseudotime
#> Pyy                             -0.7624413
#> Rbp4                            -0.7188065
#> Sst                             -0.3315519
#>      CATATGGTCACAGTACstochastic_pseudotime
#> Pyy                             -0.7623422
#> Rbp4                            -0.7186014
#> Sst                             -0.3315489
#>      TTCCCAGTCATAGCACstochastic_pseudotime
#> Pyy                             -0.7620501
#> Rbp4                            -0.7179974
#> Sst                             -0.3315401
#>      GGACGTCGTCAGAAGCstochastic_pseudotime
#> Pyy                             -0.7613065
#> Rbp4                            -0.7164676
#> Sst                             -0.3315184
#>      GGGCACTGTCTCCATCstochastic_pseudotime
#> Pyy                             -0.7603444
#> Rbp4                            -0.7145038
#> Sst                             -0.3314915
#>      CCCATACTCGTGGACCstochastic_pseudotime
#> Pyy                             -0.7593123
#> Rbp4                            -0.7124158
#> Sst                             -0.3314639
#>      TGCCCTAAGTCCGGTCstochastic_pseudotime
#> Pyy                             -0.7584403
#> Rbp4                            -0.7106663
#> Sst                             -0.3314415
#>      GCTGCGATCCCGGATGstochastic_pseudotime
#> Pyy                             -0.7574949
#> Rbp4                            -0.7087843
#> Sst                             -0.3314183
#>      TAGTGGTTCCTCAACCstochastic_pseudotime
#> Pyy                             -0.7511255
#> Rbp4                            -0.6964736
#> Sst                             -0.3312840
#>      CTACATTCACCGTTGGstochastic_pseudotime
#> Pyy                             -0.7502661
#> Rbp4                            -0.6948582
#> Sst                             -0.3312684
#>      AGGTCATGTTTAGGAAstochastic_pseudotime
#> Pyy                             -0.7496470
#> Rbp4                            -0.6937003
#> Sst                             -0.3312575
#>      CCTACCAGTTGCTCCTstochastic_pseudotime
#> Pyy                             -0.7490419
#> Rbp4                            -0.6925740
#> Sst                             -0.3312471
#>      CGTCACTTCACTTATCstochastic_pseudotime
#> Pyy                             -0.7488032
#> Rbp4                            -0.6921309
#> Sst                             -0.3312430
#>      CCTTCCCTCGTCGTTCstochastic_pseudotime
#> Pyy                             -0.7461273
#> Rbp4                            -0.6872139
#> Sst                             -0.3312001
#>      AATCCAGTCGGTCCGAstochastic_pseudotime
#> Pyy                             -0.7441367
#> Rbp4                            -0.6836132
#> Sst                             -0.3311708
#>      CAAGTTGTCTAACTGGstochastic_pseudotime
#> Pyy                             -0.7312881
#> Rbp4                            -0.6614013
#> Sst                             -0.3310232
#>      TTAGGCAGTTTGTTGGstochastic_pseudotime
#> Pyy                             -0.7293048
#> Rbp4                            -0.6581139
#> Sst                             -0.3310054
#>      AATCCAGTCAGGTTCAstochastic_pseudotime
#> Pyy                             -0.7255925
#> Rbp4                            -0.6520505
#> Sst                             -0.3309749
#>      GTACTCCCACAGGTTTstochastic_pseudotime
#> Pyy                             -0.7245596
#> Rbp4                            -0.6503834
#> Sst                             -0.3309670
#>      GCGCCAATCTTATCTGstochastic_pseudotime
#> Pyy                             -0.7224335
#> Rbp4                            -0.6469782
#> Sst                             -0.3309514
#>      TCATTTGCATTGGCGCstochastic_pseudotime
#> Pyy                             -0.7220135
#> Rbp4                            -0.6463097
#> Sst                             -0.3309485
#>      GACGTTAGTGTGAATAstochastic_pseudotime
#> Pyy                             -0.7147040
#> Rbp4                            -0.6348800
#> Sst                             -0.3309023
#>      AAAGCAAAGCTACCTAstochastic_pseudotime
#> Pyy                             -0.7122435
#> Rbp4                            -0.6311157
#> Sst                             -0.3308887
#>      ACATCAGAGAGACTTAstochastic_pseudotime
#> Pyy                             -0.7106574
#> Rbp4                            -0.6287098
#> Sst                             -0.3308805
#>      ACATACGGTCTCTCTGstochastic_pseudotime
#> Pyy                             -0.7069514
#> Rbp4                            -0.6231501
#> Sst                             -0.3308625
#>      TTAACTCGTAATCGTCstochastic_pseudotime
#> Pyy                             -0.7051758
#> Rbp4                            -0.6205159
#> Sst                             -0.3308545
#>      GTTCATTAGGCAGTCAstochastic_pseudotime
#> Pyy                             -0.7046216
#> Rbp4                            -0.6196975
#> Sst                             -0.3308521
#>      CATTATCGTCATATCGstochastic_pseudotime
#> Pyy                             -0.6986476
#> Rbp4                            -0.6109864
#> Sst                             -0.3308282
#>      TTTCCTCGTGGCGAATstochastic_pseudotime
#> Pyy                             -0.6865795
#> Rbp4                            -0.5939599
#> Sst                             -0.3307899
#>      GCAGTTACATTTCACTstochastic_pseudotime
#> Pyy                             -0.6844913
#> Rbp4                            -0.5910850
#> Sst                             -0.3307844
#>      CTTAACTTCAGCTTAGstochastic_pseudotime
#> Pyy                             -0.6802939
#> Rbp4                            -0.5853657
#> Sst                             -0.3307742
#>      AGATCTGTCCGTCAAAstochastic_pseudotime
#> Pyy                             -0.6755473
#> Rbp4                            -0.5789895
#> Sst                             -0.3307640
#>      TCTCATAAGAATGTGTstochastic_pseudotime
#> Pyy                             -0.6595530
#> Rbp4                            -0.5581600
#> Sst                             -0.3307375
#>      GTATCTTTCGTGGACCstochastic_pseudotime
#> Pyy                             -0.6569997
#> Rbp4                            -0.5549214
#> Sst                             -0.3307342
#>      GAAGCAGGTCGCGAAAstochastic_pseudotime
#> Pyy                             -0.6522815
#> Rbp4                            -0.5489952
#> Sst                             -0.3307288
#>      AACACGTGTCGACTGCstochastic_pseudotime
#> Pyy                             -0.6418770
#> Rbp4                            -0.5361802
#> Sst                             -0.3307192
#>      GTCTTCGTCTACCAGAstochastic_pseudotime
#> Pyy                             -0.6391610
#> Rbp4                            -0.5328895
#> Sst                             -0.3307172
#>      TAGTTGGCAGTAACGGstochastic_pseudotime
#> Pyy                             -0.6385462
#> Rbp4                            -0.5321476
#> Sst                             -0.3307168
#>      AATCCAGAGGGTTCCCstochastic_pseudotime
#> Pyy                             -0.6332918
#> Rbp4                            -0.5258515
#> Sst                             -0.3307134
#>      TCCCGATAGTCGTTTGstochastic_pseudotime
#> Pyy                             -0.6294911
#> Rbp4                            -0.5213454
#> Sst                             -0.3307114
#>      CGAGCACAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.6264774
#> Rbp4                            -0.5178003
#> Sst                             -0.3307100
#>      TAAGAGATCTGCGACGstochastic_pseudotime
#> Pyy                             -0.6130514
#> Rbp4                            -0.5022874
#> Sst                             -0.3307048
#>      TGGGAAGCACGTAAGGstochastic_pseudotime
#> Pyy                             -0.6111449
#> Rbp4                            -0.5001195
#> Sst                             -0.3307042
#>      CACACCTGTTATCCGAstochastic_pseudotime
#> Pyy                             -0.6079731
#> Rbp4                            -0.4965314
#> Sst                             -0.3307030
#>      GATCGCGCACCAGGCTstochastic_pseudotime
#> Pyy                             -0.6071073
#> Rbp4                            -0.4955558
#> Sst                             -0.3307027
#>      CTCGGAGGTGAAGGCTstochastic_pseudotime
#> Pyy                             -0.6059613
#> Rbp4                            -0.4942669
#> Sst                             -0.3307023
#>      AAGGTTCGTCAACTGTstochastic_pseudotime
#> Pyy                             -0.5930588
#> Rbp4                            -0.4799476
#> Sst                             -0.3306968
#>      CAGCTGGGTCTAACGTstochastic_pseudotime
#> Pyy                             -0.5800841
#> Rbp4                            -0.4658764
#> Sst                             -0.3306893
#>      CAACCAAAGCGTTCCGstochastic_pseudotime
#> Pyy                             -0.5737947
#> Rbp4                            -0.4591643
#> Sst                             -0.3306847
#>      TGGCCAGCAATGAAACstochastic_pseudotime
#> Pyy                             -0.5479548
#> Rbp4                            -0.4322565
#> Sst                             -0.3306581
#>      CTAATGGCATGGTTGTstochastic_pseudotime
#> Pyy                             -0.5419559
#> Rbp4                            -0.4261507
#> Sst                             -0.3306500
#>      CTTGGCTCAAAGGAAGstochastic_pseudotime
#> Pyy                             -0.5383523
#> Rbp4                            -0.4225065
#> Sst                             -0.3306448
#>      TCAGCAATCGCCGTGAstochastic_pseudotime
#> Pyy                             -0.5361548
#> Rbp4                            -0.4202926
#> Sst                             -0.3306415
#>      CGGTTAATCACCGGGTstochastic_pseudotime
#> Pyy                             -0.5306598
#> Rbp4                            -0.4147843
#> Sst                             -0.3306328
#>      CTACATTGTTACTGACstochastic_pseudotime
#> Pyy                             -0.5154494
#> Rbp4                            -0.3997329
#> Sst                             -0.3306055
#>      CTGCTGTAGACTGGGTstochastic_pseudotime
#> Pyy                             -0.5126073
#> Rbp4                            -0.3969509
#> Sst                             -0.3305998
#>      GCTCTGTGTCCATGATstochastic_pseudotime
#> Pyy                             -0.5069737
#> Rbp4                            -0.3914634
#> Sst                             -0.3305881
#>      CGACTTCCAGCTGCACstochastic_pseudotime
#> Pyy                             -0.5067424
#> Rbp4                            -0.3912389
#> Sst                             -0.3305876
#>      CAGCTGGCAAGTTGTCstochastic_pseudotime
#> Pyy                             -0.5017109
#> Rbp4                            -0.3863688
#> Sst                             -0.3305766
#>      CCTTACGCACAACGCCstochastic_pseudotime
#> Pyy                             -0.4969829
#> Rbp4                            -0.3818172
#> Sst                             -0.3305658
#>      TAGTGGTCATTCACTTstochastic_pseudotime
#> Pyy                             -0.4937074
#> Rbp4                            -0.3786777
#> Sst                             -0.3305581
#>      CAGGTGCCATCAGTCAstochastic_pseudotime
#> Pyy                             -0.4899906
#> Rbp4                            -0.3751284
#> Sst                             -0.3305490
#>      TACTCATAGGTGTTAAstochastic_pseudotime
#> Pyy                             -0.4841584
#> Rbp4                            -0.3695871
#> Sst                             -0.3305342
#>      AGCTCCTCAGTTTACGstochastic_pseudotime
#> Pyy                             -0.4806822
#> Rbp4                            -0.3663000
#> Sst                             -0.3305251
#>      GTCACGGCATCCGTGGstochastic_pseudotime
#> Pyy                             -0.4771400
#> Rbp4                            -0.3629625
#> Sst                             -0.3305155
#>      GTCACGGTCTTGCCGTstochastic_pseudotime
#> Pyy                             -0.4739970
#> Rbp4                            -0.3600110
#> Sst                             -0.3305069
#>      TTGACTTTCACCGGGTstochastic_pseudotime
#> Pyy                             -0.4629345
#> Rbp4                            -0.3496943
#> Sst                             -0.3304748
#>      CCTCAGTAGTGCGATGstochastic_pseudotime
#> Pyy                             -0.4583391
#> Rbp4                            -0.3454405
#> Sst                             -0.3304608
#>      TTCGGTCTCGACCAGCstochastic_pseudotime
#> Pyy                             -0.4506995
#> Rbp4                            -0.3384088
#> Sst                             -0.3304366
#>      CGTTGGGCAGGCGATAstochastic_pseudotime
#> Pyy                             -0.4284481
#> Rbp4                            -0.3181963
#> Sst                             -0.3303600
#>      CCACGGAAGTCTCCTCstochastic_pseudotime
#> Pyy                             -0.4245273
#> Rbp4                            -0.3146741
#> Sst                             -0.3303456
#>      GATGCTAGTGTTAAGAstochastic_pseudotime
#> Pyy                             -0.4030353
#> Rbp4                            -0.2955619
#> Sst                             -0.3302616
#>      GTGTTAGGTGAACCTTstochastic_pseudotime
#> Pyy                             -0.3926791
#> Rbp4                            -0.2864642
#> Sst                             -0.3302182
#>      AGCTCCTTCGCCCTTAstochastic_pseudotime
#> Pyy                             -0.3655384
#> Rbp4                            -0.2629393
#> Sst                             -0.3300959
#>      TGCGTGGGTCCCTTGTstochastic_pseudotime
#> Pyy                             -0.3629730
#> Rbp4                            -0.2607382
#> Sst                             -0.3300837
#>      GGAGCAAGTATCACCAstochastic_pseudotime
#> Pyy                             -0.3617043
#> Rbp4                            -0.2596510
#> Sst                             -0.3300776
#>      CGTAGCGTCGGAAATAstochastic_pseudotime
#> Pyy                             -0.3582023
#> Rbp4                            -0.2566548
#> Sst                             -0.3300607
#>      CTGGTCTAGGATTCGGstochastic_pseudotime
#> Pyy                             -0.3500905
#> Rbp4                            -0.2497405
#> Sst                             -0.3300208
#>      GGCGTGTTCGAGAGCAstochastic_pseudotime
#> Pyy                             -0.3230358
#> Rbp4                            -0.2269301
#> Sst                             -0.3298798
#>      AACCATGGTTGAACTCstochastic_pseudotime
#> Pyy                             -0.3099458
#> Rbp4                            -0.2160236
#> Sst                             -0.3298073
#>      ACAGCCGAGAGAACAGstochastic_pseudotime
#> Pyy                             -0.3071145
#> Rbp4                            -0.2136752
#> Sst                             -0.3297913
#>      CGAATGTGTCTGGAGAstochastic_pseudotime
#> Pyy                             -0.3046789
#> Rbp4                            -0.2116579
#> Sst                             -0.3297774
#>      GACGGCTTCAACTCTTstochastic_pseudotime
#> Pyy                             -0.2924844
#> Rbp4                            -0.2015978
#> Sst                             -0.3297064
#>      CTTTGCGAGACCACGAstochastic_pseudotime
#> Pyy                             -0.2725958
#> Rbp4                            -0.1853281
#> Sst                             -0.3295854
#>      TTCTACACATAAAGGTstochastic_pseudotime
#> Pyy                             -0.2722482
#> Rbp4                            -0.1850452
#> Sst                             -0.3295832
#>      GCATGATGTCACTGGCstochastic_pseudotime
#> Pyy                             -0.2639304
#> Rbp4                            -0.1782902
#> Sst                             -0.3295307
#>      ACACCCTTCCCAGGTGstochastic_pseudotime
#> Pyy                             -0.2476068
#> Rbp4                            -0.1651123
#> Sst                             -0.3294245
#>      CTCCTAGCACAACGTTstochastic_pseudotime
#> Pyy                             -0.2421849
#> Rbp4                            -0.1607574
#> Sst                             -0.3293883
#>      CTCTGGTCAAGTTGTCstochastic_pseudotime
#> Pyy                             -0.2399934
#> Rbp4                            -0.1590003
#> Sst                             -0.3293735
#>      CATCAGATCATCGCTCstochastic_pseudotime
#> Pyy                             -0.2192521
#> Rbp4                            -0.1424546
#> Sst                             -0.3292298
#>      CCACTACAGAGCTGGTstochastic_pseudotime
#> Pyy                             -0.2161256
#> Rbp4                            -0.1399734
#> Sst                             -0.3292075
#>      CTCGGAGGTCATTAGCstochastic_pseudotime
#> Pyy                             -0.2104512
#> Rbp4                            -0.1354786
#> Sst                             -0.3291667
#>      TGGTTAGCAGTCGTGCstochastic_pseudotime
#> Pyy                             -0.1705263
#> Rbp4                            -0.1041424
#> Sst                             -0.3288646
#>      TCGAGGCAGGCATGGTstochastic_pseudotime
#> Pyy                            -0.15329809
#> Rbp4                           -0.09076654
#> Sst                            -0.32872613
#>      GTACTTTTCTTCGGTCstochastic_pseudotime
#> Pyy                            -0.12781577
#> Rbp4                           -0.07113028
#> Sst                            -0.32851206
#>      ACGAGGATCGCCGTGAstochastic_pseudotime
#> Pyy                            -0.12712456
#> Rbp4                           -0.07059999
#> Sst                            -0.32850610
#>      ATGAGGGTCTTGACGAstochastic_pseudotime
#> Pyy                            -0.07325158
#> Rbp4                           -0.02962121
#> Sst                            -0.32801552
#>      CGATCGGAGCTGCAAGstochastic_pseudotime
#> Pyy                            -0.06961810
#> Rbp4                           -0.02688092
#> Sst                            -0.32798055
#>      GTAGGCCAGCTAACAAstochastic_pseudotime
#> Pyy                            -0.06173806
#> Rbp4                           -0.02094751
#> Sst                            -0.32790388
#>      GTCCTCACACCACGTGstochastic_pseudotime
#> Pyy                            -0.01483077
#> Rbp4                            0.01411746
#> Sst                            -0.32742311
#>      CTGAAGTTCTGCGGCAstochastic_pseudotime
#> Pyy                            0.009316151
#> Rbp4                           0.032011508
#> Sst                           -0.327158778
#>      GGACATTGTCTCCCTAstochastic_pseudotime
#> Pyy                             0.01353253
#> Rbp4                            0.03512597
#> Sst                            -0.32711141
#>      CGAGCACGTAGTAGTAstochastic_pseudotime
#> Pyy                             0.03220123
#> Rbp4                            0.04888145
#> Sst                            -0.32689722
#>      GTGCATAAGAACAATCstochastic_pseudotime
#> Pyy                             0.03842212
#> Rbp4                            0.05345309
#> Sst                            -0.32682421
#>      ATCTGCCCATGTTGACstochastic_pseudotime
#> Pyy                             0.04395866
#> Rbp4                            0.05751690
#> Sst                            -0.32675854
#>      ACTGCTCAGGTGTTAAstochastic_pseudotime
#> Pyy                             0.05747486
#> Rbp4                            0.06741892
#> Sst                            -0.32659541
#>      AGCTCTCAGACATAACstochastic_pseudotime
#> Pyy                             0.06701968
#> Rbp4                            0.07439586
#> Sst                            -0.32647779
#>      TGACGGCCATATGAGAstochastic_pseudotime
#> Pyy                             0.07600593
#> Rbp4                            0.08095309
#> Sst                            -0.32636518
#>      CAAGTTGGTGTGCCTGstochastic_pseudotime
#> Pyy                             0.09351296
#> Rbp4                            0.09369734
#> Sst                            -0.32614050
#>      GATCAGTAGGTGATATstochastic_pseudotime
#> Pyy                              0.1062553
#> Rbp4                             0.1029488
#> Sst                             -0.3259725
#>      ACGCCGAAGGGTGTTGstochastic_pseudotime
#> Pyy                              0.1696444
#> Rbp4                             0.1486989
#> Sst                             -0.3250770
#>      GGACGTCAGGGTGTTGstochastic_pseudotime
#> Pyy                              0.1781179
#> Rbp4                             0.1547834
#> Sst                             -0.3249494
#>      GGACAGACAATCCAACstochastic_pseudotime
#> Pyy                              0.1911083
#> Rbp4                             0.1640985
#> Sst                             -0.3247501
#>      TCGCGAGCACTTCTGCstochastic_pseudotime
#> Pyy                              0.2179757
#> Rbp4                             0.1833185
#> Sst                             -0.3243227
#>      TAAGTGCTCCAAATGCstochastic_pseudotime
#> Pyy                              0.2295003
#> Rbp4                             0.1915450
#> Sst                             -0.3241330
#>      TGCGTGGAGGTTCCTAstochastic_pseudotime
#> Pyy                              0.2621559
#> Rbp4                             0.2148033
#> Sst                             -0.3235736
#>      ACTATCTTCAAGATCCstochastic_pseudotime
#> Pyy                              0.2652268
#> Rbp4                             0.2169868
#> Sst                             -0.3235192
#>      TGTGTTTAGGAATGGAstochastic_pseudotime
#> Pyy                              0.2679938
#> Rbp4                             0.2189537
#> Sst                             -0.3234700
#>      GAGCAGATCAACTCTTstochastic_pseudotime
#> Pyy                              0.2764312
#> Rbp4                             0.2249485
#> Sst                             -0.3233185
#>      ACGATGTTCGCAGGCTstochastic_pseudotime
#> Pyy                              0.2813534
#> Rbp4                             0.2284437
#> Sst                             -0.3232290
#>      GTCATTTGTCGTCTTCstochastic_pseudotime
#> Pyy                              0.2823963
#> Rbp4                             0.2291842
#> Sst                             -0.3232100
#>      TTTGGTTTCCTTTCGGstochastic_pseudotime
#> Pyy                              0.2959760
#> Rbp4                             0.2388191
#> Sst                             -0.3229585
#>      AAAGCAAGTGCAGTAGstochastic_pseudotime
#> Pyy                              0.3281583
#> Rbp4                             0.2616147
#> Sst                             -0.3223372
#>      CACACCTGTACATGTCstochastic_pseudotime
#> Pyy                              0.3601428
#> Rbp4                             0.2842248
#> Sst                             -0.3216826
#>      CGATGTAAGTTTAGGAstochastic_pseudotime
#> Pyy                              0.3751247
#> Rbp4                             0.2948027
#> Sst                             -0.3213626
#>      GCGCAGTGTGCTAGCCstochastic_pseudotime
#> Pyy                              0.3769620
#> Rbp4                             0.2960993
#> Sst                             -0.3213227
#>      CTGAAACGTAGCGCAAstochastic_pseudotime
#> Pyy                              0.4042477
#> Rbp4                             0.3153456
#> Sst                             -0.3207147
#>      CGTGTCTGTAGCGTCCstochastic_pseudotime
#> Pyy                              0.4077701
#> Rbp4                             0.3178288
#> Sst                             -0.3206340
#>      AGTGGGAGTGGCAAACstochastic_pseudotime
#> Pyy                              0.4271178
#> Rbp4                             0.3314637
#> Sst                             -0.3201809
#>      AGGCCACCATGACATCstochastic_pseudotime
#> Pyy                              0.4357535
#> Rbp4                             0.3375474
#> Sst                             -0.3199734
#>      TTTGTCAAGAATGTGTstochastic_pseudotime
#> Pyy                              0.4361821
#> Rbp4                             0.3378493
#> Sst                             -0.3199630
#>      AAGGAGCAGTACATGAstochastic_pseudotime
#> Pyy                              0.4405215
#> Rbp4                             0.3409059
#> Sst                             -0.3198574
#>      ACACCGGTCATGTCTTstochastic_pseudotime
#> Pyy                              0.4412172
#> Rbp4                             0.3413959
#> Sst                             -0.3198404
#>      GAGGTGAGTCCGCTGAstochastic_pseudotime
#> Pyy                              0.4421428
#> Rbp4                             0.3420478
#> Sst                             -0.3198177
#>      CTCGAAAGTGGTGTAGstochastic_pseudotime
#> Pyy                              0.4551603
#> Rbp4                             0.3512152
#> Sst                             -0.3194947
#>      AAGGTTCAGATAGTCAstochastic_pseudotime
#> Pyy                              0.4937098
#> Rbp4                             0.3783556
#> Sst                             -0.3184905
#>      AGTGTCATCCCTTGTGstochastic_pseudotime
#> Pyy                              0.5068047
#> Rbp4                             0.3875741
#> Sst                             -0.3181323
#>      ATCACGAGTGCACCACstochastic_pseudotime
#> Pyy                              0.5254666
#> Rbp4                             0.4007127
#> Sst                             -0.3176063
#>      TACGGGCAGCCAACAGstochastic_pseudotime
#> Pyy                              0.5687299
#> Rbp4                             0.4311845
#> Sst                             -0.3163120
#>      CTCAGAACATAGGATAstochastic_pseudotime
#> Pyy                              0.5828911
#> Rbp4                             0.4411652
#> Sst                             -0.3158643
#>      TATCTCACACGCGAAAstochastic_pseudotime
#> Pyy                              0.6036208
#> Rbp4                             0.4557833
#> Sst                             -0.3151863
#>      GATTCAGTCGCCAGCAstochastic_pseudotime
#> Pyy                              0.6122917
#> Rbp4                             0.4619011
#> Sst                             -0.3148944
#>      TCACAAGCAACGATGGstochastic_pseudotime
#> Pyy                              0.6222322
#> Rbp4                             0.4689174
#> Sst                             -0.3145536
#>      TGTGTTTGTAGATTAGstochastic_pseudotime
#> Pyy                              0.6311351
#> Rbp4                             0.4752038
#> Sst                             -0.3142426
#>      GGCGTGTTCAGCTTAGstochastic_pseudotime
#> Pyy                              0.6348739
#> Rbp4                             0.4778446
#> Sst                             -0.3141104
#>      AGTAGTCGTGAGGGAGstochastic_pseudotime
#> Pyy                              0.6363129
#> Rbp4                             0.4788611
#> Sst                             -0.3140592
#>      GTACGTAGTCTCTTTAstochastic_pseudotime
#> Pyy                              0.6421269
#> Rbp4                             0.4829689
#> Sst                             -0.3138510
#>      GAATGAATCATCGGATstochastic_pseudotime
#> Pyy                              0.6421359
#> Rbp4                             0.4829754
#> Sst                             -0.3138507
#>      GCTGCGAAGAACAACTstochastic_pseudotime
#> Pyy                              0.6587535
#> Rbp4                             0.4947232
#> Sst                             -0.3132422
#>      TCTGGAAGTCATCGGCstochastic_pseudotime
#> Pyy                              0.6917780
#> Rbp4                             0.5181038
#> Sst                             -0.3119710
#>      CTGCTGTCAAAGCAATstochastic_pseudotime
#> Pyy                              0.6965099
#> Rbp4                             0.5214580
#> Sst                             -0.3117818
#>      CATCAGAGTCTCTTTAstochastic_pseudotime
#> Pyy                              0.7008619
#> Rbp4                             0.5245438
#> Sst                             -0.3116062
#>      GCAATCACAGCTGCTGstochastic_pseudotime
#> Pyy                              0.7096939
#> Rbp4                             0.5308092
#> Sst                             -0.3112450
#>      CCCAGTTTCAGCACATstochastic_pseudotime
#> Pyy                              0.7354834
#> Rbp4                             0.5491291
#> Sst                             -0.3101524
#>      CACAGTAGTGCACTTAstochastic_pseudotime
#> Pyy                              0.7453353
#> Rbp4                             0.5561381
#> Sst                             -0.3097196
#>      GAACCTATCACTTATCstochastic_pseudotime
#> Pyy                              0.7482641
#> Rbp4                             0.5582230
#> Sst                             -0.3095892
#>      GGCAATTCACGGCTACstochastic_pseudotime
#> Pyy                              0.7535708
#> Rbp4                             0.5620021
#> Sst                             -0.3093510
#>      TTGAACGGTGAGGGAGstochastic_pseudotime
#> Pyy                              0.7541134
#> Rbp4                             0.5623885
#> Sst                             -0.3093265
#>      GATTCAGAGCGATGACstochastic_pseudotime
#> Pyy                              0.7561764
#> Rbp4                             0.5638583
#> Sst                             -0.3092331
#>      GCATGCGAGGAACTGCstochastic_pseudotime
#> Pyy                              0.7602553
#> Rbp4                             0.5667650
#> Sst                             -0.3090473
#>      GAACGGAAGTCAAGCGstochastic_pseudotime
#> Pyy                              0.7705726
#> Rbp4                             0.5741229
#> Sst                             -0.3085701
#>      ACCAGTAGTTCGTTGAstochastic_pseudotime
#> Pyy                              0.7909096
#> Rbp4                             0.5886501
#> Sst                             -0.3075994
#>      GGAAAGCCAATGTAAGstochastic_pseudotime
#> Pyy                              0.7934037
#> Rbp4                             0.5904340
#> Sst                             -0.3074775
#>      CTCGTACTCTTTAGTCstochastic_pseudotime
#> Pyy                              0.8165952
#> Rbp4                             0.6070471
#> Sst                             -0.3063132
#>      GTGCTTCAGAGGGATAstochastic_pseudotime
#> Pyy                              0.8224572
#> Rbp4                             0.6112540
#> Sst                             -0.3060099
#>      AGGGAGTTCGGAGGTAstochastic_pseudotime
#> Pyy                              0.8361871
#> Rbp4                             0.6211204
#> Sst                             -0.3052847
#>      CTAACTTAGAGGTACCstochastic_pseudotime
#> Pyy                              0.8386866
#> Rbp4                             0.6229186
#> Sst                             -0.3051504
#>      CAGGTGCAGCCAGGATstochastic_pseudotime
#> Pyy                              0.8434050
#> Rbp4                             0.6263147
#> Sst                             -0.3048950
#>      GCATGATGTTTAGGAAstochastic_pseudotime
#> Pyy                              0.8584270
#> Rbp4                             0.6371428
#> Sst                             -0.3040645
#>      CCACGGAAGGGTGTTGstochastic_pseudotime
#> Pyy                              0.8877669
#> Rbp4                             0.6583644
#> Sst                             -0.3023635
#>      CATGACATCATGCATGstochastic_pseudotime
#> Pyy                              0.8878216
#> Rbp4                             0.6584041
#> Sst                             -0.3023603
#>      GACGTGCGTCACCTAAstochastic_pseudotime
#> Pyy                              0.8952590
#> Rbp4                             0.6638001
#> Sst                             -0.3019117
#>      CCACGGAAGATCCTGTstochastic_pseudotime
#> Pyy                              0.9058747
#> Rbp4                             0.6715142
#> Sst                             -0.3012588
#>      TTGGAACAGAAGATTCstochastic_pseudotime
#> Pyy                              0.9118185
#> Rbp4                             0.6758399
#> Sst                             -0.3008866
#>      CCACCTATCTAAGCCAstochastic_pseudotime
#> Pyy                              0.9154205
#> Rbp4                             0.6784636
#> Sst                             -0.3006587
#>      AGGTCATGTCTCCCTAstochastic_pseudotime
#> Pyy                              0.9299331
#> Rbp4                             0.6890528
#> Sst                             -0.2997219
#>      ATAGACCAGTGATCGGstochastic_pseudotime
#> Pyy                              0.9396138
#> Rbp4                             0.6961333
#> Sst                             -0.2990801
#>      GAACCTAAGTAGATGTstochastic_pseudotime
#> Pyy                              0.9435524
#> Rbp4                             0.6990179
#> Sst                             -0.2988150
#>      GATCGATGTTGTTTGGstochastic_pseudotime
#> Pyy                              0.9644868
#> Rbp4                             0.7143907
#> Sst                             -0.2973662
#>      TTCTCCTAGTCCATACstochastic_pseudotime
#> Pyy                              0.9649651
#> Rbp4                             0.7147427
#> Sst                             -0.2973323
#>      AACTCCCCAGTTCATGstochastic_pseudotime
#> Pyy                              0.9721728
#> Rbp4                             0.7200525
#> Sst                             -0.2968169
#>      GGGATGAGTCTGCGGTstochastic_pseudotime
#> Pyy                              0.9737081
#> Rbp4                             0.7211846
#> Sst                             -0.2967061
#>      CTCGAAATCCTGCCATstochastic_pseudotime
#> Pyy                              0.9740530
#> Rbp4                             0.7214390
#> Sst                             -0.2966811
#>      CGGGTCAAGGCGATACstochastic_pseudotime
#> Pyy                               1.001401
#> Rbp4                              0.741676
#> Sst                              -0.294638
#>      GACCTGGGTTCGTGATstochastic_pseudotime
#> Pyy                              1.0057702
#> Rbp4                             0.7449214
#> Sst                             -0.2942996
#>      CCTACACAGTACGTTCstochastic_pseudotime
#> Pyy                              1.0286499
#> Rbp4                             0.7619777
#> Sst                             -0.2924696
#>      GACCAATGTTATCCGAstochastic_pseudotime
#> Pyy                              1.0451188
#> Rbp4                             0.7743206
#> Sst                             -0.2910900
#>      CATCGGGCACGTTGGCstochastic_pseudotime
#> Pyy                              1.0505070
#> Rbp4                             0.7783715
#> Sst                             -0.2906267
#>      ACTTGTTTCACCCGAGstochastic_pseudotime
#> Pyy                              1.0677622
#> Rbp4                             0.7913876
#> Sst                             -0.2891019
#>      CCGTTCAAGCTAAACAstochastic_pseudotime
#> Pyy                              1.0748849
#> Rbp4                             0.7967807
#> Sst                             -0.2884535
#>      ATCATGGAGGAATGGAstochastic_pseudotime
#> Pyy                              1.0820688
#> Rbp4                             0.8022323
#> Sst                             -0.2877881
#>      AGGGATGGTCGACTGCstochastic_pseudotime
#> Pyy                              1.0843277
#> Rbp4                             0.8039492
#> Sst                             -0.2875764
#>      ACCAGTAAGTGATCGGstochastic_pseudotime
#> Pyy                              1.0876198
#> Rbp4                             0.8064535
#> Sst                             -0.2872657
#>      GCGACCATCACCAGGCstochastic_pseudotime
#> Pyy                              1.0892088
#> Rbp4                             0.8076633
#> Sst                             -0.2871148
#>      CGAGAAGTCGTTACAGstochastic_pseudotime
#> Pyy                              1.0964873
#> Rbp4                             0.8132128
#> Sst                             -0.2864162
#>      CTCTAATAGTTATCGCstochastic_pseudotime
#> Pyy                              1.1142710
#> Rbp4                             0.8268306
#> Sst                             -0.2846543
#>      CTACGTCAGCTAGGCAstochastic_pseudotime
#> Pyy                              1.1160173
#> Rbp4                             0.8281724
#> Sst                             -0.2844770
#>      ACTTGTTAGTGAAGTTstochastic_pseudotime
#> Pyy                              1.1292571
#> Rbp4                             0.8383737
#> Sst                             -0.2831064
#>      CGTTCTGAGGACGAAAstochastic_pseudotime
#> Pyy                              1.1387184
#> Rbp4                             0.8456947
#> Sst                             -0.2820979
#>      GTGCATAGTTCAACCAstochastic_pseudotime
#> Pyy                              1.1478173
#> Rbp4                             0.8527608
#> Sst                             -0.2811041
#>      ACAGCTAAGTCGTTTGstochastic_pseudotime
#> Pyy                              1.1559775
#> Rbp4                             0.8591198
#> Sst                             -0.2801924
#>      CACCACTAGCGATCCCstochastic_pseudotime
#> Pyy                              1.1582976
#> Rbp4                             0.8609317
#> Sst                             -0.2799296
#>      CATGCCTCACCCATTCstochastic_pseudotime
#> Pyy                              1.1750913
#> Rbp4                             0.8740995
#> Sst                             -0.2779773
#>      CAAGGCCTCAGTGCATstochastic_pseudotime
#> Pyy                              1.1772900
#> Rbp4                             0.8758306
#> Sst                             -0.2777151
#>      CCGGGATTCGACAGCCstochastic_pseudotime
#> Pyy                              1.1832029
#> Rbp4                             0.8804942
#> Sst                             -0.2770021
#>      CCACGGAAGCAATCTCstochastic_pseudotime
#> Pyy                              1.2033045
#> Rbp4                             0.8964436
#> Sst                             -0.2744888
#>      GTGGGTCTCCGAATGTstochastic_pseudotime
#> Pyy                              1.2068064
#> Rbp4                             0.8992377
#> Sst                             -0.2740363
#>      GCATGTAGTTGCGCACstochastic_pseudotime
#> Pyy                              1.2074629
#> Rbp4                             0.8997621
#> Sst                             -0.2739510
#>      TGCGTGGAGTACACCTstochastic_pseudotime
#> Pyy                              1.2133852
#> Rbp4                             0.9044998
#> Sst                             -0.2731740
#>      CGCTTCATCTCTGTCGstochastic_pseudotime
#> Pyy                              1.2263554
#> Rbp4                             0.9149251
#> Sst                             -0.2714255
#>      AGGTCATCACTTCTGCstochastic_pseudotime
#> Pyy                              1.2333707
#> Rbp4                             0.9205933
#> Sst                             -0.2704520
#>      TATGCCCAGAACAATCstochastic_pseudotime
#> Pyy                              1.2364903
#> Rbp4                             0.9231207
#> Sst                             -0.2700127
#>      ATAGACCAGCGTGTCCstochastic_pseudotime
#> Pyy                              1.2432387
#> Rbp4                             0.9286026
#> Sst                             -0.2690482
#>      CGCCAAGTCTAAGCCAstochastic_pseudotime
#> Pyy                              1.2507856
#> Rbp4                             0.9347576
#> Sst                             -0.2679464
#>      CGCGGTAAGCTGCAAGstochastic_pseudotime
#> Pyy                              1.2568899
#> Rbp4                             0.9397554
#> Sst                             -0.2670367
#>      AGAATAGCAATCACACstochastic_pseudotime
#> Pyy                              1.2615507
#> Rbp4                             0.9435834
#> Sst                             -0.2663307
#>      CGTCTACAGGGAGTAAstochastic_pseudotime
#> Pyy                              1.2627236
#> Rbp4                             0.9445483
#> Sst                             -0.2661514
#>      GGGAATGCAAAGTGCGstochastic_pseudotime
#> Pyy                              1.2643837
#> Rbp4                             0.9459153
#> Sst                             -0.2658966
#>      AACTTTCTCCAGATCAstochastic_pseudotime
#> Pyy                              1.2681036
#> Rbp4                             0.9489833
#> Sst                             -0.2653208
#>      AGAGTGGAGGCGATACstochastic_pseudotime
#> Pyy                              1.2726285
#> Rbp4                             0.9527245
#> Sst                             -0.2646114
#>      AAACGGGCAAAGAATCstochastic_pseudotime
#> Pyy                              1.2736890
#> Rbp4                             0.9536029
#> Sst                             -0.2644437
#>      GGGATGAAGGGCATGTstochastic_pseudotime
#> Pyy                              1.2832615
#> Rbp4                             0.9615578
#> Sst                             -0.2629043
#>      AGTCTTTCACGGACAAstochastic_pseudotime
#> Pyy                              1.2859342
#> Rbp4                             0.9637874
#> Sst                             -0.2624662
#>      GCAAACTCAGACAAATstochastic_pseudotime
#> Pyy                              1.2887860
#> Rbp4                             0.9661709
#> Sst                             -0.2619945
#>      ACACCAAGTATAGGTAstochastic_pseudotime
#> Pyy                              1.3044998
#> Rbp4                             0.9793856
#> Sst                             -0.2593160
#>      ACGCAGCCATTCTCATstochastic_pseudotime
#> Pyy                              1.3100604
#> Rbp4                             0.9840964
#> Sst                             -0.2583347
#>      CTAACTTCATCTGGTAstochastic_pseudotime
#> Pyy                              1.3130647
#> Rbp4                             0.9866493
#> Sst                             -0.2577969
#>      TGGGCGTGTGTGCGTCstochastic_pseudotime
#> Pyy                              1.3231851
#> Rbp4                             0.9952904
#> Sst                             -0.2559447
#>      TGAGAGGCAGGAACGTstochastic_pseudotime
#> Pyy                              1.3443160
#> Rbp4                             1.0135504
#> Sst                             -0.2518631
#>      AAGACCTAGACCTTTGstochastic_pseudotime
#> Pyy                              1.3462962
#> Rbp4                             1.0152776
#> Sst                             -0.2514649
#>      GCGCCAAAGCTCCTTCstochastic_pseudotime
#> Pyy                               1.347102
#> Rbp4                              1.015981
#> Sst                              -0.251302
#>      GCATGTATCGTGGACCstochastic_pseudotime
#> Pyy                              1.3477648
#> Rbp4                             1.0165603
#> Sst                             -0.2511677
#>      GGGAGATAGTGACATAstochastic_pseudotime
#> Pyy                              1.3526293
#> Rbp4                             1.0208205
#> Sst                             -0.2501719
#>      GTTCTCGAGAATTCCCstochastic_pseudotime
#> Pyy                              1.3563488
#> Rbp4                             1.0240899
#> Sst                             -0.2493986
#>      GGGCATCCATTTCAGGstochastic_pseudotime
#> Pyy                              1.3592086
#> Rbp4                             1.0266108
#> Sst                             -0.2487969
#>      TATTACCTCGTTGCCTstochastic_pseudotime
#> Pyy                              1.3659098
#> Rbp4                             1.0325430
#> Sst                             -0.2473617
#>      CGAATGTTCTTCATGTstochastic_pseudotime
#> Pyy                              1.3787854
#> Rbp4                             1.0440442
#> Sst                             -0.2445006
#>      CGGAGTCAGTTTAGGAstochastic_pseudotime
#> Pyy                              1.3829399
#> Rbp4                             1.0477855
#> Sst                             -0.2435468
#>      AGATTGCAGCGTGTCCstochastic_pseudotime
#> Pyy                              1.3859801
#> Rbp4                             1.0505330
#> Sst                             -0.2428391
#>      GGATGTTAGTAAGTACstochastic_pseudotime
#> Pyy                              1.3923172
#> Rbp4                             1.0562871
#> Sst                             -0.2413363
#>      GATCGATAGGATGGTCstochastic_pseudotime
#> Pyy                              1.3928430
#> Rbp4                             1.0567662
#> Sst                             -0.2412099
#>      TCAGGTACAGTCACTAstochastic_pseudotime
#> Pyy                              1.4070645
#> Rbp4                             1.0698262
#> Sst                             -0.2376872
#>      CCTTCGATCACAGGCCstochastic_pseudotime
#> Pyy                              1.4122883
#> Rbp4                             1.0746748
#> Sst                             -0.2363405
#>      GGCGACTGTAAGTTCCstochastic_pseudotime
#> Pyy                              1.4129309
#> Rbp4                             1.0752732
#> Sst                             -0.2361728
#>      TTCGGTCGTATAGTAGstochastic_pseudotime
#> Pyy                              1.4140559
#> Rbp4                             1.0763220
#> Sst                             -0.2358781
#>      CTTCTCTCAGCATACTstochastic_pseudotime
#> Pyy                              1.4146249
#> Rbp4                             1.0768529
#> Sst                             -0.2357285
#>      TTGCCGTCACCGAATTstochastic_pseudotime
#> Pyy                              1.4182340
#> Rbp4                             1.0802288
#> Sst                             -0.2347711
#>      CGTGTCTCAGAAGCACstochastic_pseudotime
#> Pyy                              1.4210961
#> Rbp4                             1.0829162
#> Sst                             -0.2340013
#>      AACTCTTCATTTGCCCstochastic_pseudotime
#> Pyy                              1.4247186
#> Rbp4                             1.0863308
#> Sst                             -0.2330132
#>      GCATGCGCAGCTCCGAstochastic_pseudotime
#> Pyy                              1.4249087
#> Rbp4                             1.0865105
#> Sst                             -0.2329609
#>      TTTGGTTAGACAGGCTstochastic_pseudotime
#> Pyy                              1.4322677
#> Rbp4                             1.0934960
#> Sst                             -0.2309028
#>      CGATGGCTCAGTGTTGstochastic_pseudotime
#> Pyy                              1.4403773
#> Rbp4                             1.1012707
#> Sst                             -0.2285547
#>      TACTCGCAGTAAGTACstochastic_pseudotime
#> Pyy                              1.4431440
#> Rbp4                             1.1039423
#> Sst                             -0.2277335
#>      TCAATCTTCGAACGGAstochastic_pseudotime
#> Pyy                              1.4467719
#> Rbp4                             1.1074607
#> Sst                             -0.2266405
#>      GCTGGGTGTCATACTGstochastic_pseudotime
#> Pyy                              1.4470518
#> Rbp4                             1.1077330
#> Sst                             -0.2265554
#>      ACGATACAGATGAGAGstochastic_pseudotime
#> Pyy                              1.4487690
#> Rbp4                             1.1094052
#> Sst                             -0.2260308
#>      GTATCTTGTACTCAACstochastic_pseudotime
#> Pyy                              1.4539273
#> Rbp4                             1.1144530
#> Sst                             -0.2244291
#>      ACGATACTCGCGCCAAstochastic_pseudotime
#> Pyy                              1.4567105
#> Rbp4                             1.1171923
#> Sst                             -0.2235482
#>      GATCGATGTAGCGATGstochastic_pseudotime
#> Pyy                              1.4617919
#> Rbp4                             1.1222228
#> Sst                             -0.2219088
#>      CCTACACTCTTTACACstochastic_pseudotime
#> Pyy                              1.4686240
#> Rbp4                             1.1290486
#> Sst                             -0.2196381
#>      AAGGTTCAGGAGCGAGstochastic_pseudotime
#> Pyy                              1.4712260
#> Rbp4                             1.1316677
#> Sst                             -0.2187525
#>      CAGTCCTCAGGGCATAstochastic_pseudotime
#> Pyy                               1.481040
#> Rbp4                              1.141648
#> Sst                              -0.215302
#>      CCCATACCATGCAACTstochastic_pseudotime
#> Pyy                              1.4849843
#> Rbp4                             1.1457072
#> Sst                             -0.2138637
#>      GCTTGAAGTCAACTGTstochastic_pseudotime
#> Pyy                              1.4851820
#> Rbp4                             1.1459113
#> Sst                             -0.2137908
#>      CTTCTCTAGTCCGTATstochastic_pseudotime
#> Pyy                              1.4902668
#> Rbp4                             1.1511880
#> Sst                             -0.2118884
#>      AGCTTGACATCTATGGstochastic_pseudotime
#> Pyy                               1.496734
#> Rbp4                              1.157971
#> Sst                              -0.209390
#>      CTCGAGGCAAGTCTGTstochastic_pseudotime
#> Pyy                              1.4980806
#> Rbp4                             1.1593930
#> Sst                             -0.2088585
#>      CAGAATCCATTGGGCCstochastic_pseudotime
#> Pyy                              1.4984190
#> Rbp4                             1.1597510
#> Sst                             -0.2087242
#>      CTTCTCTGTCAGAGGTstochastic_pseudotime
#> Pyy                              1.5064631
#> Rbp4                             1.1683319
#> Sst                             -0.2054548
#>      GATGAAACAGGCTGAAstochastic_pseudotime
#> Pyy                               1.516113
#> Rbp4                              1.178809
#> Sst                              -0.201324
#>      ATTACTCTCTGCCAGGstochastic_pseudotime
#> Pyy                              1.5320901
#> Rbp4                             1.1966355
#> Sst                             -0.1939251
#>      GGGTCTGTCCATGAACstochastic_pseudotime
#> Pyy                              1.5398800
#> Rbp4                             1.2055644
#> Sst                             -0.1900342
#>      ACCTTTATCCTGCAGGstochastic_pseudotime
#> Pyy                              1.5409138
#> Rbp4                             1.2067619
#> Sst                             -0.1895026
#>      GGGAGATAGATGTGGCstochastic_pseudotime
#> Pyy                              1.5429606
#> Rbp4                             1.2091418
#> Sst                             -0.1884392
#>      ATCCACCCAGCTCGCAstochastic_pseudotime
#> Pyy                              1.5479107
#> Rbp4                             1.2149477
#> Sst                             -0.1858055
#>      TCGAGGCAGAACAATCstochastic_pseudotime
#> Pyy                              1.5492463
#> Rbp4                             1.2165267
#> Sst                             -0.1850793
#>      GCTTCCAAGTTAAGTGstochastic_pseudotime
#> Pyy                              1.5524232
#> Rbp4                             1.2203046
#> Sst                             -0.1833247
#>      AAATGCCCAATGGAATstochastic_pseudotime
#> Pyy                              1.5527356
#> Rbp4                             1.2206777
#> Sst                             -0.1831501
#>      AGTTGGTAGAGCCTAGstochastic_pseudotime
#> Pyy                              1.5612056
#> Rbp4                             1.2309165
#> Sst                             -0.1782621
#>      GTGTGCGGTGTAAGTAstochastic_pseudotime
#> Pyy                              1.5678043
#> Rbp4                             1.2390623
#> Sst                             -0.1742377
#>      CGGAGCTAGATAGTCAstochastic_pseudotime
#> Pyy                              1.5686187
#> Rbp4                             1.2400785
#> Sst                             -0.1737269
#>      GAACGGACAAGGACACstochastic_pseudotime
#> Pyy                              1.5686453
#> Rbp4                             1.2401117
#> Sst                             -0.1737102
#>      ACAGCTATCAAGGCTTstochastic_pseudotime
#> Pyy                              1.5727029
#> Rbp4                             1.2452120
#> Sst                             -0.1711168
#>      GGGACCTCACTGTCGGstochastic_pseudotime
#> Pyy                              1.5785145
#> Rbp4                             1.2526285
#> Sst                             -0.1672554
#>      CCCAATCCAACTGGCCstochastic_pseudotime
#> Pyy                              1.5792657
#> Rbp4                             1.2535971
#> Sst                             -0.1667431
#>      AATCGGTTCGTATCAGstochastic_pseudotime
#> Pyy                              1.5853117
#> Rbp4                             1.2614797
#> Sst                             -0.1625028
#>      GACCAATCAGACGCCTstochastic_pseudotime
#> Pyy                              1.5871821
#> Rbp4                             1.2639506
#> Sst                             -0.1611471
#>      CCTCAGTCAGACACTTstochastic_pseudotime
#> Pyy                              1.5876116
#> Rbp4                             1.2645202
#> Sst                             -0.1608328
#>      GCGAGAAGTGGGTATGstochastic_pseudotime
#> Pyy                              1.5891066
#> Rbp4                             1.2665096
#> Sst                             -0.1597296
#>      ACAGCCGAGAGGTACCstochastic_pseudotime
#> Pyy                              1.5893045
#> Rbp4                             1.2667737
#> Sst                             -0.1595825
#>      ACTTGTTAGTTAAGTGstochastic_pseudotime
#> Pyy                              1.5918071
#> Rbp4                             1.2701297
#> Sst                             -0.1577006
#>      CTGAAGTGTTTGGGCCstochastic_pseudotime
#> Pyy                              1.5983987
#> Rbp4                             1.2791156
#> Sst                             -0.1525403
#>      CACTCCACATCGTCGGstochastic_pseudotime
#> Pyy                              1.6007134
#> Rbp4                             1.2823241
#> Sst                             -0.1506541
#>      CCACCTAAGACTTGAAstochastic_pseudotime
#> Pyy                              1.6066155
#> Rbp4                             1.2906381
#> Sst                             -0.1456561
#>      GTCTTCGTCCGTTGCTstochastic_pseudotime
#> Pyy                               1.611746
#> Rbp4                              1.298030
#> Sst                              -0.141075
#>      TACACGACATGTAGTCstochastic_pseudotime
#> Pyy                              1.6143522
#> Rbp4                             1.3018480
#> Sst                             -0.1386562
#>      GTGGGTCGTCCAGTATstochastic_pseudotime
#> Pyy                              1.6180482
#> Rbp4                             1.3073391
#> Sst                             -0.1351137
#>      GCGGGTTCAAGTAATGstochastic_pseudotime
#> Pyy                              1.6191397
#> Rbp4                             1.3089786
#> Sst                             -0.1340411
#>      CAGAGAGTCGGCCGATstochastic_pseudotime
#> Pyy                              1.6265573
#> Rbp4                             1.3203503
#> Sst                             -0.1264076
#>      TGCGGGTAGGGTTCCCstochastic_pseudotime
#> Pyy                              1.6293523
#> Rbp4                             1.3247462
#> Sst                             -0.1233634
#>      CATATTCTCCTGTAGAstochastic_pseudotime
#> Pyy                              1.6295314
#> Rbp4                             1.3250302
#> Sst                             -0.1231649
#>      CACAAACCACGAAAGCstochastic_pseudotime
#> Pyy                              1.6309030
#> Rbp4                             1.3272133
#> Sst                             -0.1216315
#>      TGCGTGGGTTACTGACstochastic_pseudotime
#> Pyy                              1.6310472
#> Rbp4                             1.3274437
#> Sst                             -0.1214689
#>      TCTTTCCGTGCAACTTstochastic_pseudotime
#> Pyy                              1.6310715
#> Rbp4                             1.3274825
#> Sst                             -0.1214415
#>      CCACTACCAACGCACCstochastic_pseudotime
#> Pyy                              1.6312960
#> Rbp4                             1.3278417
#> Sst                             -0.1211876
#>      TACCTTATCCGATATGstochastic_pseudotime
#> Pyy                              1.6331337
#> Rbp4                             1.3307983
#> Sst                             -0.1190843
#>      TGCGGGTTCCACTCCAstochastic_pseudotime
#> Pyy                              1.6361098
#> Rbp4                             1.3356512
#> Sst                             -0.1155774
#>      GCTTCCAAGACATAACstochastic_pseudotime
#> Pyy                              1.6369054
#> Rbp4                             1.3369626
#> Sst                             -0.1146181
#>      CAGTCCTAGGTGGGTTstochastic_pseudotime
#> Pyy                              1.6389281
#> Rbp4                             1.3403242
#> Sst                             -0.1121356
#>      CGCGGTAAGAGCTGGTstochastic_pseudotime
#> Pyy                               1.643817
#> Rbp4                              1.348622
#> Sst                              -0.105862
#>      GAGCAGAGTTTGCATGstochastic_pseudotime
#> Pyy                             1.64814271
#> Rbp4                            1.35618388
#> Sst                            -0.09995902
#>      ACGCCAGGTTACAGAAstochastic_pseudotime
#> Pyy                             1.65192992
#> Rbp4                            1.36299114
#> Sst                            -0.09448753
#>      GAAGCAGGTCTAGAGGstochastic_pseudotime
#> Pyy                             1.65310339
#> Rbp4                            1.36513834
#> Sst                            -0.09272992
#>      GTTCTCGCAAGAAAGGstochastic_pseudotime
#> Pyy                             1.65396179
#> Rbp4                            1.36672092
#> Sst                            -0.09142457
#>      ACTGATGCATGGTCTAstochastic_pseudotime
#> Pyy                             1.65514748
#> Rbp4                            1.36892391
#> Sst                            -0.08959335
#>      TACAGTGCACAGGCCTstochastic_pseudotime
#> Pyy                             1.65802448
#> Rbp4                            1.37435473
#> Sst                            -0.08500768
#>      TCTCATACAAATTGCCstochastic_pseudotime
#> Pyy                             1.65934303
#> Rbp4                            1.37688606
#> Sst                            -0.08283503
#>      TCCACACCAGGCGATAstochastic_pseudotime
#> Pyy                              1.6640296
#> Rbp4                             1.3861168
#> Sst                             -0.0747171
#>      GATTCAGAGTACGTTCstochastic_pseudotime
#> Pyy                             1.66469304
#> Rbp4                            1.38745482
#> Sst                            -0.07351438
#>      AGTCTTTGTTACGGAGstochastic_pseudotime
#> Pyy                             1.67411473
#> Rbp4                            1.40743485
#> Sst                            -0.05473337
#>      TTTGTCATCGAATGCTstochastic_pseudotime
#> Pyy                             1.67614088
#> Rbp4                            1.41200493
#> Sst                            -0.05021173
#>      GCGCAACGTTTCGCTCstochastic_pseudotime
#> Pyy                             1.67627655
#> Rbp4                            1.41231489
#> Sst                            -0.04990189
#>      ACGCCGAGTTTGCATGstochastic_pseudotime
#> Pyy                             1.68485475
#> Rbp4                            1.43308288
#> Sst                            -0.02818467
#>      CGTGAGCTCTTCTGGCstochastic_pseudotime
#> Pyy                             1.68715724
#> Rbp4                            1.43911959
#> Sst                            -0.02150204
#>      CATATGGGTTACTGACstochastic_pseudotime
#> Pyy                             1.68742678
#> Rbp4                            1.43984151
#> Sst                            -0.02069121
#>      AGGGTGAGTGCGAAACstochastic_pseudotime
#> Pyy                             1.68842071
#> Rbp4                            1.44253300
#> Sst                            -0.01764588
#>      AAACGGGTCGCATGGCstochastic_pseudotime
#> Pyy                             1.68966926
#> Rbp4                            1.44598270
#> Sst                            -0.01369059
#>      GACGTGCTCCGGCACAstochastic_pseudotime
#> Pyy                             1.69029924
#> Rbp4                            1.44775400
#> Sst                            -0.01163664
#>      TGTCCCATCTCTGCTGstochastic_pseudotime
#> Pyy                            1.690981207
#> Rbp4                           1.449695848
#> Sst                           -0.009366756
#>      TTTATGCAGGTTACCTstochastic_pseudotime
#> Pyy                             1.69688024
#> Rbp4                            1.46773701
#> Sst                             0.01266506
#>      ATCATGGGTATTCGTGstochastic_pseudotime
#> Pyy                             1.69688809
#> Rbp4                            1.46776276
#> Sst                             0.01269777
#>      GGACGTCGTCGAGATGstochastic_pseudotime
#> Pyy                              1.6994528
#> Rbp4                             1.4764806
#> Sst                              0.0239858
#>      CGACTTCTCACGCGGTstochastic_pseudotime
#> Pyy                             1.70265786
#> Rbp4                            1.48841088
#> Sst                             0.04015539
#>      TTTGTCATCTGTTTGTstochastic_pseudotime
#> Pyy                             1.70369230
#> Rbp4                            1.49257206
#> Sst                             0.04599996
#>      CGAGAAGAGGGTTCCCstochastic_pseudotime
#> Pyy                             1.70436346
#> Rbp4                            1.49536998
#> Sst                             0.04999127
#>      CACAGGCCACGTCAGCstochastic_pseudotime
#> Pyy                             1.70438296
#> Rbp4                            1.49545254
#> Sst                             0.05010982
#>      AAATGCCGTGATAAGTstochastic_pseudotime
#> Pyy                             1.70561237
#> Rbp4                            1.50081142
#> Sst                             0.05789826
#>      TTATGCTCACCCAGTGstochastic_pseudotime
#> Pyy                             1.70674291
#> Rbp4                            1.50604421
#> Sst                             0.06568644
#>      CCGTGGATCTGTCTCGstochastic_pseudotime
#> Pyy                             1.70954118
#> Rbp4                            1.52071223
#> Sst                             0.08852435
#>      TTAACTCCAGGATTGGstochastic_pseudotime
#> Pyy                             1.71059337
#> Rbp4                            1.52713561
#> Sst                             0.09901454
#>      TGAGAGGAGGTAAACTstochastic_pseudotime
#> Pyy                              1.7113843
#> Rbp4                             1.5324543
#> Sst                              0.1079367
#>      TTGAACGAGAGCAATTstochastic_pseudotime
#> Pyy                              1.7134773
#> Rbp4                             1.5500363
#> Sst                              0.1390318
#>      TGAGAGGAGCCAGGATstochastic_pseudotime
#> Pyy                              1.7143253
#> Rbp4                             1.5602834
#> Sst                              0.1583564
#>      AGTGAGGTCTCCAACCstochastic_pseudotime
#> Pyy                              1.7143915
#> Rbp4                             1.5612553
#> Sst                              0.1602376
#>      GGGAGATAGCGAAGGGstochastic_pseudotime
#> Pyy                              1.7144723
#> Rbp4                             1.5624965
#> Sst                              0.1626524
#>      CTAAGACGTATCAGTCstochastic_pseudotime
#> Pyy                              1.7146310
#> Rbp4                             1.5651381
#> Sst                              0.1678388
#>      CGGACTGTCTTGAGACstochastic_pseudotime
#> Pyy                              1.7147792
#> Rbp4                             1.5679410
#> Sst                              0.1734123
#>      GACGCGTGTCGGCTCAstochastic_pseudotime
#> Pyy                              1.7151119
#> Rbp4                             1.5768181
#> Sst                              0.1915557
#>      CCCAGTTGTCGTGGCTstochastic_pseudotime
#> Pyy                              1.7151646
#> Rbp4                             1.5791399
#> Sst                              0.1964274
#>      GTGCATATCAAGGCTTstochastic_pseudotime
#> Pyy                              1.7152365
#> Rbp4                             1.5858371
#> Sst                              0.2107826
#>      CTCGAGGGTAAGTAGTstochastic_pseudotime
#> Pyy                              1.7151552
#> Rbp4                             1.5943961
#> Sst                              0.2298035
#>      GGGAATGAGCACACAGstochastic_pseudotime
#> Pyy                              1.7132824
#> Rbp4                             1.6248783
#> Sst                              0.3042182
#>      TCGGGACGTCTGCGGTstochastic_pseudotime
#> Pyy                              1.7129991
#> Rbp4                             1.6275610
#> Sst                              0.3113038
#>      GCACTCTCAGCGTCCAstochastic_pseudotime
#> Pyy                              1.7125254
#> Rbp4                             1.6316919
#> Sst                              0.3223961
#>      CTCGTCAAGAACAATCstochastic_pseudotime
#> Pyy                              1.7119695
#> Rbp4                             1.6361031
#> Sst                              0.3344881
#>      AGCTTGAAGAGTGACCstochastic_pseudotime
#> Pyy                              1.7118262
#> Rbp4                             1.6371781
#> Sst                              0.3374741
#>      CGAATGTTCACCGGGTstochastic_pseudotime
#> Pyy                              1.7107961
#> Rbp4                             1.6443211
#> Sst                              0.3577164
#>      CCTACCATCTGCCCTAstochastic_pseudotime
#> Pyy                              1.7100841
#> Rbp4                             1.6487827
#> Sst                              0.3707207
#>      TTCTTAGTCAGTGTTGstochastic_pseudotime
#> Pyy                              1.7076439
#> Rbp4                             1.6621144
#> Sst                              0.4113065
#>      CTGGTCTAGGTGATTAstochastic_pseudotime
#> Pyy                              1.7049922
#> Rbp4                             1.6743616
#> Sst                              0.4509914
#>      GACCAATGTTGGACCCstochastic_pseudotime
#> Pyy                              1.7037655
#> Rbp4                             1.6794930
#> Sst                              0.4683391
#>      TACTCGCGTCGGCATCstochastic_pseudotime
#> Pyy                              1.7012231
#> Rbp4                             1.6893368
#> Sst                              0.5028628
#>      TACTCGCCAGATCGGAstochastic_pseudotime
#> Pyy                              1.6998672
#> Rbp4                             1.6942313
#> Sst                              0.5206586
#>      TGTCCCACATTATCTCstochastic_pseudotime
#> Pyy                              1.6997918
#> Rbp4                             1.6944972
#> Sst                              0.5216375
#>      CAGCTGGGTTCTCATTstochastic_pseudotime
#> Pyy                              1.6990547
#> Rbp4                             1.6970648
#> Sst                              0.5311576
#>      AAAGATGAGTCTTGCAstochastic_pseudotime
#> Pyy                              1.6986680
#> Rbp4                             1.6983892
#> Sst                              0.5361148
#>      CGCGTTTCAGTCTTCCstochastic_pseudotime
#> Pyy                              1.6966697
#> Rbp4                             1.7050093
#> Sst                              0.5613833
#>      GAAACTCAGAAACCTAstochastic_pseudotime
#> Pyy                               1.687360
#> Rbp4                              1.732167
#> Sst                               0.674044
#>      GGAATAAGTAAGAGAGstochastic_pseudotime
#> Pyy                              1.6791145
#> Rbp4                             1.7528954
#> Sst                              0.7706103
#>      TCTATTGTCCTATTCAstochastic_pseudotime
#> Pyy                              1.6687777
#> Rbp4                             1.7762059
#> Sst                              0.8912828
#>      TATCTCACATTGTGCAstochastic_pseudotime
#> Pyy                               1.655937
#> Rbp4                              1.802559
#> Sst                               1.044444
#>      TGGCTGGGTACCCAATstochastic_pseudotime
#> Pyy                               1.653386
#> Rbp4                              1.807560
#> Sst                               1.075621
#>      AGCATACCACAGCCCAstochastic_pseudotime
#> Pyy                               1.648753
#> Rbp4                              1.816494
#> Sst                               1.133050
#>      TTGGAACGTATATGAGstochastic_pseudotime
#> Pyy                               1.647973
#> Rbp4                              1.817984
#> Sst                               1.142843
#>      CGTCACTTCTAACTCTstochastic_pseudotime
#> Pyy                               1.640532
#> Rbp4                              1.831963
#> Sst                               1.237862
#>      CGTCCATGTCTAACGTstochastic_pseudotime
#> Pyy                               1.638600
#> Rbp4                              1.835529
#> Sst                               1.263028
#>      AGGGAGTTCCTCAATTstochastic_pseudotime
#> Pyy                               1.638020
#> Rbp4                              1.836593
#> Sst                               1.270611
#>      ACTATCTGTCAGGACAstochastic_pseudotime
#> Pyy                               1.633559
#> Rbp4                              1.844708
#> Sst                               1.329603
#>      TCACAAGAGTGACATAstochastic_pseudotime
#> Pyy                               1.624327
#> Rbp4                              1.861055
#> Sst                               1.454907
#>      GCATGTAAGCGTGTCCstochastic_pseudotime
#> Pyy                               1.623385
#> Rbp4                              1.862689
#> Sst                               1.467931
#>      GGACAGAAGGCCGAATstochastic_pseudotime
#> Pyy                               1.622070
#> Rbp4                              1.864959
#> Sst                               1.486173
#>      TCAGGTATCTGTCCGTstochastic_pseudotime
#> Pyy                               1.620911
#> Rbp4                              1.866951
#> Sst                               1.502324
#>      CTTGGCTCATCCCACTstochastic_pseudotime
#> Pyy                               1.620212
#> Rbp4                              1.868148
#> Sst                               1.512098
#>      GTCACGGAGCTATGCTstochastic_pseudotime
#> Pyy                               1.616083
#> Rbp4                              1.875145
#> Sst                               1.570265
#>      GTTTCTACAGGACCCTstochastic_pseudotime
#> Pyy                               1.614849
#> Rbp4                              1.877214
#> Sst                               1.587809
#>      GATCTAGAGCTAGCCCstochastic_pseudotime
#> Pyy                               1.607177
#> Rbp4                              1.889835
#> Sst                               1.698305
#>      TGTGGTAGTGCCTGCAstochastic_pseudotime
#> Pyy                               1.597856
#> Rbp4                              1.904619
#> Sst                               1.835772
#>      ATCTGCCCAATCCGATstochastic_pseudotime
#> Pyy                               1.594328
#> Rbp4                              1.910058
#> Sst                               1.888648
#>      GTGGGTCCATAAAGGTstochastic_pseudotime
#> Pyy                               1.593522
#> Rbp4                              1.911288
#> Sst                               1.900783
#>      AGAGCGACAGTAAGATstochastic_pseudotime
#> Pyy                               1.592834
#> Rbp4                              1.912335
#> Sst                               1.911163
#>      CGGACGTAGACACGACstochastic_pseudotime
#> Pyy                               1.586759
#> Rbp4                              1.921438
#> Sst                               2.003493
#>      AAAGATGAGTCATGCTstochastic_pseudotime
#> Pyy                               1.585908
#> Rbp4                              1.922692
#> Sst                               2.016514
#>      GATTCAGTCGTCCAGGstochastic_pseudotime
#> Pyy                               1.578043
#> Rbp4                              1.934054
#> Sst                               2.137861
#>      TACGGATAGGCACATGstochastic_pseudotime
#> Pyy                               1.575070
#> Rbp4                              1.938239
#> Sst                               2.184153
#>      AGGGATGTCGCATGATstochastic_pseudotime
#> Pyy                               1.572416
#> Rbp4                              1.941924
#> Sst                               2.225649
#>      GTATTCTTCTCTAAGGstochastic_pseudotime
#> Pyy                               1.570518
#> Rbp4                              1.944529
#> Sst                               2.255420
#>      CATCGAAGTCACACGCstochastic_pseudotime
#> Pyy                               1.569265
#> Rbp4                              1.946236
#> Sst                               2.275110
#>      GTCGTAAAGTTTAGGAstochastic_pseudotime
#> Pyy                               1.560980
#> Rbp4                              1.957254
#> Sst                               2.406047
#>      CTTAGGATCAACGGGAstochastic_pseudotime
#> Pyy                               1.553111
#> Rbp4                              1.967289
#> Sst                               2.531333
#>      GCCTCTAAGCTAGTCTstochastic_pseudotime
#> Pyy                               1.552743
#> Rbp4                              1.967748
#> Sst                               2.537207
#>      GTTCATTTCCGAACGCstochastic_pseudotime
#> Pyy                               1.552593
#> Rbp4                              1.967935
#> Sst                               2.539603
#>      CTAAGACTCTGCGGCAstochastic_pseudotime
#> Pyy                               1.538433
#> Rbp4                              1.984897
#> Sst                               2.766252
#>      AGCTCTCGTTCAGTACstochastic_pseudotime
#> Pyy                               1.534320
#> Rbp4                              1.989574
#> Sst                               2.832126
#>      GTACTTTGTACCGTATstochastic_pseudotime
#> Pyy                               1.533019
#> Rbp4                              1.991029
#> Sst                               2.852952
#>      CTGAAACGTATGGTTCstochastic_pseudotime
#> Pyy                               1.532416
#> Rbp4                              1.991700
#> Sst                               2.862593
#>      TGAGAGGAGGCGATACstochastic_pseudotime
#> Pyy                               1.528719
#> Rbp4                              1.995761
#> Sst                               2.921707
#>      ATCATCTGTGTTTGTGstochastic_pseudotime
#> Pyy                               1.528198
#> Rbp4                              1.996326
#> Sst                               2.930030
#>      TTGGCAAAGAGTAAGGstochastic_pseudotime
#> Pyy                               1.526815
#> Rbp4                              1.997817
#> Sst                               2.952106
#>      CTAGTGATCCCAGGTGstochastic_pseudotime
#> Pyy                               1.518329
#> Rbp4                              2.006690
#> Sst                               3.087099
#>      CCTAGCTGTGGCCCTAstochastic_pseudotime
#> Pyy                               1.511517
#> Rbp4                              2.013470
#> Sst                               3.194653
#>      CCTTCGAAGATAGCATstochastic_pseudotime
#> Pyy                               1.511186
#> Rbp4                              2.013791
#> Sst                               3.199853
#>      GCATGATCAAGCGAGTstochastic_pseudotime
#> Pyy                               1.507540
#> Rbp4                              2.017288
#> Sst                               3.257004
#>      AGCTTGAGTGCAGTAGstochastic_pseudotime
#> Pyy                               1.498850
#> Rbp4                              2.025273
#> Sst                               3.391875
#>      CGTCAGGTCCGCAAGCstochastic_pseudotime
#> Pyy                               1.495476
#> Rbp4                              2.028241
#> Sst                               3.443632
#>      CCGGTAGTCGGCTTGGstochastic_pseudotime
#> Pyy                               1.493247
#> Rbp4                              2.030161
#> Sst                               3.477624
#>      CAAGGCCTCTTACCTAstochastic_pseudotime
#> Pyy                               1.492230
#> Rbp4                              2.031027
#> Sst                               3.493089
#>      TTTCCTCCAGGATCGAstochastic_pseudotime
#> Pyy                               1.492134
#> Rbp4                              2.031108
#> Sst                               3.494537
#>      CAGATCAAGCTGTTCAstochastic_pseudotime
#> Pyy                               1.489144
#> Rbp4                              2.033613
#> Sst                               3.539747
#>      TCAGATGGTACACCGCstochastic_pseudotime
#> Pyy                               1.487699
#> Rbp4                              2.034803
#> Sst                               3.561478
#>      CATCGAAGTACGAAATstochastic_pseudotime
#> Pyy                               1.480389
#> Rbp4                              2.040616
#> Sst                               3.670105
#>      TACGGGCTCTTCGGTCstochastic_pseudotime
#> Pyy                               1.473823
#> Rbp4                              2.045547
#> Sst                               3.765681
#>      TGATTTCGTGAGGGTTstochastic_pseudotime
#> Pyy                               1.468746
#> Rbp4                              2.049171
#> Sst                               3.838142
#>      CCTATTACAGACAAGCstochastic_pseudotime
#> Pyy                               1.465297
#> Rbp4                              2.051540
#> Sst                               3.886593
#>      GTCGGGTCAGGCTGAAstochastic_pseudotime
#> Pyy                               1.463718
#> Rbp4                              2.052599
#> Sst                               3.908562
#>      CCCAGTTAGCATCATCstochastic_pseudotime
#> Pyy                               1.447842
#> Rbp4                              2.062381
#> Sst                               4.121215
#>      GGCTCGAAGTGACTCTstochastic_pseudotime
#> Pyy                               1.439325
#> Rbp4                              2.066981
#> Sst                               4.228535
#>      TTCGGTCAGGGCATGTstochastic_pseudotime
#> Pyy                               1.437905
#> Rbp4                              2.067704
#> Sst                               4.245923
#>      AGTGGGACAGGTCGTCstochastic_pseudotime
#> Pyy                               1.437073
#> Rbp4                              2.068122
#> Sst                               4.256046
#>      GTTTCTAAGAAGATTCstochastic_pseudotime
#> Pyy                               1.433421
#> Rbp4                              2.069906
#> Sst                               4.299876
#>      CGACCTTAGGCTAGACstochastic_pseudotime
#> Pyy                               1.430100
#> Rbp4                              2.071457
#> Sst                               4.338859
#>      TGAGAGGAGATCCCATstochastic_pseudotime
#> Pyy                               1.423483
#> Rbp4                              2.074345
#> Sst                               4.413948
#>      CTAGTGATCAAAGACAstochastic_pseudotime
#> Pyy                               1.421475
#> Rbp4                              2.075168
#> Sst                               4.436035
#>      GGACGTCGTCTAGGTTstochastic_pseudotime
#> Pyy                               1.413549
#> Rbp4                              2.078180
#> Sst                               4.519944
#>      CACATTTTCGAGGTAGstochastic_pseudotime
#> Pyy                               1.411980
#> Rbp4                              2.078731
#> Sst                               4.535914
#>      TTCTACATCGCTTAGAstochastic_pseudotime
#> Pyy                               1.390011
#> Rbp4                              2.084892
#> Sst                               4.736148
#>      TCGCGAGAGTCCGTATstochastic_pseudotime
#> Pyy                               1.388370
#> Rbp4                              2.085237
#> Sst                               4.749279
#>      TACTCATCAGATCCATstochastic_pseudotime
#> Pyy                               1.351903
#> Rbp4                              2.088834
#> Sst                               4.970799
#>      CAACCAACAACCGCCAstochastic_pseudotime
#> Pyy                               1.336348
#> Rbp4                              2.088042
#> Sst                               5.022656
#>      ACGAGGATCATCACCCstochastic_pseudotime
#> Pyy                               1.313310
#> Rbp4                              2.084378
#> Sst                               5.051346
#>      GCATGTAAGAATCTCCstochastic_pseudotime
#> Pyy                               1.263218
#> Rbp4                              2.066474
#> Sst                               4.919715
#>      AGGCCGTCACGGTGTCstochastic_pseudotime
#> Pyy                               1.259998
#> Rbp4                              2.064870
#> Sst                               4.902621
#>      CACTCCATCAACACACstochastic_pseudotime
#> Pyy                               1.258896
#> Rbp4                              2.064309
#> Sst                               4.896545
#>      CTTGGCTGTTAGAACAstochastic_pseudotime
#> Pyy                               1.256533
#> Rbp4                              2.063086
#> Sst                               4.883136
#>      CCCATACAGACTTTCGstochastic_pseudotime
#> Pyy                               1.252455
#> Rbp4                              2.060906
#> Sst                               4.858764
#>      CCGTGGACATAGAAACstochastic_pseudotime
#> Pyy                               1.249693
#> Rbp4                              2.059382
#> Sst                               4.841394
#>      TCTCTAAAGCTGGAACstochastic_pseudotime
#> Pyy                               1.249125
#> Rbp4                              2.059064
#> Sst                               4.837737
#>      GTGCTTCTCACGACTAstochastic_pseudotime
#> Pyy                               1.248195
#> Rbp4                              2.058539
#> Sst                               4.831682
#>      GGCTCGACAGCTTCGGstochastic_pseudotime
#> Pyy                               1.165609
#> Rbp4                              1.995070
#> Sst                               4.031226
#>      GTCTCGTTCTGCGTAAstochastic_pseudotime
#> Pyy                               1.146359
#> Rbp4                              1.975705
#> Sst                               3.788556
#>      GCTCCTAAGGGCATGTstochastic_pseudotime
#> Pyy                               1.088781
#> Rbp4                              1.908286
#> Sst                               3.005526
#>      GTGGGTCCATCGATTGstochastic_pseudotime
#> Pyy                               1.029660
#> Rbp4                              1.825373
#> Sst                               2.197104
#>      CCACGGAGTCTAGTCAstochastic_pseudotime
#> Pyy                               1.016820
#> Rbp4                              1.805678
#> Sst                               2.030046
#>      CTAACTTAGATCACGGstochastic_pseudotime
#> Pyy                              0.9250899
#> Rbp4                             1.6494355
#> Sst                              1.0009838
#>      TGTCCCAGTCCGTGACstochastic_pseudotime
#> Pyy                              0.9064676
#> Rbp4                             1.6147066
#> Sst                              0.8335855
#>      ATCCACCAGAGCCTAGstochastic_pseudotime
#> Pyy                              0.8903864
#> Rbp4                             1.5839889
#> Sst                              0.7009207
#>      GTAGGCCGTTACCAGTstochastic_pseudotime
#> Pyy                              0.8667432
#> Rbp4                             1.5376729
#> Sst                              0.5255735
#>      TCTGGAACACCTATCCstochastic_pseudotime
#> Pyy                              0.8331566
#> Rbp4                             1.4696971
#> Sst                              0.3150752
#>      ACGTCAAAGGTGTTAAstochastic_pseudotime
#> Pyy                             0.74236057
#> Rbp4                            1.27525730
#> Sst                            -0.06018372
#>      ACTTGTTCATGGTCATstochastic_pseudotime
#> Pyy                             0.74105373
#> Rbp4                            1.27236414
#> Sst                            -0.06387338
#>      GGACAGATCAGGCGAAstochastic_pseudotime
#> Pyy                             0.73810972
#> Rbp4                            1.26583819
#> Sst                            -0.07203649
#>      TACTTGTAGACAAAGGstochastic_pseudotime
#> Pyy                             0.72927166
#> Rbp4                            1.24617939
#> Sst                            -0.09533638
#>      AAATGCCAGACTAAGTstochastic_pseudotime
#> Pyy                              0.6764264
#> Rbp4                             1.1267807
#> Sst                             -0.2017598
#>      AAGTCTGGTTGGGACAstochastic_pseudotime
#> Pyy                              0.6456256
#> Rbp4                             1.0560188
#> Sst                             -0.2428822
#>      AGACGTTCAATGGTCTstochastic_pseudotime
#> Pyy                              0.6349080
#> Rbp4                             1.0312452
#> Sst                             -0.2543869
#>      CTACACCCATCAGTCAstochastic_pseudotime
#> Pyy                              0.5914427
#> Rbp4                             0.9302165
#> Sst                             -0.2894085
#>      CGTTCTGTCGATAGAAstochastic_pseudotime
#> Pyy                              0.5647498
#> Rbp4                             0.8679066
#> Sst                             -0.3036981
#>      GACCAATAGAATGTGTstochastic_pseudotime
#> Pyy                              0.5583705
#> Rbp4                             0.8530040
#> Sst                             -0.3064857
#>      TGACTTTAGCCAGAACstochastic_pseudotime
#> Pyy                              0.5448563
#> Rbp4                             0.8214334
#> Sst                             -0.3117123
#>      TCACAAGCAGACTCGCstochastic_pseudotime
#> Pyy                              0.5080235
#> Rbp4                             0.7355054
#> Sst                             -0.3221252
#>      CTTTGCGTCCACGTGGstochastic_pseudotime
#> Pyy                              0.4658005
#> Rbp4                             0.6375516
#> Sst                             -0.3291709
#>      CTCTGGTCAAGACACGstochastic_pseudotime
#> Pyy                              0.4110703
#> Rbp4                             0.5122033
#> Sst                             -0.3337641
#>      AGCATACGTCCGCTGAstochastic_pseudotime
#> Pyy                              0.3103629
#> Rbp4                             0.2893772
#> Sst                             -0.3364267
#>      GCTGCGATCGGCTACGstochastic_pseudotime
#> Pyy                              0.3079523
#> Rbp4                             0.2842036
#> Sst                             -0.3364500
#>      CCACCTACAGCTCGACstochastic_pseudotime
#> Pyy                              0.2842866
#> Rbp4                             0.2338768
#> Sst                             -0.3366313
#>      GGCAATTGTAGAGCTGstochastic_pseudotime
#> Pyy                              0.2522962
#> Rbp4                             0.1672673
#> Sst                             -0.3367774
#>      ACTTACTCACATGTGTstochastic_pseudotime
#> Pyy                              0.2331133
#> Rbp4                             0.1281638
#> Sst                             -0.3368300
#>      ACGAGCCCATGAGCGAstochastic_pseudotime
#> Pyy                              0.2271599
#> Rbp4                             0.1161621
#> Sst                             -0.3368427
#>      CCCAGTTGTGCGCTTGstochastic_pseudotime
#> Pyy                              0.2252602
#> Rbp4                             0.1123461
#> Sst                             -0.3368464
#>      GACGTTACATCTCGCTstochastic_pseudotime
#> Pyy                             0.14641729
#> Rbp4                           -0.03984763
#> Sst                            -0.33692027
#>      AAACGGGGTACAGTTCstochastic_pseudotime
#> Pyy                             0.05300876
#> Rbp4                           -0.20292923
#> Sst                            -0.33693356
#>      GCGCGATGTCCGAATTstochastic_pseudotime
#> Pyy                             0.05045194
#> Rbp4                           -0.20710939
#> Sst                            -0.33693366
#>      TAGTTGGTCTTGCAAGstochastic_pseudotime
#> Pyy                            0.006087735
#> Rbp4                          -0.277119555
#> Sst                           -0.336934554
#>      TGGCCAGCAAGGACTGstochastic_pseudotime
#> Pyy                            0.005014737
#> Rbp4                          -0.278753034
#> Sst                           -0.336934565
#>      TGTGGTACACATCCGGstochastic_pseudotime
#> Pyy                            0.003724758
#> Rbp4                          -0.280713056
#> Sst                           -0.336934578
#> attr(,"scaled:center")
#>       Pyy      Rbp4       Sst 
#> 62.113066 19.747098  2.682374 
#> attr(,"scaled:scale")
#>       Pyy      Rbp4       Sst 
#> 77.691894 24.451344  7.961104 
#> 
#> $cell_order
#> $cell_order$stochastic_pseudotime
#>    [1] "TCAATCTGTCAAAGATstochastic_pseudotime"
#>    [2] "CTTGGCTAGACAGACCstochastic_pseudotime"
#>    [3] "GTATTCTCACTTGGATstochastic_pseudotime"
#>    [4] "AGATTGCAGCTGTCTAstochastic_pseudotime"
#>    [5] "CCATTCGCATGCCTAAstochastic_pseudotime"
#>    [6] "GTATTCTGTCCGTCAGstochastic_pseudotime"
#>    [7] "TGGACGCTCAGGTAAAstochastic_pseudotime"
#>    [8] "ATCATCTTCTGTCAAGstochastic_pseudotime"
#>    [9] "TTCTACATCTGCAGTAstochastic_pseudotime"
#>   [10] "TATCTCAGTTCGCTAAstochastic_pseudotime"
#>   [11] "CAGCCGAAGTAGGTGCstochastic_pseudotime"
#>   [12] "AGTCTTTCACCAGGTCstochastic_pseudotime"
#>   [13] "ATTACTCAGTTGCAGGstochastic_pseudotime"
#>   [14] "CAGCAGCCACACTGCGstochastic_pseudotime"
#>   [15] "AGTGAGGGTTACGGAGstochastic_pseudotime"
#>   [16] "TGAGCATTCGGGAGTAstochastic_pseudotime"
#>   [17] "AACTCAGGTTGGAGGTstochastic_pseudotime"
#>   [18] "GGACGTCGTTTAGCTGstochastic_pseudotime"
#>   [19] "CTAGCCTGTTCGGGCTstochastic_pseudotime"
#>   [20] "TAGTGGTCAACTGCGCstochastic_pseudotime"
#>   [21] "GAATGAACAGCCTTGGstochastic_pseudotime"
#>   [22] "TGGCCAGAGCTACCTAstochastic_pseudotime"
#>   [23] "ATTACTCTCAGGCGAAstochastic_pseudotime"
#>   [24] "GATGAAAAGTTGTAGAstochastic_pseudotime"
#>   [25] "CCGTACTTCAATAAGGstochastic_pseudotime"
#>   [26] "TTAACTCCAAAGAATCstochastic_pseudotime"
#>   [27] "TATGCCCTCTTACCTAstochastic_pseudotime"
#>   [28] "CCATTCGGTAGGGACTstochastic_pseudotime"
#>   [29] "TAGGCATTCCAAAGTCstochastic_pseudotime"
#>   [30] "CTACACCGTCGCCATGstochastic_pseudotime"
#>   [31] "AACTGGTGTGCCTGTGstochastic_pseudotime"
#>   [32] "TGCGCAGGTACGCACCstochastic_pseudotime"
#>   [33] "CTCATTACAAGCGCTCstochastic_pseudotime"
#>   [34] "TGACTTTAGCGTTGCCstochastic_pseudotime"
#>   [35] "GATCGTAGTCCTAGCGstochastic_pseudotime"
#>   [36] "GGACAGATCCTGCAGGstochastic_pseudotime"
#>   [37] "TTTATGCGTCGAACAGstochastic_pseudotime"
#>   [38] "CGAACATAGCCAGTAGstochastic_pseudotime"
#>   [39] "AAGCCGCAGGCAGGTTstochastic_pseudotime"
#>   [40] "TTTATGCTCCATTCTAstochastic_pseudotime"
#>   [41] "TAGAGCTGTTGTCGCGstochastic_pseudotime"
#>   [42] "TAGAGCTGTATGAAACstochastic_pseudotime"
#>   [43] "CGTGTAATCCATGCTCstochastic_pseudotime"
#>   [44] "TGACTAGTCCACGTGGstochastic_pseudotime"
#>   [45] "TGCGCAGTCTGCCAGGstochastic_pseudotime"
#>   [46] "TGAGAGGGTAAGGGCTstochastic_pseudotime"
#>   [47] "AACTCAGTCTGGTTCCstochastic_pseudotime"
#>   [48] "CTCTGGTAGTCCCACGstochastic_pseudotime"
#>   [49] "ACAGCTATCTGGGCCAstochastic_pseudotime"
#>   [50] "CGGTTAAGTAAACACAstochastic_pseudotime"
#>   [51] "AGGCCGTGTACCAGTTstochastic_pseudotime"
#>   [52] "CGCCAAGCACCGGAAAstochastic_pseudotime"
#>   [53] "CGGTTAAAGAGCTGCAstochastic_pseudotime"
#>   [54] "ATAAGAGCATCGATTGstochastic_pseudotime"
#>   [55] "GTAACTGTCAGTGCATstochastic_pseudotime"
#>   [56] "CCTAGCTTCGATGAGGstochastic_pseudotime"
#>   [57] "GATGCTAGTAGAAGGAstochastic_pseudotime"
#>   [58] "TGTGGTACAAGACGTGstochastic_pseudotime"
#>   [59] "TCGTACCCACGCATCGstochastic_pseudotime"
#>   [60] "CTGCGGATCAGCATGTstochastic_pseudotime"
#>   [61] "CAGCCGAAGCGATATAstochastic_pseudotime"
#>   [62] "CAGCAGCAGACAAAGGstochastic_pseudotime"
#>   [63] "TCTGAGAGTCAAACTCstochastic_pseudotime"
#>   [64] "CTAGAGTCACTGTTAGstochastic_pseudotime"
#>   [65] "ACGGGCTGTGCGATAGstochastic_pseudotime"
#>   [66] "ACCGTAAGTTATCACGstochastic_pseudotime"
#>   [67] "AAGACCTCAGGTTTCAstochastic_pseudotime"
#>   [68] "TGTATTCCAATGACCTstochastic_pseudotime"
#>   [69] "CCTAGCTCATACTACGstochastic_pseudotime"
#>   [70] "CACACTCAGATCCCATstochastic_pseudotime"
#>   [71] "ACTGATGGTGAGGGTTstochastic_pseudotime"
#>   [72] "CAGCTAAAGTCCATACstochastic_pseudotime"
#>   [73] "TGCTGCTGTCCGTCAGstochastic_pseudotime"
#>   [74] "CACATAGAGCGTGTCCstochastic_pseudotime"
#>   [75] "AGGGAGTAGGGTATCGstochastic_pseudotime"
#>   [76] "CAGCGACCACAGGAGTstochastic_pseudotime"
#>   [77] "AGAGCTTGTTAAAGACstochastic_pseudotime"
#>   [78] "GCTGCTTGTCTTCAAGstochastic_pseudotime"
#>   [79] "TGGCGCATCACAACGTstochastic_pseudotime"
#>   [80] "GTCAAGTGTGGTCCGTstochastic_pseudotime"
#>   [81] "ACGTCAACAGACAAATstochastic_pseudotime"
#>   [82] "GACTAACTCAGTCCCTstochastic_pseudotime"
#>   [83] "TTGACTTAGACAAGCCstochastic_pseudotime"
#>   [84] "CGATGGCTCTGTGCAAstochastic_pseudotime"
#>   [85] "GAGCAGAAGCGTGTCCstochastic_pseudotime"
#>   [86] "ACAGCCGAGTGGTCCCstochastic_pseudotime"
#>   [87] "AGAATAGCAGTACACTstochastic_pseudotime"
#>   [88] "AGAGCGATCTCTTGATstochastic_pseudotime"
#>   [89] "ATTATCCAGTGAACATstochastic_pseudotime"
#>   [90] "CACACTCTCATGTGGTstochastic_pseudotime"
#>   [91] "CGTTAGAAGTATCTCGstochastic_pseudotime"
#>   [92] "TCTATTGCACAGGAGTstochastic_pseudotime"
#>   [93] "CTTCTCTAGGCAGGTTstochastic_pseudotime"
#>   [94] "CCTACACTCCTTTACAstochastic_pseudotime"
#>   [95] "CTGCCTAGTCTCCATCstochastic_pseudotime"
#>   [96] "GGAAAGCGTCTTGATGstochastic_pseudotime"
#>   [97] "CCGTGGATCCAAGCCGstochastic_pseudotime"
#>   [98] "ATCCACCGTACCATCAstochastic_pseudotime"
#>   [99] "AGCTTGACAGTAAGATstochastic_pseudotime"
#>  [100] "GTGCATAAGCTCCTCTstochastic_pseudotime"
#>  [101] "CAACCAATCGGGAGTAstochastic_pseudotime"
#>  [102] "CCTCTGAGTATAAACGstochastic_pseudotime"
#>  [103] "ACAGCTATCAGCATGTstochastic_pseudotime"
#>  [104] "GTGTTAGCACCATCCTstochastic_pseudotime"
#>  [105] "CATCAGACATGCATGTstochastic_pseudotime"
#>  [106] "CACACAACAGCCAGAAstochastic_pseudotime"
#>  [107] "GGGCATCCATGACATCstochastic_pseudotime"
#>  [108] "GGAATAAGTTCGGGCTstochastic_pseudotime"
#>  [109] "CATATGGCATGTTGACstochastic_pseudotime"
#>  [110] "CTACACCGTCGAATCTstochastic_pseudotime"
#>  [111] "CCCATACCACATCCAAstochastic_pseudotime"
#>  [112] "GGAGCAAAGCTTCGCGstochastic_pseudotime"
#>  [113] "CTAGCCTGTCTCTTATstochastic_pseudotime"
#>  [114] "GAGTCCGCAATTGCTGstochastic_pseudotime"
#>  [115] "CGGAGCTCATTGGGCCstochastic_pseudotime"
#>  [116] "AGGGAGTAGTGTCTCAstochastic_pseudotime"
#>  [117] "AGCGGTCGTGTATGGGstochastic_pseudotime"
#>  [118] "GCGCAACTCACAACGTstochastic_pseudotime"
#>  [119] "TTCTCAAAGCCGCCTAstochastic_pseudotime"
#>  [120] "GAACATCGTAATCGTCstochastic_pseudotime"
#>  [121] "CGAGCACAGCTGAAATstochastic_pseudotime"
#>  [122] "TGGTTCCGTCACCTAAstochastic_pseudotime"
#>  [123] "CGTCTACAGATCTGCTstochastic_pseudotime"
#>  [124] "ACGTCAACACTGCCAGstochastic_pseudotime"
#>  [125] "CTGAAACTCGCCGTGAstochastic_pseudotime"
#>  [126] "CTACCCAGTGCGCTTGstochastic_pseudotime"
#>  [127] "CTTAGGACACAGGTTTstochastic_pseudotime"
#>  [128] "CGATGGCGTTACGGAGstochastic_pseudotime"
#>  [129] "TACGGATAGAGCTGGTstochastic_pseudotime"
#>  [130] "TGGGAAGCATGACATCstochastic_pseudotime"
#>  [131] "TAGCCGGTCAACACACstochastic_pseudotime"
#>  [132] "TGACGGCAGGAATGGAstochastic_pseudotime"
#>  [133] "CGGACACAGTGTCTCAstochastic_pseudotime"
#>  [134] "TCTCTAAGTTACCGATstochastic_pseudotime"
#>  [135] "TGGCCAGGTCGGCATCstochastic_pseudotime"
#>  [136] "ACTGATGTCCCTAACCstochastic_pseudotime"
#>  [137] "TTCCCAGGTATTCTCTstochastic_pseudotime"
#>  [138] "CCTAGCTCAGACTCGCstochastic_pseudotime"
#>  [139] "AAGGCAGAGCTGAAATstochastic_pseudotime"
#>  [140] "CGGGTCAAGGCTCAGAstochastic_pseudotime"
#>  [141] "ATCATGGTCCAGAGGAstochastic_pseudotime"
#>  [142] "CTGATCCTCAGTTGACstochastic_pseudotime"
#>  [143] "AGATTGCCACCAGGCTstochastic_pseudotime"
#>  [144] "TCTGGAACACACGCTGstochastic_pseudotime"
#>  [145] "CAGCAGCGTTCACGGCstochastic_pseudotime"
#>  [146] "TTCTCCTCAATGGTCTstochastic_pseudotime"
#>  [147] "TTGACTTTCCAGTAGTstochastic_pseudotime"
#>  [148] "GGGCACTGTGGCAAACstochastic_pseudotime"
#>  [149] "CGAGCACTCTACCTGCstochastic_pseudotime"
#>  [150] "TTGACTTTCAGCAACTstochastic_pseudotime"
#>  [151] "GGTATTGGTCCCGACAstochastic_pseudotime"
#>  [152] "ACACCGGCATGTCCTCstochastic_pseudotime"
#>  [153] "CTGAAACAGACTAGGCstochastic_pseudotime"
#>  [154] "GGACAAGGTGATAAGTstochastic_pseudotime"
#>  [155] "TACCTTAAGTACATGAstochastic_pseudotime"
#>  [156] "GGGAATGAGTGTACGGstochastic_pseudotime"
#>  [157] "GAAATGACAACTGGCCstochastic_pseudotime"
#>  [158] "TTTGCGCAGACAGACCstochastic_pseudotime"
#>  [159] "CACATTTGTTTCGCTCstochastic_pseudotime"
#>  [160] "TTTGGTTCACCAGATTstochastic_pseudotime"
#>  [161] "TAAGTGCAGGGTTTCTstochastic_pseudotime"
#>  [162] "TGCTACCCAGGGCATAstochastic_pseudotime"
#>  [163] "ATTGGTGTCCAGTATGstochastic_pseudotime"
#>  [164] "ATTTCTGGTAGTGAATstochastic_pseudotime"
#>  [165] "CTAGCCTTCAGCGATTstochastic_pseudotime"
#>  [166] "CACCTTGAGAGTAATCstochastic_pseudotime"
#>  [167] "GAAATGACACATGACTstochastic_pseudotime"
#>  [168] "CGGACTGAGGTAGCTGstochastic_pseudotime"
#>  [169] "CTCAGAAGTCCTCTTGstochastic_pseudotime"
#>  [170] "AACTGGTTCTTACCGCstochastic_pseudotime"
#>  [171] "CTCGTCACAACTGCTAstochastic_pseudotime"
#>  [172] "CACAGTAGTTCCTCCAstochastic_pseudotime"
#>  [173] "GGGTTGCCAAGCTGAGstochastic_pseudotime"
#>  [174] "TACGGTACATGTTCCCstochastic_pseudotime"
#>  [175] "TGAAAGAGTTATCACGstochastic_pseudotime"
#>  [176] "TACAGTGCAATCTACGstochastic_pseudotime"
#>  [177] "GGACGTCCACGCATCGstochastic_pseudotime"
#>  [178] "CTAACTTCATTGGGCCstochastic_pseudotime"
#>  [179] "AAACCTGAGCCTTGATstochastic_pseudotime"
#>  [180] "CCTTCCCAGGCAGGTTstochastic_pseudotime"
#>  [181] "TCACGAATCCTATGTTstochastic_pseudotime"
#>  [182] "AAAGCAACACATGGGAstochastic_pseudotime"
#>  [183] "GAAGCAGTCAAGGCTTstochastic_pseudotime"
#>  [184] "AGACGTTTCTTACCTAstochastic_pseudotime"
#>  [185] "TAGACCAAGGAATCGCstochastic_pseudotime"
#>  [186] "TGAGCATTCACCCGAGstochastic_pseudotime"
#>  [187] "GGGACCTTCGAATGCTstochastic_pseudotime"
#>  [188] "GTACTCCCACGTCTCTstochastic_pseudotime"
#>  [189] "ACTGAACAGCCTCGTGstochastic_pseudotime"
#>  [190] "CTTAACTCAAGCGTAGstochastic_pseudotime"
#>  [191] "ACACCGGCATCACGTAstochastic_pseudotime"
#>  [192] "CAAGAAATCGTACCGGstochastic_pseudotime"
#>  [193] "CGAACATCAGACGTAGstochastic_pseudotime"
#>  [194] "ACGGCCACACGTGAGAstochastic_pseudotime"
#>  [195] "GAATGAACAATCGAAAstochastic_pseudotime"
#>  [196] "CTAAGACAGCCACGCTstochastic_pseudotime"
#>  [197] "TAGTTGGGTACCGTATstochastic_pseudotime"
#>  [198] "CCTTCGACAGCTGTTAstochastic_pseudotime"
#>  [199] "AGCAGCCTCTTTAGTCstochastic_pseudotime"
#>  [200] "CGTCACTCATAGTAAGstochastic_pseudotime"
#>  [201] "TACTTACCAACTGCGCstochastic_pseudotime"
#>  [202] "CATGACAGTAGCTTGTstochastic_pseudotime"
#>  [203] "CCTTCCCCAAAGTGCGstochastic_pseudotime"
#>  [204] "GGGTTGCTCGAACGGAstochastic_pseudotime"
#>  [205] "CTGGTCTCAGGAATGCstochastic_pseudotime"
#>  [206] "TGAGGGACAGACTCGCstochastic_pseudotime"
#>  [207] "GGTGTTAGTTCCTCCAstochastic_pseudotime"
#>  [208] "TGGTTCCCACCGCTAGstochastic_pseudotime"
#>  [209] "GACACGCCAAAGGAAGstochastic_pseudotime"
#>  [210] "CACAGTATCCTCCTAGstochastic_pseudotime"
#>  [211] "CGTAGGCAGCGTAGTGstochastic_pseudotime"
#>  [212] "AACGTTGGTGCCTTGGstochastic_pseudotime"
#>  [213] "CAAGGCCAGTGTGGCAstochastic_pseudotime"
#>  [214] "CTCTGGTCAGTTCATGstochastic_pseudotime"
#>  [215] "TATCTCACAACTGCTAstochastic_pseudotime"
#>  [216] "GCTGCAGGTTACCGATstochastic_pseudotime"
#>  [217] "TACGGATAGTTTAGGAstochastic_pseudotime"
#>  [218] "ACTTTCAAGGATATACstochastic_pseudotime"
#>  [219] "ATAGACCTCTAACTGGstochastic_pseudotime"
#>  [220] "GACTACAGTTCGTGATstochastic_pseudotime"
#>  [221] "AAACGGGAGATATGGTstochastic_pseudotime"
#>  [222] "CTCGTCATCGGTCTAAstochastic_pseudotime"
#>  [223] "TCTCTAAGTCGCATCGstochastic_pseudotime"
#>  [224] "GAAATGAGTCGGATCCstochastic_pseudotime"
#>  [225] "GACGGCTAGGTGATATstochastic_pseudotime"
#>  [226] "CGTCACTCACGAAATAstochastic_pseudotime"
#>  [227] "TATCAGGCACCAGATTstochastic_pseudotime"
#>  [228] "TGACTTTAGCGTTTACstochastic_pseudotime"
#>  [229] "GGGCACTCATTTCAGGstochastic_pseudotime"
#>  [230] "GAGTCCGAGTTAAGTGstochastic_pseudotime"
#>  [231] "GTTCTCGCAGACACTTstochastic_pseudotime"
#>  [232] "GTATTCTCACATTAGCstochastic_pseudotime"
#>  [233] "AACTCCCTCGAGAGCAstochastic_pseudotime"
#>  [234] "CGTCAGGAGCCTATGTstochastic_pseudotime"
#>  [235] "GTTCATTAGAAGATTCstochastic_pseudotime"
#>  [236] "TCAGCAACACCAGGCTstochastic_pseudotime"
#>  [237] "CGTCAGGCATGGGAACstochastic_pseudotime"
#>  [238] "CCTACACAGATCCGAGstochastic_pseudotime"
#>  [239] "CAGAATCGTGTAAGTAstochastic_pseudotime"
#>  [240] "TGCGTGGAGGTAAACTstochastic_pseudotime"
#>  [241] "TGATTTCAGTGACTCTstochastic_pseudotime"
#>  [242] "TCAATCTCATCTCCCAstochastic_pseudotime"
#>  [243] "CCTTACGAGATCACGGstochastic_pseudotime"
#>  [244] "TGAGCCGTCATGTCCCstochastic_pseudotime"
#>  [245] "CCACCTATCACCAGGCstochastic_pseudotime"
#>  [246] "GCAATCATCCAACCAAstochastic_pseudotime"
#>  [247] "AGAGTGGGTAGCACGAstochastic_pseudotime"
#>  [248] "AGTAGTCGTTCTGTTTstochastic_pseudotime"
#>  [249] "GCGGGTTGTTACTGACstochastic_pseudotime"
#>  [250] "GCATGATTCAAACGGGstochastic_pseudotime"
#>  [251] "TCATTACGTCAGAAGCstochastic_pseudotime"
#>  [252] "GATCTAGGTGCCTTGGstochastic_pseudotime"
#>  [253] "AGTAGTCCACATGACTstochastic_pseudotime"
#>  [254] "CGTTAGAGTGCAACGAstochastic_pseudotime"
#>  [255] "CTCGTACCAATGGATAstochastic_pseudotime"
#>  [256] "TTGCCGTTCTCTGTCGstochastic_pseudotime"
#>  [257] "ACATACGAGTAGGCCAstochastic_pseudotime"
#>  [258] "TTTACTGGTGCTAGCCstochastic_pseudotime"
#>  [259] "AGCTTGACAATAGAGTstochastic_pseudotime"
#>  [260] "CTCAGAATCAGGTTCAstochastic_pseudotime"
#>  [261] "CGTCACTAGGGTGTGTstochastic_pseudotime"
#>  [262] "CTTGGCTAGGACACCAstochastic_pseudotime"
#>  [263] "TGTCCCAAGTTTAGGAstochastic_pseudotime"
#>  [264] "TTCTACAAGATCGATAstochastic_pseudotime"
#>  [265] "AGCTCCTGTGTGCCTGstochastic_pseudotime"
#>  [266] "GCCTCTATCGCTTAGAstochastic_pseudotime"
#>  [267] "GACGTTAAGACGCAACstochastic_pseudotime"
#>  [268] "GGGTCTGCATGAAGTAstochastic_pseudotime"
#>  [269] "CTTTGCGAGGAATGGAstochastic_pseudotime"
#>  [270] "TACGGGCAGCGTCAAGstochastic_pseudotime"
#>  [271] "CCGTTCACACATGACTstochastic_pseudotime"
#>  [272] "GACGCGTCAGTCGTGCstochastic_pseudotime"
#>  [273] "ACTTACTAGCCGGTAAstochastic_pseudotime"
#>  [274] "GATGAAACATGGTTGTstochastic_pseudotime"
#>  [275] "AGGCCACAGCTGTCTAstochastic_pseudotime"
#>  [276] "CACATTTTCGCAGGCTstochastic_pseudotime"
#>  [277] "CGAGCACTCTTCGAGAstochastic_pseudotime"
#>  [278] "CAAGGCCAGCTGTTCAstochastic_pseudotime"
#>  [279] "TTGAACGGTGATGCCCstochastic_pseudotime"
#>  [280] "ACTGAGTAGTGACTCTstochastic_pseudotime"
#>  [281] "GTATCTTCAGTGGGATstochastic_pseudotime"
#>  [282] "TCTGGAATCGGTGTCGstochastic_pseudotime"
#>  [283] "CATGACACAACAACCTstochastic_pseudotime"
#>  [284] "GAGGTGAGTGACAAATstochastic_pseudotime"
#>  [285] "AACTTTCGTGACCAAGstochastic_pseudotime"
#>  [286] "TACACGACAACGATCTstochastic_pseudotime"
#>  [287] "CCACTACAGAGTACCGstochastic_pseudotime"
#>  [288] "TCTTTCCTCGTCTGCTstochastic_pseudotime"
#>  [289] "GGCCGATTCACTCTTAstochastic_pseudotime"
#>  [290] "GTAACGTGTCACTGGCstochastic_pseudotime"
#>  [291] "CGGTTAATCGTCTGAAstochastic_pseudotime"
#>  [292] "TCAGGTAAGTCAATAGstochastic_pseudotime"
#>  [293] "TGGCCAGCACGAAAGCstochastic_pseudotime"
#>  [294] "GATGAGGTCTGGGCCAstochastic_pseudotime"
#>  [295] "GCAGTTATCTTGAGACstochastic_pseudotime"
#>  [296] "TTCTACAGTACTCTCCstochastic_pseudotime"
#>  [297] "CGCTGGAGTCCATCCTstochastic_pseudotime"
#>  [298] "CGTGTAAGTCTTGATGstochastic_pseudotime"
#>  [299] "TGCGGGTTCTCAAACGstochastic_pseudotime"
#>  [300] "TACGGATGTGTTCGATstochastic_pseudotime"
#>  [301] "ATCCGAAGTGACTCATstochastic_pseudotime"
#>  [302] "CTCTGGTCAACACCTAstochastic_pseudotime"
#>  [303] "AAACGGGTCAGCTCTCstochastic_pseudotime"
#>  [304] "TGAGCCGTCCGAACGCstochastic_pseudotime"
#>  [305] "AGAGCTTAGCCGATTTstochastic_pseudotime"
#>  [306] "CGAATGTAGTTCCACAstochastic_pseudotime"
#>  [307] "CTTAGGAGTCCGAATTstochastic_pseudotime"
#>  [308] "ATTGGACAGTAATCCCstochastic_pseudotime"
#>  [309] "GGTGTTAAGGATCGCAstochastic_pseudotime"
#>  [310] "TGGACGCAGTGTTGAAstochastic_pseudotime"
#>  [311] "TCGTACCCATATACGCstochastic_pseudotime"
#>  [312] "CGACTTCCAGCCTATAstochastic_pseudotime"
#>  [313] "AAGGCAGTCCGAGCCAstochastic_pseudotime"
#>  [314] "ACGAGGAAGTTACCCAstochastic_pseudotime"
#>  [315] "TCAGCTCAGTCGTACTstochastic_pseudotime"
#>  [316] "CGTCAGGTCCTGCCATstochastic_pseudotime"
#>  [317] "TGATTTCAGTACGCCCstochastic_pseudotime"
#>  [318] "AGAGTGGGTCACCTAAstochastic_pseudotime"
#>  [319] "CGATCGGCACGCCAGTstochastic_pseudotime"
#>  [320] "TGCCCATCACTGTGTAstochastic_pseudotime"
#>  [321] "AACTTTCCACTGCCAGstochastic_pseudotime"
#>  [322] "GTCCTCATCATTTGGGstochastic_pseudotime"
#>  [323] "GTGCGGTAGTAGGCCAstochastic_pseudotime"
#>  [324] "CGACCTTAGCGCTTATstochastic_pseudotime"
#>  [325] "GTCTCGTAGGCTAGACstochastic_pseudotime"
#>  [326] "GCATGATAGCCTTGATstochastic_pseudotime"
#>  [327] "ACTTTCACAACTGCTAstochastic_pseudotime"
#>  [328] "CCAATCCAGAGTAAGGstochastic_pseudotime"
#>  [329] "ATTATCCAGAAACCTAstochastic_pseudotime"
#>  [330] "TCGGGACAGAGAACAGstochastic_pseudotime"
#>  [331] "AGCGGTCGTACCGAGAstochastic_pseudotime"
#>  [332] "TTCCCAGAGTTACCCAstochastic_pseudotime"
#>  [333] "CTCGAGGTCTGTTGAGstochastic_pseudotime"
#>  [334] "ATTATCCAGACCCACCstochastic_pseudotime"
#>  [335] "CCTCAGTGTTACCAGTstochastic_pseudotime"
#>  [336] "GCCTCTAAGTCCGGTCstochastic_pseudotime"
#>  [337] "TTCGAAGAGTACGCCCstochastic_pseudotime"
#>  [338] "GCTGCGAAGCTAGCCCstochastic_pseudotime"
#>  [339] "GCGGGTTAGACTTGAAstochastic_pseudotime"
#>  [340] "GTATCTTCATCGATGTstochastic_pseudotime"
#>  [341] "ACTTACTAGAGGGCTTstochastic_pseudotime"
#>  [342] "GGGTCTGCACAGGCCTstochastic_pseudotime"
#>  [343] "GGTATTGGTACTCGCGstochastic_pseudotime"
#>  [344] "TGTGGTATCCAGTAGTstochastic_pseudotime"
#>  [345] "CTCAGAACACTAGTACstochastic_pseudotime"
#>  [346] "CAGAATCGTCATTAGCstochastic_pseudotime"
#>  [347] "ACCAGTATCGTTTAGGstochastic_pseudotime"
#>  [348] "TTTACTGCAGCATACTstochastic_pseudotime"
#>  [349] "GCAGCCACAGCGTTCGstochastic_pseudotime"
#>  [350] "CTGTGCTCATTACCTTstochastic_pseudotime"
#>  [351] "CTCGAAAGTAAGAGGAstochastic_pseudotime"
#>  [352] "TATGCCCGTCTAGCCGstochastic_pseudotime"
#>  [353] "GGGCATCGTGGTCTCGstochastic_pseudotime"
#>  [354] "CATGACACATCGATTGstochastic_pseudotime"
#>  [355] "AACTTTCCACAGTCGCstochastic_pseudotime"
#>  [356] "GAAATGAGTCACCCAGstochastic_pseudotime"
#>  [357] "ACTTACTCAATCAGAAstochastic_pseudotime"
#>  [358] "TGCTGCTAGGATGCGTstochastic_pseudotime"
#>  [359] "AACTCAGGTGCCTGTGstochastic_pseudotime"
#>  [360] "GTATCTTGTGGCCCTAstochastic_pseudotime"
#>  [361] "GTCCTCATCTGGAGCCstochastic_pseudotime"
#>  [362] "TGACTAGTCCTATTCAstochastic_pseudotime"
#>  [363] "TGAGCATGTTAAAGACstochastic_pseudotime"
#>  [364] "TGCCAAACATGACATCstochastic_pseudotime"
#>  [365] "AAGTCTGTCCGAGCCAstochastic_pseudotime"
#>  [366] "ACATCAGTCCTTGGTCstochastic_pseudotime"
#>  [367] "GTGCGGTCAGCGTAAGstochastic_pseudotime"
#>  [368] "GTCTCGTCATAGTAAGstochastic_pseudotime"
#>  [369] "GATGAGGAGTGGACGTstochastic_pseudotime"
#>  [370] "AGAGCTTAGGTTACCTstochastic_pseudotime"
#>  [371] "GCTGGGTGTCTAGAGGstochastic_pseudotime"
#>  [372] "GGGATGAGTTGTGGAGstochastic_pseudotime"
#>  [373] "AGCGGTCGTCGAATCTstochastic_pseudotime"
#>  [374] "AGGCCGTCACGACGAAstochastic_pseudotime"
#>  [375] "CGAGAAGCACGCATCGstochastic_pseudotime"
#>  [376] "GTACTCCTCCAACCAAstochastic_pseudotime"
#>  [377] "ATGTGTGCACAAGTAAstochastic_pseudotime"
#>  [378] "TGACTAGGTGTATGGGstochastic_pseudotime"
#>  [379] "AGGGATGGTTAGGGTGstochastic_pseudotime"
#>  [380] "CGTGTAATCATTATCCstochastic_pseudotime"
#>  [381] "CTGAAACTCGATAGAAstochastic_pseudotime"
#>  [382] "TTGAACGAGGAGTTTAstochastic_pseudotime"
#>  [383] "AGCGTCGCACAGGCCTstochastic_pseudotime"
#>  [384] "CCACGGAGTTACGACTstochastic_pseudotime"
#>  [385] "CTGAAACTCTGTTTGTstochastic_pseudotime"
#>  [386] "AGATTGCAGTTCCACAstochastic_pseudotime"
#>  [387] "GTATCTTGTCGCGGTTstochastic_pseudotime"
#>  [388] "GAGTCCGAGCGATGACstochastic_pseudotime"
#>  [389] "AACCATGCACAGGAGTstochastic_pseudotime"
#>  [390] "ATGTGTGTCAAAGACAstochastic_pseudotime"
#>  [391] "CCGTACTAGGAGTACCstochastic_pseudotime"
#>  [392] "CCGTACTAGTCCGTATstochastic_pseudotime"
#>  [393] "GACGTTAAGACGCTTTstochastic_pseudotime"
#>  [394] "CCACGGAGTCAACATCstochastic_pseudotime"
#>  [395] "CATGGCGAGAGTCTGGstochastic_pseudotime"
#>  [396] "GAAATGAAGCTCCTCTstochastic_pseudotime"
#>  [397] "ACGCCAGTCCCAAGTAstochastic_pseudotime"
#>  [398] "CGGTTAAGTCTAGCCGstochastic_pseudotime"
#>  [399] "GTCTCGTTCACTTCATstochastic_pseudotime"
#>  [400] "GACCAATTCAACGGCCstochastic_pseudotime"
#>  [401] "CTCTGGTAGCAATCTCstochastic_pseudotime"
#>  [402] "TCTCTAAGTGGCAAACstochastic_pseudotime"
#>  [403] "CATCGGGTCTTCGGTCstochastic_pseudotime"
#>  [404] "CAACCAACACCTCGTTstochastic_pseudotime"
#>  [405] "CTGGTCTTCGCATGGCstochastic_pseudotime"
#>  [406] "ACTGTCCTCGCTAGCGstochastic_pseudotime"
#>  [407] "CTAAGACAGCACGCCTstochastic_pseudotime"
#>  [408] "CGGCTAGGTACTCGCGstochastic_pseudotime"
#>  [409] "TCTCATAGTGGACGATstochastic_pseudotime"
#>  [410] "GGATGTTTCCAGTATGstochastic_pseudotime"
#>  [411] "CTTAGGATCACTCTTAstochastic_pseudotime"
#>  [412] "TACGGATTCCAAGTACstochastic_pseudotime"
#>  [413] "TGAGAGGAGCCGTCGTstochastic_pseudotime"
#>  [414] "CTTGGCTCAGCCACCAstochastic_pseudotime"
#>  [415] "AATCGGTCATATGGTCstochastic_pseudotime"
#>  [416] "TACTTACCATGGGACAstochastic_pseudotime"
#>  [417] "TGAGGGAAGAGGTACCstochastic_pseudotime"
#>  [418] "TCCCGATCAACGATGGstochastic_pseudotime"
#>  [419] "GTCACAACACGTCTCTstochastic_pseudotime"
#>  [420] "TTTGGTTCACGAAGCAstochastic_pseudotime"
#>  [421] "GGGAATGTCACAAACCstochastic_pseudotime"
#>  [422] "ACAGCCGCATCGTCGGstochastic_pseudotime"
#>  [423] "CGTGAGCTCAGCGATTstochastic_pseudotime"
#>  [424] "ATCATGGGTCCGCTGAstochastic_pseudotime"
#>  [425] "GGCAATTCATCTCCCAstochastic_pseudotime"
#>  [426] "CCTATTATCATCGCTCstochastic_pseudotime"
#>  [427] "CAAGGCCCACATTTCTstochastic_pseudotime"
#>  [428] "CTCACACGTCACACGCstochastic_pseudotime"
#>  [429] "AACCATGGTTCGGCACstochastic_pseudotime"
#>  [430] "CTGCGGACAGACTCGCstochastic_pseudotime"
#>  [431] "ACGCAGCGTAATCACCstochastic_pseudotime"
#>  [432] "CCATGTCGTCGGCTCAstochastic_pseudotime"
#>  [433] "TTGCGTCCACCGATATstochastic_pseudotime"
#>  [434] "CACAAACAGTAGCGGTstochastic_pseudotime"
#>  [435] "ACTGATGCATCAGTACstochastic_pseudotime"
#>  [436] "TAGTTGGCACGTAAGGstochastic_pseudotime"
#>  [437] "CAGCTAAGTACCGAGAstochastic_pseudotime"
#>  [438] "CGTCCATTCGCGCCAAstochastic_pseudotime"
#>  [439] "TTCGGTCAGCTAGCCCstochastic_pseudotime"
#>  [440] "GTTAAGCAGGATTCGGstochastic_pseudotime"
#>  [441] "GGACAGATCCTTTCGGstochastic_pseudotime"
#>  [442] "AGTTGGTTCACATAGCstochastic_pseudotime"
#>  [443] "TTCTCAACAATAGAGTstochastic_pseudotime"
#>  [444] "GTAGTCACATCGATGTstochastic_pseudotime"
#>  [445] "GCTGCGAAGTATTGGAstochastic_pseudotime"
#>  [446] "CGTCAGGGTAGGCATGstochastic_pseudotime"
#>  [447] "AACTCAGTCCTCGCATstochastic_pseudotime"
#>  [448] "AAATGCCCACCAACCGstochastic_pseudotime"
#>  [449] "CCGTTCAAGGACCACAstochastic_pseudotime"
#>  [450] "AATCGGTTCCAAACTGstochastic_pseudotime"
#>  [451] "TCACGAAGTCCGAGTCstochastic_pseudotime"
#>  [452] "GGCTGGTAGTCACGCCstochastic_pseudotime"
#>  [453] "GGACGTCCATGGATGGstochastic_pseudotime"
#>  [454] "CAGCCGATCGTTTATCstochastic_pseudotime"
#>  [455] "TTGACTTGTTCCATGAstochastic_pseudotime"
#>  [456] "GGCTGGTTCATGCATGstochastic_pseudotime"
#>  [457] "ACACCCTAGGACAGCTstochastic_pseudotime"
#>  [458] "CGATGTAGTACTCTCCstochastic_pseudotime"
#>  [459] "TCAGGTAAGGCGTACAstochastic_pseudotime"
#>  [460] "AGCAGCCGTGTAACGGstochastic_pseudotime"
#>  [461] "CCTACCACATCGATTGstochastic_pseudotime"
#>  [462] "AAGACCTAGTAGGCCAstochastic_pseudotime"
#>  [463] "GGACAAGTCAAAGACAstochastic_pseudotime"
#>  [464] "GACTAACGTCATCCCTstochastic_pseudotime"
#>  [465] "CGCCAAGCAAAGTCAAstochastic_pseudotime"
#>  [466] "CCGGGATAGATCTGAAstochastic_pseudotime"
#>  [467] "TACACGAGTGCAACTTstochastic_pseudotime"
#>  [468] "CGTTCTGTCACAATGCstochastic_pseudotime"
#>  [469] "CAGCCGATCCTATTCAstochastic_pseudotime"
#>  [470] "AGGGATGTCGCATGGCstochastic_pseudotime"
#>  [471] "GTAGTCAGTCCGACGTstochastic_pseudotime"
#>  [472] "AGCAGCCGTCGCATCGstochastic_pseudotime"
#>  [473] "TGCTACCCATAAGACAstochastic_pseudotime"
#>  [474] "AGCTCCTTCAATACCGstochastic_pseudotime"
#>  [475] "GAATGAAAGAATGTGTstochastic_pseudotime"
#>  [476] "TTGTAGGTCAAAGTAGstochastic_pseudotime"
#>  [477] "AAACCTGGTAAGTGGCstochastic_pseudotime"
#>  [478] "AAATGCCCATGCTAGTstochastic_pseudotime"
#>  [479] "GACTGCGAGAATCTCCstochastic_pseudotime"
#>  [480] "CGATTGAAGTCAAGGCstochastic_pseudotime"
#>  [481] "CCATGTCTCCGCATAAstochastic_pseudotime"
#>  [482] "CAAGAAAAGAGTACCGstochastic_pseudotime"
#>  [483] "GGCGTGTCACCAGGTCstochastic_pseudotime"
#>  [484] "CCGTGGATCTTAGAGCstochastic_pseudotime"
#>  [485] "TGCGGGTAGACGCACAstochastic_pseudotime"
#>  [486] "AAGCCGCCAATGGACGstochastic_pseudotime"
#>  [487] "ACACTGAGTCCATGATstochastic_pseudotime"
#>  [488] "AGCTTGAAGACACGACstochastic_pseudotime"
#>  [489] "TCTATTGTCGTATCAGstochastic_pseudotime"
#>  [490] "TCGCGAGAGTACGCGAstochastic_pseudotime"
#>  [491] "TCAACGATCAACCATGstochastic_pseudotime"
#>  [492] "CACCAGGCATAAAGGTstochastic_pseudotime"
#>  [493] "TGAGCCGTCCGCTGTTstochastic_pseudotime"
#>  [494] "ACTGTCCCAGTTTACGstochastic_pseudotime"
#>  [495] "GATCGATAGAAAGTGGstochastic_pseudotime"
#>  [496] "GATCGTAAGTAGGTGCstochastic_pseudotime"
#>  [497] "CACAGGCCAAGGCTCCstochastic_pseudotime"
#>  [498] "TGTGTTTGTATGAATGstochastic_pseudotime"
#>  [499] "GGCTCGAAGATCACGGstochastic_pseudotime"
#>  [500] "GGGTCTGTCATGTGGTstochastic_pseudotime"
#>  [501] "TACAGTGGTAAACGCGstochastic_pseudotime"
#>  [502] "TCTCTAATCCCAACGGstochastic_pseudotime"
#>  [503] "AGTTGGTGTTGTCGCGstochastic_pseudotime"
#>  [504] "CCAGCGACAGACGTAGstochastic_pseudotime"
#>  [505] "TCAATCTCATCCGCGAstochastic_pseudotime"
#>  [506] "TCAGATGCACTCGACGstochastic_pseudotime"
#>  [507] "GTTACAGTCCGATATGstochastic_pseudotime"
#>  [508] "TGGACGCTCTTGTTTGstochastic_pseudotime"
#>  [509] "CTGCGGACATTAGCCAstochastic_pseudotime"
#>  [510] "CACCACTGTGTGTGCCstochastic_pseudotime"
#>  [511] "GTAACTGGTCACTGGCstochastic_pseudotime"
#>  [512] "GGCGTGTCACCGATATstochastic_pseudotime"
#>  [513] "AGAGTGGCACCTCGGAstochastic_pseudotime"
#>  [514] "CAGCTAACATGCCTAAstochastic_pseudotime"
#>  [515] "ATCCGAAGTACACCGCstochastic_pseudotime"
#>  [516] "CAAGGCCTCTCAACTTstochastic_pseudotime"
#>  [517] "ACTGATGCATGTCTCCstochastic_pseudotime"
#>  [518] "CATATGGTCACAGTACstochastic_pseudotime"
#>  [519] "TTCCCAGTCATAGCACstochastic_pseudotime"
#>  [520] "GGACGTCGTCAGAAGCstochastic_pseudotime"
#>  [521] "GGGCACTGTCTCCATCstochastic_pseudotime"
#>  [522] "CCCATACTCGTGGACCstochastic_pseudotime"
#>  [523] "TGCCCTAAGTCCGGTCstochastic_pseudotime"
#>  [524] "GCTGCGATCCCGGATGstochastic_pseudotime"
#>  [525] "TAGTGGTTCCTCAACCstochastic_pseudotime"
#>  [526] "CTACATTCACCGTTGGstochastic_pseudotime"
#>  [527] "AGGTCATGTTTAGGAAstochastic_pseudotime"
#>  [528] "CCTACCAGTTGCTCCTstochastic_pseudotime"
#>  [529] "CGTCACTTCACTTATCstochastic_pseudotime"
#>  [530] "CCTTCCCTCGTCGTTCstochastic_pseudotime"
#>  [531] "AATCCAGTCGGTCCGAstochastic_pseudotime"
#>  [532] "CAAGTTGTCTAACTGGstochastic_pseudotime"
#>  [533] "TTAGGCAGTTTGTTGGstochastic_pseudotime"
#>  [534] "AATCCAGTCAGGTTCAstochastic_pseudotime"
#>  [535] "GTACTCCCACAGGTTTstochastic_pseudotime"
#>  [536] "GCGCCAATCTTATCTGstochastic_pseudotime"
#>  [537] "TCATTTGCATTGGCGCstochastic_pseudotime"
#>  [538] "GACGTTAGTGTGAATAstochastic_pseudotime"
#>  [539] "AAAGCAAAGCTACCTAstochastic_pseudotime"
#>  [540] "ACATCAGAGAGACTTAstochastic_pseudotime"
#>  [541] "ACATACGGTCTCTCTGstochastic_pseudotime"
#>  [542] "TTAACTCGTAATCGTCstochastic_pseudotime"
#>  [543] "GTTCATTAGGCAGTCAstochastic_pseudotime"
#>  [544] "CATTATCGTCATATCGstochastic_pseudotime"
#>  [545] "TTTCCTCGTGGCGAATstochastic_pseudotime"
#>  [546] "GCAGTTACATTTCACTstochastic_pseudotime"
#>  [547] "CTTAACTTCAGCTTAGstochastic_pseudotime"
#>  [548] "AGATCTGTCCGTCAAAstochastic_pseudotime"
#>  [549] "TCTCATAAGAATGTGTstochastic_pseudotime"
#>  [550] "GTATCTTTCGTGGACCstochastic_pseudotime"
#>  [551] "GAAGCAGGTCGCGAAAstochastic_pseudotime"
#>  [552] "AACACGTGTCGACTGCstochastic_pseudotime"
#>  [553] "GTCTTCGTCTACCAGAstochastic_pseudotime"
#>  [554] "TAGTTGGCAGTAACGGstochastic_pseudotime"
#>  [555] "AATCCAGAGGGTTCCCstochastic_pseudotime"
#>  [556] "TCCCGATAGTCGTTTGstochastic_pseudotime"
#>  [557] "CGAGCACAGGCAGGTTstochastic_pseudotime"
#>  [558] "TAAGAGATCTGCGACGstochastic_pseudotime"
#>  [559] "TGGGAAGCACGTAAGGstochastic_pseudotime"
#>  [560] "CACACCTGTTATCCGAstochastic_pseudotime"
#>  [561] "GATCGCGCACCAGGCTstochastic_pseudotime"
#>  [562] "CTCGGAGGTGAAGGCTstochastic_pseudotime"
#>  [563] "AAGGTTCGTCAACTGTstochastic_pseudotime"
#>  [564] "CAGCTGGGTCTAACGTstochastic_pseudotime"
#>  [565] "CAACCAAAGCGTTCCGstochastic_pseudotime"
#>  [566] "TGGCCAGCAATGAAACstochastic_pseudotime"
#>  [567] "CTAATGGCATGGTTGTstochastic_pseudotime"
#>  [568] "CTTGGCTCAAAGGAAGstochastic_pseudotime"
#>  [569] "TCAGCAATCGCCGTGAstochastic_pseudotime"
#>  [570] "CGGTTAATCACCGGGTstochastic_pseudotime"
#>  [571] "CTACATTGTTACTGACstochastic_pseudotime"
#>  [572] "CTGCTGTAGACTGGGTstochastic_pseudotime"
#>  [573] "GCTCTGTGTCCATGATstochastic_pseudotime"
#>  [574] "CGACTTCCAGCTGCACstochastic_pseudotime"
#>  [575] "CAGCTGGCAAGTTGTCstochastic_pseudotime"
#>  [576] "CCTTACGCACAACGCCstochastic_pseudotime"
#>  [577] "TAGTGGTCATTCACTTstochastic_pseudotime"
#>  [578] "CAGGTGCCATCAGTCAstochastic_pseudotime"
#>  [579] "TACTCATAGGTGTTAAstochastic_pseudotime"
#>  [580] "AGCTCCTCAGTTTACGstochastic_pseudotime"
#>  [581] "GTCACGGCATCCGTGGstochastic_pseudotime"
#>  [582] "GTCACGGTCTTGCCGTstochastic_pseudotime"
#>  [583] "TTGACTTTCACCGGGTstochastic_pseudotime"
#>  [584] "CCTCAGTAGTGCGATGstochastic_pseudotime"
#>  [585] "TTCGGTCTCGACCAGCstochastic_pseudotime"
#>  [586] "CGTTGGGCAGGCGATAstochastic_pseudotime"
#>  [587] "CCACGGAAGTCTCCTCstochastic_pseudotime"
#>  [588] "GATGCTAGTGTTAAGAstochastic_pseudotime"
#>  [589] "GTGTTAGGTGAACCTTstochastic_pseudotime"
#>  [590] "AGCTCCTTCGCCCTTAstochastic_pseudotime"
#>  [591] "TGCGTGGGTCCCTTGTstochastic_pseudotime"
#>  [592] "GGAGCAAGTATCACCAstochastic_pseudotime"
#>  [593] "CGTAGCGTCGGAAATAstochastic_pseudotime"
#>  [594] "CTGGTCTAGGATTCGGstochastic_pseudotime"
#>  [595] "GGCGTGTTCGAGAGCAstochastic_pseudotime"
#>  [596] "AACCATGGTTGAACTCstochastic_pseudotime"
#>  [597] "ACAGCCGAGAGAACAGstochastic_pseudotime"
#>  [598] "CGAATGTGTCTGGAGAstochastic_pseudotime"
#>  [599] "GACGGCTTCAACTCTTstochastic_pseudotime"
#>  [600] "CTTTGCGAGACCACGAstochastic_pseudotime"
#>  [601] "TTCTACACATAAAGGTstochastic_pseudotime"
#>  [602] "GCATGATGTCACTGGCstochastic_pseudotime"
#>  [603] "ACACCCTTCCCAGGTGstochastic_pseudotime"
#>  [604] "CTCCTAGCACAACGTTstochastic_pseudotime"
#>  [605] "CTCTGGTCAAGTTGTCstochastic_pseudotime"
#>  [606] "CATCAGATCATCGCTCstochastic_pseudotime"
#>  [607] "CCACTACAGAGCTGGTstochastic_pseudotime"
#>  [608] "CTCGGAGGTCATTAGCstochastic_pseudotime"
#>  [609] "TGGTTAGCAGTCGTGCstochastic_pseudotime"
#>  [610] "TCGAGGCAGGCATGGTstochastic_pseudotime"
#>  [611] "GTACTTTTCTTCGGTCstochastic_pseudotime"
#>  [612] "ACGAGGATCGCCGTGAstochastic_pseudotime"
#>  [613] "ATGAGGGTCTTGACGAstochastic_pseudotime"
#>  [614] "CGATCGGAGCTGCAAGstochastic_pseudotime"
#>  [615] "GTAGGCCAGCTAACAAstochastic_pseudotime"
#>  [616] "GTCCTCACACCACGTGstochastic_pseudotime"
#>  [617] "CTGAAGTTCTGCGGCAstochastic_pseudotime"
#>  [618] "GGACATTGTCTCCCTAstochastic_pseudotime"
#>  [619] "CGAGCACGTAGTAGTAstochastic_pseudotime"
#>  [620] "GTGCATAAGAACAATCstochastic_pseudotime"
#>  [621] "ATCTGCCCATGTTGACstochastic_pseudotime"
#>  [622] "ACTGCTCAGGTGTTAAstochastic_pseudotime"
#>  [623] "AGCTCTCAGACATAACstochastic_pseudotime"
#>  [624] "TGACGGCCATATGAGAstochastic_pseudotime"
#>  [625] "CAAGTTGGTGTGCCTGstochastic_pseudotime"
#>  [626] "GATCAGTAGGTGATATstochastic_pseudotime"
#>  [627] "ACGCCGAAGGGTGTTGstochastic_pseudotime"
#>  [628] "GGACGTCAGGGTGTTGstochastic_pseudotime"
#>  [629] "GGACAGACAATCCAACstochastic_pseudotime"
#>  [630] "TCGCGAGCACTTCTGCstochastic_pseudotime"
#>  [631] "TAAGTGCTCCAAATGCstochastic_pseudotime"
#>  [632] "TGCGTGGAGGTTCCTAstochastic_pseudotime"
#>  [633] "ACTATCTTCAAGATCCstochastic_pseudotime"
#>  [634] "TGTGTTTAGGAATGGAstochastic_pseudotime"
#>  [635] "GAGCAGATCAACTCTTstochastic_pseudotime"
#>  [636] "ACGATGTTCGCAGGCTstochastic_pseudotime"
#>  [637] "GTCATTTGTCGTCTTCstochastic_pseudotime"
#>  [638] "TTTGGTTTCCTTTCGGstochastic_pseudotime"
#>  [639] "AAAGCAAGTGCAGTAGstochastic_pseudotime"
#>  [640] "CACACCTGTACATGTCstochastic_pseudotime"
#>  [641] "CGATGTAAGTTTAGGAstochastic_pseudotime"
#>  [642] "GCGCAGTGTGCTAGCCstochastic_pseudotime"
#>  [643] "CTGAAACGTAGCGCAAstochastic_pseudotime"
#>  [644] "CGTGTCTGTAGCGTCCstochastic_pseudotime"
#>  [645] "AGTGGGAGTGGCAAACstochastic_pseudotime"
#>  [646] "AGGCCACCATGACATCstochastic_pseudotime"
#>  [647] "TTTGTCAAGAATGTGTstochastic_pseudotime"
#>  [648] "AAGGAGCAGTACATGAstochastic_pseudotime"
#>  [649] "ACACCGGTCATGTCTTstochastic_pseudotime"
#>  [650] "GAGGTGAGTCCGCTGAstochastic_pseudotime"
#>  [651] "CTCGAAAGTGGTGTAGstochastic_pseudotime"
#>  [652] "AAGGTTCAGATAGTCAstochastic_pseudotime"
#>  [653] "AGTGTCATCCCTTGTGstochastic_pseudotime"
#>  [654] "ATCACGAGTGCACCACstochastic_pseudotime"
#>  [655] "TACGGGCAGCCAACAGstochastic_pseudotime"
#>  [656] "CTCAGAACATAGGATAstochastic_pseudotime"
#>  [657] "TATCTCACACGCGAAAstochastic_pseudotime"
#>  [658] "GATTCAGTCGCCAGCAstochastic_pseudotime"
#>  [659] "TCACAAGCAACGATGGstochastic_pseudotime"
#>  [660] "TGTGTTTGTAGATTAGstochastic_pseudotime"
#>  [661] "GGCGTGTTCAGCTTAGstochastic_pseudotime"
#>  [662] "AGTAGTCGTGAGGGAGstochastic_pseudotime"
#>  [663] "GTACGTAGTCTCTTTAstochastic_pseudotime"
#>  [664] "GAATGAATCATCGGATstochastic_pseudotime"
#>  [665] "GCTGCGAAGAACAACTstochastic_pseudotime"
#>  [666] "TCTGGAAGTCATCGGCstochastic_pseudotime"
#>  [667] "CTGCTGTCAAAGCAATstochastic_pseudotime"
#>  [668] "CATCAGAGTCTCTTTAstochastic_pseudotime"
#>  [669] "GCAATCACAGCTGCTGstochastic_pseudotime"
#>  [670] "CCCAGTTTCAGCACATstochastic_pseudotime"
#>  [671] "CACAGTAGTGCACTTAstochastic_pseudotime"
#>  [672] "GAACCTATCACTTATCstochastic_pseudotime"
#>  [673] "GGCAATTCACGGCTACstochastic_pseudotime"
#>  [674] "TTGAACGGTGAGGGAGstochastic_pseudotime"
#>  [675] "GATTCAGAGCGATGACstochastic_pseudotime"
#>  [676] "GCATGCGAGGAACTGCstochastic_pseudotime"
#>  [677] "GAACGGAAGTCAAGCGstochastic_pseudotime"
#>  [678] "ACCAGTAGTTCGTTGAstochastic_pseudotime"
#>  [679] "GGAAAGCCAATGTAAGstochastic_pseudotime"
#>  [680] "CTCGTACTCTTTAGTCstochastic_pseudotime"
#>  [681] "GTGCTTCAGAGGGATAstochastic_pseudotime"
#>  [682] "AGGGAGTTCGGAGGTAstochastic_pseudotime"
#>  [683] "CTAACTTAGAGGTACCstochastic_pseudotime"
#>  [684] "CAGGTGCAGCCAGGATstochastic_pseudotime"
#>  [685] "GCATGATGTTTAGGAAstochastic_pseudotime"
#>  [686] "CCACGGAAGGGTGTTGstochastic_pseudotime"
#>  [687] "CATGACATCATGCATGstochastic_pseudotime"
#>  [688] "GACGTGCGTCACCTAAstochastic_pseudotime"
#>  [689] "CCACGGAAGATCCTGTstochastic_pseudotime"
#>  [690] "TTGGAACAGAAGATTCstochastic_pseudotime"
#>  [691] "CCACCTATCTAAGCCAstochastic_pseudotime"
#>  [692] "AGGTCATGTCTCCCTAstochastic_pseudotime"
#>  [693] "ATAGACCAGTGATCGGstochastic_pseudotime"
#>  [694] "GAACCTAAGTAGATGTstochastic_pseudotime"
#>  [695] "GATCGATGTTGTTTGGstochastic_pseudotime"
#>  [696] "TTCTCCTAGTCCATACstochastic_pseudotime"
#>  [697] "AACTCCCCAGTTCATGstochastic_pseudotime"
#>  [698] "GGGATGAGTCTGCGGTstochastic_pseudotime"
#>  [699] "CTCGAAATCCTGCCATstochastic_pseudotime"
#>  [700] "CGGGTCAAGGCGATACstochastic_pseudotime"
#>  [701] "GACCTGGGTTCGTGATstochastic_pseudotime"
#>  [702] "CCTACACAGTACGTTCstochastic_pseudotime"
#>  [703] "GACCAATGTTATCCGAstochastic_pseudotime"
#>  [704] "CATCGGGCACGTTGGCstochastic_pseudotime"
#>  [705] "ACTTGTTTCACCCGAGstochastic_pseudotime"
#>  [706] "CCGTTCAAGCTAAACAstochastic_pseudotime"
#>  [707] "ATCATGGAGGAATGGAstochastic_pseudotime"
#>  [708] "AGGGATGGTCGACTGCstochastic_pseudotime"
#>  [709] "ACCAGTAAGTGATCGGstochastic_pseudotime"
#>  [710] "GCGACCATCACCAGGCstochastic_pseudotime"
#>  [711] "CGAGAAGTCGTTACAGstochastic_pseudotime"
#>  [712] "CTCTAATAGTTATCGCstochastic_pseudotime"
#>  [713] "CTACGTCAGCTAGGCAstochastic_pseudotime"
#>  [714] "ACTTGTTAGTGAAGTTstochastic_pseudotime"
#>  [715] "CGTTCTGAGGACGAAAstochastic_pseudotime"
#>  [716] "GTGCATAGTTCAACCAstochastic_pseudotime"
#>  [717] "ACAGCTAAGTCGTTTGstochastic_pseudotime"
#>  [718] "CACCACTAGCGATCCCstochastic_pseudotime"
#>  [719] "CATGCCTCACCCATTCstochastic_pseudotime"
#>  [720] "CAAGGCCTCAGTGCATstochastic_pseudotime"
#>  [721] "CCGGGATTCGACAGCCstochastic_pseudotime"
#>  [722] "CCACGGAAGCAATCTCstochastic_pseudotime"
#>  [723] "GTGGGTCTCCGAATGTstochastic_pseudotime"
#>  [724] "GCATGTAGTTGCGCACstochastic_pseudotime"
#>  [725] "TGCGTGGAGTACACCTstochastic_pseudotime"
#>  [726] "CGCTTCATCTCTGTCGstochastic_pseudotime"
#>  [727] "AGGTCATCACTTCTGCstochastic_pseudotime"
#>  [728] "TATGCCCAGAACAATCstochastic_pseudotime"
#>  [729] "ATAGACCAGCGTGTCCstochastic_pseudotime"
#>  [730] "CGCCAAGTCTAAGCCAstochastic_pseudotime"
#>  [731] "CGCGGTAAGCTGCAAGstochastic_pseudotime"
#>  [732] "AGAATAGCAATCACACstochastic_pseudotime"
#>  [733] "CGTCTACAGGGAGTAAstochastic_pseudotime"
#>  [734] "GGGAATGCAAAGTGCGstochastic_pseudotime"
#>  [735] "AACTTTCTCCAGATCAstochastic_pseudotime"
#>  [736] "AGAGTGGAGGCGATACstochastic_pseudotime"
#>  [737] "AAACGGGCAAAGAATCstochastic_pseudotime"
#>  [738] "GGGATGAAGGGCATGTstochastic_pseudotime"
#>  [739] "AGTCTTTCACGGACAAstochastic_pseudotime"
#>  [740] "GCAAACTCAGACAAATstochastic_pseudotime"
#>  [741] "ACACCAAGTATAGGTAstochastic_pseudotime"
#>  [742] "ACGCAGCCATTCTCATstochastic_pseudotime"
#>  [743] "CTAACTTCATCTGGTAstochastic_pseudotime"
#>  [744] "TGGGCGTGTGTGCGTCstochastic_pseudotime"
#>  [745] "TGAGAGGCAGGAACGTstochastic_pseudotime"
#>  [746] "AAGACCTAGACCTTTGstochastic_pseudotime"
#>  [747] "GCGCCAAAGCTCCTTCstochastic_pseudotime"
#>  [748] "GCATGTATCGTGGACCstochastic_pseudotime"
#>  [749] "GGGAGATAGTGACATAstochastic_pseudotime"
#>  [750] "GTTCTCGAGAATTCCCstochastic_pseudotime"
#>  [751] "GGGCATCCATTTCAGGstochastic_pseudotime"
#>  [752] "TATTACCTCGTTGCCTstochastic_pseudotime"
#>  [753] "CGAATGTTCTTCATGTstochastic_pseudotime"
#>  [754] "CGGAGTCAGTTTAGGAstochastic_pseudotime"
#>  [755] "AGATTGCAGCGTGTCCstochastic_pseudotime"
#>  [756] "GGATGTTAGTAAGTACstochastic_pseudotime"
#>  [757] "GATCGATAGGATGGTCstochastic_pseudotime"
#>  [758] "TCAGGTACAGTCACTAstochastic_pseudotime"
#>  [759] "CCTTCGATCACAGGCCstochastic_pseudotime"
#>  [760] "GGCGACTGTAAGTTCCstochastic_pseudotime"
#>  [761] "TTCGGTCGTATAGTAGstochastic_pseudotime"
#>  [762] "CTTCTCTCAGCATACTstochastic_pseudotime"
#>  [763] "TTGCCGTCACCGAATTstochastic_pseudotime"
#>  [764] "CGTGTCTCAGAAGCACstochastic_pseudotime"
#>  [765] "AACTCTTCATTTGCCCstochastic_pseudotime"
#>  [766] "GCATGCGCAGCTCCGAstochastic_pseudotime"
#>  [767] "TTTGGTTAGACAGGCTstochastic_pseudotime"
#>  [768] "CGATGGCTCAGTGTTGstochastic_pseudotime"
#>  [769] "TACTCGCAGTAAGTACstochastic_pseudotime"
#>  [770] "TCAATCTTCGAACGGAstochastic_pseudotime"
#>  [771] "GCTGGGTGTCATACTGstochastic_pseudotime"
#>  [772] "ACGATACAGATGAGAGstochastic_pseudotime"
#>  [773] "GTATCTTGTACTCAACstochastic_pseudotime"
#>  [774] "ACGATACTCGCGCCAAstochastic_pseudotime"
#>  [775] "GATCGATGTAGCGATGstochastic_pseudotime"
#>  [776] "CCTACACTCTTTACACstochastic_pseudotime"
#>  [777] "AAGGTTCAGGAGCGAGstochastic_pseudotime"
#>  [778] "CAGTCCTCAGGGCATAstochastic_pseudotime"
#>  [779] "CCCATACCATGCAACTstochastic_pseudotime"
#>  [780] "GCTTGAAGTCAACTGTstochastic_pseudotime"
#>  [781] "CTTCTCTAGTCCGTATstochastic_pseudotime"
#>  [782] "AGCTTGACATCTATGGstochastic_pseudotime"
#>  [783] "CTCGAGGCAAGTCTGTstochastic_pseudotime"
#>  [784] "CAGAATCCATTGGGCCstochastic_pseudotime"
#>  [785] "CTTCTCTGTCAGAGGTstochastic_pseudotime"
#>  [786] "GATGAAACAGGCTGAAstochastic_pseudotime"
#>  [787] "ATTACTCTCTGCCAGGstochastic_pseudotime"
#>  [788] "GGGTCTGTCCATGAACstochastic_pseudotime"
#>  [789] "ACCTTTATCCTGCAGGstochastic_pseudotime"
#>  [790] "GGGAGATAGATGTGGCstochastic_pseudotime"
#>  [791] "ATCCACCCAGCTCGCAstochastic_pseudotime"
#>  [792] "TCGAGGCAGAACAATCstochastic_pseudotime"
#>  [793] "GCTTCCAAGTTAAGTGstochastic_pseudotime"
#>  [794] "AAATGCCCAATGGAATstochastic_pseudotime"
#>  [795] "AGTTGGTAGAGCCTAGstochastic_pseudotime"
#>  [796] "GTGTGCGGTGTAAGTAstochastic_pseudotime"
#>  [797] "CGGAGCTAGATAGTCAstochastic_pseudotime"
#>  [798] "GAACGGACAAGGACACstochastic_pseudotime"
#>  [799] "ACAGCTATCAAGGCTTstochastic_pseudotime"
#>  [800] "GGGACCTCACTGTCGGstochastic_pseudotime"
#>  [801] "CCCAATCCAACTGGCCstochastic_pseudotime"
#>  [802] "AATCGGTTCGTATCAGstochastic_pseudotime"
#>  [803] "GACCAATCAGACGCCTstochastic_pseudotime"
#>  [804] "CCTCAGTCAGACACTTstochastic_pseudotime"
#>  [805] "GCGAGAAGTGGGTATGstochastic_pseudotime"
#>  [806] "ACAGCCGAGAGGTACCstochastic_pseudotime"
#>  [807] "ACTTGTTAGTTAAGTGstochastic_pseudotime"
#>  [808] "CTGAAGTGTTTGGGCCstochastic_pseudotime"
#>  [809] "CACTCCACATCGTCGGstochastic_pseudotime"
#>  [810] "CCACCTAAGACTTGAAstochastic_pseudotime"
#>  [811] "GTCTTCGTCCGTTGCTstochastic_pseudotime"
#>  [812] "TACACGACATGTAGTCstochastic_pseudotime"
#>  [813] "GTGGGTCGTCCAGTATstochastic_pseudotime"
#>  [814] "GCGGGTTCAAGTAATGstochastic_pseudotime"
#>  [815] "CAGAGAGTCGGCCGATstochastic_pseudotime"
#>  [816] "TGCGGGTAGGGTTCCCstochastic_pseudotime"
#>  [817] "CATATTCTCCTGTAGAstochastic_pseudotime"
#>  [818] "CACAAACCACGAAAGCstochastic_pseudotime"
#>  [819] "TGCGTGGGTTACTGACstochastic_pseudotime"
#>  [820] "TCTTTCCGTGCAACTTstochastic_pseudotime"
#>  [821] "CCACTACCAACGCACCstochastic_pseudotime"
#>  [822] "TACCTTATCCGATATGstochastic_pseudotime"
#>  [823] "TGCGGGTTCCACTCCAstochastic_pseudotime"
#>  [824] "GCTTCCAAGACATAACstochastic_pseudotime"
#>  [825] "CAGTCCTAGGTGGGTTstochastic_pseudotime"
#>  [826] "CGCGGTAAGAGCTGGTstochastic_pseudotime"
#>  [827] "GAGCAGAGTTTGCATGstochastic_pseudotime"
#>  [828] "ACGCCAGGTTACAGAAstochastic_pseudotime"
#>  [829] "GAAGCAGGTCTAGAGGstochastic_pseudotime"
#>  [830] "GTTCTCGCAAGAAAGGstochastic_pseudotime"
#>  [831] "ACTGATGCATGGTCTAstochastic_pseudotime"
#>  [832] "TACAGTGCACAGGCCTstochastic_pseudotime"
#>  [833] "TCTCATACAAATTGCCstochastic_pseudotime"
#>  [834] "TCCACACCAGGCGATAstochastic_pseudotime"
#>  [835] "GATTCAGAGTACGTTCstochastic_pseudotime"
#>  [836] "AGTCTTTGTTACGGAGstochastic_pseudotime"
#>  [837] "TTTGTCATCGAATGCTstochastic_pseudotime"
#>  [838] "GCGCAACGTTTCGCTCstochastic_pseudotime"
#>  [839] "ACGCCGAGTTTGCATGstochastic_pseudotime"
#>  [840] "CGTGAGCTCTTCTGGCstochastic_pseudotime"
#>  [841] "CATATGGGTTACTGACstochastic_pseudotime"
#>  [842] "AGGGTGAGTGCGAAACstochastic_pseudotime"
#>  [843] "AAACGGGTCGCATGGCstochastic_pseudotime"
#>  [844] "GACGTGCTCCGGCACAstochastic_pseudotime"
#>  [845] "TGTCCCATCTCTGCTGstochastic_pseudotime"
#>  [846] "TTTATGCAGGTTACCTstochastic_pseudotime"
#>  [847] "ATCATGGGTATTCGTGstochastic_pseudotime"
#>  [848] "GGACGTCGTCGAGATGstochastic_pseudotime"
#>  [849] "CGACTTCTCACGCGGTstochastic_pseudotime"
#>  [850] "TTTGTCATCTGTTTGTstochastic_pseudotime"
#>  [851] "CGAGAAGAGGGTTCCCstochastic_pseudotime"
#>  [852] "CACAGGCCACGTCAGCstochastic_pseudotime"
#>  [853] "AAATGCCGTGATAAGTstochastic_pseudotime"
#>  [854] "TTATGCTCACCCAGTGstochastic_pseudotime"
#>  [855] "CCGTGGATCTGTCTCGstochastic_pseudotime"
#>  [856] "TTAACTCCAGGATTGGstochastic_pseudotime"
#>  [857] "TGAGAGGAGGTAAACTstochastic_pseudotime"
#>  [858] "TTGAACGAGAGCAATTstochastic_pseudotime"
#>  [859] "TGAGAGGAGCCAGGATstochastic_pseudotime"
#>  [860] "AGTGAGGTCTCCAACCstochastic_pseudotime"
#>  [861] "GGGAGATAGCGAAGGGstochastic_pseudotime"
#>  [862] "CTAAGACGTATCAGTCstochastic_pseudotime"
#>  [863] "CGGACTGTCTTGAGACstochastic_pseudotime"
#>  [864] "GACGCGTGTCGGCTCAstochastic_pseudotime"
#>  [865] "CCCAGTTGTCGTGGCTstochastic_pseudotime"
#>  [866] "GTGCATATCAAGGCTTstochastic_pseudotime"
#>  [867] "CTCGAGGGTAAGTAGTstochastic_pseudotime"
#>  [868] "GGGAATGAGCACACAGstochastic_pseudotime"
#>  [869] "TCGGGACGTCTGCGGTstochastic_pseudotime"
#>  [870] "GCACTCTCAGCGTCCAstochastic_pseudotime"
#>  [871] "CTCGTCAAGAACAATCstochastic_pseudotime"
#>  [872] "AGCTTGAAGAGTGACCstochastic_pseudotime"
#>  [873] "CGAATGTTCACCGGGTstochastic_pseudotime"
#>  [874] "CCTACCATCTGCCCTAstochastic_pseudotime"
#>  [875] "TTCTTAGTCAGTGTTGstochastic_pseudotime"
#>  [876] "CTGGTCTAGGTGATTAstochastic_pseudotime"
#>  [877] "GACCAATGTTGGACCCstochastic_pseudotime"
#>  [878] "TACTCGCGTCGGCATCstochastic_pseudotime"
#>  [879] "TACTCGCCAGATCGGAstochastic_pseudotime"
#>  [880] "TGTCCCACATTATCTCstochastic_pseudotime"
#>  [881] "CAGCTGGGTTCTCATTstochastic_pseudotime"
#>  [882] "AAAGATGAGTCTTGCAstochastic_pseudotime"
#>  [883] "CGCGTTTCAGTCTTCCstochastic_pseudotime"
#>  [884] "GAAACTCAGAAACCTAstochastic_pseudotime"
#>  [885] "GGAATAAGTAAGAGAGstochastic_pseudotime"
#>  [886] "TCTATTGTCCTATTCAstochastic_pseudotime"
#>  [887] "TATCTCACATTGTGCAstochastic_pseudotime"
#>  [888] "TGGCTGGGTACCCAATstochastic_pseudotime"
#>  [889] "AGCATACCACAGCCCAstochastic_pseudotime"
#>  [890] "TTGGAACGTATATGAGstochastic_pseudotime"
#>  [891] "CGTCACTTCTAACTCTstochastic_pseudotime"
#>  [892] "CGTCCATGTCTAACGTstochastic_pseudotime"
#>  [893] "AGGGAGTTCCTCAATTstochastic_pseudotime"
#>  [894] "ACTATCTGTCAGGACAstochastic_pseudotime"
#>  [895] "TCACAAGAGTGACATAstochastic_pseudotime"
#>  [896] "GCATGTAAGCGTGTCCstochastic_pseudotime"
#>  [897] "GGACAGAAGGCCGAATstochastic_pseudotime"
#>  [898] "TCAGGTATCTGTCCGTstochastic_pseudotime"
#>  [899] "CTTGGCTCATCCCACTstochastic_pseudotime"
#>  [900] "GTCACGGAGCTATGCTstochastic_pseudotime"
#>  [901] "GTTTCTACAGGACCCTstochastic_pseudotime"
#>  [902] "GATCTAGAGCTAGCCCstochastic_pseudotime"
#>  [903] "TGTGGTAGTGCCTGCAstochastic_pseudotime"
#>  [904] "ATCTGCCCAATCCGATstochastic_pseudotime"
#>  [905] "GTGGGTCCATAAAGGTstochastic_pseudotime"
#>  [906] "AGAGCGACAGTAAGATstochastic_pseudotime"
#>  [907] "CGGACGTAGACACGACstochastic_pseudotime"
#>  [908] "AAAGATGAGTCATGCTstochastic_pseudotime"
#>  [909] "GATTCAGTCGTCCAGGstochastic_pseudotime"
#>  [910] "TACGGATAGGCACATGstochastic_pseudotime"
#>  [911] "AGGGATGTCGCATGATstochastic_pseudotime"
#>  [912] "GTATTCTTCTCTAAGGstochastic_pseudotime"
#>  [913] "CATCGAAGTCACACGCstochastic_pseudotime"
#>  [914] "GTCGTAAAGTTTAGGAstochastic_pseudotime"
#>  [915] "CTTAGGATCAACGGGAstochastic_pseudotime"
#>  [916] "GCCTCTAAGCTAGTCTstochastic_pseudotime"
#>  [917] "GTTCATTTCCGAACGCstochastic_pseudotime"
#>  [918] "CTAAGACTCTGCGGCAstochastic_pseudotime"
#>  [919] "AGCTCTCGTTCAGTACstochastic_pseudotime"
#>  [920] "GTACTTTGTACCGTATstochastic_pseudotime"
#>  [921] "CTGAAACGTATGGTTCstochastic_pseudotime"
#>  [922] "TGAGAGGAGGCGATACstochastic_pseudotime"
#>  [923] "ATCATCTGTGTTTGTGstochastic_pseudotime"
#>  [924] "TTGGCAAAGAGTAAGGstochastic_pseudotime"
#>  [925] "CTAGTGATCCCAGGTGstochastic_pseudotime"
#>  [926] "CCTAGCTGTGGCCCTAstochastic_pseudotime"
#>  [927] "CCTTCGAAGATAGCATstochastic_pseudotime"
#>  [928] "GCATGATCAAGCGAGTstochastic_pseudotime"
#>  [929] "AGCTTGAGTGCAGTAGstochastic_pseudotime"
#>  [930] "CGTCAGGTCCGCAAGCstochastic_pseudotime"
#>  [931] "CCGGTAGTCGGCTTGGstochastic_pseudotime"
#>  [932] "CAAGGCCTCTTACCTAstochastic_pseudotime"
#>  [933] "TTTCCTCCAGGATCGAstochastic_pseudotime"
#>  [934] "CAGATCAAGCTGTTCAstochastic_pseudotime"
#>  [935] "TCAGATGGTACACCGCstochastic_pseudotime"
#>  [936] "CATCGAAGTACGAAATstochastic_pseudotime"
#>  [937] "TACGGGCTCTTCGGTCstochastic_pseudotime"
#>  [938] "TGATTTCGTGAGGGTTstochastic_pseudotime"
#>  [939] "CCTATTACAGACAAGCstochastic_pseudotime"
#>  [940] "GTCGGGTCAGGCTGAAstochastic_pseudotime"
#>  [941] "CCCAGTTAGCATCATCstochastic_pseudotime"
#>  [942] "GGCTCGAAGTGACTCTstochastic_pseudotime"
#>  [943] "TTCGGTCAGGGCATGTstochastic_pseudotime"
#>  [944] "AGTGGGACAGGTCGTCstochastic_pseudotime"
#>  [945] "GTTTCTAAGAAGATTCstochastic_pseudotime"
#>  [946] "CGACCTTAGGCTAGACstochastic_pseudotime"
#>  [947] "TGAGAGGAGATCCCATstochastic_pseudotime"
#>  [948] "CTAGTGATCAAAGACAstochastic_pseudotime"
#>  [949] "GGACGTCGTCTAGGTTstochastic_pseudotime"
#>  [950] "CACATTTTCGAGGTAGstochastic_pseudotime"
#>  [951] "TTCTACATCGCTTAGAstochastic_pseudotime"
#>  [952] "TCGCGAGAGTCCGTATstochastic_pseudotime"
#>  [953] "TACTCATCAGATCCATstochastic_pseudotime"
#>  [954] "CAACCAACAACCGCCAstochastic_pseudotime"
#>  [955] "ACGAGGATCATCACCCstochastic_pseudotime"
#>  [956] "GCATGTAAGAATCTCCstochastic_pseudotime"
#>  [957] "AGGCCGTCACGGTGTCstochastic_pseudotime"
#>  [958] "CACTCCATCAACACACstochastic_pseudotime"
#>  [959] "CTTGGCTGTTAGAACAstochastic_pseudotime"
#>  [960] "CCCATACAGACTTTCGstochastic_pseudotime"
#>  [961] "CCGTGGACATAGAAACstochastic_pseudotime"
#>  [962] "TCTCTAAAGCTGGAACstochastic_pseudotime"
#>  [963] "GTGCTTCTCACGACTAstochastic_pseudotime"
#>  [964] "GGCTCGACAGCTTCGGstochastic_pseudotime"
#>  [965] "GTCTCGTTCTGCGTAAstochastic_pseudotime"
#>  [966] "GCTCCTAAGGGCATGTstochastic_pseudotime"
#>  [967] "GTGGGTCCATCGATTGstochastic_pseudotime"
#>  [968] "CCACGGAGTCTAGTCAstochastic_pseudotime"
#>  [969] "CTAACTTAGATCACGGstochastic_pseudotime"
#>  [970] "TGTCCCAGTCCGTGACstochastic_pseudotime"
#>  [971] "ATCCACCAGAGCCTAGstochastic_pseudotime"
#>  [972] "GTAGGCCGTTACCAGTstochastic_pseudotime"
#>  [973] "TCTGGAACACCTATCCstochastic_pseudotime"
#>  [974] "ACGTCAAAGGTGTTAAstochastic_pseudotime"
#>  [975] "ACTTGTTCATGGTCATstochastic_pseudotime"
#>  [976] "GGACAGATCAGGCGAAstochastic_pseudotime"
#>  [977] "TACTTGTAGACAAAGGstochastic_pseudotime"
#>  [978] "AAATGCCAGACTAAGTstochastic_pseudotime"
#>  [979] "AAGTCTGGTTGGGACAstochastic_pseudotime"
#>  [980] "AGACGTTCAATGGTCTstochastic_pseudotime"
#>  [981] "CTACACCCATCAGTCAstochastic_pseudotime"
#>  [982] "CGTTCTGTCGATAGAAstochastic_pseudotime"
#>  [983] "GACCAATAGAATGTGTstochastic_pseudotime"
#>  [984] "TGACTTTAGCCAGAACstochastic_pseudotime"
#>  [985] "TCACAAGCAGACTCGCstochastic_pseudotime"
#>  [986] "CTTTGCGTCCACGTGGstochastic_pseudotime"
#>  [987] "CTCTGGTCAAGACACGstochastic_pseudotime"
#>  [988] "AGCATACGTCCGCTGAstochastic_pseudotime"
#>  [989] "GCTGCGATCGGCTACGstochastic_pseudotime"
#>  [990] "CCACCTACAGCTCGACstochastic_pseudotime"
#>  [991] "GGCAATTGTAGAGCTGstochastic_pseudotime"
#>  [992] "ACTTACTCACATGTGTstochastic_pseudotime"
#>  [993] "ACGAGCCCATGAGCGAstochastic_pseudotime"
#>  [994] "CCCAGTTGTGCGCTTGstochastic_pseudotime"
#>  [995] "GACGTTACATCTCGCTstochastic_pseudotime"
#>  [996] "AAACGGGGTACAGTTCstochastic_pseudotime"
#>  [997] "GCGCGATGTCCGAATTstochastic_pseudotime"
#>  [998] "TAGTTGGTCTTGCAAGstochastic_pseudotime"
#>  [999] "TGGCCAGCAAGGACTGstochastic_pseudotime"
#> [1000] "TGTGGTACACATCCGGstochastic_pseudotime"
#> 
#> 
#> $feature_split
#> NULL
#> 
#> $cell_metadata
#>                             cells Pseudotime_assign stochastic_pseudotime
#> AAACCTGAGCCTTGAT AAACCTGAGCCTTGAT       0.193109800           0.193109800
#> AAACCTGGTAAGTGGC AAACCTGGTAAGTGGC       0.672484152           0.672484152
#> AAACGGGAGATATGGT AAACGGGAGATATGGT       0.199749366           0.199749366
#> AAACGGGCAAAGAATC AAACGGGCAAAGAATC       0.822802141           0.822802141
#> AAACGGGGTACAGTTC AAACGGGGTACAGTTC       0.993088334           0.993088334
#> AAACGGGTCAGCTCTC AAACGGGTCAGCTCTC       0.217492023           0.217492023
#> AAACGGGTCGCATGGC AAACGGGTCGCATGGC       0.834355274           0.834355274
#> AAAGATGAGTCATGCT AAAGATGAGTCATGCT       0.852785605           0.852785605
#> AAAGATGAGTCTTGCA AAAGATGAGTCTTGCA       0.842462818           0.842462818
#> AAAGCAAAGCTACCTA AAAGCAAAGCTACCTA       0.752758865           0.752758865
#> AAAGCAACACATGGGA AAAGCAACACATGGGA       0.193369909           0.193369909
#> AAAGCAAGTGCAGTAG AAAGCAAGTGCAGTAG       0.806681518           0.806681518
#> AAATGCCAGACTAAGT AAATGCCAGACTAAGT       0.926839197           0.926839197
#> AAATGCCCAATGGAAT AAATGCCCAATGGAAT       0.828659058           0.828659058
#> AAATGCCCACCAACCG AAATGCCCACCAACCG       0.584617356           0.584617356
#> AAATGCCCATGCTAGT AAATGCCCATGCTAGT       0.672874689           0.672874689
#> AAATGCCGTGATAAGT AAATGCCGTGATAAGT       0.835917534           0.835917534
#> AACACGTGTCGACTGC AACACGTGTCGACTGC       0.766611002           0.766611002
#> AACCATGCACAGGAGT AACCATGCACAGGAGT       0.453953386           0.453953386
#> AACCATGGTTCGGCAC AACCATGGTTCGGCAC       0.566112628           0.566112628
#> AACCATGGTTGAACTC AACCATGGTTGAACTC       0.789448085           0.789448085
#> AACGTTGGTGCCTTGG AACGTTGGTGCCTTGG       0.197966982           0.197966982
#> AACTCAGGTGCCTGTG AACTCAGGTGCCTGTG       0.369067636           0.369067636
#> AACTCAGGTTGGAGGT AACTCAGGTTGGAGGT       0.025372223           0.025372223
#> AACTCAGTCCTCGCAT AACTCAGTCCTCGCAT       0.583781259           0.583781259
#> AACTCAGTCTGGTTCC AACTCAGTCTGGTTCC       0.055718604           0.055718604
#> AACTCCCCAGTTCATG AACTCCCCAGTTCATG       0.817829897           0.817829897
#> AACTCCCTCGAGAGCA AACTCCCTCGAGAGCA       0.201304618           0.201304618
#> AACTCTTCATTTGCCC AACTCTTCATTTGCCC       0.825641721           0.825641721
#> AACTGGTGTGCCTGTG AACTGGTGTGCCTGTG       0.040091568           0.040091568
#> AACTGGTTCTTACCGC AACTGGTTCTTACCGC       0.190505970           0.190505970
#> AACTTTCCACAGTCGC AACTTTCCACAGTCGC       0.361048276           0.361048276
#> AACTTTCCACTGCCAG AACTTTCCACTGCCAG       0.234487320           0.234487320
#> AACTTTCGTGACCAAG AACTTTCGTGACCAAG       0.212534168           0.212534168
#> AACTTTCTCCAGATCA AACTTTCTCCAGATCA       0.822704390           0.822704390
#> AAGACCTAGACCTTTG AAGACCTAGACCTTTG       0.824112047           0.824112047
#> AAGACCTAGTAGGCCA AAGACCTAGTAGGCCA       0.628515230           0.628515230
#> AAGACCTCAGGTTTCA AAGACCTCAGGTTTCA       0.070582897           0.070582897
#> AAGCCGCAGGCAGGTT AAGCCGCAGGCAGGTT       0.048070035           0.048070035
#> AAGCCGCCAATGGACG AAGCCGCCAATGGACG       0.681090206           0.681090206
#> AAGGAGCAGTACATGA AAGGAGCAGTACATGA       0.808843799           0.808843799
#> AAGGCAGAGCTGAAAT AAGGCAGAGCTGAAAT       0.174685862           0.174685862
#> AAGGCAGTCCGAGCCA AAGGCAGTCCGAGCCA       0.221495884           0.221495884
#> AAGGTTCAGATAGTCA AAGGTTCAGATAGTCA       0.809822875           0.809822875
#> AAGGTTCAGGAGCGAG AAGGTTCAGGAGCGAG       0.826637205           0.826637205
#> AAGGTTCGTCAACTGT AAGGTTCGTCAACTGT       0.772310943           0.772310943
#> AAGTCTGGTTGGGACA AAGTCTGGTTGGGACA       0.929562362           0.929562362
#> AAGTCTGTCCGAGCCA AAGTCTGTCCGAGCCA       0.403139286           0.403139286
#> AATCCAGAGGGTTCCC AATCCAGAGGGTTCCC       0.767757193           0.767757193
#> AATCCAGTCAGGTTCA AATCCAGTCAGGTTCA       0.748447330           0.748447330
#> AATCCAGTCGGTCCGA AATCCAGTCGGTCCGA       0.740378589           0.740378589
#> AATCGGTCATATGGTC AATCGGTCATATGGTC       0.529235912           0.529235912
#> AATCGGTTCCAAACTG AATCGGTTCCAAACTG       0.588842677           0.588842677
#> AATCGGTTCGTATCAG AATCGGTTCGTATCAG       0.829621964           0.829621964
#> ACACCAAGTATAGGTA ACACCAAGTATAGGTA       0.823348452           0.823348452
#> ACACCCTAGGACAGCT ACACCCTAGGACAGCT       0.610480715           0.610480715
#> ACACCCTTCCCAGGTG ACACCCTTCCCAGGTG       0.791814426           0.791814426
#> ACACCGGCATCACGTA ACACCGGCATCACGTA       0.195044576           0.195044576
#> ACACCGGCATGTCCTC ACACCGGCATGTCCTC       0.183745900           0.183745900
#> ACACCGGTCATGTCTT ACACCGGTCATGTCTT       0.808856775           0.808856775
#> ACACTGAGTCCATGAT ACACTGAGTCCATGAT       0.684773147           0.684773147
#> ACAGCCGAGAGAACAG ACAGCCGAGAGAACAG       0.789561590           0.789561590
#> ACAGCCGAGAGGTACC ACAGCCGAGAGGTACC       0.829748774           0.829748774
#> ACAGCCGAGTGGTCCC ACAGCCGAGTGGTCCC       0.086771749           0.086771749
#> ACAGCCGCATCGTCGG ACAGCCGCATCGTCGG       0.551080671           0.551080671
#> ACAGCTAAGTCGTTTG ACAGCTAAGTCGTTTG       0.820804668           0.820804668
#> ACAGCTATCAAGGCTT ACAGCTATCAAGGCTT       0.829235052           0.829235052
#> ACAGCTATCAGCATGT ACAGCTATCAGCATGT       0.129647842           0.129647842
#> ACAGCTATCTGGGCCA ACAGCTATCTGGGCCA       0.057342104           0.057342104
#> ACATACGAGTAGGCCA ACATACGAGTAGGCCA       0.206595911           0.206595911
#> ACATACGGTCTCTCTG ACATACGGTCTCTCTG       0.754237922           0.754237922
#> ACATCAGAGAGACTTA ACATCAGAGAGACTTA       0.753214030           0.753214030
#> ACATCAGTCCTTGGTC ACATCAGTCCTTGGTC       0.407011189           0.407011189
#> ACCAGTAAGTGATCGG ACCAGTAAGTGATCGG       0.819686580           0.819686580
#> ACCAGTAGTTCGTTGA ACCAGTAGTTCGTTGA       0.814911233           0.814911233
#> ACCAGTATCGTTTAGG ACCAGTATCGTTTAGG       0.317036006           0.317036006
#> ACCGTAAGTTATCACG ACCGTAAGTTATCACG       0.069983952           0.069983952
#> ACCTTTATCCTGCAGG ACCTTTATCCTGCAGG       0.828336191           0.828336191
#> ACGAGCCCATGAGCGA ACGAGCCCATGAGCGA       0.971368609           0.971368609
#> ACGAGGAAGTTACCCA ACGAGGAAGTTACCCA       0.222219161           0.222219161
#> ACGAGGATCATCACCC ACGAGGATCATCACCC       0.874968240           0.874968240
#> ACGAGGATCGCCGTGA ACGAGGATCGCCGTGA       0.795761165           0.795761165
#> ACGATACAGATGAGAG ACGATACAGATGAGAG       0.826146010           0.826146010
#> ACGATACTCGCGCCAA ACGATACTCGCGCCAA       0.826317265           0.826317265
#> ACGATGTTCGCAGGCT ACGATGTTCGCAGGCT       0.805737133           0.805737133
#> ACGCAGCCATTCTCAT ACGCAGCCATTCTCAT       0.823448438           0.823448438
#> ACGCAGCGTAATCACC ACGCAGCGTAATCACC       0.567099989           0.567099989
#> ACGCCAGGTTACAGAA ACGCCAGGTTACAGAA       0.832138922           0.832138922
#> ACGCCAGTCCCAAGTA ACGCCAGTCCCAAGTA       0.471753369           0.471753369
#> ACGCCGAAGGGTGTTG ACGCCGAAGGGTGTTG       0.803356733           0.803356733
#> ACGCCGAGTTTGCATG ACGCCGAGTTTGCATG       0.833999814           0.833999814
#> ACGGCCACACGTGAGA ACGGCCACACGTGAGA       0.195806143           0.195806143
#> ACGGGCTGTGCGATAG ACGGGCTGTGCGATAG       0.069462274           0.069462274
#> ACGTCAAAGGTGTTAA ACGTCAAAGGTGTTAA       0.921117564           0.921117564
#> ACGTCAACACTGCCAG ACGTCAACACTGCCAG       0.160837557           0.160837557
#> ACGTCAACAGACAAAT ACGTCAACAGACAAAT       0.078366646           0.078366646
#> ACTATCTGTCAGGACA ACTATCTGTCAGGACA       0.848700580           0.848700580
#> ACTATCTTCAAGATCC ACTATCTTCAAGATCC       0.805405010           0.805405010
#> ACTGAACAGCCTCGTG ACTGAACAGCCTCGTG       0.194666436           0.194666436
#> ACTGAGTAGTGACTCT ACTGAGTAGTGACTCT       0.212036994           0.212036994
#> ACTGATGCATCAGTAC ACTGATGCATCAGTAC       0.571499166           0.571499166
#> ACTGATGCATGGTCTA ACTGATGCATGGTCTA       0.832292148           0.832292148
#> ACTGATGCATGTCTCC ACTGATGCATGTCTCC       0.727630150           0.727630150
#> ACTGATGGTGAGGGTT ACTGATGGTGAGGGTT       0.072933103           0.072933103
#> ACTGATGTCCCTAACC ACTGATGTCCCTAACC       0.173048170           0.173048170
#> ACTGCTCAGGTGTTAA ACTGCTCAGGTGTTAA       0.800745202           0.800745202
#> ACTGTCCCAGTTTACG ACTGTCCCAGTTTACG       0.700944895           0.700944895
#> ACTGTCCTCGCTAGCG ACTGTCCTCGCTAGCG       0.506457243           0.506457243
#> ACTTACTAGAGGGCTT ACTTACTAGAGGGCTT       0.302984864           0.302984864
#> ACTTACTAGCCGGTAA ACTTACTAGCCGGTAA       0.210992072           0.210992072
#> ACTTACTCAATCAGAA ACTTACTCAATCAGAA       0.365119055           0.365119055
#> ACTTACTCACATGTGT ACTTACTCACATGTGT       0.970686456           0.970686456
#> ACTTGTTAGTGAAGTT ACTTGTTAGTGAAGTT       0.820365054           0.820365054
#> ACTTGTTAGTTAAGTG ACTTGTTAGTTAAGTG       0.829829393           0.829829393
#> ACTTGTTCATGGTCAT ACTTGTTCATGGTCAT       0.921229647           0.921229647
#> ACTTGTTTCACCCGAG ACTTGTTTCACCCGAG       0.819365217           0.819365217
#> ACTTTCAAGGATATAC ACTTTCAAGGATATAC       0.199365043           0.199365043
#> ACTTTCACAACTGCTA ACTTTCACAACTGCTA       0.250226007           0.250226007
#> AGAATAGCAATCACAC AGAATAGCAATCACAC       0.822590165           0.822590165
#> AGAATAGCAGTACACT AGAATAGCAGTACACT       0.088540569           0.088540569
#> AGACGTTCAATGGTCT AGACGTTCAATGGTCT       0.930518092           0.930518092
#> AGACGTTTCTTACCTA AGACGTTTCTTACCTA       0.193660447           0.193660447
#> AGAGCGACAGTAAGAT AGAGCGACAGTAAGAT       0.852197508           0.852197508
#> AGAGCGATCTCTTGAT AGAGCGATCTCTTGAT       0.090744366           0.090744366
#> AGAGCTTAGCCGATTT AGAGCTTAGCCGATTT       0.219369860           0.219369860
#> AGAGCTTAGGTTACCT AGAGCTTAGGTTACCT       0.415186896           0.415186896
#> AGAGCTTGTTAAAGAC AGAGCTTGTTAAAGAC       0.075239020           0.075239020
#> AGAGTGGAGGCGATAC AGAGTGGAGGCGATAC       0.822783552           0.822783552
#> AGAGTGGCACCTCGGA AGAGTGGCACCTCGGA       0.726633971           0.726633971
#> AGAGTGGGTAGCACGA AGAGTGGGTAGCACGA       0.203902489           0.203902489
#> AGAGTGGGTCACCTAA AGAGTGGGTCACCTAA       0.229264788           0.229264788
#> AGATCTGTCCGTCAAA AGATCTGTCCGTCAAA       0.761233706           0.761233706
#> AGATTGCAGCGTGTCC AGATTGCAGCGTGTCC       0.824867337           0.824867337
#> AGATTGCAGCTGTCTA AGATTGCAGCTGTCTA       0.010182326           0.010182326
#> AGATTGCAGTTCCACA AGATTGCAGTTCCACA       0.451024851           0.451024851
#> AGATTGCCACCAGGCT AGATTGCCACCAGGCT       0.178250973           0.178250973
#> AGCAGCCGTCGCATCG AGCAGCCGTCGCATCG       0.652205253           0.652205253
#> AGCAGCCGTGTAACGG AGCAGCCGTGTAACGG       0.626401469           0.626401469
#> AGCAGCCTCTTTAGTC AGCAGCCTCTTTAGTC       0.196529127           0.196529127
#> AGCATACCACAGCCCA AGCATACCACAGCCCA       0.847377157           0.847377157
#> AGCATACGTCCGCTGA AGCATACGTCCGCTGA       0.962124525           0.962124525
#> AGCGGTCGTACCGAGA AGCGGTCGTACCGAGA       0.258112889           0.258112889
#> AGCGGTCGTCGAATCT AGCGGTCGTCGAATCT       0.418068919           0.418068919
#> AGCGGTCGTGTATGGG AGCGGTCGTGTATGGG       0.152865591           0.152865591
#> AGCGTCGCACAGGCCT AGCGTCGCACAGGCCT       0.447052297           0.447052297
#> AGCTCCTCAGTTTACG AGCTCCTCAGTTTACG       0.781031071           0.781031071
#> AGCTCCTGTGTGCCTG AGCTCCTGTGTGCCTG       0.209598845           0.209598845
#> AGCTCCTTCAATACCG AGCTCCTTCAATACCG       0.663681073           0.663681073
#> AGCTCCTTCGCCCTTA AGCTCCTTCGCCCTTA       0.787080892           0.787080892
#> AGCTCTCAGACATAAC AGCTCTCAGACATAAC       0.800977639           0.800977639
#> AGCTCTCGTTCAGTAC AGCTCTCGTTCAGTAC       0.857113580           0.857113580
#> AGCTTGAAGACACGAC AGCTTGAAGACACGAC       0.688852198           0.688852198
#> AGCTTGAAGAGTGACC AGCTTGAAGAGTGACC       0.840249528           0.840249528
#> AGCTTGACAATAGAGT AGCTTGACAATAGAGT       0.207138613           0.207138613
#> AGCTTGACAGTAAGAT AGCTTGACAGTAAGAT       0.121445498           0.121445498
#> AGCTTGACATCTATGG AGCTTGACATCTATGG       0.827224451           0.827224451
#> AGCTTGAGTGCAGTAG AGCTTGAGTGCAGTAG       0.860042937           0.860042937
#> AGGCCACAGCTGTCTA AGGCCACAGCTGTCTA       0.211561901           0.211561901
#> AGGCCACCATGACATC AGGCCACCATGACATC       0.808754740           0.808754740
#> AGGCCGTCACGACGAA AGGCCGTCACGACGAA       0.427030181           0.427030181
#> AGGCCGTCACGGTGTC AGGCCGTCACGGTGTC       0.879187591           0.879187591
#> AGGCCGTGTACCAGTT AGGCCGTGTACCAGTT       0.057542427           0.057542427
#> AGGGAGTAGGGTATCG AGGGAGTAGGGTATCG       0.074512563           0.074512563
#> AGGGAGTAGTGTCTCA AGGGAGTAGTGTCTCA       0.151211294           0.151211294
#> AGGGAGTTCCTCAATT AGGGAGTTCCTCAATT       0.848313141           0.848313141
#> AGGGAGTTCGGAGGTA AGGGAGTTCGGAGGTA       0.815646432           0.815646432
#> AGGGATGGTCGACTGC AGGGATGGTCGACTGC       0.819633217           0.819633217
#> AGGGATGGTTAGGGTG AGGGATGGTTAGGGTG       0.436155506           0.436155506
#> AGGGATGTCGCATGAT AGGGATGTCGCATGAT       0.853926163           0.853926163
#> AGGGATGTCGCATGGC AGGGATGTCGCATGGC       0.648289045           0.648289045
#> AGGGTGAGTGCGAAAC AGGGTGAGTGCGAAAC       0.834259790           0.834259790
#> AGGTCATCACTTCTGC AGGTCATCACTTCTGC       0.822104196           0.822104196
#> AGGTCATGTCTCCCTA AGGTCATGTCTCCCTA       0.817153674           0.817153674
#> AGGTCATGTTTAGGAA AGGTCATGTTTAGGAA       0.737242866           0.737242866
#> AGTAGTCCACATGACT AGTAGTCCACATGACT       0.204868994           0.204868994
#> AGTAGTCGTGAGGGAG AGTAGTCGTGAGGGAG       0.812334053           0.812334053
#> AGTAGTCGTTCTGTTT AGTAGTCGTTCTGTTT       0.203908820           0.203908820
#> AGTCTTTCACCAGGTC AGTCTTTCACCAGGTC       0.021267612           0.021267612
#> AGTCTTTCACGGACAA AGTCTTTCACGGACAA       0.823017775           0.823017775
#> AGTCTTTGTTACGGAG AGTCTTTGTTACGGAG       0.833305565           0.833305565
#> AGTGAGGGTTACGGAG AGTGAGGGTTACGGAG       0.023426298           0.023426298
#> AGTGAGGTCTCCAACC AGTGAGGTCTCCAACC       0.837749619           0.837749619
#> AGTGGGACAGGTCGTC AGTGGGACAGGTCGTC       0.865073663           0.865073663
#> AGTGGGAGTGGCAAAC AGTGGGAGTGGCAAAC       0.808592870           0.808592870
#> AGTGTCATCCCTTGTG AGTGTCATCCCTTGTG       0.810060061           0.810060061
#> AGTTGGTAGAGCCTAG AGTTGGTAGAGCCTAG       0.828898362           0.828898362
#> AGTTGGTGTTGTCGCG AGTTGGTGTTGTCGCG       0.710593010           0.710593010
#> AGTTGGTTCACATAGC AGTTGGTTCACATAGC       0.580112046           0.580112046
#> ATAAGAGCATCGATTG ATAAGAGCATCGATTG       0.061202682           0.061202682
#> ATAGACCAGCGTGTCC ATAGACCAGCGTGTCC       0.822273452           0.822273452
#> ATAGACCAGTGATCGG ATAGACCAGTGATCGG       0.817308694           0.817308694
#> ATAGACCTCTAACTGG ATAGACCTCTAACTGG       0.199467202           0.199467202
#> ATCACGAGTGCACCAC ATCACGAGTGCACCAC       0.810395597           0.810395597
#> ATCATCTGTGTTTGTG ATCATCTGTGTTTGTG       0.857621626           0.857621626
#> ATCATCTTCTGTCAAG ATCATCTTCTGTCAAG       0.017257097           0.017257097
#> ATCATGGAGGAATGGA ATCATGGAGGAATGGA       0.819596622           0.819596622
#> ATCATGGGTATTCGTG ATCATGGGTATTCGTG       0.834965551           0.834965551
#> ATCATGGGTCCGCTGA ATCATGGGTCCGCTGA       0.556066906           0.556066906
#> ATCATGGTCCAGAGGA ATCATGGTCCAGAGGA       0.177783799           0.177783799
#> ATCCACCAGAGCCTAG ATCCACCAGAGCCTAG       0.908711310           0.908711310
#> ATCCACCCAGCTCGCA ATCCACCCAGCTCGCA       0.828525786           0.828525786
#> ATCCACCGTACCATCA ATCCACCGTACCATCA       0.118444318           0.118444318
#> ATCCGAAGTACACCGC ATCCGAAGTACACCGC       0.726853482           0.726853482
#> ATCCGAAGTGACTCAT ATCCGAAGTGACTCAT       0.217032938           0.217032938
#> ATCTGCCCAATCCGAT ATCTGCCCAATCCGAT       0.852070443           0.852070443
#> ATCTGCCCATGTTGAC ATCTGCCCATGTTGAC       0.800412446           0.800412446
#> ATGAGGGTCTTGACGA ATGAGGGTCTTGACGA       0.797325075           0.797325075
#> ATGTGTGCACAAGTAA ATGTGTGCACAAGTAA       0.428294823           0.428294823
#> ATGTGTGTCAAAGACA ATGTGTGTCAAAGACA       0.454539445           0.454539445
#> ATTACTCAGTTGCAGG ATTACTCAGTTGCAGG       0.022044234           0.022044234
#> ATTACTCTCAGGCGAA ATTACTCTCAGGCGAA       0.029460475           0.029460475
#> ATTACTCTCTGCCAGG ATTACTCTCTGCCAGG       0.828102911           0.828102911
#> ATTATCCAGAAACCTA ATTATCCAGAAACCTA       0.253109144           0.253109144
#> ATTATCCAGACCCACC ATTATCCAGACCCACC       0.262607798           0.262607798
#> ATTATCCAGTGAACAT ATTATCCAGTGAACAT       0.095041495           0.095041495
#> ATTGGACAGTAATCCC ATTGGACAGTAATCCC       0.220647047           0.220647047
#> ATTGGTGTCCAGTATG ATTGGTGTCCAGTATG       0.187598052           0.187598052
#> ATTTCTGGTAGTGAAT ATTTCTGGTAGTGAAT       0.187897346           0.187897346
#> CAACCAAAGCGTTCCG CAACCAAAGCGTTCCG       0.774135836           0.774135836
#> CAACCAACAACCGCCA CAACCAACAACCGCCA       0.873139249           0.873139249
#> CAACCAACACCTCGTT CAACCAACACCTCGTT       0.502415557           0.502415557
#> CAACCAATCGGGAGTA CAACCAATCGGGAGTA       0.126725915           0.126725915
#> CAAGAAAAGAGTACCG CAAGAAAAGAGTACCG       0.674623524           0.674623524
#> CAAGAAATCGTACCGG CAAGAAATCGTACCGG       0.195166976           0.195166976
#> CAAGGCCAGCTGTTCA CAAGGCCAGCTGTTCA       0.211913584           0.211913584
#> CAAGGCCAGTGTGGCA CAAGGCCAGTGTGGCA       0.198578959           0.198578959
#> CAAGGCCCACATTTCT CAAGGCCCACATTTCT       0.565313016           0.565313016
#> CAAGGCCTCAGTGCAT CAAGGCCTCAGTGCAT       0.821158232           0.821158232
#> CAAGGCCTCTCAACTT CAAGGCCTCTCAACTT       0.726867844           0.726867844
#> CAAGGCCTCTTACCTA CAAGGCCTCTTACCTA       0.860586048           0.860586048
#> CAAGTTGGTGTGCCTG CAAGTTGGTGTGCCTG       0.801612218           0.801612218
#> CAAGTTGTCTAACTGG CAAGTTGTCTAACTGG       0.746284481           0.746284481
#> CACAAACAGTAGCGGT CACAAACAGTAGCGGT       0.567547289           0.567547289
#> CACAAACCACGAAAGC CACAAACCACGAAAGC       0.831230089           0.831230089
#> CACACAACAGCCAGAA CACACAACAGCCAGAA       0.132635218           0.132635218
#> CACACCTGTACATGTC CACACCTGTACATGTC       0.807311254           0.807311254
#> CACACCTGTTATCCGA CACACCTGTTATCCGA       0.770756978           0.770756978
#> CACACTCAGATCCCAT CACACTCAGATCCCAT       0.072081902           0.072081902
#> CACACTCTCATGTGGT CACACTCTCATGTGGT       0.096907048           0.096907048
#> CACAGGCCAAGGCTCC CACAGGCCAAGGCTCC       0.706094404           0.706094404
#> CACAGGCCACGTCAGC CACAGGCCACGTCAGC       0.835760937           0.835760937
#> CACAGTAGTGCACTTA CACAGTAGTGCACTTA       0.814163804           0.814163804
#> CACAGTAGTTCCTCCA CACAGTAGTTCCTCCA       0.190850316           0.190850316
#> CACAGTATCCTCCTAG CACAGTATCCTCCTAG       0.197936436           0.197936436
#> CACATAGAGCGTGTCC CACATAGAGCGTGTCC       0.074194965           0.074194965
#> CACATTTGTTTCGCTC CACATTTGTTTCGCTC       0.186276653           0.186276653
#> CACATTTTCGAGGTAG CACATTTTCGAGGTAG       0.867096208           0.867096208
#> CACATTTTCGCAGGCT CACATTTTCGCAGGCT       0.211625741           0.211625741
#> CACCACTAGCGATCCC CACCACTAGCGATCCC       0.820843023           0.820843023
#> CACCACTGTGTGTGCC CACCACTGTGTGTGCC       0.722490313           0.722490313
#> CACCAGGCATAAAGGT CACCAGGCATAAAGGT       0.697634315           0.697634315
#> CACCTTGAGAGTAATC CACCTTGAGAGTAATC       0.188848077           0.188848077
#> CACTCCACATCGTCGG CACTCCACATCGTCGG       0.830123887           0.830123887
#> CACTCCATCAACACAC CACTCCATCAACACAC       0.879274691           0.879274691
#> CAGAATCCATTGGGCC CAGAATCCATTGGGCC       0.827264487           0.827264487
#> CAGAATCGTCATTAGC CAGAATCGTCATTAGC       0.316604073           0.316604073
#> CAGAATCGTGTAAGTA CAGAATCGTGTAAGTA       0.202379191           0.202379191
#> CAGAGAGTCGGCCGAT CAGAGAGTCGGCCGAT       0.831058606           0.831058606
#> CAGATCAAGCTGTTCA CAGATCAAGCTGTTCA       0.860838833           0.860838833
#> CAGCAGCAGACAAAGG CAGCAGCAGACAAAGG       0.064777090           0.064777090
#> CAGCAGCCACACTGCG CAGCAGCCACACTGCG       0.022216769           0.022216769
#> CAGCAGCGTTCACGGC CAGCAGCGTTCACGGC       0.179779873           0.179779873
#> CAGCCGAAGCGATATA CAGCCGAAGCGATATA       0.064612407           0.064612407
#> CAGCCGAAGTAGGTGC CAGCCGAAGTAGGTGC       0.020515596           0.020515596
#> CAGCCGATCCTATTCA CAGCCGATCCTATTCA       0.643951799           0.643951799
#> CAGCCGATCGTTTATC CAGCCGATCGTTTATC       0.604912498           0.604912498
#> CAGCGACCACAGGAGT CAGCGACCACAGGAGT       0.074825028           0.074825028
#> CAGCTAAAGTCCATAC CAGCTAAAGTCCATAC       0.074068167           0.074068167
#> CAGCTAACATGCCTAA CAGCTAACATGCCTAA       0.726813561           0.726813561
#> CAGCTAAGTACCGAGA CAGCTAAGTACCGAGA       0.572940442           0.572940442
#> CAGCTGGCAAGTTGTC CAGCTGGCAAGTTGTC       0.779684214           0.779684214
#> CAGCTGGGTCTAACGT CAGCTGGGTCTAACGT       0.773560306           0.773560306
#> CAGCTGGGTTCTCATT CAGCTGGGTTCTCATT       0.842412789           0.842412789
#> CAGGTGCAGCCAGGAT CAGGTGCAGCCAGGAT       0.815763085           0.815763085
#> CAGGTGCCATCAGTCA CAGGTGCCATCAGTCA       0.780446826           0.780446826
#> CAGTCCTAGGTGGGTT CAGTCCTAGGTGGGTT       0.831560187           0.831560187
#> CAGTCCTCAGGGCATA CAGTCCTCAGGGCATA       0.826859099           0.826859099
#> CATATGGCATGTTGAC CATATGGCATGTTGAC       0.136762778           0.136762778
#> CATATGGGTTACTGAC CATATGGGTTACTGAC       0.834185510           0.834185510
#> CATATGGTCACAGTAC CATATGGTCACAGTAC       0.727722276           0.727722276
#> CATATTCTCCTGTAGA CATATTCTCCTGTAGA       0.831175446           0.831175446
#> CATCAGACATGCATGT CATCAGACATGCATGT       0.131132914           0.131132914
#> CATCAGAGTCTCTTTA CATCAGAGTCTCTTTA       0.813425364           0.813425364
#> CATCAGATCATCGCTC CATCAGATCATCGCTC       0.792809103           0.792809103
#> CATCGAAGTACGAAAT CATCGAAGTACGAAAT       0.861554887           0.861554887
#> CATCGAAGTCACACGC CATCGAAGTCACACGC       0.854191585           0.854191585
#> CATCGGGCACGTTGGC CATCGGGCACGTTGGC       0.819086876           0.819086876
#> CATCGGGTCTTCGGTC CATCGGGTCTTCGGTC       0.496334917           0.496334917
#> CATGACACAACAACCT CATGACACAACAACCT       0.212349437           0.212349437
#> CATGACACATCGATTG CATGACACATCGATTG       0.354201888           0.354201888
#> CATGACAGTAGCTTGT CATGACAGTAGCTTGT       0.197123715           0.197123715
#> CATGACATCATGCATG CATGACATCATGCATG       0.816478369           0.816478369
#> CATGCCTCACCCATTC CATGCCTCACCCATTC       0.821121623           0.821121623
#> CATGGCGAGAGTCTGG CATGGCGAGAGTCTGG       0.467835807           0.467835807
#> CATTATCGTCATATCG CATTATCGTCATATCG       0.756352656           0.756352656
#> CCAATCCAGAGTAAGG CCAATCCAGAGTAAGG       0.250422389           0.250422389
#> CCACCTAAGACTTGAA CCACCTAAGACTTGAA       0.830326092           0.830326092
#> CCACCTACAGCTCGAC CCACCTACAGCTCGAC       0.964957297           0.964957297
#> CCACCTATCACCAGGC CCACCTATCACCAGGC       0.203883985           0.203883985
#> CCACCTATCTAAGCCA CCACCTATCTAAGCCA       0.816921160           0.816921160
#> CCACGGAAGATCCTGT CCACGGAAGATCCTGT       0.816768115           0.816768115
#> CCACGGAAGCAATCTC CCACGGAAGCAATCTC       0.821593942           0.821593942
#> CCACGGAAGGGTGTTG CCACGGAAGGGTGTTG       0.816477489           0.816477489
#> CCACGGAAGTCTCCTC CCACGGAAGTCTCCTC       0.784217729           0.784217729
#> CCACGGAGTCAACATC CCACGGAGTCAACATC       0.467095752           0.467095752
#> CCACGGAGTCTAGTCA CCACGGAGTCTAGTCA       0.898464524           0.898464524
#> CCACGGAGTTACGACT CCACGGAGTTACGACT       0.449594603           0.449594603
#> CCACTACAGAGCTGGT CCACTACAGAGCTGGT       0.792916028           0.792916028
#> CCACTACAGAGTACCG CCACTACAGAGTACCG       0.213116452           0.213116452
#> CCACTACCAACGCACC CCACTACCAACGCACC       0.831245837           0.831245837
#> CCAGCGACAGACGTAG CCAGCGACAGACGTAG       0.711189513           0.711189513
#> CCATGTCGTCGGCTCA CCATGTCGTCGGCTCA       0.567440167           0.567440167
#> CCATGTCTCCGCATAA CCATGTCTCCGCATAA       0.674129458           0.674129458
#> CCATTCGCATGCCTAA CCATTCGCATGCCTAA       0.010620761           0.010620761
#> CCATTCGGTAGGGACT CCATTCGGTAGGGACT       0.031485015           0.031485015
#> CCCAATCCAACTGGCC CCCAATCCAACTGGCC       0.829433972           0.829433972
#> CCCAGTTAGCATCATC CCCAGTTAGCATCATC       0.864202319           0.864202319
#> CCCAGTTGTCGTGGCT CCCAGTTGTCGTGGCT       0.838316957           0.838316957
#> CCCAGTTGTGCGCTTG CCCAGTTGTGCGCTTG       0.971587003           0.971587003
#> CCCAGTTTCAGCACAT CCCAGTTTCAGCACAT       0.814001062           0.814001062
#> CCCATACAGACTTTCG CCCATACAGACTTTCG       0.879783598           0.879783598
#> CCCATACCACATCCAA CCCATACCACATCCAA       0.142518980           0.142518980
#> CCCATACCATGCAACT CCCATACCATGCAACT       0.826949656           0.826949656
#> CCCATACTCGTGGACC CCCATACTCGTGGACC       0.730378466           0.730378466
#> CCGGGATAGATCTGAA CCGGGATAGATCTGAA       0.636534297           0.636534297
#> CCGGGATTCGACAGCC CCGGGATTCGACAGCC       0.821256841           0.821256841
#> CCGGTAGTCGGCTTGG CCGGTAGTCGGCTTGG       0.860502618           0.860502618
#> CCGTACTAGGAGTACC CCGTACTAGGAGTACC       0.456227033           0.456227033
#> CCGTACTAGTCCGTAT CCGTACTAGTCCGTAT       0.456590238           0.456590238
#> CCGTACTTCAATAAGG CCGTACTTCAATAAGG       0.030546057           0.030546057
#> CCGTGGACATAGAAAC CCGTGGACATAGAAAC       0.880001772           0.880001772
#> CCGTGGATCCAAGCCG CCGTGGATCCAAGCCG       0.116709345           0.116709345
#> CCGTGGATCTGTCTCG CCGTGGATCTGTCTCG       0.836507048           0.836507048
#> CCGTGGATCTTAGAGC CCGTGGATCTTAGAGC       0.677495379           0.677495379
#> CCGTTCAAGCTAAACA CCGTTCAAGCTAAACA       0.819480348           0.819480348
#> CCGTTCAAGGACCACA CCGTTCAAGGACCACA       0.586999192           0.586999192
#> CCGTTCACACATGACT CCGTTCACACATGACT       0.210670572           0.210670572
#> CCTACACAGATCCGAG CCTACACAGATCCGAG       0.202365511           0.202365511
#> CCTACACAGTACGTTC CCTACACAGTACGTTC       0.818735269           0.818735269
#> CCTACACTCCTTTACA CCTACACTCCTTTACA       0.112423434           0.112423434
#> CCTACACTCTTTACAC CCTACACTCTTTACAC       0.826579157           0.826579157
#> CCTACCACATCGATTG CCTACCACATCGATTG       0.626445302           0.626445302
#> CCTACCAGTTGCTCCT CCTACCAGTTGCTCCT       0.737609955           0.737609955
#> CCTACCATCTGCCCTA CCTACCATCTGCCCTA       0.840654317           0.840654317
#> CCTAGCTCAGACTCGC CCTAGCTCAGACTCGC       0.174537960           0.174537960
#> CCTAGCTCATACTACG CCTAGCTCATACTACG       0.071668140           0.071668140
#> CCTAGCTGTGGCCCTA CCTAGCTGTGGCCCTA       0.859000636           0.859000636
#> CCTAGCTTCGATGAGG CCTAGCTTCGATGAGG       0.063140488           0.063140488
#> CCTATTACAGACAAGC CCTATTACAGACAAGC       0.862785230           0.862785230
#> CCTATTATCATCGCTC CCTATTATCATCGCTC       0.560722184           0.560722184
#> CCTCAGTAGTGCGATG CCTCAGTAGTGCGATG       0.782363581           0.782363581
#> CCTCAGTCAGACACTT CCTCAGTCAGACACTT       0.829694740           0.829694740
#> CCTCAGTGTTACCAGT CCTCAGTGTTACCAGT       0.266080436           0.266080436
#> CCTCTGAGTATAAACG CCTCTGAGTATAAACG       0.128468766           0.128468766
#> CCTTACGAGATCACGG CCTTACGAGATCACGG       0.203181726           0.203181726
#> CCTTACGCACAACGCC CCTTACGCACAACGCC       0.779995621           0.779995621
#> CCTTCCCAGGCAGGTT CCTTCCCAGGCAGGTT       0.193209821           0.193209821
#> CCTTCCCCAAAGTGCG CCTTCCCCAAAGTGCG       0.197276708           0.197276708
#> CCTTCCCTCGTCGTTC CCTTCCCTCGTCGTTC       0.739297136           0.739297136
#> CCTTCGAAGATAGCAT CCTTCGAAGATAGCAT       0.859027906           0.859027906
#> CCTTCGACAGCTGTTA CCTTCGACAGCTGTTA       0.196286929           0.196286929
#> CCTTCGATCACAGGCC CCTTCGATCACAGGCC       0.825388605           0.825388605
#> CGAACATAGCCAGTAG CGAACATAGCCAGTAG       0.047450984           0.047450984
#> CGAACATCAGACGTAG CGAACATCAGACGTAG       0.195233549           0.195233549
#> CGAATGTAGTTCCACA CGAATGTAGTTCCACA       0.219430120           0.219430120
#> CGAATGTGTCTGGAGA CGAATGTGTCTGGAGA       0.789658733           0.789658733
#> CGAATGTTCACCGGGT CGAATGTTCACCGGGT       0.840497916           0.840497916
#> CGAATGTTCTTCATGT CGAATGTTCTTCATGT       0.824727843           0.824727843
#> CGACCTTAGCGCTTAT CGACCTTAGCGCTTAT       0.244326762           0.244326762
#> CGACCTTAGGCTAGAC CGACCTTAGGCTAGAC       0.865636797           0.865636797
#> CGACTTCCAGCCTATA CGACTTCCAGCCTATA       0.221026909           0.221026909
#> CGACTTCCAGCTGCAC CGACTTCCAGCTGCAC       0.779346978           0.779346978
#> CGACTTCTCACGCGGT CGACTTCTCACGCGGT       0.835556506           0.835556506
#> CGAGAAGAGGGTTCCC CGAGAAGAGGGTTCCC       0.835758531           0.835758531
#> CGAGAAGCACGCATCG CGAGAAGCACGCATCG       0.427094005           0.427094005
#> CGAGAAGTCGTTACAG CGAGAAGTCGTTACAG       0.819830504           0.819830504
#> CGAGCACAGCTGAAAT CGAGCACAGCTGAAAT       0.157872317           0.157872317
#> CGAGCACAGGCAGGTT CGAGCACAGGCAGGTT       0.768616445           0.768616445
#> CGAGCACGTAGTAGTA CGAGCACGTAGTAGTA       0.800119443           0.800119443
#> CGAGCACTCTACCTGC CGAGCACTCTACCTGC       0.182711435           0.182711435
#> CGAGCACTCTTCGAGA CGAGCACTCTTCGAGA       0.211688027           0.211688027
#> CGATCGGAGCTGCAAG CGATCGGAGCTGCAAG       0.797426897           0.797426897
#> CGATCGGCACGCCAGT CGATCGGCACGCCAGT       0.231056120           0.231056120
#> CGATGGCGTTACGGAG CGATGGCGTTACGGAG       0.166614819           0.166614819
#> CGATGGCTCAGTGTTG CGATGGCTCAGTGTTG       0.825967710           0.825967710
#> CGATGGCTCTGTGCAA CGATGGCTCTGTGCAA       0.084399561           0.084399561
#> CGATGTAAGTTTAGGA CGATGTAAGTTTAGGA       0.807602153           0.807602153
#> CGATGTAGTACTCTCC CGATGTAGTACTCTCC       0.610712070           0.610712070
#> CGATTGAAGTCAAGGC CGATTGAAGTCAAGGC       0.673993252           0.673993252
#> CGCCAAGCAAAGTCAA CGCCAAGCAAAGTCAA       0.632421842           0.632421842
#> CGCCAAGCACCGGAAA CGCCAAGCACCGGAAA       0.058991011           0.058991011
#> CGCCAAGTCTAAGCCA CGCCAAGTCTAAGCCA       0.822403550           0.822403550
#> CGCGGTAAGAGCTGGT CGCGGTAAGAGCTGGT       0.831770849           0.831770849
#> CGCGGTAAGCTGCAAG CGCGGTAAGCTGCAAG       0.822509215           0.822509215
#> CGCGTTTCAGTCTTCC CGCGTTTCAGTCTTCC       0.842714430           0.842714430
#> CGCTGGAGTCCATCCT CGCTGGAGTCCATCCT       0.214872304           0.214872304
#> CGCTTCATCTCTGTCG CGCTTCATCTCTGTCG       0.821984431           0.821984431
#> CGGACACAGTGTCTCA CGGACACAGTGTCTCA       0.169876145           0.169876145
#> CGGACGTAGACACGAC CGGACGTAGACACGAC       0.852713465           0.852713465
#> CGGACTGAGGTAGCTG CGGACTGAGGTAGCTG       0.190166665           0.190166665
#> CGGACTGTCTTGAGAC CGGACTGTCTTGAGAC       0.837960260           0.837960260
#> CGGAGCTAGATAGTCA CGGAGCTAGATAGTCA       0.829113797           0.829113797
#> CGGAGCTCATTGGGCC CGGAGCTCATTGGGCC       0.148883556           0.148883556
#> CGGAGTCAGTTTAGGA CGGAGTCAGTTTAGGA       0.824808242           0.824808242
#> CGGCTAGGTACTCGCG CGGCTAGGTACTCGCG       0.513888557           0.513888557
#> CGGGTCAAGGCGATAC CGGGTCAAGGCGATAC       0.818298049           0.818298049
#> CGGGTCAAGGCTCAGA CGGGTCAAGGCTCAGA       0.175148448           0.175148448
#> CGGTTAAAGAGCTGCA CGGTTAAAGAGCTGCA       0.059616628           0.059616628
#> CGGTTAAGTAAACACA CGGTTAAGTAAACACA       0.057519815           0.057519815
#> CGGTTAAGTCTAGCCG CGGTTAAGTCTAGCCG       0.473222909           0.473222909
#> CGGTTAATCACCGGGT CGGTTAATCACCGGGT       0.777653917           0.777653917
#> CGGTTAATCGTCTGAA CGGTTAATCGTCTGAA       0.213729119           0.213729119
#> CGTAGCGTCGGAAATA CGTAGCGTCGGAAATA       0.787409632           0.787409632
#> CGTAGGCAGCGTAGTG CGTAGGCAGCGTAGTG       0.197961949           0.197961949
#> CGTCACTAGGGTGTGT CGTCACTAGGGTGTGT       0.208252981           0.208252981
#> CGTCACTCACGAAATA CGTCACTCACGAAATA       0.200490908           0.200490908
#> CGTCACTCATAGTAAG CGTCACTCATAGTAAG       0.196560414           0.196560414
#> CGTCACTTCACTTATC CGTCACTTCACTTATC       0.737753140           0.737753140
#> CGTCACTTCTAACTCT CGTCACTTCTAACTCT       0.848094620           0.848094620
#> CGTCAGGAGCCTATGT CGTCAGGAGCCTATGT       0.201460235           0.201460235
#> CGTCAGGCATGGGAAC CGTCAGGCATGGGAAC       0.202357593           0.202357593
#> CGTCAGGGTAGGCATG CGTCAGGGTAGGCATG       0.583327551           0.583327551
#> CGTCAGGTCCGCAAGC CGTCAGGTCCGCAAGC       0.860319824           0.860319824
#> CGTCAGGTCCTGCCAT CGTCAGGTCCTGCCAT       0.228462203           0.228462203
#> CGTCCATGTCTAACGT CGTCCATGTCTAACGT       0.848262768           0.848262768
#> CGTCCATTCGCGCCAA CGTCCATTCGCGCCAA       0.573553436           0.573553436
#> CGTCTACAGATCTGCT CGTCTACAGATCTGCT       0.158673489           0.158673489
#> CGTCTACAGGGAGTAA CGTCTACAGGGAGTAA       0.822610574           0.822610574
#> CGTGAGCTCAGCGATT CGTGAGCTCAGCGATT       0.555298443           0.555298443
#> CGTGAGCTCTTCTGGC CGTGAGCTCTTCTGGC       0.834165619           0.834165619
#> CGTGTAAGTCTTGATG CGTGTAAGTCTTGATG       0.215184401           0.215184401
#> CGTGTAATCATTATCC CGTGTAATCATTATCC       0.442673046           0.442673046
#> CGTGTAATCCATGCTC CGTGTAATCCATGCTC       0.052102056           0.052102056
#> CGTGTCTCAGAAGCAC CGTGTCTCAGAAGCAC       0.825567465           0.825567465
#> CGTGTCTGTAGCGTCC CGTGTCTGTAGCGTCC       0.808227504           0.808227504
#> CGTTAGAAGTATCTCG CGTTAGAAGTATCTCG       0.097681810           0.097681810
#> CGTTAGAGTGCAACGA CGTTAGAGTGCAACGA       0.205536405           0.205536405
#> CGTTCTGAGGACGAAA CGTTCTGAGGACGAAA       0.820520288           0.820520288
#> CGTTCTGTCACAATGC CGTTCTGTCACAATGC       0.641511674           0.641511674
#> CGTTCTGTCGATAGAA CGTTCTGTCGATAGAA       0.936887770           0.936887770
#> CGTTGGGCAGGCGATA CGTTGGGCAGGCGATA       0.784011671           0.784011671
#> CTAACTTAGAGGTACC CTAACTTAGAGGTACC       0.815686843           0.815686843
#> CTAACTTAGATCACGG CTAACTTAGATCACGG       0.905872786           0.905872786
#> CTAACTTCATCTGGTA CTAACTTCATCTGGTA       0.823502649           0.823502649
#> CTAACTTCATTGGGCC CTAACTTCATTGGGCC       0.192997551           0.192997551
#> CTAAGACAGCACGCCT CTAAGACAGCACGCCT       0.510419568           0.510419568
#> CTAAGACAGCCACGCT CTAAGACAGCCACGCT       0.196001990           0.196001990
#> CTAAGACGTATCAGTC CTAAGACGTATCAGTC       0.837871744           0.837871744
#> CTAAGACTCTGCGGCA CTAAGACTCTGCGGCA       0.856771643           0.856771643
#> CTAATGGCATGGTTGT CTAATGGCATGGTTGT       0.776796486           0.776796486
#> CTACACCCATCAGTCA CTACACCCATCAGTCA       0.934440179           0.934440179
#> CTACACCGTCGAATCT CTACACCGTCGAATCT       0.141270714           0.141270714
#> CTACACCGTCGCCATG CTACACCGTCGCCATG       0.038221965           0.038221965
#> CTACATTCACCGTTGG CTACATTCACCGTTGG       0.736860737           0.736860737
#> CTACATTGTTACTGAC CTACATTGTTACTGAC       0.778748482           0.778748482
#> CTACCCAGTGCGCTTG CTACCCAGTGCGCTTG       0.162349090           0.162349090
#> CTACGTCAGCTAGGCA CTACGTCAGCTAGGCA       0.820148534           0.820148534
#> CTAGAGTCACTGTTAG CTAGAGTCACTGTTAG       0.067680482           0.067680482
#> CTAGCCTGTCTCTTAT CTAGCCTGTCTCTTAT       0.146956574           0.146956574
#> CTAGCCTGTTCGGGCT CTAGCCTGTTCGGGCT       0.025797520           0.025797520
#> CTAGCCTTCAGCGATT CTAGCCTTCAGCGATT       0.187897721           0.187897721
#> CTAGTGATCAAAGACA CTAGTGATCAAAGACA       0.866332108           0.866332108
#> CTAGTGATCCCAGGTG CTAGTGATCCCAGGTG       0.858438395           0.858438395
#> CTCACACGTCACACGC CTCACACGTCACACGC       0.565650973           0.565650973
#> CTCAGAACACTAGTAC CTCAGAACACTAGTAC       0.315677890           0.315677890
#> CTCAGAACATAGGATA CTCAGAACATAGGATA       0.811410864           0.811410864
#> CTCAGAAGTCCTCTTG CTCAGAAGTCCTCTTG       0.190253562           0.190253562
#> CTCAGAATCAGGTTCA CTCAGAATCAGGTTCA       0.207772573           0.207772573
#> CTCATTACAAGCGCTC CTCATTACAAGCGCTC       0.043183874           0.043183874
#> CTCCTAGCACAACGTT CTCCTAGCACAACGTT       0.792008223           0.792008223
#> CTCGAAAGTAAGAGGA CTCGAAAGTAAGAGGA       0.329489132           0.329489132
#> CTCGAAAGTGGTGTAG CTCGAAAGTGGTGTAG       0.809115867           0.809115867
#> CTCGAAATCCTGCCAT CTCGAAATCCTGCCAT       0.817859998           0.817859998
#> CTCGAGGCAAGTCTGT CTCGAGGCAAGTCTGT       0.827256432           0.827256432
#> CTCGAGGGTAAGTAGT CTCGAGGGTAAGTAGT       0.838810947           0.838810947
#> CTCGAGGTCTGTTGAG CTCGAGGTCTGTTGAG       0.262000208           0.262000208
#> CTCGGAGGTCATTAGC CTCGGAGGTCATTAGC       0.793108754           0.793108754
#> CTCGGAGGTGAAGGCT CTCGGAGGTGAAGGCT       0.770974741           0.770974741
#> CTCGTACCAATGGATA CTCGTACCAATGGATA       0.205588710           0.205588710
#> CTCGTACTCTTTAGTC CTCGTACTCTTTAGTC       0.815329082           0.815329082
#> CTCGTCAAGAACAATC CTCGTCAAGAACAATC       0.840212360           0.840212360
#> CTCGTCACAACTGCTA CTCGTCACAACTGCTA       0.190777365           0.190777365
#> CTCGTCATCGGTCTAA CTCGTCATCGGTCTAA       0.199761795           0.199761795
#> CTCTAATAGTTATCGC CTCTAATAGTTATCGC       0.820120034           0.820120034
#> CTCTGGTAGCAATCTC CTCTGGTAGCAATCTC       0.492618649           0.492618649
#> CTCTGGTAGTCCCACG CTCTGGTAGTCCCACG       0.057080119           0.057080119
#> CTCTGGTCAACACCTA CTCTGGTCAACACCTA       0.217452547           0.217452547
#> CTCTGGTCAAGACACG CTCTGGTCAAGACACG       0.951662508           0.951662508
#> CTCTGGTCAAGTTGTC CTCTGGTCAAGTTGTC       0.792086057           0.792086057
#> CTCTGGTCAGTTCATG CTCTGGTCAGTTCATG       0.198914928           0.198914928
#> CTGAAACAGACTAGGC CTGAAACAGACTAGGC       0.183824479           0.183824479
#> CTGAAACGTAGCGCAA CTGAAACGTAGCGCAA       0.808160573           0.808160573
#> CTGAAACGTATGGTTC CTGAAACGTATGGTTC       0.857271665           0.857271665
#> CTGAAACTCGATAGAA CTGAAACTCGATAGAA       0.446004958           0.446004958
#> CTGAAACTCGCCGTGA CTGAAACTCGCCGTGA       0.162024398           0.162024398
#> CTGAAACTCTGTTTGT CTGAAACTCTGTTTGT       0.449680335           0.449680335
#> CTGAAGTGTTTGGGCC CTGAAGTGTTTGGGCC       0.830046164           0.830046164
#> CTGAAGTTCTGCGGCA CTGAAGTTCTGCGGCA       0.799539263           0.799539263
#> CTGATCCTCAGTTGAC CTGATCCTCAGTTGAC       0.178130367           0.178130367
#> CTGCCTAGTCTCCATC CTGCCTAGTCTCCATC       0.115692047           0.115692047
#> CTGCGGACAGACTCGC CTGCGGACAGACTCGC       0.566408815           0.566408815
#> CTGCGGACATTAGCCA CTGCGGACATTAGCCA       0.718554394           0.718554394
#> CTGCGGATCAGCATGT CTGCGGATCAGCATGT       0.064542472           0.064542472
#> CTGCTGTAGACTGGGT CTGCTGTAGACTGGGT       0.778945977           0.778945977
#> CTGCTGTCAAAGCAAT CTGCTGTCAAAGCAAT       0.813352551           0.813352551
#> CTGGTCTAGGATTCGG CTGGTCTAGGATTCGG       0.787766964           0.787766964
#> CTGGTCTAGGTGATTA CTGGTCTAGGTGATTA       0.841570421           0.841570421
#> CTGGTCTCAGGAATGC CTGGTCTCAGGAATGC       0.197369403           0.197369403
#> CTGGTCTTCGCATGGC CTGGTCTTCGCATGGC       0.505509545           0.505509545
#> CTGTGCTCATTACCTT CTGTGCTCATTACCTT       0.323052280           0.323052280
#> CTTAACTCAAGCGTAG CTTAACTCAAGCGTAG       0.194689590           0.194689590
#> CTTAACTTCAGCTTAG CTTAACTTCAGCTTAG       0.760329582           0.760329582
#> CTTAGGACACAGGTTT CTTAGGACACAGGTTT       0.162497323           0.162497323
#> CTTAGGAGTCCGAATT CTTAGGAGTCCGAATT       0.219614785           0.219614785
#> CTTAGGATCAACGGGA CTTAGGATCAACGGGA       0.855547166           0.855547166
#> CTTAGGATCACTCTTA CTTAGGATCACTCTTA       0.522006379           0.522006379
#> CTTCTCTAGGCAGGTT CTTCTCTAGGCAGGTT       0.107299719           0.107299719
#> CTTCTCTAGTCCGTAT CTTCTCTAGTCCGTAT       0.827072246           0.827072246
#> CTTCTCTCAGCATACT CTTCTCTCAGCATACT       0.825435827           0.825435827
#> CTTCTCTGTCAGAGGT CTTCTCTGTCAGAGGT       0.827458016           0.827458016
#> CTTGGCTAGACAGACC CTTGGCTAGACAGACC       0.003347141           0.003347141
#> CTTGGCTAGGACACCA CTTGGCTAGGACACCA       0.208296954           0.208296954
#> CTTGGCTCAAAGGAAG CTTGGCTCAAAGGAAG       0.777074417           0.777074417
#> CTTGGCTCAGCCACCA CTTGGCTCAGCCACCA       0.527067136           0.527067136
#> CTTGGCTCATCCCACT CTTGGCTCATCCCACT       0.849854411           0.849854411
#> CTTGGCTGTTAGAACA CTTGGCTGTTAGAACA       0.879461355           0.879461355
#> CTTTGCGAGACCACGA CTTTGCGAGACCACGA       0.790897668           0.790897668
#> CTTTGCGAGGAATGGA CTTTGCGAGGAATGGA       0.210555244           0.210555244
#> CTTTGCGTCCACGTGG CTTTGCGTCCACGTGG       0.946254520           0.946254520
#> GAAACTCAGAAACCTA GAAACTCAGAAACCTA       0.843774087           0.843774087
#> GAAATGAAGCTCCTCT GAAATGAAGCTCCTCT       0.470383442           0.470383442
#> GAAATGACAACTGGCC GAAATGACAACTGGCC       0.185799008           0.185799008
#> GAAATGACACATGACT GAAATGACACATGACT       0.190162030           0.190162030
#> GAAATGAGTCACCCAG GAAATGAGTCACCCAG       0.361327874           0.361327874
#> GAAATGAGTCGGATCC GAAATGAGTCGGATCC       0.200002219           0.200002219
#> GAACATCGTAATCGTC GAACATCGTAATCGTC       0.157484330           0.157484330
#> GAACCTAAGTAGATGT GAACCTAAGTAGATGT       0.817371752           0.817371752
#> GAACCTATCACTTATC GAACCTATCACTTATC       0.814212099           0.814212099
#> GAACGGAAGTCAAGCG GAACGGAAGTCAAGCG       0.814578738           0.814578738
#> GAACGGACAAGGACAC GAACGGACAAGGACAC       0.829114579           0.829114579
#> GAAGCAGGTCGCGAAA GAAGCAGGTCGCGAAA       0.765114287           0.765114287
#> GAAGCAGGTCTAGAGG GAAGCAGGTCTAGAGG       0.832194292           0.832194292
#> GAAGCAGTCAAGGCTT GAAGCAGTCAAGGCTT       0.193543818           0.193543818
#> GAATGAAAGAATGTGT GAATGAAAGAATGTGT       0.670103465           0.670103465
#> GAATGAACAATCGAAA GAATGAACAATCGAAA       0.195983200           0.195983200
#> GAATGAACAGCCTTGG GAATGAACAGCCTTGG       0.027394431           0.027394431
#> GAATGAATCATCGGAT GAATGAATCATCGGAT       0.812433540           0.812433540
#> GACACGCCAAAGGAAG GACACGCCAAAGGAAG       0.197834592           0.197834592
#> GACCAATAGAATGTGT GACCAATAGAATGTGT       0.937477362           0.937477362
#> GACCAATCAGACGCCT GACCAATCAGACGCCT       0.829681097           0.829681097
#> GACCAATGTTATCCGA GACCAATGTTATCCGA       0.819000106           0.819000106
#> GACCAATGTTGGACCC GACCAATGTTGGACCC       0.841758334           0.841758334
#> GACCAATTCAACGGCC GACCAATTCAACGGCC       0.480309570           0.480309570
#> GACCTGGGTTCGTGAT GACCTGGGTTCGTGAT       0.818368087           0.818368087
#> GACGCGTCAGTCGTGC GACGCGTCAGTCGTGC       0.210698205           0.210698205
#> GACGCGTGTCGGCTCA GACGCGTGTCGGCTCA       0.838242600           0.838242600
#> GACGGCTAGGTGATAT GACGGCTAGGTGATAT       0.200165045           0.200165045
#> GACGGCTTCAACTCTT GACGGCTTCAACTCTT       0.790138369           0.790138369
#> GACGTGCGTCACCTAA GACGTGCGTCACCTAA       0.816597790           0.816597790
#> GACGTGCTCCGGCACA GACGTGCTCCGGCACA       0.834404424           0.834404424
#> GACGTTAAGACGCAAC GACGTTAAGACGCAAC       0.210008862           0.210008862
#> GACGTTAAGACGCTTT GACGTTAAGACGCTTT       0.461743219           0.461743219
#> GACGTTACATCTCGCT GACGTTACATCTCGCT       0.980981521           0.980981521
#> GACGTTAGTGTGAATA GACGTTAGTGTGAATA       0.752031271           0.752031271
#> GACTAACGTCATCCCT GACTAACGTCATCCCT       0.631847767           0.631847767
#> GACTAACTCAGTCCCT GACTAACTCAGTCCCT       0.078481834           0.078481834
#> GACTACAGTTCGTGAT GACTACAGTTCGTGAT       0.199631682           0.199631682
#> GACTGCGAGAATCTCC GACTGCGAGAATCTCC       0.673218409           0.673218409
#> GAGCAGAAGCGTGTCC GAGCAGAAGCGTGTCC       0.085083800           0.085083800
#> GAGCAGAGTTTGCATG GAGCAGAGTTTGCATG       0.831964019           0.831964019
#> GAGCAGATCAACTCTT GAGCAGATCAACTCTT       0.805636142           0.805636142
#> GAGGTGAGTCCGCTGA GAGGTGAGTCCGCTGA       0.808874032           0.808874032
#> GAGGTGAGTGACAAAT GAGGTGAGTGACAAAT       0.212401881           0.212401881
#> GAGTCCGAGCGATGAC GAGTCCGAGCGATGAC       0.452461546           0.452461546
#> GAGTCCGAGTTAAGTG GAGTCCGAGTTAAGTG       0.200923651           0.200923651
#> GAGTCCGCAATTGCTG GAGTCCGCAATTGCTG       0.148372394           0.148372394
#> GATCAGTAGGTGATAT GATCAGTAGGTGATAT       0.801912106           0.801912106
#> GATCGATAGAAAGTGG GATCGATAGAAAGTGG       0.702546237           0.702546237
#> GATCGATAGGATGGTC GATCGATAGGATGGTC       0.825001567           0.825001567
#> GATCGATGTAGCGATG GATCGATGTAGCGATG       0.826428210           0.826428210
#> GATCGATGTTGTTTGG GATCGATGTTGTTTGG       0.817706859           0.817706859
#> GATCGCGCACCAGGCT GATCGCGCACCAGGCT       0.770851021           0.770851021
#> GATCGTAAGTAGGTGC GATCGTAAGTAGGTGC       0.703706497           0.703706497
#> GATCGTAGTCCTAGCG GATCGTAGTCCTAGCG       0.043970628           0.043970628
#> GATCTAGAGCTAGCCC GATCTAGAGCTAGCCC       0.850973753           0.850973753
#> GATCTAGGTGCCTTGG GATCTAGGTGCCTTGG       0.204761440           0.204761440
#> GATGAAAAGTTGTAGA GATGAAAAGTTGTAGA       0.030479691           0.030479691
#> GATGAAACAGGCTGAA GATGAAACAGGCTGAA       0.827695538           0.827695538
#> GATGAAACATGGTTGT GATGAAACATGGTTGT       0.211448618           0.211448618
#> GATGAGGAGTGGACGT GATGAGGAGTGGACGT       0.413533668           0.413533668
#> GATGAGGTCTGGGCCA GATGAGGTCTGGGCCA       0.214485499           0.214485499
#> GATGCTAGTAGAAGGA GATGCTAGTAGAAGGA       0.063548947           0.063548947
#> GATGCTAGTGTTAAGA GATGCTAGTGTTAAGA       0.785309907           0.785309907
#> GATTCAGAGCGATGAC GATTCAGAGCGATGAC       0.814342380           0.814342380
#> GATTCAGAGTACGTTC GATTCAGAGTACGTTC       0.832775629           0.832775629
#> GATTCAGTCGCCAGCA GATTCAGTCGCCAGCA       0.811921328           0.811921328
#> GATTCAGTCGTCCAGG GATTCAGTCGTCCAGG       0.853451267           0.853451267
#> GCAAACTCAGACAAAT GCAAACTCAGACAAAT       0.823068269           0.823068269
#> GCAATCACAGCTGCTG GCAATCACAGCTGCTG       0.813572816           0.813572816
#> GCAATCATCCAACCAA GCAATCATCCAACCAA       0.203884102           0.203884102
#> GCACTCTCAGCGTCCA GCACTCTCAGCGTCCA       0.840060411           0.840060411
#> GCAGCCACAGCGTTCG GCAGCCACAGCGTTCG       0.322739169           0.322739169
#> GCAGTTACATTTCACT GCAGTTACATTTCACT       0.759491577           0.759491577
#> GCAGTTATCTTGAGAC GCAGTTATCTTGAGAC       0.214511750           0.214511750
#> GCATGATAGCCTTGAT GCATGATAGCCTTGAT       0.248604289           0.248604289
#> GCATGATCAAGCGAGT GCATGATCAAGCGAGT       0.859328294           0.859328294
#> GCATGATGTCACTGGC GCATGATGTCACTGGC       0.791220089           0.791220089
#> GCATGATGTTTAGGAA GCATGATGTTTAGGAA       0.816005456           0.816005456
#> GCATGATTCAAACGGG GCATGATTCAAACGGG       0.204064073           0.204064073
#> GCATGCGAGGAACTGC GCATGCGAGGAACTGC       0.814409434           0.814409434
#> GCATGCGCAGCTCCGA GCATGCGCAGCTCCGA       0.825645631           0.825645631
#> GCATGTAAGAATCTCC GCATGTAAGAATCTCC       0.878933088           0.878933088
#> GCATGTAAGCGTGTCC GCATGTAAGCGTGTCC       0.849580847           0.849580847
#> GCATGTAGTTGCGCAC GCATGTAGTTGCGCAC       0.821664060           0.821664060
#> GCATGTATCGTGGACC GCATGTATCGTGGACC       0.824139425           0.824139425
#> GCCTCTAAGCTAGTCT GCCTCTAAGCTAGTCT       0.855577941           0.855577941
#> GCCTCTAAGTCCGGTC GCCTCTAAGTCCGGTC       0.275752807           0.275752807
#> GCCTCTATCGCTTAGA GCCTCTATCGCTTAGA       0.209843711           0.209843711
#> GCGACCATCACCAGGC GCGACCATCACCAGGC       0.819712350           0.819712350
#> GCGAGAAGTGGGTATG GCGAGAAGTGGGTATG       0.829742438           0.829742438
#> GCGCAACGTTTCGCTC GCGCAACGTTTCGCTC       0.833436420           0.833436420
#> GCGCAACTCACAACGT GCGCAACTCACAACGT       0.154247686           0.154247686
#> GCGCAGTGTGCTAGCC GCGCAGTGTGCTAGCC       0.807637653           0.807637653
#> GCGCCAAAGCTCCTTC GCGCCAAAGCTCCTTC       0.824127064           0.824127064
#> GCGCCAATCTTATCTG GCGCCAATCTTATCTG       0.749554681           0.749554681
#> GCGCGATGTCCGAATT GCGCGATGTCCGAATT       0.993437165           0.993437165
#> GCGGGTTAGACTTGAA GCGGGTTAGACTTGAA       0.291653603           0.291653603
#> GCGGGTTCAAGTAATG GCGGGTTCAAGTAATG       0.830776382           0.830776382
#> GCGGGTTGTTACTGAC GCGGGTTGTTACTGAC       0.203945635           0.203945635
#> GCTCCTAAGGGCATGT GCTCCTAAGGGCATGT       0.892725541           0.892725541
#> GCTCTGTGTCCATGAT GCTCTGTGTCCATGAT       0.779331331           0.779331331
#> GCTGCAGGTTACCGAT GCTGCAGGTTACCGAT       0.199051792           0.199051792
#> GCTGCGAAGAACAACT GCTGCGAAGAACAACT       0.812716285           0.812716285
#> GCTGCGAAGCTAGCCC GCTGCGAAGCTAGCCC       0.281019378           0.281019378
#> GCTGCGAAGTATTGGA GCTGCGAAGTATTGGA       0.582813572           0.582813572
#> GCTGCGATCCCGGATG GCTGCGATCCCGGATG       0.731837454           0.731837454
#> GCTGCGATCGGCTACG GCTGCGATCGGCTACG       0.962384087           0.962384087
#> GCTGCTTGTCTTCAAG GCTGCTTGTCTTCAAG       0.075838386           0.075838386
#> GCTGGGTGTCATACTG GCTGGGTGTCATACTG       0.826109309           0.826109309
#> GCTGGGTGTCTAGAGG GCTGGGTGTCTAGAGG       0.415702282           0.415702282
#> GCTTCCAAGACATAAC GCTTCCAAGACATAAC       0.831475233           0.831475233
#> GCTTCCAAGTTAAGTG GCTTCCAAGTTAAGTG       0.828650365           0.828650365
#> GCTTGAAGTCAACTGT GCTTGAAGTCAACTGT       0.826954217           0.826954217
#> GGAAAGCCAATGTAAG GGAAAGCCAATGTAAG       0.814951905           0.814951905
#> GGAAAGCGTCTTGATG GGAAAGCGTCTTGATG       0.116117713           0.116117713
#> GGAATAAGTAAGAGAG GGAATAAGTAAGAGAG       0.844614146           0.844614146
#> GGAATAAGTTCGGGCT GGAATAAGTTCGGGCT       0.133954470           0.133954470
#> GGACAAGGTGATAAGT GGACAAGGTGATAAGT       0.184359523           0.184359523
#> GGACAAGTCAAAGACA GGACAAGTCAAAGACA       0.630306486           0.630306486
#> GGACAGAAGGCCGAAT GGACAGAAGGCCGAAT       0.849694246           0.849694246
#> GGACAGACAATCCAAC GGACAGACAATCCAAC       0.803829413           0.803829413
#> GGACAGATCAGGCGAA GGACAGATCAGGCGAA       0.921482331           0.921482331
#> GGACAGATCCTGCAGG GGACAGATCCTGCAGG       0.044200433           0.044200433
#> GGACAGATCCTTTCGG GGACAGATCCTTTCGG       0.576407555           0.576407555
#> GGACATTGTCTCCCTA GGACATTGTCTCCCTA       0.799647163           0.799647163
#> GGACGTCAGGGTGTTG GGACGTCAGGGTGTTG       0.803544278           0.803544278
#> GGACGTCCACGCATCG GGACGTCCACGCATCG       0.192957672           0.192957672
#> GGACGTCCATGGATGG GGACGTCCATGGATGG       0.604376826           0.604376826
#> GGACGTCGTCAGAAGC GGACGTCGTCAGAAGC       0.728664483           0.728664483
#> GGACGTCGTCGAGATG GGACGTCGTCGAGATG       0.835213543           0.835213543
#> GGACGTCGTCTAGGTT GGACGTCGTCTAGGTT       0.866970060           0.866970060
#> GGACGTCGTTTAGCTG GGACGTCGTTTAGCTG       0.025448610           0.025448610
#> GGAGCAAAGCTTCGCG GGAGCAAAGCTTCGCG       0.145408838           0.145408838
#> GGAGCAAGTATCACCA GGAGCAAGTATCACCA       0.787253377           0.787253377
#> GGATGTTAGTAAGTAC GGATGTTAGTAAGTAC       0.824991241           0.824991241
#> GGATGTTTCCAGTATG GGATGTTTCCAGTATG       0.520488822           0.520488822
#> GGCAATTCACGGCTAC GGCAATTCACGGCTAC       0.814299507           0.814299507
#> GGCAATTCATCTCCCA GGCAATTCATCTCCCA       0.559265129           0.559265129
#> GGCAATTGTAGAGCTG GGCAATTGTAGAGCTG       0.968511170           0.968511170
#> GGCCGATTCACTCTTA GGCCGATTCACTCTTA       0.213404735           0.213404735
#> GGCGACTGTAAGTTCC GGCGACTGTAAGTTCC       0.825401575           0.825401575
#> GGCGTGTCACCAGGTC GGCGTGTCACCAGGTC       0.675959523           0.675959523
#> GGCGTGTCACCGATAT GGCGTGTCACCGATAT       0.725287192           0.725287192
#> GGCGTGTTCAGCTTAG GGCGTGTTCAGCTTAG       0.812309434           0.812309434
#> GGCGTGTTCGAGAGCA GGCGTGTTCGAGAGCA       0.788915010           0.788915010
#> GGCTCGAAGATCACGG GGCTCGAAGATCACGG       0.707434444           0.707434444
#> GGCTCGAAGTGACTCT GGCTCGAAGTGACTCT       0.864891657           0.864891657
#> GGCTCGACAGCTTCGG GGCTCGACAGCTTCGG       0.886641603           0.886641603
#> GGCTGGTAGTCACGCC GGCTGGTAGTCACGCC       0.598068596           0.598068596
#> GGCTGGTTCATGCATG GGCTGGTTCATGCATG       0.609780034           0.609780034
#> GGGAATGAGCACACAG GGGAATGAGCACACAG       0.839827495           0.839827495
#> GGGAATGAGTGTACGG GGGAATGAGTGTACGG       0.185291218           0.185291218
#> GGGAATGCAAAGTGCG GGGAATGCAAAGTGCG       0.822639487           0.822639487
#> GGGAATGTCACAAACC GGGAATGTCACAAACC       0.549968981           0.549968981
#> GGGACCTCACTGTCGG GGGACCTCACTGTCGG       0.829410938           0.829410938
#> GGGACCTTCGAATGCT GGGACCTTCGAATGCT       0.194037641           0.194037641
#> GGGAGATAGATGTGGC GGGAGATAGATGTGGC       0.828391217           0.828391217
#> GGGAGATAGCGAAGGG GGGAGATAGCGAAGGG       0.837788595           0.837788595
#> GGGAGATAGTGACATA GGGAGATAGTGACATA       0.824230400           0.824230400
#> GGGATGAAGGGCATGT GGGATGAAGGGCATGT       0.822970550           0.822970550
#> GGGATGAGTCTGCGGT GGGATGAGTCTGCGGT       0.817854475           0.817854475
#> GGGATGAGTTGTGGAG GGGATGAGTTGTGGAG       0.415913489           0.415913489
#> GGGCACTCATTTCAGG GGGCACTCATTTCAGG       0.200881642           0.200881642
#> GGGCACTGTCTCCATC GGGCACTGTCTCCATC       0.729507089           0.729507089
#> GGGCACTGTGGCAAAC GGGCACTGTGGCAAAC       0.181709199           0.181709199
#> GGGCATCCATGACATC GGGCATCCATGACATC       0.133111667           0.133111667
#> GGGCATCCATTTCAGG GGGCATCCATTTCAGG       0.824354181           0.824354181
#> GGGCATCGTGGTCTCG GGGCATCGTGGTCTCG       0.353855710           0.353855710
#> GGGTCTGCACAGGCCT GGGTCTGCACAGGCCT       0.303789502           0.303789502
#> GGGTCTGCATGAAGTA GGGTCTGCATGAAGTA       0.210133162           0.210133162
#> GGGTCTGTCATGTGGT GGGTCTGTCATGTGGT       0.709135447           0.709135447
#> GGGTCTGTCCATGAAC GGGTCTGTCCATGAAC       0.828308532           0.828308532
#> GGGTTGCCAAGCTGAG GGGTTGCCAAGCTGAG       0.191083672           0.191083672
#> GGGTTGCTCGAACGGA GGGTTGCTCGAACGGA       0.197352133           0.197352133
#> GGTATTGGTACTCGCG GGTATTGGTACTCGCG       0.308643288           0.308643288
#> GGTATTGGTCCCGACA GGTATTGGTCCCGACA       0.183717671           0.183717671
#> GGTGTTAAGGATCGCA GGTGTTAAGGATCGCA       0.220700057           0.220700057
#> GGTGTTAGTTCCTCCA GGTGTTAGTTCCTCCA       0.197506173           0.197506173
#> GTAACGTGTCACTGGC GTAACGTGTCACTGGC       0.213553393           0.213553393
#> GTAACTGGTCACTGGC GTAACTGGTCACTGGC       0.724395370           0.724395370
#> GTAACTGTCAGTGCAT GTAACTGTCAGTGCAT       0.062248752           0.062248752
#> GTACGTAGTCTCTTTA GTACGTAGTCTCTTTA       0.812433385           0.812433385
#> GTACTCCCACAGGTTT GTACTCCCACAGGTTT       0.748816105           0.748816105
#> GTACTCCCACGTCTCT GTACTCCCACGTCTCT       0.194483613           0.194483613
#> GTACTCCTCCAACCAA GTACTCCTCCAACCAA       0.428005504           0.428005504
#> GTACTTTGTACCGTAT GTACTTTGTACCGTAT       0.857221641           0.857221641
#> GTACTTTTCTTCGGTC GTACTTTTCTTCGGTC       0.795740398           0.795740398
#> GTAGGCCAGCTAACAA GTAGGCCAGCTAACAA       0.797646238           0.797646238
#> GTAGGCCGTTACCAGT GTAGGCCGTTACCAGT       0.910658552           0.910658552
#> GTAGTCACATCGATGT GTAGTCACATCGATGT       0.582153395           0.582153395
#> GTAGTCAGTCCGACGT GTAGTCAGTCCGACGT       0.650972869           0.650972869
#> GTATCTTCAGTGGGAT GTATCTTCAGTGGGAT       0.212279673           0.212279673
#> GTATCTTCATCGATGT GTATCTTCATCGATGT       0.296309912           0.296309912
#> GTATCTTGTACTCAAC GTATCTTGTACTCAAC       0.826256957           0.826256957
#> GTATCTTGTCGCGGTT GTATCTTGTCGCGGTT       0.452418441           0.452418441
#> GTATCTTGTGGCCCTA GTATCTTGTGGCCCTA       0.371538774           0.371538774
#> GTATCTTTCGTGGACC GTATCTTTCGTGGACC       0.764391816           0.764391816
#> GTATTCTCACATTAGC GTATTCTCACATTAGC       0.201164447           0.201164447
#> GTATTCTCACTTGGAT GTATTCTCACTTGGAT       0.009020684           0.009020684
#> GTATTCTGTCCGTCAG GTATTCTGTCCGTCAG       0.014072783           0.014072783
#> GTATTCTTCTCTAAGG GTATTCTTCTCTAAGG       0.854086097           0.854086097
#> GTCAAGTGTGGTCCGT GTCAAGTGTGGTCCGT       0.077110665           0.077110665
#> GTCACAACACGTCTCT GTCACAACACGTCTCT       0.547833321           0.547833321
#> GTCACGGAGCTATGCT GTCACGGAGCTATGCT       0.850209714           0.850209714
#> GTCACGGCATCCGTGG GTCACGGCATCCGTGG       0.781248705           0.781248705
#> GTCACGGTCTTGCCGT GTCACGGTCTTGCCGT       0.781439722           0.781439722
#> GTCATTTGTCGTCTTC GTCATTTGTCGTCTTC       0.805758490           0.805758490
#> GTCCTCACACCACGTG GTCCTCACACCACGTG       0.798912179           0.798912179
#> GTCCTCATCATTTGGG GTCCTCATCATTTGGG       0.237734072           0.237734072
#> GTCCTCATCTGGAGCC GTCCTCATCTGGAGCC       0.379061056           0.379061056
#> GTCGGGTCAGGCTGAA GTCGGGTCAGGCTGAA       0.862913679           0.862913679
#> GTCGTAAAGTTTAGGA GTCGTAAAGTTTAGGA       0.854887886           0.854887886
#> GTCTCGTAGGCTAGAC GTCTCGTAGGCTAGAC       0.245226646           0.245226646
#> GTCTCGTCATAGTAAG GTCTCGTCATAGTAAG       0.413160838           0.413160838
#> GTCTCGTTCACTTCAT GTCTCGTTCACTTCAT       0.476170903           0.476170903
#> GTCTCGTTCTGCGTAA GTCTCGTTCTGCGTAA       0.888163243           0.888163243
#> GTCTTCGTCCGTTGCT GTCTTCGTCCGTTGCT       0.830506858           0.830506858
#> GTCTTCGTCTACCAGA GTCTTCGTCTACCAGA       0.766981744           0.766981744
#> GTGCATAAGAACAATC GTGCATAAGAACAATC       0.800274890           0.800274890
#> GTGCATAAGCTCCTCT GTGCATAAGCTCCTCT       0.125036206           0.125036206
#> GTGCATAGTTCAACCA GTGCATAGTTCAACCA       0.820670010           0.820670010
#> GTGCATATCAAGGCTT GTGCATATCAAGGCTT       0.838532644           0.838532644
#> GTGCGGTAGTAGGCCA GTGCGGTAGTAGGCCA       0.238496927           0.238496927
#> GTGCGGTCAGCGTAAG GTGCGGTCAGCGTAAG       0.410200458           0.410200458
#> GTGCTTCAGAGGGATA GTGCTTCAGAGGGATA       0.815424150           0.815424150
#> GTGCTTCTCACGACTA GTGCTTCTCACGACTA       0.880120124           0.880120124
#> GTGGGTCCATAAAGGT GTGGGTCCATAAAGGT       0.852138994           0.852138994
#> GTGGGTCCATCGATTG GTGGGTCCATCGATTG       0.897436611           0.897436611
#> GTGGGTCGTCCAGTAT GTGGGTCGTCCAGTAT       0.830735888           0.830735888
#> GTGGGTCTCCGAATGT GTGGGTCTCCGAATGT       0.821652981           0.821652981
#> GTGTGCGGTGTAAGTA GTGTGCGGTGTAAGTA       0.829089840           0.829089840
#> GTGTTAGCACCATCCT GTGTTAGCACCATCCT       0.131030061           0.131030061
#> GTGTTAGGTGAACCTT GTGTTAGGTGAACCTT       0.785815137           0.785815137
#> GTTAAGCAGGATTCGG GTTAAGCAGGATTCGG       0.576176002           0.576176002
#> GTTACAGTCCGATATG GTTACAGTCCGATATG       0.714704756           0.714704756
#> GTTCATTAGAAGATTC GTTCATTAGAAGATTC       0.201901852           0.201901852
#> GTTCATTAGGCAGTCA GTTCATTAGGCAGTCA       0.754855029           0.754855029
#> GTTCATTTCCGAACGC GTTCATTTCCGAACGC       0.855590491           0.855590491
#> GTTCTCGAGAATTCCC GTTCTCGAGAATTCCC       0.824300271           0.824300271
#> GTTCTCGCAAGAAAGG GTTCTCGCAAGAAAGG       0.832235164           0.832235164
#> GTTCTCGCAGACACTT GTTCTCGCAGACACTT       0.201018528           0.201018528
#> GTTTCTAAGAAGATTC GTTTCTAAGAAGATTC       0.865368698           0.865368698
#> GTTTCTACAGGACCCT GTTTCTACAGGACCCT       0.850315821           0.850315821
#> TAAGAGATCTGCGACG TAAGAGATCTGCGACG       0.770195024           0.770195024
#> TAAGTGCAGGGTTTCT TAAGTGCAGGGTTTCT       0.186890362           0.186890362
#> TAAGTGCTCCAAATGC TAAGTGCTCCAAATGC       0.804656016           0.804656016
#> TACACGACAACGATCT TACACGACAACGATCT       0.213015394           0.213015394
#> TACACGACATGTAGTC TACACGACATGTAGTC       0.830600612           0.830600612
#> TACACGAGTGCAACTT TACACGAGTGCAACTT       0.637962372           0.637962372
#> TACAGTGCAATCTACG TACAGTGCAATCTACG       0.192323943           0.192323943
#> TACAGTGCACAGGCCT TACAGTGCACAGGCCT       0.832433066           0.832433066
#> TACAGTGGTAAACGCG TACAGTGGTAAACGCG       0.709330688           0.709330688
#> TACCTTAAGTACATGA TACCTTAAGTACATGA       0.184728989           0.184728989
#> TACCTTATCCGATATG TACCTTATCCGATATG       0.831320021           0.831320021
#> TACGGATAGAGCTGGT TACGGATAGAGCTGGT       0.166654060           0.166654060
#> TACGGATAGGCACATG TACGGATAGGCACATG       0.853702325           0.853702325
#> TACGGATAGTTTAGGA TACGGATAGTTTAGGA       0.199343637           0.199343637
#> TACGGATGTGTTCGAT TACGGATGTGTTCGAT       0.216004017           0.216004017
#> TACGGATTCCAAGTAC TACGGATTCCAAGTAC       0.524754663           0.524754663
#> TACGGGCAGCCAACAG TACGGGCAGCCAACAG       0.811162805           0.811162805
#> TACGGGCAGCGTCAAG TACGGGCAGCGTCAAG       0.210609367           0.210609367
#> TACGGGCTCTTCGGTC TACGGGCTCTTCGGTC       0.862090807           0.862090807
#> TACGGTACATGTTCCC TACGGTACATGTTCCC       0.191738650           0.191738650
#> TACTCATAGGTGTTAA TACTCATAGGTGTTAA       0.780815012           0.780815012
#> TACTCATCAGATCCAT TACTCATCAGATCCAT       0.871901740           0.871901740
#> TACTCGCAGTAAGTAC TACTCGCAGTAAGTAC       0.826026204           0.826026204
#> TACTCGCCAGATCGGA TACTCGCCAGATCGGA       0.842306082           0.842306082
#> TACTCGCGTCGGCATC TACTCGCGTCGGCATC       0.842122834           0.842122834
#> TACTTACCAACTGCGC TACTTACCAACTGCGC       0.196793226           0.196793226
#> TACTTACCATGGGACA TACTTACCATGGGACA       0.530375014           0.530375014
#> TACTTGTAGACAAAGG TACTTGTAGACAAAGG       0.922242466           0.922242466
#> TAGACCAAGGAATCGC TAGACCAAGGAATCGC       0.193950012           0.193950012
#> TAGAGCTGTATGAAAC TAGAGCTGTATGAAAC       0.051208945           0.051208945
#> TAGAGCTGTTGTCGCG TAGAGCTGTTGTCGCG       0.049321242           0.049321242
#> TAGCCGGTCAACACAC TAGCCGGTCAACACAC       0.168875071           0.168875071
#> TAGGCATTCCAAAGTC TAGGCATTCCAAAGTC       0.037977141           0.037977141
#> TAGTGGTCAACTGCGC TAGTGGTCAACTGCGC       0.025848352           0.025848352
#> TAGTGGTCATTCACTT TAGTGGTCATTCACTT       0.780208343           0.780208343
#> TAGTGGTTCCTCAACC TAGTGGTTCCTCAACC       0.736319252           0.736319252
#> TAGTTGGCACGTAAGG TAGTTGGCACGTAAGG       0.572181731           0.572181731
#> TAGTTGGCAGTAACGG TAGTTGGCAGTAACGG       0.767064602           0.767064602
#> TAGTTGGGTACCGTAT TAGTTGGGTACCGTAT       0.196190018           0.196190018
#> TAGTTGGTCTTGCAAG TAGTTGGTCTTGCAAG       0.999659154           0.999659154
#> TATCAGGCACCAGATT TATCAGGCACCAGATT       0.200674637           0.200674637
#> TATCTCACAACTGCTA TATCTCACAACTGCTA       0.198934040           0.198934040
#> TATCTCACACGCGAAA TATCTCACACGCGAAA       0.811771402           0.811771402
#> TATCTCACATTGTGCA TATCTCACATTGTGCA       0.846746610           0.846746610
#> TATCTCAGTTCGCTAA TATCTCAGTTCGCTAA       0.019292876           0.019292876
#> TATGCCCAGAACAATC TATGCCCAGAACAATC       0.822157601           0.822157601
#> TATGCCCGTCTAGCCG TATGCCCGTCTAGCCG       0.337016568           0.337016568
#> TATGCCCTCTTACCTA TATGCCCTCTTACCTA       0.030966241           0.030966241
#> TATTACCTCGTTGCCT TATTACCTCGTTGCCT       0.824481158           0.824481158
#> TCAACGATCAACCATG TCAACGATCAACCATG       0.697595703           0.697595703
#> TCAATCTCATCCGCGA TCAATCTCATCCGCGA       0.712090257           0.712090257
#> TCAATCTCATCTCCCA TCAATCTCATCTCCCA       0.202919339           0.202919339
#> TCAATCTGTCAAAGAT TCAATCTGTCAAAGAT       0.000000000           0.000000000
#> TCAATCTTCGAACGGA TCAATCTTCGAACGGA       0.826103336           0.826103336
#> TCACAAGAGTGACATA TCACAAGAGTGACATA       0.849499521           0.849499521
#> TCACAAGCAACGATGG TCACAAGCAACGATGG       0.812092583           0.812092583
#> TCACAAGCAGACTCGC TCACAAGCAGACTCGC       0.942197487           0.942197487
#> TCACGAAGTCCGAGTC TCACGAAGTCCGAGTC       0.597311392           0.597311392
#> TCACGAATCCTATGTT TCACGAATCCTATGTT       0.193285546           0.193285546
#> TCAGATGCACTCGACG TCAGATGCACTCGACG       0.713302027           0.713302027
#> TCAGATGGTACACCGC TCAGATGGTACACCGC       0.860957157           0.860957157
#> TCAGCAACACCAGGCT TCAGCAACACCAGGCT       0.202024080           0.202024080
#> TCAGCAATCGCCGTGA TCAGCAATCGCCGTGA       0.777241847           0.777241847
#> TCAGCTCAGTCGTACT TCAGCTCAGTCGTACT       0.223010571           0.223010571
#> TCAGGTAAGGCGTACA TCAGGTAAGGCGTACA       0.619475967           0.619475967
#> TCAGGTAAGTCAATAG TCAGGTAAGTCAATAG       0.213930079           0.213930079
#> TCAGGTACAGTCACTA TCAGGTACAGTCACTA       0.825283605           0.825283605
#> TCAGGTATCTGTCCGT TCAGGTATCTGTCCGT       0.849794159           0.849794159
#> TCATTACGTCAGAAGC TCATTACGTCAGAAGC       0.204618423           0.204618423
#> TCATTTGCATTGGCGC TCATTTGCATTGGCGC       0.749697422           0.749697422
#> TCCACACCAGGCGATA TCCACACCAGGCGATA       0.832740466           0.832740466
#> TCCCGATAGTCGTTTG TCCCGATAGTCGTTTG       0.768241639           0.768241639
#> TCCCGATCAACGATGG TCCCGATCAACGATGG       0.540321885           0.540321885
#> TCGAGGCAGAACAATC TCGAGGCAGAACAATC       0.828562465           0.828562465
#> TCGAGGCAGGCATGGT TCGAGGCAGGCATGGT       0.794961309           0.794961309
#> TCGCGAGAGTACGCGA TCGCGAGAGTACGCGA       0.690803027           0.690803027
#> TCGCGAGAGTCCGTAT TCGCGAGAGTCCGTAT       0.868990298           0.868990298
#> TCGCGAGCACTTCTGC TCGCGAGCACTTCTGC       0.804410334           0.804410334
#> TCGGGACAGAGAACAG TCGGGACAGAGAACAG       0.256299790           0.256299790
#> TCGGGACGTCTGCGGT TCGGGACGTCTGCGGT       0.839918942           0.839918942
#> TCGTACCCACGCATCG TCGTACCCACGCATCG       0.064000405           0.064000405
#> TCGTACCCATATACGC TCGTACCCATATACGC       0.220951002           0.220951002
#> TCTATTGCACAGGAGT TCTATTGCACAGGAGT       0.105017634           0.105017634
#> TCTATTGTCCTATTCA TCTATTGTCCTATTCA       0.845593292           0.845593292
#> TCTATTGTCGTATCAG TCTATTGTCGTATCAG       0.689210156           0.689210156
#> TCTCATAAGAATGTGT TCTCATAAGAATGTGT       0.763988374           0.763988374
#> TCTCATACAAATTGCC TCTCATACAAATTGCC       0.832498966           0.832498966
#> TCTCATAGTGGACGAT TCTCATAGTGGACGAT       0.517250822           0.517250822
#> TCTCTAAAGCTGGAAC TCTCTAAAGCTGGAAC       0.880046632           0.880046632
#> TCTCTAAGTCGCATCG TCTCTAAGTCGCATCG       0.199815201           0.199815201
#> TCTCTAAGTGGCAAAC TCTCTAAGTGGCAAAC       0.494136263           0.494136263
#> TCTCTAAGTTACCGAT TCTCTAAGTTACCGAT       0.170275658           0.170275658
#> TCTCTAATCCCAACGG TCTCTAATCCCAACGG       0.709892144           0.709892144
#> TCTGAGAGTCAAACTC TCTGAGAGTCAAACTC       0.067549735           0.067549735
#> TCTGGAACACACGCTG TCTGGAACACACGCTG       0.179757093           0.179757093
#> TCTGGAACACCTATCC TCTGGAACACCTATCC       0.913445240           0.913445240
#> TCTGGAAGTCATCGGC TCTGGAAGTCATCGGC       0.813273260           0.813273260
#> TCTGGAATCGGTGTCG TCTGGAATCGGTGTCG       0.212337644           0.212337644
#> TCTTTCCGTGCAACTT TCTTTCCGTGCAACTT       0.831236834           0.831236834
#> TCTTTCCTCGTCTGCT TCTTTCCTCGTCTGCT       0.213394528           0.213394528
#> TGAAAGAGTTATCACG TGAAAGAGTTATCACG       0.191933630           0.191933630
#> TGACGGCAGGAATGGA TGACGGCAGGAATGGA       0.169065942           0.169065942
#> TGACGGCCATATGAGA TGACGGCCATATGAGA       0.801194600           0.801194600
#> TGACTAGGTGTATGGG TGACTAGGTGTATGGG       0.435601396           0.435601396
#> TGACTAGTCCACGTGG TGACTAGTCCACGTGG       0.052407709           0.052407709
#> TGACTAGTCCTATTCA TGACTAGTCCTATTCA       0.386221485           0.386221485
#> TGACTTTAGCCAGAAC TGACTTTAGCCAGAAC       0.938732515           0.938732515
#> TGACTTTAGCGTTGCC TGACTTTAGCGTTGCC       0.043354617           0.043354617
#> TGACTTTAGCGTTTAC TGACTTTAGCGTTTAC       0.200825544           0.200825544
#> TGAGAGGAGATCCCAT TGAGAGGAGATCCCAT       0.866170355           0.866170355
#> TGAGAGGAGCCAGGAT TGAGAGGAGCCAGGAT       0.837719139           0.837719139
#> TGAGAGGAGCCGTCGT TGAGAGGAGCCGTCGT       0.526156786           0.526156786
#> TGAGAGGAGGCGATAC TGAGAGGAGGCGATAC       0.857578421           0.857578421
#> TGAGAGGAGGTAAACT TGAGAGGAGGTAAACT       0.836860998           0.836860998
#> TGAGAGGCAGGAACGT TGAGAGGCAGGAACGT       0.824075196           0.824075196
#> TGAGAGGGTAAGGGCT TGAGAGGGTAAGGGCT       0.054871155           0.054871155
#> TGAGCATGTTAAAGAC TGAGCATGTTAAAGAC       0.394295664           0.394295664
#> TGAGCATTCACCCGAG TGAGCATTCACCCGAG       0.193956264           0.193956264
#> TGAGCATTCGGGAGTA TGAGCATTCGGGAGTA       0.024688077           0.024688077
#> TGAGCCGTCATGTCCC TGAGCCGTCATGTCCC       0.203400549           0.203400549
#> TGAGCCGTCCGAACGC TGAGCCGTCCGAACGC       0.218440872           0.218440872
#> TGAGCCGTCCGCTGTT TGAGCCGTCCGCTGTT       0.700462459           0.700462459
#> TGAGGGAAGAGGTACC TGAGGGAAGAGGTACC       0.533653793           0.533653793
#> TGAGGGACAGACTCGC TGAGGGACAGACTCGC       0.197503792           0.197503792
#> TGATTTCAGTACGCCC TGATTTCAGTACGCCC       0.228951497           0.228951497
#> TGATTTCAGTGACTCT TGATTTCAGTGACTCT       0.202867377           0.202867377
#> TGATTTCGTGAGGGTT TGATTTCGTGAGGGTT       0.862504529           0.862504529
#> TGCCAAACATGACATC TGCCAAACATGACATC       0.399635042           0.399635042
#> TGCCCATCACTGTGTA TGCCCATCACTGTGTA       0.232475093           0.232475093
#> TGCCCTAAGTCCGGTC TGCCCTAAGTCCGGTC       0.731090015           0.731090015
#> TGCGCAGGTACGCACC TGCGCAGGTACGCACC       0.041919601           0.041919601
#> TGCGCAGTCTGCCAGG TGCGCAGTCTGCCAGG       0.052642637           0.052642637
#> TGCGGGTAGACGCACA TGCGGGTAGACGCACA       0.678925814           0.678925814
#> TGCGGGTAGGGTTCCC TGCGGGTAGGGTTCCC       0.831168344           0.831168344
#> TGCGGGTTCCACTCCA TGCGGGTTCCACTCCA       0.831442153           0.831442153
#> TGCGGGTTCTCAAACG TGCGGGTTCTCAAACG       0.215768877           0.215768877
#> TGCGTGGAGGTAAACT TGCGTGGAGGTAAACT       0.202778301           0.202778301
#> TGCGTGGAGGTTCCTA TGCGTGGAGGTTCCTA       0.805341356           0.805341356
#> TGCGTGGAGTACACCT TGCGTGGAGTACACCT       0.821764163           0.821764163
#> TGCGTGGGTCCCTTGT TGCGTGGGTCCCTTGT       0.787196466           0.787196466
#> TGCGTGGGTTACTGAC TGCGTGGGTTACTGAC       0.831235862           0.831235862
#> TGCTACCCAGGGCATA TGCTACCCAGGGCATA       0.187462850           0.187462850
#> TGCTACCCATAAGACA TGCTACCCATAAGACA       0.663392293           0.663392293
#> TGCTGCTAGGATGCGT TGCTGCTAGGATGCGT       0.366138875           0.366138875
#> TGCTGCTGTCCGTCAG TGCTGCTGTCCGTCAG       0.074103058           0.074103058
#> TGGACGCAGTGTTGAA TGGACGCAGTGTTGAA       0.220845749           0.220845749
#> TGGACGCTCAGGTAAA TGGACGCTCAGGTAAA       0.016878213           0.016878213
#> TGGACGCTCTTGTTTG TGGACGCTCTTGTTTG       0.715041196           0.715041196
#> TGGCCAGAGCTACCTA TGGCCAGAGCTACCTA       0.028022796           0.028022796
#> TGGCCAGCAAGGACTG TGGCCAGCAAGGACTG       0.999813804           0.999813804
#> TGGCCAGCAATGAAAC TGGCCAGCAATGAAAC       0.776324193           0.776324193
#> TGGCCAGCACGAAAGC TGGCCAGCACGAAAGC       0.214173484           0.214173484
#> TGGCCAGGTCGGCATC TGGCCAGGTCGGCATC       0.172934248           0.172934248
#> TGGCGCATCACAACGT TGGCGCATCACAACGT       0.077049936           0.077049936
#> TGGCTGGGTACCCAAT TGGCTGGGTACCCAAT       0.846971209           0.846971209
#> TGGGAAGCACGTAAGG TGGGAAGCACGTAAGG       0.770408107           0.770408107
#> TGGGAAGCATGACATC TGGGAAGCATGACATC       0.168084971           0.168084971
#> TGGGCGTGTGTGCGTC TGGGCGTGTGTGCGTC       0.823686290           0.823686290
#> TGGTTAGCAGTCGTGC TGGTTAGCAGTCGTGC       0.794418946           0.794418946
#> TGGTTCCCACCGCTAG TGGTTCCCACCGCTAG       0.197833847           0.197833847
#> TGGTTCCGTCACCTAA TGGTTCCGTCACCTAA       0.158133590           0.158133590
#> TGTATTCCAATGACCT TGTATTCCAATGACCT       0.071387511           0.071387511
#> TGTCCCAAGTTTAGGA TGTCCCAAGTTTAGGA       0.208431017           0.208431017
#> TGTCCCACATTATCTC TGTCCCACATTATCTC       0.842316076           0.842316076
#> TGTCCCAGTCCGTGAC TGTCCCAGTCCGTGAC       0.907393185           0.907393185
#> TGTCCCATCTCTGCTG TGTCCCATCTCTGCTG       0.834458403           0.834458403
#> TGTGGTACAAGACGTG TGTGGTACAAGACGTG       0.063845356           0.063845356
#> TGTGGTACACATCCGG TGTGGTACACATCCGG       1.000000000           1.000000000
#> TGTGGTAGTGCCTGCA TGTGGTAGTGCCTGCA       0.851769956           0.851769956
#> TGTGGTATCCAGTAGT TGTGGTATCCAGTAGT       0.309206159           0.309206159
#> TGTGTTTAGGAATGGA TGTGTTTAGGAATGGA       0.805462252           0.805462252
#> TGTGTTTGTAGATTAG TGTGTTTGTAGATTAG       0.812245408           0.812245408
#> TGTGTTTGTATGAATG TGTGTTTGTATGAATG       0.706510183           0.706510183
#> TTAACTCCAAAGAATC TTAACTCCAAAGAATC       0.030751379           0.030751379
#> TTAACTCCAGGATTGG TTAACTCCAGGATTGG       0.836700094           0.836700094
#> TTAACTCGTAATCGTC TTAACTCGTAATCGTC       0.754710003           0.754710003
#> TTAGGCAGTTTGTTGG TTAGGCAGTTTGTTGG       0.747063824           0.747063824
#> TTATGCTCACCCAGTG TTATGCTCACCCAGTG       0.836071312           0.836071312
#> TTCCCAGAGTTACCCA TTCCCAGAGTTACCCA       0.258589791           0.258589791
#> TTCCCAGGTATTCTCT TTCCCAGGTATTCTCT       0.174496549           0.174496549
#> TTCCCAGTCATAGCAC TTCCCAGTCATAGCAC       0.727991894           0.727991894
#> TTCGAAGAGTACGCCC TTCGAAGAGTACGCCC       0.280874789           0.280874789
#> TTCGGTCAGCTAGCCC TTCGGTCAGCTAGCCC       0.576021892           0.576021892
#> TTCGGTCAGGGCATGT TTCGGTCAGGGCATGT       0.865006418           0.865006418
#> TTCGGTCGTATAGTAG TTCGGTCGTATAGTAG       0.825424313           0.825424313
#> TTCGGTCTCGACCAGC TTCGGTCTCGACCAGC       0.782798506           0.782798506
#> TTCTACAAGATCGATA TTCTACAAGATCGATA       0.209004502           0.209004502
#> TTCTACACATAAAGGT TTCTACACATAAAGGT       0.790910696           0.790910696
#> TTCTACAGTACTCTCC TTCTACAGTACTCTCC       0.214774608           0.214774608
#> TTCTACATCGCTTAGA TTCTACATCGCTTAGA       0.868858910           0.868858910
#> TTCTACATCTGCAGTA TTCTACATCTGCAGTA       0.018174713           0.018174713
#> TTCTCAAAGCCGCCTA TTCTCAAAGCCGCCTA       0.156957165           0.156957165
#> TTCTCAACAATAGAGT TTCTCAACAATAGAGT       0.581354279           0.581354279
#> TTCTCCTAGTCCATAC TTCTCCTAGTCCATAC       0.817714516           0.817714516
#> TTCTCCTCAATGGTCT TTCTCCTCAATGGTCT       0.180058919           0.180058919
#> TTCTTAGTCAGTGTTG TTCTTAGTCAGTGTTG       0.841127589           0.841127589
#> TTGAACGAGAGCAATT TTGAACGAGAGCAATT       0.837399910           0.837399910
#> TTGAACGAGGAGTTTA TTGAACGAGGAGTTTA       0.446488496           0.446488496
#> TTGAACGGTGAGGGAG TTGAACGGTGAGGGAG       0.814308437           0.814308437
#> TTGAACGGTGATGCCC TTGAACGGTGATGCCC       0.212015692           0.212015692
#> TTGACTTAGACAAGCC TTGACTTAGACAAGCC       0.083756281           0.083756281
#> TTGACTTGTTCCATGA TTGACTTGTTCCATGA       0.606829782           0.606829782
#> TTGACTTTCACCGGGT TTGACTTTCACCGGGT       0.782097080           0.782097080
#> TTGACTTTCAGCAACT TTGACTTTCAGCAACT       0.183020164           0.183020164
#> TTGACTTTCCAGTAGT TTGACTTTCCAGTAGT       0.181631233           0.181631233
#> TTGCCGTCACCGAATT TTGCCGTCACCGAATT       0.825509085           0.825509085
#> TTGCCGTTCTCTGTCG TTGCCGTTCTCTGTCG       0.206223194           0.206223194
#> TTGCGTCCACCGATAT TTGCGTCCACCGATAT       0.567491356           0.567491356
#> TTGGAACAGAAGATTC TTGGAACAGAAGATTC       0.816863422           0.816863422
#> TTGGAACGTATATGAG TTGGAACGTATATGAG       0.847445435           0.847445435
#> TTGGCAAAGAGTAAGG TTGGCAAAGAGTAAGG       0.857736237           0.857736237
#> TTGTAGGTCAAAGTAG TTGTAGGTCAAAGTAG       0.672216233           0.672216233
#> TTTACTGCAGCATACT TTTACTGCAGCATACT       0.319484923           0.319484923
#> TTTACTGGTGCTAGCC TTTACTGGTGCTAGCC       0.206954883           0.206954883
#> TTTATGCAGGTTACCT TTTATGCAGGTTACCT       0.834964822           0.834964822
#> TTTATGCGTCGAACAG TTTATGCGTCGAACAG       0.047343165           0.047343165
#> TTTATGCTCCATTCTA TTTATGCTCCATTCTA       0.048675264           0.048675264
#> TTTCCTCCAGGATCGA TTTCCTCCAGGATCGA       0.860593870           0.860593870
#> TTTCCTCGTGGCGAAT TTTCCTCGTGGCGAAT       0.759060198           0.759060198
#> TTTGCGCAGACAGACC TTTGCGCAGACAGACC       0.186044789           0.186044789
#> TTTGGTTAGACAGGCT TTTGGTTAGACAGGCT       0.825797837           0.825797837
#> TTTGGTTCACCAGATT TTTGGTTCACCAGATT       0.186416522           0.186416522
#> TTTGGTTCACGAAGCA TTTGGTTCACGAAGCA       0.547896950           0.547896950
#> TTTGGTTTCCTTTCGG TTTGGTTTCCTTTCGG       0.806035228           0.806035228
#> TTTGTCAAGAATGTGT TTTGTCAAGAATGTGT       0.808762754           0.808762754
#> TTTGTCATCGAATGCT TTTGTCATCGAATGCT       0.833428092           0.833428092
#> TTTGTCATCTGTTTGT TTTGTCATCTGTTTGT       0.835677129           0.835677129
#>                    SubCellType
#> AAACCTGAGCCTTGAT        Ductal
#> AAACCTGGTAAGTGGC  Ngn3-high-EP
#> AAACGGGAGATATGGT        Ductal
#> AAACGGGCAAAGAATC          Beta
#> AAACGGGGTACAGTTC          Beta
#> AAACGGGTCAGCTCTC   Ngn3-low-EP
#> AAACGGGTCGCATGGC          Beta
#> AAAGATGAGTCATGCT          Beta
#> AAAGATGAGTCTTGCA          Beta
#> AAAGCAAAGCTACCTA Pre-endocrine
#> AAAGCAACACATGGGA        Ductal
#> AAAGCAAGTGCAGTAG          Beta
#> AAATGCCAGACTAAGT          Beta
#> AAATGCCCAATGGAAT         Alpha
#> AAATGCCCACCAACCG  Ngn3-high-EP
#> AAATGCCCATGCTAGT  Ngn3-high-EP
#> AAATGCCGTGATAAGT         Alpha
#> AACACGTGTCGACTGC Pre-endocrine
#> AACCATGCACAGGAGT  Ngn3-high-EP
#> AACCATGGTTCGGCAC  Ngn3-high-EP
#> AACCATGGTTGAACTC Pre-endocrine
#> AACGTTGGTGCCTTGG        Ductal
#> AACTCAGGTGCCTGTG  Ngn3-high-EP
#> AACTCAGGTTGGAGGT        Ductal
#> AACTCAGTCCTCGCAT  Ngn3-high-EP
#> AACTCAGTCTGGTTCC        Ductal
#> AACTCCCCAGTTCATG         Alpha
#> AACTCCCTCGAGAGCA        Ductal
#> AACTCTTCATTTGCCC          Beta
#> AACTGGTGTGCCTGTG        Ductal
#> AACTGGTTCTTACCGC        Ductal
#> AACTTTCCACAGTCGC  Ngn3-high-EP
#> AACTTTCCACTGCCAG        Ductal
#> AACTTTCGTGACCAAG   Ngn3-low-EP
#> AACTTTCTCCAGATCA         Alpha
#> AAGACCTAGACCTTTG         Alpha
#> AAGACCTAGTAGGCCA  Ngn3-high-EP
#> AAGACCTCAGGTTTCA        Ductal
#> AAGCCGCAGGCAGGTT        Ductal
#> AAGCCGCCAATGGACG  Ngn3-high-EP
#> AAGGAGCAGTACATGA         Alpha
#> AAGGCAGAGCTGAAAT        Ductal
#> AAGGCAGTCCGAGCCA   Ngn3-low-EP
#> AAGGTTCAGATAGTCA         Alpha
#> AAGGTTCAGGAGCGAG         Alpha
#> AAGGTTCGTCAACTGT Pre-endocrine
#> AAGTCTGGTTGGGACA          Beta
#> AAGTCTGTCCGAGCCA  Ngn3-high-EP
#> AATCCAGAGGGTTCCC Pre-endocrine
#> AATCCAGTCAGGTTCA Pre-endocrine
#> AATCCAGTCGGTCCGA Pre-endocrine
#> AATCGGTCATATGGTC  Ngn3-high-EP
#> AATCGGTTCCAAACTG  Ngn3-high-EP
#> AATCGGTTCGTATCAG         Alpha
#> ACACCAAGTATAGGTA          Beta
#> ACACCCTAGGACAGCT  Ngn3-high-EP
#> ACACCCTTCCCAGGTG Pre-endocrine
#> ACACCGGCATCACGTA        Ductal
#> ACACCGGCATGTCCTC        Ductal
#> ACACCGGTCATGTCTT          Beta
#> ACACTGAGTCCATGAT Pre-endocrine
#> ACAGCCGAGAGAACAG Pre-endocrine
#> ACAGCCGAGAGGTACC         Alpha
#> ACAGCCGAGTGGTCCC        Ductal
#> ACAGCCGCATCGTCGG  Ngn3-high-EP
#> ACAGCTAAGTCGTTTG         Alpha
#> ACAGCTATCAAGGCTT         Alpha
#> ACAGCTATCAGCATGT        Ductal
#> ACAGCTATCTGGGCCA        Ductal
#> ACATACGAGTAGGCCA  Ngn3-high-EP
#> ACATACGGTCTCTCTG Pre-endocrine
#> ACATCAGAGAGACTTA Pre-endocrine
#> ACATCAGTCCTTGGTC  Ngn3-high-EP
#> ACCAGTAAGTGATCGG         Alpha
#> ACCAGTAGTTCGTTGA          Beta
#> ACCAGTATCGTTTAGG  Ngn3-high-EP
#> ACCGTAAGTTATCACG        Ductal
#> ACCTTTATCCTGCAGG         Alpha
#> ACGAGCCCATGAGCGA          Beta
#> ACGAGGAAGTTACCCA   Ngn3-low-EP
#> ACGAGGATCATCACCC          Beta
#> ACGAGGATCGCCGTGA Pre-endocrine
#> ACGATACAGATGAGAG         Alpha
#> ACGATACTCGCGCCAA         Alpha
#> ACGATGTTCGCAGGCT Pre-endocrine
#> ACGCAGCCATTCTCAT Pre-endocrine
#> ACGCAGCGTAATCACC  Ngn3-high-EP
#> ACGCCAGGTTACAGAA         Alpha
#> ACGCCAGTCCCAAGTA  Ngn3-high-EP
#> ACGCCGAAGGGTGTTG          Beta
#> ACGCCGAGTTTGCATG          Beta
#> ACGGCCACACGTGAGA        Ductal
#> ACGGGCTGTGCGATAG        Ductal
#> ACGTCAAAGGTGTTAA          Beta
#> ACGTCAACACTGCCAG        Ductal
#> ACGTCAACAGACAAAT        Ductal
#> ACTATCTGTCAGGACA          Beta
#> ACTATCTTCAAGATCC Pre-endocrine
#> ACTGAACAGCCTCGTG        Ductal
#> ACTGAGTAGTGACTCT        Ductal
#> ACTGATGCATCAGTAC  Ngn3-high-EP
#> ACTGATGCATGGTCTA         Alpha
#> ACTGATGCATGTCTCC Pre-endocrine
#> ACTGATGGTGAGGGTT        Ductal
#> ACTGATGTCCCTAACC   Ngn3-low-EP
#> ACTGCTCAGGTGTTAA Pre-endocrine
#> ACTGTCCCAGTTTACG Pre-endocrine
#> ACTGTCCTCGCTAGCG  Ngn3-high-EP
#> ACTTACTAGAGGGCTT  Ngn3-high-EP
#> ACTTACTAGCCGGTAA   Ngn3-low-EP
#> ACTTACTCAATCAGAA  Ngn3-high-EP
#> ACTTACTCACATGTGT          Beta
#> ACTTGTTAGTGAAGTT          Beta
#> ACTTGTTAGTTAAGTG       Epsilon
#> ACTTGTTCATGGTCAT          Beta
#> ACTTGTTTCACCCGAG          Beta
#> ACTTTCAAGGATATAC        Ductal
#> ACTTTCACAACTGCTA  Ngn3-high-EP
#> AGAATAGCAATCACAC       Epsilon
#> AGAATAGCAGTACACT        Ductal
#> AGACGTTCAATGGTCT          Beta
#> AGACGTTTCTTACCTA   Ngn3-low-EP
#> AGAGCGACAGTAAGAT       Epsilon
#> AGAGCGATCTCTTGAT        Ductal
#> AGAGCTTAGCCGATTT        Ductal
#> AGAGCTTAGGTTACCT  Ngn3-high-EP
#> AGAGCTTGTTAAAGAC        Ductal
#> AGAGTGGAGGCGATAC          Beta
#> AGAGTGGCACCTCGGA Pre-endocrine
#> AGAGTGGGTAGCACGA        Ductal
#> AGAGTGGGTCACCTAA   Ngn3-low-EP
#> AGATCTGTCCGTCAAA Pre-endocrine
#> AGATTGCAGCGTGTCC         Alpha
#> AGATTGCAGCTGTCTA        Ductal
#> AGATTGCAGTTCCACA  Ngn3-high-EP
#> AGATTGCCACCAGGCT        Ductal
#> AGCAGCCGTCGCATCG  Ngn3-high-EP
#> AGCAGCCGTGTAACGG  Ngn3-high-EP
#> AGCAGCCTCTTTAGTC   Ngn3-low-EP
#> AGCATACCACAGCCCA       Epsilon
#> AGCATACGTCCGCTGA          Beta
#> AGCGGTCGTACCGAGA  Ngn3-high-EP
#> AGCGGTCGTCGAATCT  Ngn3-high-EP
#> AGCGGTCGTGTATGGG        Ductal
#> AGCGTCGCACAGGCCT  Ngn3-high-EP
#> AGCTCCTCAGTTTACG Pre-endocrine
#> AGCTCCTGTGTGCCTG        Ductal
#> AGCTCCTTCAATACCG  Ngn3-high-EP
#> AGCTCCTTCGCCCTTA Pre-endocrine
#> AGCTCTCAGACATAAC Pre-endocrine
#> AGCTCTCGTTCAGTAC         Delta
#> AGCTTGAAGACACGAC  Ngn3-high-EP
#> AGCTTGAAGAGTGACC          Beta
#> AGCTTGACAATAGAGT   Ngn3-low-EP
#> AGCTTGACAGTAAGAT        Ductal
#> AGCTTGACATCTATGG          Beta
#> AGCTTGAGTGCAGTAG          Beta
#> AGGCCACAGCTGTCTA   Ngn3-low-EP
#> AGGCCACCATGACATC          Beta
#> AGGCCGTCACGACGAA  Ngn3-high-EP
#> AGGCCGTCACGGTGTC          Beta
#> AGGCCGTGTACCAGTT   Ngn3-low-EP
#> AGGGAGTAGGGTATCG  Ngn3-high-EP
#> AGGGAGTAGTGTCTCA        Ductal
#> AGGGAGTTCCTCAATT          Beta
#> AGGGAGTTCGGAGGTA         Alpha
#> AGGGATGGTCGACTGC         Alpha
#> AGGGATGGTTAGGGTG  Ngn3-high-EP
#> AGGGATGTCGCATGAT          Beta
#> AGGGATGTCGCATGGC  Ngn3-high-EP
#> AGGGTGAGTGCGAAAC       Epsilon
#> AGGTCATCACTTCTGC         Alpha
#> AGGTCATGTCTCCCTA         Alpha
#> AGGTCATGTTTAGGAA Pre-endocrine
#> AGTAGTCCACATGACT   Ngn3-low-EP
#> AGTAGTCGTGAGGGAG          Beta
#> AGTAGTCGTTCTGTTT        Ductal
#> AGTCTTTCACCAGGTC        Ductal
#> AGTCTTTCACGGACAA         Alpha
#> AGTCTTTGTTACGGAG          Beta
#> AGTGAGGGTTACGGAG        Ductal
#> AGTGAGGTCTCCAACC          Beta
#> AGTGGGACAGGTCGTC          Beta
#> AGTGGGAGTGGCAAAC Pre-endocrine
#> AGTGTCATCCCTTGTG Pre-endocrine
#> AGTTGGTAGAGCCTAG         Alpha
#> AGTTGGTGTTGTCGCG Pre-endocrine
#> AGTTGGTTCACATAGC  Ngn3-high-EP
#> ATAAGAGCATCGATTG        Ductal
#> ATAGACCAGCGTGTCC         Alpha
#> ATAGACCAGTGATCGG         Alpha
#> ATAGACCTCTAACTGG        Ductal
#> ATCACGAGTGCACCAC          Beta
#> ATCATCTGTGTTTGTG         Delta
#> ATCATCTTCTGTCAAG        Ductal
#> ATCATGGAGGAATGGA         Alpha
#> ATCATGGGTATTCGTG         Alpha
#> ATCATGGGTCCGCTGA  Ngn3-high-EP
#> ATCATGGTCCAGAGGA        Ductal
#> ATCCACCAGAGCCTAG          Beta
#> ATCCACCCAGCTCGCA          Beta
#> ATCCACCGTACCATCA        Ductal
#> ATCCGAAGTACACCGC Pre-endocrine
#> ATCCGAAGTGACTCAT        Ductal
#> ATCTGCCCAATCCGAT         Delta
#> ATCTGCCCATGTTGAC          Beta
#> ATGAGGGTCTTGACGA Pre-endocrine
#> ATGTGTGCACAAGTAA  Ngn3-high-EP
#> ATGTGTGTCAAAGACA  Ngn3-high-EP
#> ATTACTCAGTTGCAGG        Ductal
#> ATTACTCTCAGGCGAA        Ductal
#> ATTACTCTCTGCCAGG         Alpha
#> ATTATCCAGAAACCTA  Ngn3-high-EP
#> ATTATCCAGACCCACC  Ngn3-high-EP
#> ATTATCCAGTGAACAT        Ductal
#> ATTGGACAGTAATCCC   Ngn3-low-EP
#> ATTGGTGTCCAGTATG        Ductal
#> ATTTCTGGTAGTGAAT   Ngn3-low-EP
#> CAACCAAAGCGTTCCG Pre-endocrine
#> CAACCAACAACCGCCA          Beta
#> CAACCAACACCTCGTT  Ngn3-high-EP
#> CAACCAATCGGGAGTA        Ductal
#> CAAGAAAAGAGTACCG  Ngn3-high-EP
#> CAAGAAATCGTACCGG        Ductal
#> CAAGGCCAGCTGTTCA        Ductal
#> CAAGGCCAGTGTGGCA        Ductal
#> CAAGGCCCACATTTCT  Ngn3-high-EP
#> CAAGGCCTCAGTGCAT         Alpha
#> CAAGGCCTCTCAACTT Pre-endocrine
#> CAAGGCCTCTTACCTA          Beta
#> CAAGTTGGTGTGCCTG Pre-endocrine
#> CAAGTTGTCTAACTGG Pre-endocrine
#> CACAAACAGTAGCGGT  Ngn3-high-EP
#> CACAAACCACGAAAGC          Beta
#> CACACAACAGCCAGAA        Ductal
#> CACACCTGTACATGTC          Beta
#> CACACCTGTTATCCGA Pre-endocrine
#> CACACTCAGATCCCAT        Ductal
#> CACACTCTCATGTGGT   Ngn3-low-EP
#> CACAGGCCAAGGCTCC Pre-endocrine
#> CACAGGCCACGTCAGC         Alpha
#> CACAGTAGTGCACTTA         Alpha
#> CACAGTAGTTCCTCCA        Ductal
#> CACAGTATCCTCCTAG        Ductal
#> CACATAGAGCGTGTCC        Ductal
#> CACATTTGTTTCGCTC        Ductal
#> CACATTTTCGAGGTAG          Beta
#> CACATTTTCGCAGGCT        Ductal
#> CACCACTAGCGATCCC         Alpha
#> CACCACTGTGTGTGCC Pre-endocrine
#> CACCAGGCATAAAGGT Pre-endocrine
#> CACCTTGAGAGTAATC        Ductal
#> CACTCCACATCGTCGG         Alpha
#> CACTCCATCAACACAC          Beta
#> CAGAATCCATTGGGCC         Alpha
#> CAGAATCGTCATTAGC  Ngn3-high-EP
#> CAGAATCGTGTAAGTA   Ngn3-low-EP
#> CAGAGAGTCGGCCGAT          Beta
#> CAGATCAAGCTGTTCA          Beta
#> CAGCAGCAGACAAAGG        Ductal
#> CAGCAGCCACACTGCG        Ductal
#> CAGCAGCGTTCACGGC        Ductal
#> CAGCCGAAGCGATATA        Ductal
#> CAGCCGAAGTAGGTGC        Ductal
#> CAGCCGATCCTATTCA  Ngn3-high-EP
#> CAGCCGATCGTTTATC  Ngn3-high-EP
#> CAGCGACCACAGGAGT        Ductal
#> CAGCTAAAGTCCATAC        Ductal
#> CAGCTAACATGCCTAA Pre-endocrine
#> CAGCTAAGTACCGAGA  Ngn3-high-EP
#> CAGCTGGCAAGTTGTC Pre-endocrine
#> CAGCTGGGTCTAACGT Pre-endocrine
#> CAGCTGGGTTCTCATT          Beta
#> CAGGTGCAGCCAGGAT          Beta
#> CAGGTGCCATCAGTCA Pre-endocrine
#> CAGTCCTAGGTGGGTT         Alpha
#> CAGTCCTCAGGGCATA         Alpha
#> CATATGGCATGTTGAC        Ductal
#> CATATGGGTTACTGAC          Beta
#> CATATGGTCACAGTAC Pre-endocrine
#> CATATTCTCCTGTAGA         Alpha
#> CATCAGACATGCATGT        Ductal
#> CATCAGAGTCTCTTTA Pre-endocrine
#> CATCAGATCATCGCTC Pre-endocrine
#> CATCGAAGTACGAAAT          Beta
#> CATCGAAGTCACACGC         Delta
#> CATCGGGCACGTTGGC       Epsilon
#> CATCGGGTCTTCGGTC  Ngn3-high-EP
#> CATGACACAACAACCT   Ngn3-low-EP
#> CATGACACATCGATTG  Ngn3-high-EP
#> CATGACAGTAGCTTGT   Ngn3-low-EP
#> CATGACATCATGCATG         Alpha
#> CATGCCTCACCCATTC         Alpha
#> CATGGCGAGAGTCTGG  Ngn3-high-EP
#> CATTATCGTCATATCG Pre-endocrine
#> CCAATCCAGAGTAAGG        Ductal
#> CCACCTAAGACTTGAA         Alpha
#> CCACCTACAGCTCGAC          Beta
#> CCACCTATCACCAGGC        Ductal
#> CCACCTATCTAAGCCA         Alpha
#> CCACGGAAGATCCTGT         Alpha
#> CCACGGAAGCAATCTC          Beta
#> CCACGGAAGGGTGTTG          Beta
#> CCACGGAAGTCTCCTC Pre-endocrine
#> CCACGGAGTCAACATC  Ngn3-high-EP
#> CCACGGAGTCTAGTCA          Beta
#> CCACGGAGTTACGACT  Ngn3-high-EP
#> CCACTACAGAGCTGGT Pre-endocrine
#> CCACTACAGAGTACCG   Ngn3-low-EP
#> CCACTACCAACGCACC         Alpha
#> CCAGCGACAGACGTAG Pre-endocrine
#> CCATGTCGTCGGCTCA  Ngn3-high-EP
#> CCATGTCTCCGCATAA  Ngn3-high-EP
#> CCATTCGCATGCCTAA        Ductal
#> CCATTCGGTAGGGACT        Ductal
#> CCCAATCCAACTGGCC         Alpha
#> CCCAGTTAGCATCATC          Beta
#> CCCAGTTGTCGTGGCT       Epsilon
#> CCCAGTTGTGCGCTTG          Beta
#> CCCAGTTTCAGCACAT          Beta
#> CCCATACAGACTTTCG          Beta
#> CCCATACCACATCCAA        Ductal
#> CCCATACCATGCAACT          Beta
#> CCCATACTCGTGGACC Pre-endocrine
#> CCGGGATAGATCTGAA  Ngn3-high-EP
#> CCGGGATTCGACAGCC         Alpha
#> CCGGTAGTCGGCTTGG          Beta
#> CCGTACTAGGAGTACC  Ngn3-high-EP
#> CCGTACTAGTCCGTAT  Ngn3-high-EP
#> CCGTACTTCAATAAGG   Ngn3-low-EP
#> CCGTGGACATAGAAAC          Beta
#> CCGTGGATCCAAGCCG        Ductal
#> CCGTGGATCTGTCTCG          Beta
#> CCGTGGATCTTAGAGC  Ngn3-high-EP
#> CCGTTCAAGCTAAACA         Alpha
#> CCGTTCAAGGACCACA  Ngn3-high-EP
#> CCGTTCACACATGACT   Ngn3-low-EP
#> CCTACACAGATCCGAG   Ngn3-low-EP
#> CCTACACAGTACGTTC         Alpha
#> CCTACACTCCTTTACA        Ductal
#> CCTACACTCTTTACAC         Alpha
#> CCTACCACATCGATTG  Ngn3-high-EP
#> CCTACCAGTTGCTCCT Pre-endocrine
#> CCTACCATCTGCCCTA          Beta
#> CCTAGCTCAGACTCGC        Ductal
#> CCTAGCTCATACTACG        Ductal
#> CCTAGCTGTGGCCCTA          Beta
#> CCTAGCTTCGATGAGG        Ductal
#> CCTATTACAGACAAGC         Delta
#> CCTATTATCATCGCTC  Ngn3-high-EP
#> CCTCAGTAGTGCGATG Pre-endocrine
#> CCTCAGTCAGACACTT       Epsilon
#> CCTCAGTGTTACCAGT  Ngn3-high-EP
#> CCTCTGAGTATAAACG        Ductal
#> CCTTACGAGATCACGG        Ductal
#> CCTTACGCACAACGCC Pre-endocrine
#> CCTTCCCAGGCAGGTT        Ductal
#> CCTTCCCCAAAGTGCG        Ductal
#> CCTTCCCTCGTCGTTC Pre-endocrine
#> CCTTCGAAGATAGCAT       Epsilon
#> CCTTCGACAGCTGTTA   Ngn3-low-EP
#> CCTTCGATCACAGGCC          Beta
#> CGAACATAGCCAGTAG        Ductal
#> CGAACATCAGACGTAG        Ductal
#> CGAATGTAGTTCCACA        Ductal
#> CGAATGTGTCTGGAGA Pre-endocrine
#> CGAATGTTCACCGGGT Pre-endocrine
#> CGAATGTTCTTCATGT          Beta
#> CGACCTTAGCGCTTAT   Ngn3-low-EP
#> CGACCTTAGGCTAGAC          Beta
#> CGACTTCCAGCCTATA        Ductal
#> CGACTTCCAGCTGCAC Pre-endocrine
#> CGACTTCTCACGCGGT         Alpha
#> CGAGAAGAGGGTTCCC          Beta
#> CGAGAAGCACGCATCG  Ngn3-high-EP
#> CGAGAAGTCGTTACAG         Alpha
#> CGAGCACAGCTGAAAT        Ductal
#> CGAGCACAGGCAGGTT Pre-endocrine
#> CGAGCACGTAGTAGTA Pre-endocrine
#> CGAGCACTCTACCTGC        Ductal
#> CGAGCACTCTTCGAGA   Ngn3-low-EP
#> CGATCGGAGCTGCAAG Pre-endocrine
#> CGATCGGCACGCCAGT  Ngn3-high-EP
#> CGATGGCGTTACGGAG        Ductal
#> CGATGGCTCAGTGTTG         Alpha
#> CGATGGCTCTGTGCAA        Ductal
#> CGATGTAAGTTTAGGA         Alpha
#> CGATGTAGTACTCTCC  Ngn3-high-EP
#> CGATTGAAGTCAAGGC Pre-endocrine
#> CGCCAAGCAAAGTCAA  Ngn3-high-EP
#> CGCCAAGCACCGGAAA        Ductal
#> CGCCAAGTCTAAGCCA          Beta
#> CGCGGTAAGAGCTGGT       Epsilon
#> CGCGGTAAGCTGCAAG          Beta
#> CGCGTTTCAGTCTTCC          Beta
#> CGCTGGAGTCCATCCT        Ductal
#> CGCTTCATCTCTGTCG         Alpha
#> CGGACACAGTGTCTCA        Ductal
#> CGGACGTAGACACGAC         Delta
#> CGGACTGAGGTAGCTG        Ductal
#> CGGACTGTCTTGAGAC          Beta
#> CGGAGCTAGATAGTCA          Beta
#> CGGAGCTCATTGGGCC        Ductal
#> CGGAGTCAGTTTAGGA         Alpha
#> CGGCTAGGTACTCGCG  Ngn3-high-EP
#> CGGGTCAAGGCGATAC         Alpha
#> CGGGTCAAGGCTCAGA        Ductal
#> CGGTTAAAGAGCTGCA        Ductal
#> CGGTTAAGTAAACACA   Ngn3-low-EP
#> CGGTTAAGTCTAGCCG  Ngn3-high-EP
#> CGGTTAATCACCGGGT Pre-endocrine
#> CGGTTAATCGTCTGAA        Ductal
#> CGTAGCGTCGGAAATA Pre-endocrine
#> CGTAGGCAGCGTAGTG        Ductal
#> CGTCACTAGGGTGTGT        Ductal
#> CGTCACTCACGAAATA   Ngn3-low-EP
#> CGTCACTCATAGTAAG   Ngn3-low-EP
#> CGTCACTTCACTTATC Pre-endocrine
#> CGTCACTTCTAACTCT       Epsilon
#> CGTCAGGAGCCTATGT   Ngn3-low-EP
#> CGTCAGGCATGGGAAC   Ngn3-low-EP
#> CGTCAGGGTAGGCATG  Ngn3-high-EP
#> CGTCAGGTCCGCAAGC         Delta
#> CGTCAGGTCCTGCCAT        Ductal
#> CGTCCATGTCTAACGT          Beta
#> CGTCCATTCGCGCCAA  Ngn3-high-EP
#> CGTCTACAGATCTGCT        Ductal
#> CGTCTACAGGGAGTAA         Alpha
#> CGTGAGCTCAGCGATT  Ngn3-high-EP
#> CGTGAGCTCTTCTGGC         Alpha
#> CGTGTAAGTCTTGATG   Ngn3-low-EP
#> CGTGTAATCATTATCC  Ngn3-high-EP
#> CGTGTAATCCATGCTC        Ductal
#> CGTGTCTCAGAAGCAC         Alpha
#> CGTGTCTGTAGCGTCC         Alpha
#> CGTTAGAAGTATCTCG        Ductal
#> CGTTAGAGTGCAACGA        Ductal
#> CGTTCTGAGGACGAAA         Alpha
#> CGTTCTGTCACAATGC  Ngn3-high-EP
#> CGTTCTGTCGATAGAA          Beta
#> CGTTGGGCAGGCGATA Pre-endocrine
#> CTAACTTAGAGGTACC         Alpha
#> CTAACTTAGATCACGG          Beta
#> CTAACTTCATCTGGTA         Alpha
#> CTAACTTCATTGGGCC        Ductal
#> CTAAGACAGCACGCCT  Ngn3-high-EP
#> CTAAGACAGCCACGCT   Ngn3-low-EP
#> CTAAGACGTATCAGTC         Alpha
#> CTAAGACTCTGCGGCA       Epsilon
#> CTAATGGCATGGTTGT Pre-endocrine
#> CTACACCCATCAGTCA          Beta
#> CTACACCGTCGAATCT        Ductal
#> CTACACCGTCGCCATG        Ductal
#> CTACATTCACCGTTGG Pre-endocrine
#> CTACATTGTTACTGAC Pre-endocrine
#> CTACCCAGTGCGCTTG        Ductal
#> CTACGTCAGCTAGGCA Pre-endocrine
#> CTAGAGTCACTGTTAG        Ductal
#> CTAGCCTGTCTCTTAT        Ductal
#> CTAGCCTGTTCGGGCT        Ductal
#> CTAGCCTTCAGCGATT        Ductal
#> CTAGTGATCAAAGACA          Beta
#> CTAGTGATCCCAGGTG          Beta
#> CTCACACGTCACACGC  Ngn3-high-EP
#> CTCAGAACACTAGTAC  Ngn3-high-EP
#> CTCAGAACATAGGATA          Beta
#> CTCAGAAGTCCTCTTG        Ductal
#> CTCAGAATCAGGTTCA   Ngn3-low-EP
#> CTCATTACAAGCGCTC        Ductal
#> CTCCTAGCACAACGTT Pre-endocrine
#> CTCGAAAGTAAGAGGA  Ngn3-high-EP
#> CTCGAAAGTGGTGTAG         Alpha
#> CTCGAAATCCTGCCAT          Beta
#> CTCGAGGCAAGTCTGT         Alpha
#> CTCGAGGGTAAGTAGT       Epsilon
#> CTCGAGGTCTGTTGAG  Ngn3-high-EP
#> CTCGGAGGTCATTAGC Pre-endocrine
#> CTCGGAGGTGAAGGCT Pre-endocrine
#> CTCGTACCAATGGATA        Ductal
#> CTCGTACTCTTTAGTC         Alpha
#> CTCGTCAAGAACAATC       Epsilon
#> CTCGTCACAACTGCTA        Ductal
#> CTCGTCATCGGTCTAA        Ductal
#> CTCTAATAGTTATCGC         Alpha
#> CTCTGGTAGCAATCTC  Ngn3-high-EP
#> CTCTGGTAGTCCCACG        Ductal
#> CTCTGGTCAACACCTA        Ductal
#> CTCTGGTCAAGACACG          Beta
#> CTCTGGTCAAGTTGTC Pre-endocrine
#> CTCTGGTCAGTTCATG   Ngn3-low-EP
#> CTGAAACAGACTAGGC        Ductal
#> CTGAAACGTAGCGCAA          Beta
#> CTGAAACGTATGGTTC         Delta
#> CTGAAACTCGATAGAA  Ngn3-high-EP
#> CTGAAACTCGCCGTGA        Ductal
#> CTGAAACTCTGTTTGT  Ngn3-high-EP
#> CTGAAGTGTTTGGGCC         Alpha
#> CTGAAGTTCTGCGGCA Pre-endocrine
#> CTGATCCTCAGTTGAC        Ductal
#> CTGCCTAGTCTCCATC        Ductal
#> CTGCGGACAGACTCGC  Ngn3-high-EP
#> CTGCGGACATTAGCCA Pre-endocrine
#> CTGCGGATCAGCATGT        Ductal
#> CTGCTGTAGACTGGGT Pre-endocrine
#> CTGCTGTCAAAGCAAT         Alpha
#> CTGGTCTAGGATTCGG Pre-endocrine
#> CTGGTCTAGGTGATTA       Epsilon
#> CTGGTCTCAGGAATGC        Ductal
#> CTGGTCTTCGCATGGC  Ngn3-high-EP
#> CTGTGCTCATTACCTT  Ngn3-high-EP
#> CTTAACTCAAGCGTAG        Ductal
#> CTTAACTTCAGCTTAG Pre-endocrine
#> CTTAGGACACAGGTTT        Ductal
#> CTTAGGAGTCCGAATT        Ductal
#> CTTAGGATCAACGGGA         Delta
#> CTTAGGATCACTCTTA  Ngn3-high-EP
#> CTTCTCTAGGCAGGTT        Ductal
#> CTTCTCTAGTCCGTAT Pre-endocrine
#> CTTCTCTCAGCATACT          Beta
#> CTTCTCTGTCAGAGGT       Epsilon
#> CTTGGCTAGACAGACC        Ductal
#> CTTGGCTAGGACACCA        Ductal
#> CTTGGCTCAAAGGAAG Pre-endocrine
#> CTTGGCTCAGCCACCA  Ngn3-high-EP
#> CTTGGCTCATCCCACT          Beta
#> CTTGGCTGTTAGAACA          Beta
#> CTTTGCGAGACCACGA Pre-endocrine
#> CTTTGCGAGGAATGGA   Ngn3-low-EP
#> CTTTGCGTCCACGTGG          Beta
#> GAAACTCAGAAACCTA          Beta
#> GAAATGAAGCTCCTCT  Ngn3-high-EP
#> GAAATGACAACTGGCC   Ngn3-low-EP
#> GAAATGACACATGACT        Ductal
#> GAAATGAGTCACCCAG  Ngn3-high-EP
#> GAAATGAGTCGGATCC        Ductal
#> GAACATCGTAATCGTC        Ductal
#> GAACCTAAGTAGATGT         Alpha
#> GAACCTATCACTTATC          Beta
#> GAACGGAAGTCAAGCG         Alpha
#> GAACGGACAAGGACAC          Beta
#> GAAGCAGGTCGCGAAA Pre-endocrine
#> GAAGCAGGTCTAGAGG       Epsilon
#> GAAGCAGTCAAGGCTT        Ductal
#> GAATGAAAGAATGTGT  Ngn3-high-EP
#> GAATGAACAATCGAAA        Ductal
#> GAATGAACAGCCTTGG        Ductal
#> GAATGAATCATCGGAT          Beta
#> GACACGCCAAAGGAAG        Ductal
#> GACCAATAGAATGTGT          Beta
#> GACCAATCAGACGCCT Pre-endocrine
#> GACCAATGTTATCCGA         Alpha
#> GACCAATGTTGGACCC          Beta
#> GACCAATTCAACGGCC  Ngn3-high-EP
#> GACCTGGGTTCGTGAT         Alpha
#> GACGCGTCAGTCGTGC        Ductal
#> GACGCGTGTCGGCTCA          Beta
#> GACGGCTAGGTGATAT        Ductal
#> GACGGCTTCAACTCTT Pre-endocrine
#> GACGTGCGTCACCTAA Pre-endocrine
#> GACGTGCTCCGGCACA         Alpha
#> GACGTTAAGACGCAAC        Ductal
#> GACGTTAAGACGCTTT  Ngn3-high-EP
#> GACGTTACATCTCGCT          Beta
#> GACGTTAGTGTGAATA Pre-endocrine
#> GACTAACGTCATCCCT  Ngn3-high-EP
#> GACTAACTCAGTCCCT        Ductal
#> GACTACAGTTCGTGAT        Ductal
#> GACTGCGAGAATCTCC  Ngn3-high-EP
#> GAGCAGAAGCGTGTCC        Ductal
#> GAGCAGAGTTTGCATG          Beta
#> GAGCAGATCAACTCTT Pre-endocrine
#> GAGGTGAGTCCGCTGA Pre-endocrine
#> GAGGTGAGTGACAAAT        Ductal
#> GAGTCCGAGCGATGAC  Ngn3-high-EP
#> GAGTCCGAGTTAAGTG        Ductal
#> GAGTCCGCAATTGCTG        Ductal
#> GATCAGTAGGTGATAT          Beta
#> GATCGATAGAAAGTGG Pre-endocrine
#> GATCGATAGGATGGTC         Alpha
#> GATCGATGTAGCGATG          Beta
#> GATCGATGTTGTTTGG         Alpha
#> GATCGCGCACCAGGCT Pre-endocrine
#> GATCGTAAGTAGGTGC Pre-endocrine
#> GATCGTAGTCCTAGCG        Ductal
#> GATCTAGAGCTAGCCC          Beta
#> GATCTAGGTGCCTTGG        Ductal
#> GATGAAAAGTTGTAGA   Ngn3-low-EP
#> GATGAAACAGGCTGAA         Alpha
#> GATGAAACATGGTTGT        Ductal
#> GATGAGGAGTGGACGT  Ngn3-high-EP
#> GATGAGGTCTGGGCCA        Ductal
#> GATGCTAGTAGAAGGA        Ductal
#> GATGCTAGTGTTAAGA Pre-endocrine
#> GATTCAGAGCGATGAC Pre-endocrine
#> GATTCAGAGTACGTTC         Alpha
#> GATTCAGTCGCCAGCA Pre-endocrine
#> GATTCAGTCGTCCAGG         Delta
#> GCAAACTCAGACAAAT          Beta
#> GCAATCACAGCTGCTG         Alpha
#> GCAATCATCCAACCAA        Ductal
#> GCACTCTCAGCGTCCA       Epsilon
#> GCAGCCACAGCGTTCG  Ngn3-high-EP
#> GCAGTTACATTTCACT Pre-endocrine
#> GCAGTTATCTTGAGAC        Ductal
#> GCATGATAGCCTTGAT  Ngn3-high-EP
#> GCATGATCAAGCGAGT          Beta
#> GCATGATGTCACTGGC         Alpha
#> GCATGATGTTTAGGAA          Beta
#> GCATGATTCAAACGGG        Ductal
#> GCATGCGAGGAACTGC Pre-endocrine
#> GCATGCGCAGCTCCGA          Beta
#> GCATGTAAGAATCTCC          Beta
#> GCATGTAAGCGTGTCC       Epsilon
#> GCATGTAGTTGCGCAC         Alpha
#> GCATGTATCGTGGACC         Alpha
#> GCCTCTAAGCTAGTCT         Delta
#> GCCTCTAAGTCCGGTC  Ngn3-high-EP
#> GCCTCTATCGCTTAGA   Ngn3-low-EP
#> GCGACCATCACCAGGC         Alpha
#> GCGAGAAGTGGGTATG         Alpha
#> GCGCAACGTTTCGCTC          Beta
#> GCGCAACTCACAACGT        Ductal
#> GCGCAGTGTGCTAGCC Pre-endocrine
#> GCGCCAAAGCTCCTTC         Alpha
#> GCGCCAATCTTATCTG Pre-endocrine
#> GCGCGATGTCCGAATT          Beta
#> GCGGGTTAGACTTGAA  Ngn3-high-EP
#> GCGGGTTCAAGTAATG         Alpha
#> GCGGGTTGTTACTGAC        Ductal
#> GCTCCTAAGGGCATGT          Beta
#> GCTCTGTGTCCATGAT Pre-endocrine
#> GCTGCAGGTTACCGAT        Ductal
#> GCTGCGAAGAACAACT         Alpha
#> GCTGCGAAGCTAGCCC   Ngn3-low-EP
#> GCTGCGAAGTATTGGA  Ngn3-high-EP
#> GCTGCGATCCCGGATG Pre-endocrine
#> GCTGCGATCGGCTACG          Beta
#> GCTGCTTGTCTTCAAG        Ductal
#> GCTGGGTGTCATACTG          Beta
#> GCTGGGTGTCTAGAGG  Ngn3-high-EP
#> GCTTCCAAGACATAAC         Alpha
#> GCTTCCAAGTTAAGTG         Alpha
#> GCTTGAAGTCAACTGT          Beta
#> GGAAAGCCAATGTAAG Pre-endocrine
#> GGAAAGCGTCTTGATG        Ductal
#> GGAATAAGTAAGAGAG          Beta
#> GGAATAAGTTCGGGCT        Ductal
#> GGACAAGGTGATAAGT        Ductal
#> GGACAAGTCAAAGACA  Ngn3-high-EP
#> GGACAGAAGGCCGAAT          Beta
#> GGACAGACAATCCAAC Pre-endocrine
#> GGACAGATCAGGCGAA          Beta
#> GGACAGATCCTGCAGG        Ductal
#> GGACAGATCCTTTCGG  Ngn3-high-EP
#> GGACATTGTCTCCCTA       Epsilon
#> GGACGTCAGGGTGTTG Pre-endocrine
#> GGACGTCCACGCATCG   Ngn3-low-EP
#> GGACGTCCATGGATGG  Ngn3-high-EP
#> GGACGTCGTCAGAAGC Pre-endocrine
#> GGACGTCGTCGAGATG       Epsilon
#> GGACGTCGTCTAGGTT          Beta
#> GGACGTCGTTTAGCTG        Ductal
#> GGAGCAAAGCTTCGCG        Ductal
#> GGAGCAAGTATCACCA Pre-endocrine
#> GGATGTTAGTAAGTAC         Alpha
#> GGATGTTTCCAGTATG  Ngn3-high-EP
#> GGCAATTCACGGCTAC         Alpha
#> GGCAATTCATCTCCCA  Ngn3-high-EP
#> GGCAATTGTAGAGCTG          Beta
#> GGCCGATTCACTCTTA   Ngn3-low-EP
#> GGCGACTGTAAGTTCC          Beta
#> GGCGTGTCACCAGGTC  Ngn3-high-EP
#> GGCGTGTCACCGATAT Pre-endocrine
#> GGCGTGTTCAGCTTAG         Alpha
#> GGCGTGTTCGAGAGCA Pre-endocrine
#> GGCTCGAAGATCACGG Pre-endocrine
#> GGCTCGAAGTGACTCT         Delta
#> GGCTCGACAGCTTCGG          Beta
#> GGCTGGTAGTCACGCC  Ngn3-high-EP
#> GGCTGGTTCATGCATG  Ngn3-high-EP
#> GGGAATGAGCACACAG       Epsilon
#> GGGAATGAGTGTACGG        Ductal
#> GGGAATGCAAAGTGCG          Beta
#> GGGAATGTCACAAACC  Ngn3-high-EP
#> GGGACCTCACTGTCGG Pre-endocrine
#> GGGACCTTCGAATGCT        Ductal
#> GGGAGATAGATGTGGC          Beta
#> GGGAGATAGCGAAGGG          Beta
#> GGGAGATAGTGACATA          Beta
#> GGGATGAAGGGCATGT         Alpha
#> GGGATGAGTCTGCGGT          Beta
#> GGGATGAGTTGTGGAG  Ngn3-high-EP
#> GGGCACTCATTTCAGG        Ductal
#> GGGCACTGTCTCCATC Pre-endocrine
#> GGGCACTGTGGCAAAC        Ductal
#> GGGCATCCATGACATC        Ductal
#> GGGCATCCATTTCAGG         Alpha
#> GGGCATCGTGGTCTCG  Ngn3-high-EP
#> GGGTCTGCACAGGCCT  Ngn3-high-EP
#> GGGTCTGCATGAAGTA  Ngn3-high-EP
#> GGGTCTGTCATGTGGT Pre-endocrine
#> GGGTCTGTCCATGAAC          Beta
#> GGGTTGCCAAGCTGAG   Ngn3-low-EP
#> GGGTTGCTCGAACGGA        Ductal
#> GGTATTGGTACTCGCG  Ngn3-high-EP
#> GGTATTGGTCCCGACA        Ductal
#> GGTGTTAAGGATCGCA        Ductal
#> GGTGTTAGTTCCTCCA        Ductal
#> GTAACGTGTCACTGGC        Ductal
#> GTAACTGGTCACTGGC Pre-endocrine
#> GTAACTGTCAGTGCAT        Ductal
#> GTACGTAGTCTCTTTA          Beta
#> GTACTCCCACAGGTTT Pre-endocrine
#> GTACTCCCACGTCTCT   Ngn3-low-EP
#> GTACTCCTCCAACCAA  Ngn3-high-EP
#> GTACTTTGTACCGTAT         Delta
#> GTACTTTTCTTCGGTC          Beta
#> GTAGGCCAGCTAACAA Pre-endocrine
#> GTAGGCCGTTACCAGT          Beta
#> GTAGTCACATCGATGT  Ngn3-high-EP
#> GTAGTCAGTCCGACGT  Ngn3-high-EP
#> GTATCTTCAGTGGGAT        Ductal
#> GTATCTTCATCGATGT  Ngn3-high-EP
#> GTATCTTGTACTCAAC         Alpha
#> GTATCTTGTCGCGGTT  Ngn3-high-EP
#> GTATCTTGTGGCCCTA  Ngn3-high-EP
#> GTATCTTTCGTGGACC Pre-endocrine
#> GTATTCTCACATTAGC        Ductal
#> GTATTCTCACTTGGAT        Ductal
#> GTATTCTGTCCGTCAG   Ngn3-low-EP
#> GTATTCTTCTCTAAGG         Delta
#> GTCAAGTGTGGTCCGT        Ductal
#> GTCACAACACGTCTCT  Ngn3-high-EP
#> GTCACGGAGCTATGCT       Epsilon
#> GTCACGGCATCCGTGG Pre-endocrine
#> GTCACGGTCTTGCCGT Pre-endocrine
#> GTCATTTGTCGTCTTC          Beta
#> GTCCTCACACCACGTG          Beta
#> GTCCTCATCATTTGGG  Ngn3-high-EP
#> GTCCTCATCTGGAGCC  Ngn3-high-EP
#> GTCGGGTCAGGCTGAA          Beta
#> GTCGTAAAGTTTAGGA       Epsilon
#> GTCTCGTAGGCTAGAC  Ngn3-high-EP
#> GTCTCGTCATAGTAAG  Ngn3-high-EP
#> GTCTCGTTCACTTCAT  Ngn3-high-EP
#> GTCTCGTTCTGCGTAA          Beta
#> GTCTTCGTCCGTTGCT         Alpha
#> GTCTTCGTCTACCAGA Pre-endocrine
#> GTGCATAAGAACAATC Pre-endocrine
#> GTGCATAAGCTCCTCT        Ductal
#> GTGCATAGTTCAACCA         Alpha
#> GTGCATATCAAGGCTT       Epsilon
#> GTGCGGTAGTAGGCCA        Ductal
#> GTGCGGTCAGCGTAAG  Ngn3-high-EP
#> GTGCTTCAGAGGGATA         Alpha
#> GTGCTTCTCACGACTA          Beta
#> GTGGGTCCATAAAGGT         Delta
#> GTGGGTCCATCGATTG          Beta
#> GTGGGTCGTCCAGTAT         Alpha
#> GTGGGTCTCCGAATGT          Beta
#> GTGTGCGGTGTAAGTA         Alpha
#> GTGTTAGCACCATCCT        Ductal
#> GTGTTAGGTGAACCTT Pre-endocrine
#> GTTAAGCAGGATTCGG  Ngn3-high-EP
#> GTTACAGTCCGATATG  Ngn3-high-EP
#> GTTCATTAGAAGATTC   Ngn3-low-EP
#> GTTCATTAGGCAGTCA Pre-endocrine
#> GTTCATTTCCGAACGC         Delta
#> GTTCTCGAGAATTCCC         Alpha
#> GTTCTCGCAAGAAAGG         Alpha
#> GTTCTCGCAGACACTT        Ductal
#> GTTTCTAAGAAGATTC         Delta
#> GTTTCTACAGGACCCT       Epsilon
#> TAAGAGATCTGCGACG Pre-endocrine
#> TAAGTGCAGGGTTTCT        Ductal
#> TAAGTGCTCCAAATGC          Beta
#> TACACGACAACGATCT        Ductal
#> TACACGACATGTAGTC         Alpha
#> TACACGAGTGCAACTT  Ngn3-high-EP
#> TACAGTGCAATCTACG   Ngn3-low-EP
#> TACAGTGCACAGGCCT         Alpha
#> TACAGTGGTAAACGCG Pre-endocrine
#> TACCTTAAGTACATGA        Ductal
#> TACCTTATCCGATATG          Beta
#> TACGGATAGAGCTGGT        Ductal
#> TACGGATAGGCACATG          Beta
#> TACGGATAGTTTAGGA   Ngn3-low-EP
#> TACGGATGTGTTCGAT        Ductal
#> TACGGATTCCAAGTAC  Ngn3-high-EP
#> TACGGGCAGCCAACAG          Beta
#> TACGGGCAGCGTCAAG        Ductal
#> TACGGGCTCTTCGGTC          Beta
#> TACGGTACATGTTCCC        Ductal
#> TACTCATAGGTGTTAA Pre-endocrine
#> TACTCATCAGATCCAT          Beta
#> TACTCGCAGTAAGTAC         Alpha
#> TACTCGCCAGATCGGA          Beta
#> TACTCGCGTCGGCATC          Beta
#> TACTTACCAACTGCGC        Ductal
#> TACTTACCATGGGACA  Ngn3-high-EP
#> TACTTGTAGACAAAGG          Beta
#> TAGACCAAGGAATCGC        Ductal
#> TAGAGCTGTATGAAAC        Ductal
#> TAGAGCTGTTGTCGCG        Ductal
#> TAGCCGGTCAACACAC        Ductal
#> TAGGCATTCCAAAGTC   Ngn3-low-EP
#> TAGTGGTCAACTGCGC        Ductal
#> TAGTGGTCATTCACTT Pre-endocrine
#> TAGTGGTTCCTCAACC Pre-endocrine
#> TAGTTGGCACGTAAGG  Ngn3-high-EP
#> TAGTTGGCAGTAACGG Pre-endocrine
#> TAGTTGGGTACCGTAT        Ductal
#> TAGTTGGTCTTGCAAG          Beta
#> TATCAGGCACCAGATT        Ductal
#> TATCTCACAACTGCTA        Ductal
#> TATCTCACACGCGAAA         Alpha
#> TATCTCACATTGTGCA          Beta
#> TATCTCAGTTCGCTAA        Ductal
#> TATGCCCAGAACAATC Pre-endocrine
#> TATGCCCGTCTAGCCG  Ngn3-high-EP
#> TATGCCCTCTTACCTA   Ngn3-low-EP
#> TATTACCTCGTTGCCT         Alpha
#> TCAACGATCAACCATG Pre-endocrine
#> TCAATCTCATCCGCGA  Ngn3-high-EP
#> TCAATCTCATCTCCCA   Ngn3-low-EP
#> TCAATCTGTCAAAGAT   Ngn3-low-EP
#> TCAATCTTCGAACGGA         Alpha
#> TCACAAGAGTGACATA          Beta
#> TCACAAGCAACGATGG         Alpha
#> TCACAAGCAGACTCGC          Beta
#> TCACGAAGTCCGAGTC  Ngn3-high-EP
#> TCACGAATCCTATGTT   Ngn3-low-EP
#> TCAGATGCACTCGACG Pre-endocrine
#> TCAGATGGTACACCGC          Beta
#> TCAGCAACACCAGGCT        Ductal
#> TCAGCAATCGCCGTGA Pre-endocrine
#> TCAGCTCAGTCGTACT   Ngn3-low-EP
#> TCAGGTAAGGCGTACA  Ngn3-high-EP
#> TCAGGTAAGTCAATAG        Ductal
#> TCAGGTACAGTCACTA         Alpha
#> TCAGGTATCTGTCCGT          Beta
#> TCATTACGTCAGAAGC        Ductal
#> TCATTTGCATTGGCGC Pre-endocrine
#> TCCACACCAGGCGATA         Alpha
#> TCCCGATAGTCGTTTG Pre-endocrine
#> TCCCGATCAACGATGG  Ngn3-high-EP
#> TCGAGGCAGAACAATC         Alpha
#> TCGAGGCAGGCATGGT Pre-endocrine
#> TCGCGAGAGTACGCGA Pre-endocrine
#> TCGCGAGAGTCCGTAT          Beta
#> TCGCGAGCACTTCTGC Pre-endocrine
#> TCGGGACAGAGAACAG  Ngn3-high-EP
#> TCGGGACGTCTGCGGT       Epsilon
#> TCGTACCCACGCATCG   Ngn3-low-EP
#> TCGTACCCATATACGC        Ductal
#> TCTATTGCACAGGAGT        Ductal
#> TCTATTGTCCTATTCA       Epsilon
#> TCTATTGTCGTATCAG Pre-endocrine
#> TCTCATAAGAATGTGT Pre-endocrine
#> TCTCATACAAATTGCC          Beta
#> TCTCATAGTGGACGAT  Ngn3-high-EP
#> TCTCTAAAGCTGGAAC          Beta
#> TCTCTAAGTCGCATCG        Ductal
#> TCTCTAAGTGGCAAAC  Ngn3-high-EP
#> TCTCTAAGTTACCGAT        Ductal
#> TCTCTAATCCCAACGG Pre-endocrine
#> TCTGAGAGTCAAACTC        Ductal
#> TCTGGAACACACGCTG        Ductal
#> TCTGGAACACCTATCC          Beta
#> TCTGGAAGTCATCGGC          Beta
#> TCTGGAATCGGTGTCG        Ductal
#> TCTTTCCGTGCAACTT          Beta
#> TCTTTCCTCGTCTGCT   Ngn3-low-EP
#> TGAAAGAGTTATCACG   Ngn3-low-EP
#> TGACGGCAGGAATGGA        Ductal
#> TGACGGCCATATGAGA          Beta
#> TGACTAGGTGTATGGG  Ngn3-high-EP
#> TGACTAGTCCACGTGG   Ngn3-low-EP
#> TGACTAGTCCTATTCA  Ngn3-high-EP
#> TGACTTTAGCCAGAAC          Beta
#> TGACTTTAGCGTTGCC        Ductal
#> TGACTTTAGCGTTTAC        Ductal
#> TGAGAGGAGATCCCAT          Beta
#> TGAGAGGAGCCAGGAT       Epsilon
#> TGAGAGGAGCCGTCGT  Ngn3-high-EP
#> TGAGAGGAGGCGATAC       Epsilon
#> TGAGAGGAGGTAAACT       Epsilon
#> TGAGAGGCAGGAACGT          Beta
#> TGAGAGGGTAAGGGCT        Ductal
#> TGAGCATGTTAAAGAC  Ngn3-high-EP
#> TGAGCATTCACCCGAG        Ductal
#> TGAGCATTCGGGAGTA        Ductal
#> TGAGCCGTCATGTCCC        Ductal
#> TGAGCCGTCCGAACGC        Ductal
#> TGAGCCGTCCGCTGTT  Ngn3-high-EP
#> TGAGGGAAGAGGTACC  Ngn3-high-EP
#> TGAGGGACAGACTCGC        Ductal
#> TGATTTCAGTACGCCC  Ngn3-high-EP
#> TGATTTCAGTGACTCT        Ductal
#> TGATTTCGTGAGGGTT         Delta
#> TGCCAAACATGACATC  Ngn3-high-EP
#> TGCCCATCACTGTGTA   Ngn3-low-EP
#> TGCCCTAAGTCCGGTC Pre-endocrine
#> TGCGCAGGTACGCACC        Ductal
#> TGCGCAGTCTGCCAGG   Ngn3-low-EP
#> TGCGGGTAGACGCACA  Ngn3-high-EP
#> TGCGGGTAGGGTTCCC         Alpha
#> TGCGGGTTCCACTCCA       Epsilon
#> TGCGGGTTCTCAAACG        Ductal
#> TGCGTGGAGGTAAACT   Ngn3-low-EP
#> TGCGTGGAGGTTCCTA Pre-endocrine
#> TGCGTGGAGTACACCT          Beta
#> TGCGTGGGTCCCTTGT Pre-endocrine
#> TGCGTGGGTTACTGAC Pre-endocrine
#> TGCTACCCAGGGCATA   Ngn3-low-EP
#> TGCTACCCATAAGACA  Ngn3-high-EP
#> TGCTGCTAGGATGCGT  Ngn3-high-EP
#> TGCTGCTGTCCGTCAG        Ductal
#> TGGACGCAGTGTTGAA   Ngn3-low-EP
#> TGGACGCTCAGGTAAA        Ductal
#> TGGACGCTCTTGTTTG Pre-endocrine
#> TGGCCAGAGCTACCTA        Ductal
#> TGGCCAGCAAGGACTG          Beta
#> TGGCCAGCAATGAAAC Pre-endocrine
#> TGGCCAGCACGAAAGC        Ductal
#> TGGCCAGGTCGGCATC        Ductal
#> TGGCGCATCACAACGT        Ductal
#> TGGCTGGGTACCCAAT          Beta
#> TGGGAAGCACGTAAGG Pre-endocrine
#> TGGGAAGCATGACATC        Ductal
#> TGGGCGTGTGTGCGTC         Alpha
#> TGGTTAGCAGTCGTGC Pre-endocrine
#> TGGTTCCCACCGCTAG        Ductal
#> TGGTTCCGTCACCTAA        Ductal
#> TGTATTCCAATGACCT        Ductal
#> TGTCCCAAGTTTAGGA   Ngn3-low-EP
#> TGTCCCACATTATCTC       Epsilon
#> TGTCCCAGTCCGTGAC          Beta
#> TGTCCCATCTCTGCTG       Epsilon
#> TGTGGTACAAGACGTG        Ductal
#> TGTGGTACACATCCGG          Beta
#> TGTGGTAGTGCCTGCA       Epsilon
#> TGTGGTATCCAGTAGT  Ngn3-high-EP
#> TGTGTTTAGGAATGGA          Beta
#> TGTGTTTGTAGATTAG          Beta
#> TGTGTTTGTATGAATG Pre-endocrine
#> TTAACTCCAAAGAATC        Ductal
#> TTAACTCCAGGATTGG       Epsilon
#> TTAACTCGTAATCGTC Pre-endocrine
#> TTAGGCAGTTTGTTGG Pre-endocrine
#> TTATGCTCACCCAGTG          Beta
#> TTCCCAGAGTTACCCA  Ngn3-high-EP
#> TTCCCAGGTATTCTCT   Ngn3-low-EP
#> TTCCCAGTCATAGCAC Pre-endocrine
#> TTCGAAGAGTACGCCC  Ngn3-high-EP
#> TTCGGTCAGCTAGCCC  Ngn3-high-EP
#> TTCGGTCAGGGCATGT         Delta
#> TTCGGTCGTATAGTAG       Epsilon
#> TTCGGTCTCGACCAGC Pre-endocrine
#> TTCTACAAGATCGATA   Ngn3-low-EP
#> TTCTACACATAAAGGT Pre-endocrine
#> TTCTACAGTACTCTCC   Ngn3-low-EP
#> TTCTACATCGCTTAGA          Beta
#> TTCTACATCTGCAGTA        Ductal
#> TTCTCAAAGCCGCCTA        Ductal
#> TTCTCAACAATAGAGT  Ngn3-high-EP
#> TTCTCCTAGTCCATAC          Beta
#> TTCTCCTCAATGGTCT        Ductal
#> TTCTTAGTCAGTGTTG       Epsilon
#> TTGAACGAGAGCAATT          Beta
#> TTGAACGAGGAGTTTA  Ngn3-high-EP
#> TTGAACGGTGAGGGAG          Beta
#> TTGAACGGTGATGCCC        Ductal
#> TTGACTTAGACAAGCC        Ductal
#> TTGACTTGTTCCATGA  Ngn3-high-EP
#> TTGACTTTCACCGGGT Pre-endocrine
#> TTGACTTTCAGCAACT        Ductal
#> TTGACTTTCCAGTAGT        Ductal
#> TTGCCGTCACCGAATT         Alpha
#> TTGCCGTTCTCTGTCG        Ductal
#> TTGCGTCCACCGATAT  Ngn3-high-EP
#> TTGGAACAGAAGATTC         Alpha
#> TTGGAACGTATATGAG          Beta
#> TTGGCAAAGAGTAAGG       Epsilon
#> TTGTAGGTCAAAGTAG  Ngn3-high-EP
#> TTTACTGCAGCATACT  Ngn3-high-EP
#> TTTACTGGTGCTAGCC        Ductal
#> TTTATGCAGGTTACCT         Alpha
#> TTTATGCGTCGAACAG        Ductal
#> TTTATGCTCCATTCTA        Ductal
#> TTTCCTCCAGGATCGA         Delta
#> TTTCCTCGTGGCGAAT Pre-endocrine
#> TTTGCGCAGACAGACC        Ductal
#> TTTGGTTAGACAGGCT         Alpha
#> TTTGGTTCACCAGATT        Ductal
#> TTTGGTTCACGAAGCA  Ngn3-high-EP
#> TTTGGTTTCCTTTCGG          Beta
#> TTTGTCAAGAATGTGT Pre-endocrine
#> TTTGTCATCGAATGCT         Alpha
#> TTTGTCATCTGTTTGT       Epsilon
#> 
#> $feature_metadata
#>      features stochastic_pseudotimepeaktime stochastic_pseudotimeexp_ncells
#> Pyy       Pyy                     0.8379160                             665
#> Rbp4     Rbp4                     0.8679776                             647
#> Sst       Sst                     0.8792311                              64
#>      stochastic_pseudotimer.sq stochastic_pseudotimedev.expl
#> Pyy                 0.34620653                     0.6992168
#> Rbp4                0.43812557                     0.7368411
#> Sst                 0.06459176                     0.5798457
#>      stochastic_pseudotimepadjust  peaktime feature_split index
#> Pyy                             0 0.8379160            NA     1
#> Rbp4                            0 0.8679776            NA     2
#> Sst                             0 0.8792311            NA     3
#> 
#> $enrichment
#> NULL
#> 

DynamicPlot(
  pancreas_sub,
  lineages = "stochastic_pseudotime",
  features = c("Rbp4", "Sst", "Pyy"),
  assay = "RNA",
  group.by = "SubCellType"
)
```
