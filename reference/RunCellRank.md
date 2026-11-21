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
#> ℹ [2026-11-21 11:01:49] Start standard processing workflow...
#> ℹ [2026-11-21 11:01:49] Checking a list of <Seurat>...
#> ! [2026-11-21 11:01:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-21 11:01:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:01:49] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:01:49] Use the separate HVF from `srt_list`
#> ℹ [2026-11-21 11:01:49] Number of available HVF: 2000
#> ℹ [2026-11-21 11:01:49] Finished check
#> ℹ [2026-11-21 11:01:49] Perform `ScaleData()`
#> ℹ [2026-11-21 11:01:49] Perform pca linear dimension reduction
#> ℹ [2026-11-21 11:01:50] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-21 11:01:50] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-21 11:01:50] Reorder clusters...
#> ℹ [2026-11-21 11:01:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-21 11:01:50] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-21 11:01:58] Standard processing workflow completed
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
#> ℹ [2026-11-21 11:01:58] Running scanpy-compatible preprocessing (15998 features -> filter + normalize)...
#> ℹ [2026-11-21 11:02:03] Running scVelo "stochastic" mode with `backend = 'cpp'` (9699 features)
#> ✔ [2026-11-21 11:02:07] scVelo "stochastic" mode completed
#> ✔ [2026-11-21 11:02:07] scVelo cpp backend completed
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
#> ! [2026-11-21 11:02:07] RunCellRank(backend = "cpp") uses an approximate C++ implementation; results may not be completely identical to the reference python workflow.
#> ✔ [2026-11-21 11:02:09] CellRank cpp backend ("gpcca" + "pseudotime_connectivity_combined" kernel) completed

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
#> ℹ [2026-11-21 11:02:10] Start find dynamic features
#> ℹ [2026-11-21 11:02:10] Data type is raw counts
#> ℹ [2026-11-21 11:02:10] Number of candidate features (union): 3
#> ℹ [2026-11-21 11:02:11] Data type is raw counts
#> ℹ [2026-11-21 11:02:11] Calculating dynamic features for "stochastic_pseudotime"...
#> ℹ [2026-11-21 11:02:11] Using 1 core
#> ⠙ [2026-11-21 11:02:11] Running for 1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-11-21 11:02:11] Completed 3 tasks in 254ms
#> 
#> ℹ [2026-11-21 11:02:11] Building results
#> ✔ [2026-11-21 11:02:11] Find dynamic features done

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
#>      GCATGCGAGGAACTGCstochastic_pseudotime
#> Pyy                             -0.6430859
#> Sst                             -0.3143440
#> Rbp4                             3.4701665
#>      CATCGGGCACGTTGGCstochastic_pseudotime
#> Pyy                             -0.8318066
#> Sst                             -0.3143440
#> Rbp4                             0.5910230
#>      CTCGGAGGTCATTAGCstochastic_pseudotime
#> Pyy                             -0.8533196
#> Sst                             -0.3143440
#> Rbp4                             0.3431818
#>      GGACATTGTCTCCCTAstochastic_pseudotime
#> Pyy                           -0.886108304
#> Sst                           -0.314343934
#> Rbp4                          -0.009974129
#>      AGCTCCTTCGCCCTTAstochastic_pseudotime
#> Pyy                            -0.88700350
#> Sst                            -0.31434393
#> Rbp4                           -0.01926689
#>      CAGGTGCCATCAGTCAstochastic_pseudotime
#> Pyy                             -0.9088828
#> Sst                             -0.3143438
#> Rbp4                            -0.2429373
#>      AGCTCCTCAGTTTACGstochastic_pseudotime
#> Pyy                             -0.9092423
#> Sst                             -0.3143438
#> Rbp4                            -0.2465871
#>      GACGGCTTCAACTCTTstochastic_pseudotime
#> Pyy                             -0.9345826
#> Sst                             -0.3143427
#> Rbp4                            -0.5164182
#>      GTCTTCGTCTACCAGAstochastic_pseudotime
#> Pyy                             -0.9428299
#> Sst                             -0.3143372
#> Rbp4                            -0.6540374
#>      ACGAGGATCGCCGTGAstochastic_pseudotime
#> Pyy                             -0.9429956
#> Sst                             -0.3143324
#> Rbp4                            -0.6866626
#>      TCCCGATAGTCGTTTGstochastic_pseudotime
#> Pyy                             -0.9394630
#> Sst                             -0.3143126
#> Rbp4                            -0.7371833
#>      GTTCATTAGGCAGTCAstochastic_pseudotime
#> Pyy                             -0.9391759
#> Sst                             -0.3143112
#> Rbp4                            -0.7390683
#>      TCATTTGCATTGGCGCstochastic_pseudotime
#> Pyy                             -0.9332575
#> Sst                             -0.3142800
#> Rbp4                            -0.7641347
#>      ACATCAGAGAGACTTAstochastic_pseudotime
#> Pyy                             -0.9296658
#> Sst                             -0.3142586
#> Rbp4                            -0.7730127
#>      GTACTCCCACAGGTTTstochastic_pseudotime
#> Pyy                             -0.9286318
#> Sst                             -0.3142520
#> Rbp4                            -0.7750953
#>      CAGCTGGGTCTAACGTstochastic_pseudotime
#> Pyy                             -0.9210915
#> Sst                             -0.3141983
#> Rbp4                            -0.7861799
#>      GCGCCAATCTTATCTGstochastic_pseudotime
#> Pyy                             -0.9187033
#> Sst                             -0.3141791
#> Rbp4                            -0.7886115
#>      AATCCAGTCGGTCCGAstochastic_pseudotime
#> Pyy                             -0.9134520
#> Sst                             -0.3141332
#> Rbp4                            -0.7927275
#>      AATCCAGTCAGGTTCAstochastic_pseudotime
#> Pyy                             -0.9089249
#> Sst                             -0.3140895
#> Rbp4                            -0.7952223
#>      CCTACCAGTTGCTCCTstochastic_pseudotime
#> Pyy                             -0.9031817
#> Sst                             -0.3140285
#> Rbp4                            -0.7973325
#>      CAAGTTGTCTAACTGGstochastic_pseudotime
#> Pyy                             -0.8989578
#> Sst                             -0.3139797
#> Rbp4                            -0.7982893
#>      GAAGCAGGTCTAGAGGstochastic_pseudotime
#> Pyy                             -0.8978829
#> Sst                             -0.3139667
#> Rbp4                            -0.7984645
#>      CTGCTGTAGACTGGGTstochastic_pseudotime
#> Pyy                             -0.8781630
#> Sst                             -0.3136903
#> Rbp4                            -0.7980583
#>      CTTAACTTCAGCTTAGstochastic_pseudotime
#> Pyy                             -0.8612381
#> Sst                             -0.3133963
#> Rbp4                            -0.7940160
#>      CTAATGGCATGGTTGTstochastic_pseudotime
#> Pyy                             -0.8568632
#> Sst                             -0.3133120
#> Rbp4                            -0.7925997
#>      CATCAGATCATCGCTCstochastic_pseudotime
#> Pyy                             -0.8561447
#> Sst                             -0.3132979
#> Rbp4                            -0.7923550
#>      AGGTCATGTTTAGGAAstochastic_pseudotime
#> Pyy                             -0.8523749
#> Sst                             -0.3132221
#> Rbp4                            -0.7910190
#>      AATCCAGAGGGTTCCCstochastic_pseudotime
#> Pyy                             -0.8480249
#> Sst                             -0.3131318
#> Rbp4                            -0.7893753
#>      TCTCATAAGAATGTGTstochastic_pseudotime
#> Pyy                             -0.8417794
#> Sst                             -0.3129965
#> Rbp4                            -0.7868426
#>      GTCACGGCATCCGTGGstochastic_pseudotime
#> Pyy                             -0.8397175
#> Sst                             -0.3129505
#> Rbp4                            -0.7859656
#>      AGATCTGTCCGTCAAAstochastic_pseudotime
#> Pyy                             -0.8380234
#> Sst                             -0.3129121
#> Rbp4                            -0.7852308
#>      TAGTTGGCAGTAACGGstochastic_pseudotime
#> Pyy                             -0.8379676
#> Sst                             -0.3129109
#> Rbp4                            -0.7852063
#>      CGACTTCCAGCTGCACstochastic_pseudotime
#> Pyy                             -0.8371879
#> Sst                             -0.3128931
#> Rbp4                            -0.7848638
#>      TTAACTCGTAATCGTCstochastic_pseudotime
#> Pyy                             -0.8366541
#> Sst                             -0.3128808
#> Rbp4                            -0.7846278
#>      GTAACTGGTCACTGGCstochastic_pseudotime
#> Pyy                             -0.8366172
#> Sst                             -0.3128800
#> Rbp4                            -0.7846114
#>      AAAGCAAAGCTACCTAstochastic_pseudotime
#> Pyy                             -0.8338556
#> Sst                             -0.3128159
#> Rbp4                            -0.7833712
#>      CTTGGCTCAAAGGAAGstochastic_pseudotime
#> Pyy                             -0.8310580
#> Sst                             -0.3127497
#> Rbp4                            -0.7820839
#>      AGGGTGAGTGCGAAACstochastic_pseudotime
#> Pyy                             -0.8308390
#> Sst                             -0.3127445
#> Rbp4                            -0.7819818
#>      CATTATCGTCATATCGstochastic_pseudotime
#> Pyy                             -0.8305965
#> Sst                             -0.3127387
#> Rbp4                            -0.7818686
#>      TACTCATAGGTGTTAAstochastic_pseudotime
#> Pyy                             -0.8269007
#> Sst                             -0.3126491
#> Rbp4                            -0.7801167
#>      TTAGGCAGTTTGTTGGstochastic_pseudotime
#> Pyy                             -0.8220001
#> Sst                             -0.3125272
#> Rbp4                            -0.7777204
#>      TCAGCAATCGCCGTGAstochastic_pseudotime
#> Pyy                             -0.8175592
#> Sst                             -0.3124137
#> Rbp4                            -0.7754822
#>      CGTCACTTCACTTATCstochastic_pseudotime
#> Pyy                             -0.8148844
#> Sst                             -0.3123439
#> Rbp4                            -0.7741054
#>      CGTTGGGCAGGCGATAstochastic_pseudotime
#> Pyy                             -0.8129489
#> Sst                             -0.3122927
#> Rbp4                            -0.7730964
#>      AACACGTGTCGACTGCstochastic_pseudotime
#> Pyy                             -0.8119731
#> Sst                             -0.3122667
#> Rbp4                            -0.7725837
#>      CCTTCCCTCGTCGTTCstochastic_pseudotime
#> Pyy                             -0.8119647
#> Sst                             -0.3122665
#> Rbp4                            -0.7725794
#>      CTACATTGTTACTGACstochastic_pseudotime
#> Pyy                             -0.8073563
#> Sst                             -0.3121420
#> Rbp4                            -0.7701239
#>      GAAGCAGGTCGCGAAAstochastic_pseudotime
#> Pyy                             -0.7976514
#> Sst                             -0.3118703
#> Rbp4                            -0.7647822
#>      GCAGTTACATTTCACTstochastic_pseudotime
#> Pyy                             -0.7974501
#> Sst                             -0.3118645
#> Rbp4                            -0.7646691
#>      GATCGCGCACCAGGCTstochastic_pseudotime
#> Pyy                             -0.7942238
#> Sst                             -0.3117714
#> Rbp4                            -0.7628450
#>      TATGCCCAGAACAATCstochastic_pseudotime
#> Pyy                             -0.7935893
#> Sst                             -0.3117529
#> Rbp4                            -0.7624837
#>      CAGCTGGCAAGTTGTCstochastic_pseudotime
#> Pyy                             -0.7933683
#> Sst                             -0.3117464
#> Rbp4                            -0.7623577
#>      CAACCAAAGCGTTCCGstochastic_pseudotime
#> Pyy                             -0.7916688
#> Sst                             -0.3116966
#> Rbp4                            -0.7613851
#>      CTACATTCACCGTTGGstochastic_pseudotime
#> Pyy                             -0.7828252
#> Sst                             -0.3114317
#> Rbp4                            -0.7562342
#>      TTGACTTTCACCGGGTstochastic_pseudotime
#> Pyy                             -0.7797719
#> Sst                             -0.3113381
#> Rbp4                            -0.7544229
#>      CCTTACGCACAACGCCstochastic_pseudotime
#> Pyy                             -0.7794050
#> Sst                             -0.3113267
#> Rbp4                            -0.7542042
#>      GTATCTTTCGTGGACCstochastic_pseudotime
#> Pyy                             -0.7782659
#> Sst                             -0.3112915
#> Rbp4                            -0.7535237
#>      CGAGCACAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.7782124
#> Sst                             -0.3112898
#> Rbp4                            -0.7534917
#>      CGGTTAATCACCGGGTstochastic_pseudotime
#> Pyy                             -0.7764164
#> Sst                             -0.3112339
#> Rbp4                            -0.7524142
#>      GGGCACTGTCTCCATCstochastic_pseudotime
#> Pyy                             -0.7750363
#> Sst                             -0.3111906
#> Rbp4                            -0.7515827
#>      ACATACGGTCTCTCTGstochastic_pseudotime
#> Pyy                             -0.7742791
#> Sst                             -0.3111668
#> Rbp4                            -0.7511252
#>      CCCATACTCGTGGACCstochastic_pseudotime
#> Pyy                             -0.7742465
#> Sst                             -0.3111658
#> Rbp4                            -0.7511055
#>      TTCCCAGTCATAGCACstochastic_pseudotime
#> Pyy                             -0.7734517
#> Sst                             -0.3111407
#> Rbp4                            -0.7506243
#>      TGGCCAGCAATGAAACstochastic_pseudotime
#> Pyy                             -0.7731112
#> Sst                             -0.3111300
#> Rbp4                            -0.7504178
#>      CTCGGAGGTGAAGGCTstochastic_pseudotime
#> Pyy                             -0.7703543
#> Sst                             -0.3110424
#> Rbp4                            -0.7487395
#>      GACGTTAGTGTGAATAstochastic_pseudotime
#> Pyy                             -0.7665027
#> Sst                             -0.3109186
#> Rbp4                            -0.7463759
#>      CGTAGCGTCGGAAATAstochastic_pseudotime
#> Pyy                             -0.7655964
#> Sst                             -0.3108892
#> Rbp4                            -0.7458167
#>      GGAGCAAGTATCACCAstochastic_pseudotime
#> Pyy                             -0.7641375
#> Sst                             -0.3108418
#> Rbp4                            -0.7449140
#>      CCTCAGTAGTGCGATGstochastic_pseudotime
#> Pyy                             -0.7621492
#> Sst                             -0.3107768
#> Rbp4                            -0.7436790
#>      TAGTGGTTCCTCAACCstochastic_pseudotime
#> Pyy                             -0.7615342
#> Sst                             -0.3107566
#> Rbp4                            -0.7432959
#>      TGGGAAGCACGTAAGGstochastic_pseudotime
#> Pyy                             -0.7591126
#> Sst                             -0.3106767
#> Rbp4                            -0.7417826
#>      AGAGTGGCACCTCGGAstochastic_pseudotime
#> Pyy                             -0.7558019
#> Sst                             -0.3105664
#> Rbp4                            -0.7397014
#>      TAGTGGTCATTCACTTstochastic_pseudotime
#> Pyy                             -0.7536330
#> Sst                             -0.3104936
#> Rbp4                            -0.7383304
#>      CAAGGCCTCTCAACTTstochastic_pseudotime
#> Pyy                             -0.7532993
#> Sst                             -0.3104824
#> Rbp4                            -0.7381190
#>      TTTCCTCGTGGCGAATstochastic_pseudotime
#> Pyy                             -0.7501893
#> Sst                             -0.3103771
#> Rbp4                            -0.7361420
#>      CTTTGCGAGACCACGAstochastic_pseudotime
#> Pyy                             -0.7498596
#> Sst                             -0.3103659
#> Rbp4                            -0.7359318
#>      ACTGATGCATGTCTCCstochastic_pseudotime
#> Pyy                             -0.7478056
#> Sst                             -0.3102957
#> Rbp4                            -0.7346190
#>      CCACTACAGAGCTGGTstochastic_pseudotime
#> Pyy                             -0.7474975
#> Sst                             -0.3102852
#> Rbp4                            -0.7344217
#>      CATATGGTCACAGTACstochastic_pseudotime
#> Pyy                             -0.7427915
#> Sst                             -0.3101229
#> Rbp4                            -0.7313946
#>      CTGGTCTAGGATTCGGstochastic_pseudotime
#> Pyy                             -0.7406643
#> Sst                             -0.3100489
#> Rbp4                            -0.7300184
#>      ATCCGAAGTACACCGCstochastic_pseudotime
#> Pyy                             -0.7394826
#> Sst                             -0.3100075
#> Rbp4                            -0.7292518
#>      GTTACAGTCCGATATGstochastic_pseudotime
#> Pyy                             -0.7391292
#> Sst                             -0.3099952
#> Rbp4                            -0.7290223
#>      GAGGTGAGTCCGCTGAstochastic_pseudotime
#> Pyy                             -0.7365872
#> Sst                             -0.3099058
#> Rbp4                            -0.7273675
#>      AAGGTTCGTCAACTGTstochastic_pseudotime
#> Pyy                             -0.7343294
#> Sst                             -0.3098260
#> Rbp4                            -0.7258922
#>      CAGCTAACATGCCTAAstochastic_pseudotime
#> Pyy                             -0.7294976
#> Sst                             -0.3096537
#> Rbp4                            -0.7227185
#>      TGCCCTAAGTCCGGTCstochastic_pseudotime
#> Pyy                             -0.7250162
#> Sst                             -0.3094921
#> Rbp4                            -0.7197555
#>      GGACGTCGTCAGAAGCstochastic_pseudotime
#> Pyy                             -0.7249543
#> Sst                             -0.3094899
#> Rbp4                            -0.7197144
#>      TTCGGTCTCGACCAGCstochastic_pseudotime
#> Pyy                             -0.7211710
#> Sst                             -0.3093523
#> Rbp4                            -0.7171989
#>      ACACCCTTCCCAGGTGstochastic_pseudotime
#> Pyy                             -0.7194566
#> Sst                             -0.3092895
#> Rbp4                            -0.7160549
#>      CACCACTGTGTGTGCCstochastic_pseudotime
#> Pyy                             -0.7169558
#> Sst                             -0.3091976
#> Rbp4                            -0.7143818
#>      GCTGCGATCCCGGATGstochastic_pseudotime
#> Pyy                             -0.7079458
#> Sst                             -0.3088626
#> Rbp4                            -0.7083123
#>      ATGAGGGTCTTGACGAstochastic_pseudotime
#> Pyy                             -0.7018786
#> Sst                             -0.3086337
#> Rbp4                            -0.7041906
#>      AACCATGGTTGAACTCstochastic_pseudotime
#> Pyy                             -0.7003534
#> Sst                             -0.3085758
#> Rbp4                            -0.7031504
#>      AGTGGGAGTGGCAAACstochastic_pseudotime
#> Pyy                             -0.7002950
#> Sst                             -0.3085736
#> Rbp4                            -0.7031106
#>      CCAGCGACAGACGTAGstochastic_pseudotime
#> Pyy                             -0.6879154
#> Sst                             -0.3080981
#> Rbp4                            -0.6946088
#>      CCACGGAAGTCTCCTCstochastic_pseudotime
#> Pyy                             -0.6821273
#> Sst                             -0.3078726
#> Rbp4                            -0.6906006
#>      ACTGCTCAGGTGTTAAstochastic_pseudotime
#> Pyy                             -0.6768249
#> Sst                             -0.3076644
#> Rbp4                            -0.6869114
#>      CACACCTGTTATCCGAstochastic_pseudotime
#> Pyy                             -0.6767249
#> Sst                             -0.3076605
#> Rbp4                            -0.6868417
#>      GGCGTGTCACCGATATstochastic_pseudotime
#> Pyy                             -0.6713224
#> Sst                             -0.3074468
#> Rbp4                            -0.6830663
#>      TTCTACACATAAAGGTstochastic_pseudotime
#> Pyy                             -0.6646918
#> Sst                             -0.3071825
#> Rbp4                            -0.6784114
#>      AGTTGGTGTTGTCGCGstochastic_pseudotime
#> Pyy                             -0.6566121
#> Sst                             -0.3068577
#> Rbp4                            -0.6727093
#>      TAAGAGATCTGCGACGstochastic_pseudotime
#> Pyy                             -0.6531965
#> Sst                             -0.3067196
#> Rbp4                            -0.6702894
#>      GTGTTAGGTGAACCTTstochastic_pseudotime
#> Pyy                             -0.6419354
#> Sst                             -0.3062607
#> Rbp4                            -0.6622737
#>      CTGCGGACATTAGCCAstochastic_pseudotime
#> Pyy                             -0.6330499
#> Sst                             -0.3058954
#> Rbp4                            -0.6559110
#>      GTCACGGTCTTGCCGTstochastic_pseudotime
#> Pyy                             -0.6299434
#> Sst                             -0.3057670
#> Rbp4                            -0.6536791
#>      TGGACGCTCTTGTTTGstochastic_pseudotime
#> Pyy                             -0.6231340
#> Sst                             -0.3054846
#> Rbp4                            -0.6487738
#>      CTCCTAGCACAACGTTstochastic_pseudotime
#> Pyy                             -0.6082346
#> Sst                             -0.3048621
#> Rbp4                            -0.6379829
#>      TCAATCTCATCCGCGAstochastic_pseudotime
#> Pyy                             -0.6040204
#> Sst                             -0.3046850
#> Rbp4                            -0.6349172
#>      TCTCTAATCCCAACGGstochastic_pseudotime
#> Pyy                             -0.5529457
#> Sst                             -0.3025151
#> Rbp4                            -0.5973541
#>      GCTCTGTGTCCATGATstochastic_pseudotime
#> Pyy                             -0.5458128
#> Sst                             -0.3022098
#> Rbp4                            -0.5920558
#>      TACAGTGGTAAACGCGstochastic_pseudotime
#> Pyy                             -0.5454212
#> Sst                             -0.3021930
#> Rbp4                            -0.5917645
#>      TCGAGGCAGGCATGGTstochastic_pseudotime
#> Pyy                             -0.5416792
#> Sst                             -0.3020328
#> Rbp4                            -0.5889801
#>      GGCGTGTTCGAGAGCAstochastic_pseudotime
#> Pyy                             -0.5359879
#> Sst                             -0.3017889
#> Rbp4                            -0.5847395
#>      GGGTCTGTCATGTGGTstochastic_pseudotime
#> Pyy                             -0.5346856
#> Sst                             -0.3017331
#> Rbp4                            -0.5837682
#>      GATGCTAGTGTTAAGAstochastic_pseudotime
#> Pyy                             -0.5303109
#> Sst                             -0.3015456
#> Rbp4                            -0.5805030
#>      CTACGTCAGCTAGGCAstochastic_pseudotime
#> Pyy                             -0.5240634
#> Sst                             -0.3012777
#> Rbp4                            -0.5758336
#>      TCAGATGCACTCGACGstochastic_pseudotime
#> Pyy                             -0.4911884
#> Sst                             -0.2998686
#> Rbp4                            -0.5511557
#>      GTGCATATCAAGGCTTstochastic_pseudotime
#> Pyy                             -0.4863432
#> Sst                             -0.2996612
#> Rbp4                            -0.5475054
#>      CTCTGGTCAAGTTGTCstochastic_pseudotime
#> Pyy                             -0.4648170
#> Sst                             -0.2987412
#> Rbp4                            -0.5312535
#>      GAGCAGATCAACTCTTstochastic_pseudotime
#> Pyy                             -0.4547421
#> Sst                             -0.2983117
#> Rbp4                            -0.5236299
#>      CACAGGCCAAGGCTCCstochastic_pseudotime
#> Pyy                             -0.4459091
#> Sst                             -0.2979360
#> Rbp4                            -0.5169380
#>      TGCGTGGGTCCCTTGTstochastic_pseudotime
#> Pyy                             -0.4452906
#> Sst                             -0.2979098
#> Rbp4                            -0.5164692
#>      GGCTCGAAGATCACGGstochastic_pseudotime
#> Pyy                             -0.4386756
#> Sst                             -0.2976290
#> Rbp4                            -0.5114527
#>      TGTGTTTGTATGAATGstochastic_pseudotime
#> Pyy                             -0.4314510
#> Sst                             -0.2973229
#> Rbp4                            -0.5059700
#>      GATCGTAAGTAGGTGCstochastic_pseudotime
#> Pyy                             -0.4287373
#> Sst                             -0.2972081
#> Rbp4                            -0.5039094
#>      GATCGATAGAAAGTGGstochastic_pseudotime
#> Pyy                             -0.4056602
#> Sst                             -0.2962358
#> Rbp4                            -0.4863653
#>      GTGCATAAGAACAATCstochastic_pseudotime
#> Pyy                             -0.3612911
#> Sst                             -0.2943910
#> Rbp4                            -0.4525452
#>      TGGTTAGCAGTCGTGCstochastic_pseudotime
#> Pyy                             -0.3108397
#> Sst                             -0.2923428
#> Rbp4                            -0.4139852
#>      AAGGTTCAGATAGTCAstochastic_pseudotime
#> Pyy                             -0.3019337
#> Sst                             -0.2919877
#> Rbp4                            -0.4071698
#>      ACTGTCCCAGTTTACGstochastic_pseudotime
#> Pyy                             -0.2935131
#> Sst                             -0.2916538
#> Rbp4                            -0.4007238
#>      CACCAGGCATAAAGGTstochastic_pseudotime
#> Pyy                             -0.2700343
#> Sst                             -0.2907330
#> Rbp4                            -0.3827418
#>      CTGAAGTTCTGCGGCAstochastic_pseudotime
#> Pyy                             -0.2141790
#> Sst                             -0.2886085
#> Rbp4                            -0.3399205
#>      GCATGATGTCACTGGCstochastic_pseudotime
#> Pyy                             -0.1934821
#> Sst                             -0.2878467
#> Rbp4                            -0.3240416
#>      TCAACGATCAACCATGstochastic_pseudotime
#> Pyy                             -0.1804990
#> Sst                             -0.2873763
#> Rbp4                            -0.3140783
#>      ACAGCCGAGAGAACAGstochastic_pseudotime
#> Pyy                             -0.1787710
#> Sst                             -0.2873141
#> Rbp4                            -0.3127521
#>      GACGTGCGTCACCTAAstochastic_pseudotime
#> Pyy                             -0.1719360
#> Sst                             -0.2870693
#> Rbp4                            -0.3075061
#>      TCGGGACGTCTGCGGTstochastic_pseudotime
#> Pyy                             -0.1464502
#> Sst                             -0.2861710
#> Rbp4                            -0.2879413
#>      TGAGCCGTCCGCTGTTstochastic_pseudotime
#> Pyy                             -0.1402739
#> Sst                             -0.2859569
#> Rbp4                            -0.2831992
#>      TCGCGAGAGTACGCGAstochastic_pseudotime
#> Pyy                            -0.09152524
#> Sst                            -0.28431780
#> Rbp4                           -0.24576094
#>      GATTCAGAGCGATGACstochastic_pseudotime
#> Pyy                            -0.07976931
#> Sst                            -0.28393639
#> Rbp4                           -0.23673055
#>      TCTATTGTCGTATCAGstochastic_pseudotime
#> Pyy                            -0.03163259
#> Sst                            -0.28243309
#> Rbp4                           -0.19974763
#>      AGCTTGAAGACACGACstochastic_pseudotime
#> Pyy                            -0.02978192
#> Sst                            -0.28237721
#> Rbp4                           -0.19832560
#>      CGAATGTGTCTGGAGAstochastic_pseudotime
#> Pyy                            -0.02532447
#> Sst                            -0.28224319
#> Rbp4                           -0.19490050
#>      CGTTCTGAGGACGAAAstochastic_pseudotime
#> Pyy                           -0.007676362
#> Sst                           -0.281720839
#> Rbp4                          -0.181339072
#>      ACTATCTTCAAGATCCstochastic_pseudotime
#> Pyy                             0.08658458
#> Sst                            -0.27915962
#> Rbp4                           -0.10888981
#>      ACACTGAGTCCATGATstochastic_pseudotime
#> Pyy                             0.08863592
#> Sst                            -0.27910826
#> Rbp4                           -0.10731287
#>      GCAATCACAGCTGCTGstochastic_pseudotime
#> Pyy                             0.09406518
#> Sst                            -0.27897326
#> Rbp4                           -0.10313912
#>      AAGCCGCCAATGGACGstochastic_pseudotime
#> Pyy                             0.09984490
#> Sst                            -0.27883100
#> Rbp4                           -0.09869587
#>      GTAGGCCAGCTAACAAstochastic_pseudotime
#> Pyy                             0.11350800
#> Sst                            -0.27850074
#> Rbp4                           -0.08819178
#>      GATCAGTAGGTGATATstochastic_pseudotime
#> Pyy                              0.1205933
#> Sst                             -0.2783328
#> Rbp4                            -0.0827444
#>      TGCGTGGAGGTTCCTAstochastic_pseudotime
#> Pyy                              0.1222350
#> Sst                             -0.2782942
#> Rbp4                            -0.0814822
#>      GGACAGACAATCCAACstochastic_pseudotime
#> Pyy                             0.13295509
#> Sst                            -0.27804536
#> Rbp4                           -0.07324004
#>      CATGACATCATGCATGstochastic_pseudotime
#> Pyy                             0.13973162
#> Sst                            -0.27789076
#> Rbp4                           -0.06802968
#>      CATGCCTCACCCATTCstochastic_pseudotime
#> Pyy                             0.20027046
#> Sst                            -0.27660392
#> Rbp4                           -0.02147566
#>      TGCGGGTTCCACTCCAstochastic_pseudotime
#> Pyy                             0.20505231
#> Sst                            -0.27650956
#> Rbp4                           -0.01779788
#>      TGTCCCACATTATCTCstochastic_pseudotime
#> Pyy                            0.222082476
#> Sst                           -0.276182224
#> Rbp4                          -0.004699027
#>      CAAGTTGGTGTGCCTGstochastic_pseudotime
#> Pyy                             0.24658899
#> Sst                            -0.27573516
#> Rbp4                            0.01415246
#>      CTCGAAAGTGGTGTAGstochastic_pseudotime
#> Pyy                             0.27382222
#> Sst                            -0.27527170
#> Rbp4                            0.03510472
#>      TCGCGAGCACTTCTGCstochastic_pseudotime
#> Pyy                             0.28294172
#> Sst                            -0.27512438
#> Rbp4                            0.04212176
#>      CGAGCACGTAGTAGTAstochastic_pseudotime
#> Pyy                             0.28943403
#> Sst                            -0.27502190
#> Rbp4                            0.04711757
#>      TTGGAACAGAAGATTCstochastic_pseudotime
#> Pyy                             0.29840671
#> Sst                            -0.27488359
#> Rbp4                            0.05402239
#>      ACGATGTTCGCAGGCTstochastic_pseudotime
#> Pyy                             0.31554148
#> Sst                            -0.27463011
#> Rbp4                            0.06720952
#>      GGACGTCAGGGTGTTGstochastic_pseudotime
#> Pyy                             0.32305784
#> Sst                            -0.27452335
#> Rbp4                            0.07299476
#>      TGCGGGTAGACGCACAstochastic_pseudotime
#> Pyy                              0.3671001
#> Sst                             -0.2739521
#> Rbp4                             0.1069008
#>      CCATGTCTCCGCATAAstochastic_pseudotime
#> Pyy                              0.3677049
#> Sst                             -0.2739449
#> Rbp4                             0.1073665
#>      CAAGAAAAGAGTACCGstochastic_pseudotime
#> Pyy                              0.3794884
#> Sst                             -0.2738081
#> Rbp4                             0.1164404
#>      CCGTGGATCTTAGAGCstochastic_pseudotime
#> Pyy                              0.3927588
#> Sst                             -0.2736621
#> Rbp4                             0.1266605
#>      GGCGTGTCACCAGGTCstochastic_pseudotime
#> Pyy                              0.4025492
#> Sst                             -0.2735598
#> Rbp4                             0.1342014
#>      GAACCTAAGTAGATGTstochastic_pseudotime
#> Pyy                              0.4636151
#> Sst                             -0.2730255
#> Rbp4                             0.1812552
#>      CGAGAAGTCGTTACAGstochastic_pseudotime
#> Pyy                              0.4826216
#> Sst                             -0.2728957
#> Rbp4                             0.1959078
#>      AAGGAGCAGTACATGAstochastic_pseudotime
#> Pyy                              0.4867808
#> Sst                             -0.2728696
#> Rbp4                             0.1991147
#>      CGTCTACAGGGAGTAAstochastic_pseudotime
#> Pyy                              0.5124411
#> Sst                             -0.2727271
#> Rbp4                             0.2189041
#>      TTGCCGTCACCGAATTstochastic_pseudotime
#> Pyy                              0.5129877
#> Sst                             -0.2727244
#> Rbp4                             0.2193257
#>      CGATTGAAGTCAAGGCstochastic_pseudotime
#> Pyy                              0.5144504
#> Sst                             -0.2727173
#> Rbp4                             0.2204540
#>      AGCTCTCAGACATAACstochastic_pseudotime
#> Pyy                              0.5443702
#> Sst                             -0.2725938
#> Rbp4                             0.2435386
#>      AAATGCCCATGCTAGTstochastic_pseudotime
#> Pyy                              0.5524555
#> Sst                             -0.2725678
#> Rbp4                             0.2497787
#>      GACTGCGAGAATCTCCstochastic_pseudotime
#> Pyy                              0.5735164
#> Sst                             -0.2725147
#> Rbp4                             0.2660371
#>      ACGCCGAAGGGTGTTGstochastic_pseudotime
#> Pyy                              0.5977234
#> Sst                             -0.2724800
#> Rbp4                             0.2847318
#>      TTGTAGGTCAAAGTAGstochastic_pseudotime
#> Pyy                              0.6208190
#> Sst                             -0.2724729
#> Rbp4                             0.3025761
#>      CGATGTAAGTTTAGGAstochastic_pseudotime
#> Pyy                              0.6434242
#> Sst                             -0.2724906
#> Rbp4                             0.3200495
#>      GGAAAGCCAATGTAAGstochastic_pseudotime
#> Pyy                              0.6485352
#> Sst                             -0.2724980
#> Rbp4                             0.3240014
#>      AAACCTGGTAAGTGGCstochastic_pseudotime
#> Pyy                              0.6582563
#> Sst                             -0.2725155
#> Rbp4                             0.3315190
#>      AGGGAGTTCGGAGGTAstochastic_pseudotime
#> Pyy                              0.6667454
#> Sst                             -0.2725345
#> Rbp4                             0.3380853
#>      CTCGAGGGTAAGTAGTstochastic_pseudotime
#> Pyy                              0.6939703
#> Sst                             -0.2726183
#> Rbp4                             0.3591520
#>      CTTCTCTAGTCCGTATstochastic_pseudotime
#> Pyy                              0.7234825
#> Sst                             -0.2727491
#> Rbp4                             0.3820046
#>      GAATGAAAGAATGTGTstochastic_pseudotime
#> Pyy                              0.7285666
#> Sst                             -0.2727758
#> Rbp4                             0.3859431
#>      ATCATGGAGGAATGGAstochastic_pseudotime
#> Pyy                              0.7419177
#> Sst                             -0.2728517
#> Rbp4                             0.3962885
#>      GATTCAGTCGCCAGCAstochastic_pseudotime
#> Pyy                              0.7449336
#> Sst                             -0.2728701
#> Rbp4                             0.3986260
#>      CAGGTGCAGCCAGGATstochastic_pseudotime
#> Pyy                              0.7801171
#> Sst                             -0.2731157
#> Rbp4                             0.4259096
#>      ACGCAGCCATTCTCATstochastic_pseudotime
#> Pyy                              0.7947217
#> Sst                             -0.2732350
#> Rbp4                             0.4372432
#>      GCGCAGTGTGCTAGCCstochastic_pseudotime
#> Pyy                              0.8056712
#> Sst                             -0.2733310
#> Rbp4                             0.4457439
#>      GACCAATCAGACGCCTstochastic_pseudotime
#> Pyy                              0.8165522
#> Sst                             -0.2734321
#> Rbp4                             0.4541942
#>      CGGGTCAAGGCGATACstochastic_pseudotime
#> Pyy                              0.8425979
#> Sst                             -0.2736968
#> Rbp4                             0.4744341
#>      CTGCTGTCAAAGCAATstochastic_pseudotime
#> Pyy                              0.8579253
#> Sst                             -0.2738676
#> Rbp4                             0.4863536
#>      CGTGTCTCAGAAGCACstochastic_pseudotime
#> Pyy                              0.8777448
#> Sst                             -0.2741051
#> Rbp4                             0.5017762
#>      CCTACACAGTACGTTCstochastic_pseudotime
#> Pyy                              0.9013291
#> Sst                             -0.2744120
#> Rbp4                             0.5201437
#>      AACTCCCCAGTTCATGstochastic_pseudotime
#> Pyy                              0.9079702
#> Sst                             -0.2745032
#> Rbp4                             0.5253190
#>      AGCTCCTTCAATACCGstochastic_pseudotime
#> Pyy                              0.9134627
#> Sst                             -0.2745802
#> Rbp4                             0.5296003
#>      TGCGTGGGTTACTGACstochastic_pseudotime
#> Pyy                              0.9190847
#> Sst                             -0.2746606
#> Rbp4                             0.5339834
#>      AAGGTTCAGGAGCGAGstochastic_pseudotime
#> Pyy                              0.9378171
#> Sst                             -0.2749392
#> Rbp4                             0.5485959
#>      TGTCCCATCTCTGCTGstochastic_pseudotime
#> Pyy                              0.9737042
#> Sst                             -0.2755201
#> Rbp4                             0.5766247
#>      TGCTACCCATAAGACAstochastic_pseudotime
#> Pyy                              0.9855099
#> Sst                             -0.2757248
#> Rbp4                             0.5858559
#>      ACCAGTAGTTCGTTGAstochastic_pseudotime
#> Pyy                              1.0305810
#> Sst                             -0.2765687
#> Rbp4                             0.6211504
#>      CGATCGGAGCTGCAAGstochastic_pseudotime
#> Pyy                              1.0594777
#> Sst                             -0.2771625
#> Rbp4                             0.6438262
#>      CCACTACCAACGCACCstochastic_pseudotime
#> Pyy                              1.0742114
#> Sst                             -0.2774813
#> Rbp4                             0.6554034
#>      CATCAGAGTCTCTTTAstochastic_pseudotime
#> Pyy                              1.0905001
#> Sst                             -0.2778466
#> Rbp4                             0.6682155
#>      AGTCTTTCACGGACAAstochastic_pseudotime
#> Pyy                              1.1675728
#> Sst                             -0.2797615
#> Rbp4                             0.7290436
#>      TTTGTCATCTGTTTGTstochastic_pseudotime
#> Pyy                              1.1746254
#> Sst                             -0.2799526
#> Rbp4                             0.7346287
#>      AACTTTCTCCAGATCAstochastic_pseudotime
#> Pyy                              1.1808198
#> Sst                             -0.2801227
#> Rbp4                             0.7395371
#>      AAAGCAAGTGCAGTAGstochastic_pseudotime
#> Pyy                              1.1830118
#> Sst                             -0.2801834
#> Rbp4                             0.7412746
#>      CGTGTCTGTAGCGTCCstochastic_pseudotime
#> Pyy                              1.2018319
#> Sst                             -0.2807157
#> Rbp4                             0.7562076
#>      CAGTCCTCAGGGCATAstochastic_pseudotime
#> Pyy                              1.2062738
#> Sst                             -0.2808443
#> Rbp4                             0.7597361
#>      AGGGATGGTCGACTGCstochastic_pseudotime
#> Pyy                              1.2246679
#> Sst                             -0.2813889
#> Rbp4                             0.7743645
#>      CCACGGAAGATCCTGTstochastic_pseudotime
#> Pyy                              1.2426077
#> Sst                             -0.2819395
#> Rbp4                             0.7886593
#>      TTGAACGGTGAGGGAGstochastic_pseudotime
#> Pyy                              1.2586929
#> Sst                             -0.2824499
#> Rbp4                             0.8015014
#>      GTGCTTCAGAGGGATAstochastic_pseudotime
#> Pyy                              1.2632853
#> Sst                             -0.2825986
#> Rbp4                             0.8051724
#>      TCACAAGCAACGATGGstochastic_pseudotime
#> Pyy                              1.2676347
#> Sst                             -0.2827407
#> Rbp4                             0.8086512
#>      TTCGGTCGTATAGTAGstochastic_pseudotime
#> Pyy                              1.3200813
#> Sst                             -0.2845544
#> Rbp4                             0.8507634
#>      CCTTCGATCACAGGCCstochastic_pseudotime
#> Pyy                              1.3255744
#> Sst                             -0.2847558
#> Rbp4                             0.8551936
#>      AGGCCACCATGACATCstochastic_pseudotime
#> Pyy                              1.3330606
#> Sst                             -0.2850338
#> Rbp4                             0.8612377
#>      TTTGTCAAGAATGTGTstochastic_pseudotime
#> Pyy                              1.3394082
#> Sst                             -0.2852730
#> Rbp4                             0.8663687
#>      CCTCAGTCAGACACTTstochastic_pseudotime
#> Pyy                              1.3504908
#> Sst                             -0.2856984
#> Rbp4                             0.8753416
#>      CGCTTCATCTCTGTCGstochastic_pseudotime
#> Pyy                              1.3524382
#> Sst                             -0.2857742
#> Rbp4                             0.8769202
#>      GTGCATAGTTCAACCAstochastic_pseudotime
#> Pyy                              1.3573696
#> Sst                             -0.2859676
#> Rbp4                             0.8809204
#>      ATAGACCAGCGTGTCCstochastic_pseudotime
#> Pyy                              1.3676799
#> Sst                             -0.2863785
#> Rbp4                             0.8892969
#>      GGGAATGAGCACACAGstochastic_pseudotime
#> Pyy                              1.3714818
#> Sst                             -0.2865323
#> Rbp4                             0.8923904
#>      GTACTTTTCTTCGGTCstochastic_pseudotime
#> Pyy                              1.3788608
#> Sst                             -0.2868344
#> Rbp4                             0.8984019
#>      AGCAGCCGTCGCATCGstochastic_pseudotime
#> Pyy                              1.3848827
#> Sst                             -0.2870845
#> Rbp4                             0.9033154
#>      GTAGTCAGTCCGACGTstochastic_pseudotime
#> Pyy                              1.3865368
#> Sst                             -0.2871537
#> Rbp4                             0.9046663
#>      CTTCTCTGTCAGAGGTstochastic_pseudotime
#> Pyy                              1.4049883
#> Sst                             -0.2879425
#> Rbp4                             0.9197745
#>      GTGGGTCGTCCAGTATstochastic_pseudotime
#> Pyy                              1.4136687
#> Sst                             -0.2883242
#> Rbp4                             0.9269085
#>      GATGAAACAGGCTGAAstochastic_pseudotime
#> Pyy                              1.4349312
#> Sst                             -0.2892889
#> Rbp4                             0.9444645
#>      ATCTGCCCATGTTGACstochastic_pseudotime
#> Pyy                              1.4412099
#> Sst                             -0.2895821
#> Rbp4                             0.9496731
#>      TTAACTCCAGGATTGGstochastic_pseudotime
#> Pyy                              1.4851246
#> Sst                             -0.2917469
#> Rbp4                             0.9864859
#>      AGAATAGCAATCACACstochastic_pseudotime
#> Pyy                              1.4910887
#> Sst                             -0.2920576
#> Rbp4                             0.9915472
#>      GGGACCTCACTGTCGGstochastic_pseudotime
#> Pyy                              1.5139446
#> Sst                             -0.2932896
#> Rbp4                             1.0111179
#>      AGGGATGTCGCATGGCstochastic_pseudotime
#> Pyy                              1.5143689
#> Sst                             -0.2933132
#> Rbp4                             1.0114841
#>      GGACGTCGTCGAGATGstochastic_pseudotime
#> Pyy                              1.5444401
#> Sst                             -0.2950519
#> Rbp4                             1.0377791
#>      TCAATCTTCGAACGGAstochastic_pseudotime
#> Pyy                              1.5464850
#> Sst                             -0.2951757
#> Rbp4                             1.0395953
#>      CTCGTCAAGAACAATCstochastic_pseudotime
#> Pyy                              1.5526191
#> Sst                             -0.2955521
#> Rbp4                             1.0450694
#>      CACAGTAGTGCACTTAstochastic_pseudotime
#> Pyy                              1.5684472
#> Sst                             -0.2965604
#> Rbp4                             1.0593988
#>      AGTGTCATCCCTTGTGstochastic_pseudotime
#> Pyy                              1.5735846
#> Sst                             -0.2969008
#> Rbp4                             1.0641239
#>      GACCTGGGTTCGTGATstochastic_pseudotime
#> Pyy                              1.5753479
#> Sst                             -0.2970193
#> Rbp4                             1.0657555
#>      GCTGCGAAGAACAACTstochastic_pseudotime
#> Pyy                              1.5773422
#> Sst                             -0.2971544
#> Rbp4                             1.0676071
#>      CAGCCGATCCTATTCAstochastic_pseudotime
#> Pyy                              1.5831913
#> Sst                             -0.2975575
#> Rbp4                             1.0730790
#>      CACCACTAGCGATCCCstochastic_pseudotime
#> Pyy                              1.5833810
#> Sst                             -0.2975708
#> Rbp4                             1.0732575
#>      GCTTCCAAGTTAAGTGstochastic_pseudotime
#> Pyy                              1.5838701
#> Sst                             -0.2976051
#> Rbp4                             1.0737184
#>      ACACCGGTCATGTCTTstochastic_pseudotime
#> Pyy                               1.590741
#> Sst                              -0.298095
#> Rbp4                              1.080245
#>      TGACGGCCATATGAGAstochastic_pseudotime
#> Pyy                              1.5980957
#> Sst                             -0.2986401
#> Rbp4                             1.0873576
#>      GCATGTAGTTGCGCACstochastic_pseudotime
#> Pyy                              1.6047954
#> Sst                             -0.2991585
#> Rbp4                             1.0939752
#>      CATATTCTCCTGTAGAstochastic_pseudotime
#> Pyy                              1.6213118
#> Sst                             -0.3005603
#> Rbp4                             1.1110854
#>      CGTTCTGTCACAATGCstochastic_pseudotime
#> Pyy                              1.6252523
#> Sst                             -0.3009313
#> Rbp4                             1.1154080
#>      AATCGGTTCGTATCAGstochastic_pseudotime
#> Pyy                              1.6258992
#> Sst                             -0.3009941
#> Rbp4                             1.1161298
#>      GAACGGAAGTCAAGCGstochastic_pseudotime
#> Pyy                              1.6294068
#> Sst                             -0.3013448
#> Rbp4                             1.1201136
#>      GGCGTGTTCAGCTTAGstochastic_pseudotime
#> Pyy                              1.6344088
#> Sst                             -0.3018825
#> Rbp4                             1.1260465
#>      CTAACTTAGAGGTACCstochastic_pseudotime
#> Pyy                              1.6393934
#> Sst                             -0.3024812
#> Rbp4                             1.1323825
#>      CTCTAATAGTTATCGCstochastic_pseudotime
#> Pyy                              1.6414929
#> Sst                             -0.3027611
#> Rbp4                             1.1352391
#>      TACACGAGTGCAACTTstochastic_pseudotime
#> Pyy                              1.6416333
#> Sst                             -0.3027806
#> Rbp4                             1.1354353
#>      CTCGTACTCTTTAGTCstochastic_pseudotime
#> Pyy                               1.642233
#> Sst                              -0.302865
#> Rbp4                              1.136282
#>      CAAGGCCTCAGTGCATstochastic_pseudotime
#> Pyy                              1.6436773
#> Sst                             -0.3030779
#> Rbp4                             1.1383849
#>      AGGTCATCACTTCTGCstochastic_pseudotime
#> Pyy                              1.6450222
#> Sst                             -0.3032904
#> Rbp4                             1.1404410
#>      ATCACGAGTGCACCACstochastic_pseudotime
#> Pyy                               1.645025
#> Sst                              -0.303291
#> Rbp4                              1.140446
#>      ACAGCTAAGTCGTTTGstochastic_pseudotime
#> Pyy                              1.6492563
#> Sst                             -0.3041231
#> Rbp4                             1.1480293
#>      GTGGGTCTCCGAATGTstochastic_pseudotime
#> Pyy                              1.6500586
#> Sst                             -0.3043394
#> Rbp4                             1.1498666
#>      CCGGGATAGATCTGAAstochastic_pseudotime
#> Pyy                              1.6505336
#> Sst                             -0.3044886
#> Rbp4                             1.1510991
#>      GGATGTTAGTAAGTACstochastic_pseudotime
#> Pyy                              1.6518332
#> Sst                             -0.3052077
#> Rbp4                             1.1565985
#>      TTCTTAGTCAGTGTTGstochastic_pseudotime
#> Pyy                              1.6518778
#> Sst                             -0.3053623
#> Rbp4                             1.1576761
#>      GCGACCATCACCAGGCstochastic_pseudotime
#> Pyy                              1.6518771
#> Sst                             -0.3053813
#> Rbp4                             1.1578063
#>      ACCAGTAAGTGATCGGstochastic_pseudotime
#> Pyy                              1.6518230
#> Sst                             -0.3055296
#> Rbp4                             1.1587970
#>      TAAGTGCTCCAAATGCstochastic_pseudotime
#> Pyy                               1.650575
#> Sst                              -0.306155
#> Rbp4                              1.162525
#>      CGCGGTAAGAGCTGGTstochastic_pseudotime
#> Pyy                              1.6490253
#> Sst                             -0.3065129
#> Rbp4                             1.1642884
#>      CCGTTCAAGCTAAACAstochastic_pseudotime
#> Pyy                              1.6482662
#> Sst                             -0.3066473
#> Rbp4                             1.1648727
#>      ATAGACCAGTGATCGGstochastic_pseudotime
#> Pyy                              1.6481682
#> Sst                             -0.3066635
#> Rbp4                             1.1649398
#>      GACTAACGTCATCCCTstochastic_pseudotime
#> Pyy                              1.6481473
#> Sst                             -0.3066669
#> Rbp4                             1.1649538
#>      GGCAATTCACGGCTACstochastic_pseudotime
#> Pyy                               1.647402
#> Sst                              -0.306782
#> Rbp4                              1.165412
#>      CCTACACTCTTTACACstochastic_pseudotime
#> Pyy                              1.6468545
#> Sst                             -0.3068595
#> Rbp4                             1.1657010
#>      GCATGTATCGTGGACCstochastic_pseudotime
#> Pyy                              1.6464675
#> Sst                             -0.3069114
#> Rbp4                             1.1658853
#>      CGCCAAGCAAAGTCAAstochastic_pseudotime
#> Pyy                              1.6449590
#> Sst                             -0.3070951
#> Rbp4                             1.1664763
#>      CGAATGTTCACCGGGTstochastic_pseudotime
#> Pyy                              1.6446194
#> Sst                             -0.3071331
#> Rbp4                             1.1665861
#>      ACGATACTCGCGCCAAstochastic_pseudotime
#> Pyy                              1.6438208
#> Sst                             -0.3072183
#> Rbp4                             1.1668168
#>      GGACAAGTCAAAGACAstochastic_pseudotime
#> Pyy                              1.6409357
#> Sst                             -0.3074889
#> Rbp4                             1.1673959
#>      AAGACCTAGTAGGCCAstochastic_pseudotime
#> Pyy                              1.6378446
#> Sst                             -0.3077326
#> Rbp4                             1.1677016
#>      CCACCTATCTAAGCCAstochastic_pseudotime
#> Pyy                              1.6342131
#> Sst                             -0.3079778
#> Rbp4                             1.1677804
#>      AGCATACCACAGCCCAstochastic_pseudotime
#> Pyy                              1.6309411
#> Sst                             -0.3081713
#> Rbp4                             1.1676644
#>      CCTACCACATCGATTGstochastic_pseudotime
#> Pyy                              1.6212591
#> Sst                             -0.3086418
#> Rbp4                             1.1666309
#>      ACGATACAGATGAGAGstochastic_pseudotime
#> Pyy                              1.6190006
#> Sst                             -0.3087356
#> Rbp4                             1.1662817
#>      AGCAGCCGTGTAACGGstochastic_pseudotime
#> Pyy                              1.6142440
#> Sst                             -0.3089183
#> Rbp4                             1.1654471
#>      CCGGGATTCGACAGCCstochastic_pseudotime
#> Pyy                               1.613509
#> Sst                              -0.308945
#> Rbp4                              1.165308
#>      ACTTGTTAGTTAAGTGstochastic_pseudotime
#> Pyy                              1.6066307
#> Sst                             -0.3091769
#> Rbp4                             1.1638839
#>      TGAGAGGAGGTAAACTstochastic_pseudotime
#> Pyy                              1.6041231
#> Sst                             -0.3092546
#> Rbp4                             1.1633188
#>      AAGACCTAGACCTTTGstochastic_pseudotime
#> Pyy                              1.6012676
#> Sst                             -0.3093391
#> Rbp4                             1.1626493
#>      TGCGGGTAGGGTTCCCstochastic_pseudotime
#> Pyy                              1.5979952
#> Sst                             -0.3094312
#> Rbp4                             1.1618511
#>      GATCGATAGGATGGTCstochastic_pseudotime
#> Pyy                              1.5974695
#> Sst                             -0.3094456
#> Rbp4                             1.1617200
#>      CCCAATCCAACTGGCCstochastic_pseudotime
#> Pyy                              1.5961858
#> Sst                             -0.3094802
#> Rbp4                             1.1613968
#>      TGTGTTTAGGAATGGAstochastic_pseudotime
#> Pyy                              1.5947051
#> Sst                             -0.3095193
#> Rbp4                             1.1610185
#>      GATCGATGTTGTTTGGstochastic_pseudotime
#> Pyy                              1.5908952
#> Sst                             -0.3096161
#> Rbp4                             1.1600204
#>      GTCCTCACACCACGTGstochastic_pseudotime
#> Pyy                              1.5890265
#> Sst                             -0.3096617
#> Rbp4                             1.1595186
#>      TATCTCACACGCGAAAstochastic_pseudotime
#> Pyy                              1.5640531
#> Sst                             -0.3101776
#> Rbp4                             1.1522105
#>      CTGGTCTAGGTGATTAstochastic_pseudotime
#> Pyy                              1.5612080
#> Sst                             -0.3102274
#> Rbp4                             1.1513217
#>      TTCTCCTAGTCCATACstochastic_pseudotime
#> Pyy                              1.5543415
#> Sst                             -0.3103417
#> Rbp4                             1.1491400
#>      TCGAGGCAGAACAATCstochastic_pseudotime
#> Pyy                              1.5470475
#> Sst                             -0.3104547
#> Rbp4                             1.1467719
#>      CCCAGTTGTCGTGGCTstochastic_pseudotime
#> Pyy                              1.5410028
#> Sst                             -0.3105424
#> Rbp4                             1.1447748
#>      TCAGGTAAGGCGTACAstochastic_pseudotime
#> Pyy                              1.5283784
#> Sst                             -0.3107106
#> Rbp4                             1.1405179
#>      TACTCGCAGTAAGTACstochastic_pseudotime
#> Pyy                              1.5162841
#> Sst                             -0.3108552
#> Rbp4                             1.1363490
#>      GACCAATGTTATCCGAstochastic_pseudotime
#> Pyy                              1.5059951
#> Sst                             -0.3109672
#> Rbp4                             1.1327460
#>      TCTGGAAGTCATCGGCstochastic_pseudotime
#> Pyy                              1.4793298
#> Sst                             -0.3112185
#> Rbp4                             1.1232278
#>      GTACGTAGTCTCTTTAstochastic_pseudotime
#> Pyy                              1.4784030
#> Sst                             -0.3112264
#> Rbp4                             1.1228933
#>      CGTCACTTCTAACTCTstochastic_pseudotime
#> Pyy                              1.4714350
#> Sst                             -0.3112838
#> Rbp4                             1.1203723
#>      GGGCATCCATTTCAGGstochastic_pseudotime
#> Pyy                              1.4686167
#> Sst                             -0.3113062
#> Rbp4                             1.1193497
#>      GCATGATGTTTAGGAAstochastic_pseudotime
#> Pyy                              1.4627297
#> Sst                             -0.3113515
#> Rbp4                             1.1172090
#>      GTATCTTGTACTCAACstochastic_pseudotime
#> Pyy                               1.447998
#> Sst                              -0.311457
#> Rbp4                              1.111829
#>      CTAACTTCATCTGGTAstochastic_pseudotime
#> Pyy                              1.4437697
#> Sst                             -0.3114853
#> Rbp4                             1.1102797
#>      GTCACGGAGCTATGCTstochastic_pseudotime
#> Pyy                              1.4340636
#> Sst                             -0.3115472
#> Rbp4                             1.1067178
#>      GTTCTCGAGAATTCCCstochastic_pseudotime
#> Pyy                              1.3979696
#> Sst                             -0.3117442
#> Rbp4                             1.0934418
#>      TCAGGTACAGTCACTAstochastic_pseudotime
#> Pyy                              1.3930726
#> Sst                             -0.3117674
#> Rbp4                             1.0916415
#>      GAATGAATCATCGGATstochastic_pseudotime
#> Pyy                              1.3873404
#> Sst                             -0.3117935
#> Rbp4                             1.0895358
#>      CGATGTAGTACTCTCCstochastic_pseudotime
#> Pyy                              1.3816770
#> Sst                             -0.3118183
#> Rbp4                             1.0874573
#>      ACACCCTAGGACAGCTstochastic_pseudotime
#> Pyy                              1.3619270
#> Sst                             -0.3118973
#> Rbp4                             1.0802301
#>      TGAGAGGAGCCAGGATstochastic_pseudotime
#> Pyy                              1.3608055
#> Sst                             -0.3119015
#> Rbp4                             1.0798210
#>      CGCGGTAAGCTGCAAGstochastic_pseudotime
#> Pyy                              1.3600897
#> Sst                             -0.3119041
#> Rbp4                             1.0795599
#>      GGCTGGTTCATGCATGstochastic_pseudotime
#> Pyy                               1.357094
#> Sst                              -0.311915
#> Rbp4                              1.078468
#>      TGCGTGGAGTACACCTstochastic_pseudotime
#> Pyy                              1.3510003
#> Sst                             -0.3119364
#> Rbp4                             1.0762503
#>      TGAGAGGCAGGAACGTstochastic_pseudotime
#> Pyy                              1.3491179
#> Sst                             -0.3119428
#> Rbp4                             1.0755663
#>      GCACTCTCAGCGTCCAstochastic_pseudotime
#> Pyy                              1.3479375
#> Sst                             -0.3119468
#> Rbp4                             1.0751376
#>      GGGATGAAGGGCATGTstochastic_pseudotime
#> Pyy                              1.3345404
#> Sst                             -0.3119895
#> Rbp4                             1.0702873
#>      CCCAGTTTCAGCACATstochastic_pseudotime
#> Pyy                              1.3088390
#> Sst                             -0.3120594
#> Rbp4                             1.0610728
#>      TTGACTTGTTCCATGAstochastic_pseudotime
#> Pyy                              1.3024084
#> Sst                             -0.3120746
#> Rbp4                             1.0587888
#>      TGTGTTTGTAGATTAGstochastic_pseudotime
#> Pyy                              1.3008512
#> Sst                             -0.3120781
#> Rbp4                             1.0582372
#>      TGGGCGTGTGTGCGTCstochastic_pseudotime
#> Pyy                              1.2912600
#> Sst                             -0.3120987
#> Rbp4                             1.0548524
#>      CAGCCGATCGTTTATCstochastic_pseudotime
#> Pyy                              1.2611579
#> Sst                             -0.3121509
#> Rbp4                             1.0443892
#>      GGACGTCCATGGATGGstochastic_pseudotime
#> Pyy                              1.2607046
#> Sst                             -0.3121516
#> Rbp4                             1.0442336
#>      AGGTCATGTCTCCCTAstochastic_pseudotime
#> Pyy                               1.245398
#> Sst                              -0.312171
#> Rbp4                              1.039021
#>      AGATTGCAGCGTGTCCstochastic_pseudotime
#> Pyy                              1.2204213
#> Sst                             -0.3121926
#> Rbp4                             1.0306907
#>      CACTCCACATCGTCGGstochastic_pseudotime
#> Pyy                              1.2021379
#> Sst                             -0.3122006
#> Rbp4                             1.0247477
#>      TTTGGTTAGACAGGCTstochastic_pseudotime
#> Pyy                              1.1952445
#> Sst                             -0.3122019
#> Rbp4                             1.0225440
#>      GCGCCAAAGCTCCTTCstochastic_pseudotime
#> Pyy                              1.1892109
#> Sst                             -0.3122022
#> Rbp4                             1.0206320
#>      CTGAAGTGTTTGGGCCstochastic_pseudotime
#> Pyy                              1.1876022
#> Sst                             -0.3122021
#> Rbp4                             1.0201248
#>      CGGAGTCAGTTTAGGAstochastic_pseudotime
#> Pyy                              1.1863873
#> Sst                             -0.3122021
#> Rbp4                             1.0197426
#>      CAGTCCTAGGTGGGTTstochastic_pseudotime
#> Pyy                               1.185231
#> Sst                              -0.312202
#> Rbp4                              1.019379
#>      AGTTGGTAGAGCCTAGstochastic_pseudotime
#> Pyy                              1.1836962
#> Sst                             -0.3122019
#> Rbp4                             1.0188980
#>      CAGAATCCATTGGGCCstochastic_pseudotime
#> Pyy                              1.1805660
#> Sst                             -0.3122014
#> Rbp4                             1.0179196
#>      GGGTCTGTCCATGAACstochastic_pseudotime
#> Pyy                              1.1382438
#> Sst                             -0.3121762
#> Rbp4                             1.0050919
#>      GTGTGCGGTGTAAGTAstochastic_pseudotime
#> Pyy                              1.1362444
#> Sst                             -0.3121741
#> Rbp4                             1.0045041
#>      GCATGTAAGCGTGTCCstochastic_pseudotime
#> Pyy                              1.1354675
#> Sst                             -0.3121733
#> Rbp4                             1.0042761
#>      ACCTTTATCCTGCAGGstochastic_pseudotime
#> Pyy                              1.1351429
#> Sst                             -0.3121729
#> Rbp4                             1.0041810
#>      GGCTGGTAGTCACGCCstochastic_pseudotime
#> Pyy                              1.1254629
#> Sst                             -0.3121615
#> Rbp4                             1.0013622
#>      ACAGCTATCAAGGCTTstochastic_pseudotime
#> Pyy                              1.1231738
#> Sst                             -0.3121585
#> Rbp4                             1.0007011
#>      TCACGAAGTCCGAGTCstochastic_pseudotime
#> Pyy                              1.1168664
#> Sst                             -0.3121498
#> Rbp4                             0.9988904
#>      ACAGCCGAGAGGTACCstochastic_pseudotime
#> Pyy                              1.1131146
#> Sst                             -0.3121442
#> Rbp4                             0.9978209
#>      TATTACCTCGTTGCCTstochastic_pseudotime
#> Pyy                              1.1091817
#> Sst                             -0.3121380
#> Rbp4                             0.9967058
#>      GAACCTATCACTTATCstochastic_pseudotime
#> Pyy                              1.0595111
#> Sst                             -0.3120312
#> Rbp4                             0.9831467
#>      CGATGGCTCAGTGTTGstochastic_pseudotime
#> Pyy                              1.0333209
#> Sst                             -0.3119516
#> Rbp4                             0.9763855
#>      GTCATTTGTCGTCTTCstochastic_pseudotime
#> Pyy                              1.0239599
#> Sst                             -0.3119189
#> Rbp4                             0.9740337
#>      TTTATGCAGGTTACCTstochastic_pseudotime
#> Pyy                              1.0237080
#> Sst                             -0.3119180
#> Rbp4                             0.9739709
#>      CACACCTGTACATGTCstochastic_pseudotime
#> Pyy                              1.0181676
#> Sst                             -0.3118974
#> Rbp4                             0.9725956
#>      AAATGCCCAATGGAATstochastic_pseudotime
#> Pyy                              1.0049277
#> Sst                             -0.3118448
#> Rbp4                             0.9693578
#>      ATTACTCTCTGCCAGGstochastic_pseudotime
#> Pyy                              1.0018898
#> Sst                             -0.3118320
#> Rbp4                             0.9686245
#>      TGTGGTAGTGCCTGCAstochastic_pseudotime
#> Pyy                              0.9875298
#> Sst                             -0.3117678
#> Rbp4                             0.9652079
#>      GCGGGTTCAAGTAATGstochastic_pseudotime
#> Pyy                              0.9839645
#> Sst                             -0.3117508
#> Rbp4                             0.9643723
#>      TCTATTGTCCTATTCAstochastic_pseudotime
#> Pyy                              0.9749321
#> Sst                             -0.3117059
#> Rbp4                             0.9622778
#>      TTTGGTTTCCTTTCGGstochastic_pseudotime
#> Pyy                              0.9687909
#> Sst                             -0.3116737
#> Rbp4                             0.9608725
#>      GGGAATGCAAAGTGCGstochastic_pseudotime
#> Pyy                              0.9586111
#> Sst                             -0.3116175
#> Rbp4                             0.9585763
#>      GCGAGAAGTGGGTATGstochastic_pseudotime
#> Pyy                              0.9487400
#> Sst                             -0.3115593
#> Rbp4                             0.9563899
#>      TACACGACATGTAGTCstochastic_pseudotime
#> Pyy                              0.9317254
#> Sst                             -0.3114496
#> Rbp4                             0.9527151
#>      CCGTTCAAGGACCACAstochastic_pseudotime
#> Pyy                              0.9300019
#> Sst                             -0.3114378
#> Rbp4                             0.9523495
#>      AATCGGTTCCAAACTGstochastic_pseudotime
#> Pyy                              0.9280480
#> Sst                             -0.3114243
#> Rbp4                             0.9519366
#>      CTGAAACGTAGCGCAAstochastic_pseudotime
#> Pyy                              0.9277674
#> Sst                             -0.3114223
#> Rbp4                             0.9518774
#>      GTCTTCGTCCGTTGCTstochastic_pseudotime
#> Pyy                              0.9254205
#> Sst                             -0.3114058
#> Rbp4                             0.9513839
#>      CCACCTAAGACTTGAAstochastic_pseudotime
#> Pyy                              0.9213699
#> Sst                             -0.3113767
#> Rbp4                             0.9505375
#>      ACTTGTTTCACCCGAGstochastic_pseudotime
#> Pyy                              0.9200946
#> Sst                             -0.3113674
#> Rbp4                             0.9502724
#>      CGCCAAGTCTAAGCCAstochastic_pseudotime
#> Pyy                              0.9180931
#> Sst                             -0.3113526
#> Rbp4                             0.9498578
#>      CGTGAGCTCTTCTGGCstochastic_pseudotime
#> Pyy                              0.9055405
#> Sst                             -0.3112554
#> Rbp4                             0.9472968
#>      AAATGCCCACCAACCGstochastic_pseudotime
#> Pyy                              0.8921646
#> Sst                             -0.3111428
#> Rbp4                             0.9446428
#>      CGTCAGGGTAGGCATGstochastic_pseudotime
#> Pyy                              0.8748498
#> Sst                             -0.3109819
#> Rbp4                             0.9413252
#>      CGACTTCTCACGCGGTstochastic_pseudotime
#> Pyy                              0.8747388
#> Sst                             -0.3109808
#> Rbp4                             0.9413044
#>      CCACGGAAGCAATCTCstochastic_pseudotime
#> Pyy                              0.8735618
#> Sst                             -0.3109692
#> Rbp4                             0.9410839
#>      GCTGCGAAGTATTGGAstochastic_pseudotime
#> Pyy                              0.8732114
#> Sst                             -0.3109657
#> Rbp4                             0.9410184
#>      CTCGAGGCAAGTCTGTstochastic_pseudotime
#> Pyy                              0.8671604
#> Sst                             -0.3109044
#> Rbp4                             0.9398956
#>      AACTCAGTCCTCGCATstochastic_pseudotime
#> Pyy                              0.8552190
#> Sst                             -0.3107762
#> Rbp4                             0.9377300
#>      GGGATGAGTCTGCGGTstochastic_pseudotime
#> Pyy                              0.8547586
#> Sst                             -0.3107710
#> Rbp4                             0.9376478
#>      GATTCAGAGTACGTTCstochastic_pseudotime
#> Pyy                              0.8487144
#> Sst                             -0.3107019
#> Rbp4                             0.9365788
#>      AGTAGTCGTGAGGGAGstochastic_pseudotime
#> Pyy                              0.8469535
#> Sst                             -0.3106813
#> Rbp4                             0.9362706
#>      TTCTCAACAATAGAGTstochastic_pseudotime
#> Pyy                              0.8407816
#> Sst                             -0.3106070
#> Rbp4                             0.9352025
#>      GTAGTCACATCGATGTstochastic_pseudotime
#> Pyy                              0.8333731
#> Sst                             -0.3105136
#> Rbp4                             0.9339452
#>      TACAGTGCACAGGCCTstochastic_pseudotime
#> Pyy                              0.8160890
#> Sst                             -0.3102770
#> Rbp4                             0.9311199
#>      AGTTGGTTCACATAGCstochastic_pseudotime
#> Pyy                              0.8133400
#> Sst                             -0.3102367
#> Rbp4                             0.9306848
#>      ACGCCAGGTTACAGAAstochastic_pseudotime
#> Pyy                              0.8107805
#> Sst                             -0.3101985
#> Rbp4                             0.9302833
#>      TTTGTCATCGAATGCTstochastic_pseudotime
#> Pyy                              0.8079624
#> Sst                             -0.3101557
#> Rbp4                             0.9298453
#>      GCTTCCAAGACATAACstochastic_pseudotime
#> Pyy                              0.8031241
#> Sst                             -0.3100801
#> Rbp4                             0.9291033
#>      TCCACACCAGGCGATAstochastic_pseudotime
#> Pyy                              0.8030044
#> Sst                             -0.3100782
#> Rbp4                             0.9290851
#>      ACTGATGCATGGTCTAstochastic_pseudotime
#> Pyy                              0.7816519
#> Sst                             -0.3097124
#> Rbp4                             0.9259670
#>      AACTCTTCATTTGCCCstochastic_pseudotime
#> Pyy                              0.7715430
#> Sst                             -0.3095188
#> Rbp4                             0.9245827
#>      CTCAGAACATAGGATAstochastic_pseudotime
#> Pyy                              0.7579232
#> Sst                             -0.3092343
#> Rbp4                             0.9228160
#>      GGACAGATCCTTTCGGstochastic_pseudotime
#> Pyy                              0.7560658
#> Sst                             -0.3091933
#> Rbp4                             0.9225840
#>      TACGGGCAGCCAACAGstochastic_pseudotime
#> Pyy                              0.7493197
#> Sst                             -0.3090393
#> Rbp4                             0.9217604
#>      GTTAAGCAGGATTCGGstochastic_pseudotime
#> Pyy                              0.7490603
#> Sst                             -0.3090332
#> Rbp4                             0.9217294
#>      GTTCTCGCAAGAAAGGstochastic_pseudotime
#> Pyy                              0.7444560
#> Sst                             -0.3089232
#> Rbp4                             0.9211853
#>      CTCGAAATCCTGCCATstochastic_pseudotime
#> Pyy                              0.7414940
#> Sst                             -0.3088505
#> Rbp4                             0.9208427
#>      TTCGGTCAGCTAGCCCstochastic_pseudotime
#> Pyy                              0.7346261
#> Sst                             -0.3086752
#> Rbp4                             0.9200716
#>      AAATGCCGTGATAAGTstochastic_pseudotime
#> Pyy                              0.7235122
#> Sst                             -0.3083709
#> Rbp4                             0.9188943
#>      CGTCCATTCGCGCCAAstochastic_pseudotime
#> Pyy                              0.7128816
#> Sst                             -0.3080535
#> Rbp4                             0.9178533
#>      CAGCTAAGTACCGAGAstochastic_pseudotime
#> Pyy                              0.7059079
#> Sst                             -0.3078300
#> Rbp4                             0.9172177
#>      GCATGCGCAGCTCCGAstochastic_pseudotime
#> Pyy                              0.7047838
#> Sst                             -0.3077928
#> Rbp4                             0.9171189
#>      GCTTGAAGTCAACTGTstochastic_pseudotime
#> Pyy                              0.7041556
#> Sst                             -0.3077718
#> Rbp4                             0.9170641
#>      GGGAGATAGTGACATAstochastic_pseudotime
#> Pyy                              0.7025405
#> Sst                             -0.3077174
#> Rbp4                             0.9169247
#>      TAGTTGGCACGTAAGGstochastic_pseudotime
#> Pyy                              0.6864750
#> Sst                             -0.3071341
#> Rbp4                             0.9156569
#>      ACTGATGCATCAGTACstochastic_pseudotime
#> Pyy                              0.6821130
#> Sst                             -0.3069613
#> Rbp4                             0.9153514
#>      GTTTCTACAGGACCCTstochastic_pseudotime
#> Pyy                              0.6798439
#> Sst                             -0.3068689
#> Rbp4                             0.9151994
#>      AGAGTGGAGGCGATACstochastic_pseudotime
#> Pyy                              0.6722896
#> Sst                             -0.3065473
#> Rbp4                             0.9147278
#>      TTGCGTCCACCGATATstochastic_pseudotime
#> Pyy                              0.6372986
#> Sst                             -0.3047285
#> Rbp4                             0.9133066
#>      ACGCAGCGTAATCACCstochastic_pseudotime
#> Pyy                              0.6332895
#> Sst                             -0.3044783
#> Rbp4                             0.9132327
#>      CCATGTCGTCGGCTCAstochastic_pseudotime
#> Pyy                              0.6317470
#> Sst                             -0.3043793
#> Rbp4                             0.9132096
#>      CTGCGGACAGACTCGCstochastic_pseudotime
#> Pyy                              0.6305597
#> Sst                             -0.3043021
#> Rbp4                             0.9131938
#>      CACAAACAGTAGCGGTstochastic_pseudotime
#> Pyy                              0.6264708
#> Sst                             -0.3040289
#> Rbp4                             0.9131536
#>      CTCACACGTCACACGCstochastic_pseudotime
#> Pyy                              0.6203395
#> Sst                             -0.3035971
#> Rbp4                             0.9131351
#>      AACCATGGTTCGGCACstochastic_pseudotime
#> Pyy                              0.6193023
#> Sst                             -0.3035213
#> Rbp4                             0.9131371
#>      CAGAGAGTCGGCCGATstochastic_pseudotime
#> Pyy                              0.6012420
#> Sst                             -0.3020565
#> Rbp4                             0.9134257
#>      CAAGGCCCACATTTCTstochastic_pseudotime
#> Pyy                              0.5998540
#> Sst                             -0.3019314
#> Rbp4                             0.9134690
#>      GACGTGCTCCGGCACAstochastic_pseudotime
#> Pyy                              0.5945889
#> Sst                             -0.3014386
#> Rbp4                             0.9136630
#>      CCCATACCATGCAACTstochastic_pseudotime
#> Pyy                              0.5913916
#> Sst                             -0.3011244
#> Rbp4                             0.9138043
#>      CTAAGACGTATCAGTCstochastic_pseudotime
#> Pyy                              0.5869880
#> Sst                             -0.3006719
#> Rbp4                             0.9140294
#>      ATCATGGGTATTCGTGstochastic_pseudotime
#> Pyy                              0.5813060
#> Sst                             -0.3000518
#> Rbp4                             0.9143744
#>      AGAGCGACAGTAAGATstochastic_pseudotime
#> Pyy                              0.5765951
#> Sst                             -0.2995041
#> Rbp4                             0.9147096
#>      CTTCTCTCAGCATACTstochastic_pseudotime
#> Pyy                              0.5733545
#> Sst                             -0.2991081
#> Rbp4                             0.9149673
#>      CCACGGAAGGGTGTTGstochastic_pseudotime
#> Pyy                              0.5689103
#> Sst                             -0.2985380
#> Rbp4                             0.9153585
#>      CCTATTATCATCGCTCstochastic_pseudotime
#> Pyy                              0.5659847
#> Sst                             -0.2981443
#> Rbp4                             0.9156409
#>      GGCAATTCATCTCCCAstochastic_pseudotime
#> Pyy                              0.5525888
#> Sst                             -0.2961305
#> Rbp4                             0.9172090
#>      GCTGGGTGTCATACTGstochastic_pseudotime
#> Pyy                              0.5508115
#> Sst                             -0.2958341
#> Rbp4                             0.9174538
#>      GAGCAGAGTTTGCATGstochastic_pseudotime
#> Pyy                              0.5365461
#> Sst                             -0.2931519
#> Rbp4                             0.9197768
#>      GCAAACTCAGACAAATstochastic_pseudotime
#> Pyy                              0.5313641
#> Sst                             -0.2920218
#> Rbp4                             0.9207970
#>      ATCATGGGTCCGCTGAstochastic_pseudotime
#> Pyy                              0.5211037
#> Sst                             -0.2894765
#> Rbp4                             0.9231451
#>      GGGAGATAGATGTGGCstochastic_pseudotime
#> Pyy                              0.5176025
#> Sst                             -0.2884984
#> Rbp4                             0.9240587
#>      CGTGAGCTCAGCGATTstochastic_pseudotime
#> Pyy                              0.5131103
#> Sst                             -0.2871476
#> Rbp4                             0.9253257
#>      GAACGGACAAGGACACstochastic_pseudotime
#> Pyy                              0.5089995
#> Sst                             -0.2858055
#> Rbp4                             0.9265872
#>      CTAAGACTCTGCGGCAstochastic_pseudotime
#> Pyy                              0.5026640
#> Sst                             -0.2835067
#> Rbp4                             0.9287467
#>      ACAGCCGCATCGTCGGstochastic_pseudotime
#> Pyy                              0.4902679
#> Sst                             -0.2779595
#> Rbp4                             0.9339042
#>      GTCGTAAAGTTTAGGAstochastic_pseudotime
#> Pyy                              0.4902220
#> Sst                             -0.2779357
#> Rbp4                             0.9339260
#>      GGGAATGTCACAAACCstochastic_pseudotime
#> Pyy                              0.4892534
#> Sst                             -0.2774275
#> Rbp4                             0.9343928
#>      GTCACAACACGTCTCTstochastic_pseudotime
#> Pyy                              0.4777987
#> Sst                             -0.2702536
#> Rbp4                             0.9408551
#>      TTTGGTTCACGAAGCAstochastic_pseudotime
#> Pyy                              0.4762996
#> Sst                             -0.2691148
#> Rbp4                             0.9418581
#>      GCGCAACGTTTCGCTCstochastic_pseudotime
#> Pyy                              0.4681080
#> Sst                             -0.2616634
#> Rbp4                             0.9482632
#>      TACCTTATCCGATATGstochastic_pseudotime
#> Pyy                              0.4652271
#> Sst                             -0.2583746
#> Rbp4                             0.9510055
#>      TGAGAGGAGGCGATACstochastic_pseudotime
#> Pyy                              0.4553950
#> Sst                             -0.2418819
#> Rbp4                             0.9640523
#>      TCCCGATCAACGATGGstochastic_pseudotime
#> Pyy                              0.4498453
#> Sst                             -0.2209560
#> Rbp4                             0.9792049
#>      TTGGCAAAGAGTAAGGstochastic_pseudotime
#> Pyy                              0.4490892
#> Sst                             -0.2141890
#> Rbp4                             0.9838357
#>      CCTTCGAAGATAGCATstochastic_pseudotime
#> Pyy                              0.4488182
#> Sst                             -0.2101413
#> Rbp4                             0.9865515
#>      GGCGACTGTAAGTTCCstochastic_pseudotime
#> Pyy                              0.4507531
#> Sst                             -0.1719296
#> Rbp4                             1.0103931
#>      TGAGGGAAGAGGTACCstochastic_pseudotime
#> Pyy                              0.4508157
#> Sst                             -0.1714351
#> Rbp4                             1.0106822
#>      ATCTGCCCAATCCGATstochastic_pseudotime
#> Pyy                              0.4520096
#> Sst                             -0.1628403
#> Rbp4                             1.0156365
#>      CACAGGCCACGTCAGCstochastic_pseudotime
#> Pyy                              0.4574206
#> Sst                             -0.1326483
#> Rbp4                             1.0320643
#>      TACTTACCATGGGACAstochastic_pseudotime
#> Pyy                              0.4589593
#> Sst                             -0.1252484
#> Rbp4                             1.0358816
#>      AATCGGTCATATGGTCstochastic_pseudotime
#> Pyy                              0.4621492
#> Sst                             -0.1107375
#> Rbp4                             1.0431535
#>      TTGAACGAGAGCAATTstochastic_pseudotime
#> Pyy                              0.4628953
#> Sst                             -0.1074657
#> Rbp4                             1.0447560
#>      CTTGGCTCAGCCACCAstochastic_pseudotime
#> Pyy                             0.46652228
#> Sst                            -0.09204027
#> Rbp4                            1.05213925
#>      CGAATGTTCTTCATGTstochastic_pseudotime
#> Pyy                              0.4693848
#> Sst                             -0.0802930
#> Rbp4                             1.0575824
#>      TGAGAGGAGCCGTCGTstochastic_pseudotime
#> Pyy                             0.47725839
#> Sst                            -0.04919162
#> Rbp4                            1.07131993
#>      AGCTTGACATCTATGGstochastic_pseudotime
#> Pyy                             0.47739353
#> Sst                            -0.04866877
#> Rbp4                            1.07154314
#>      ACTTGTTAGTGAAGTTstochastic_pseudotime
#> Pyy                             0.47982055
#> Sst                            -0.03932434
#> Rbp4                            1.07549225
#>      TACGGATTCCAAGTACstochastic_pseudotime
#> Pyy                             0.47989594
#> Sst                            -0.03903534
#> Rbp4                            1.07561320
#>      CTTAGGATCACTCTTAstochastic_pseudotime
#> Pyy                            0.492385004
#> Sst                            0.008145009
#> Rbp4                           1.094480652
#>      GGATGTTTCCAGTATGstochastic_pseudotime
#> Pyy                             0.49905886
#> Sst                             0.03305492
#> Rbp4                            1.10381148
#>      ACGCCGAGTTTGCATGstochastic_pseudotime
#> Pyy                             0.49907478
#> Sst                             0.03311421
#> Rbp4                            1.10383323
#>      TCTCATAGTGGACGATstochastic_pseudotime
#> Pyy                              0.5245617
#> Sst                              0.1284037
#> Rbp4                             1.1363183
#>      CGGCTAGGTACTCGCGstochastic_pseudotime
#> Pyy                              0.5427252
#> Sst                              0.1975243
#> Rbp4                             1.1573448
#>      CACAAACCACGAAAGCstochastic_pseudotime
#> Pyy                              0.5460296
#> Sst                              0.2102564
#> Rbp4                             1.1610269
#>      CTAAGACAGCACGCCTstochastic_pseudotime
#> Pyy                              0.5784550
#> Sst                              0.3381917
#> Rbp4                             1.1953157
#>      GTGGGTCCATAAAGGTstochastic_pseudotime
#> Pyy                              0.5809097
#> Sst                              0.3481110
#> Rbp4                             1.1977936
#>      ACACCAAGTATAGGTAstochastic_pseudotime
#> Pyy                              0.5953690
#> Sst                              0.4072402
#> Rbp4                             1.2121043
#>      CTGGTCTTCGCATGGCstochastic_pseudotime
#> Pyy                              0.6043675
#> Sst                              0.4446465
#> Rbp4                             1.2207815
#>      ACTGTCCTCGCTAGCGstochastic_pseudotime
#> Pyy                              0.6155619
#> Sst                              0.4918370
#> Rbp4                             1.2313526
#>      TCTCATACAAATTGCCstochastic_pseudotime
#> Pyy                              0.6236010
#> Sst                              0.5261754
#> Rbp4                             1.2388012
#>      GATCGATGTAGCGATGstochastic_pseudotime
#> Pyy                              0.6538193
#> Sst                              0.6585761
#> Rbp4                             1.2658363
#>      CAACCAACACCTCGTTstochastic_pseudotime
#> Pyy                              0.6590034
#> Sst                              0.6818110
#> Rbp4                             1.2703335
#>      CATATGGGTTACTGACstochastic_pseudotime
#> Pyy                              0.6708317
#> Sst                              0.7353863
#> Rbp4                             1.2804525
#>      AAACGGGCAAAGAATCstochastic_pseudotime
#> Pyy                              0.6891928
#> Sst                              0.8200670
#> Rbp4                             1.2957889
#>      CATCGGGTCTTCGGTCstochastic_pseudotime
#> Pyy                              0.7310498
#> Sst                              1.0196936
#> Rbp4                             1.3292149
#>      TCTCTAAGTGGCAAACstochastic_pseudotime
#> Pyy                              0.7446423
#> Sst                              1.0863810
#> Rbp4                             1.3396477
#>      CTCTGGTAGCAATCTCstochastic_pseudotime
#> Pyy                               0.766246
#> Sst                               1.194127
#> Rbp4                              1.355835
#>      CCGTGGATCTGTCTCGstochastic_pseudotime
#> Pyy                              0.8333465
#> Sst                              1.5410677
#> Rbp4                             1.4032492
#>      CGGAGCTAGATAGTCAstochastic_pseudotime
#> Pyy                              0.8445969
#> Sst                              1.6008334
#> Rbp4                             1.4108006
#>      GACCAATTCAACGGCCstochastic_pseudotime
#> Pyy                              0.9437851
#> Sst                              2.1421507
#> Rbp4                             1.4727755
#>      CATCGAAGTCACACGCstochastic_pseudotime
#> Pyy                              0.9440937
#> Sst                              2.1438649
#> Rbp4                             1.4729558
#>      GATTCAGTCGTCCAGGstochastic_pseudotime
#> Pyy                              0.9795565
#> Sst                              2.3415562
#> Rbp4                             1.4931640
#>      TACTCGCGTCGGCATCstochastic_pseudotime
#> Pyy                              0.9827291
#> Sst                              2.3592927
#> Rbp4                             1.4949227
#>      GTCTCGTTCACTTCATstochastic_pseudotime
#> Pyy                               0.990031
#> Sst                               2.400134
#> Rbp4                              1.498940
#>      GCCTCTAAGCTAGTCTstochastic_pseudotime
#> Pyy                               0.998469
#> Sst                               2.447356
#> Rbp4                              1.503529
#>      ATCCACCCAGCTCGCAstochastic_pseudotime
#> Pyy                               1.010850
#> Sst                               2.516664
#> Rbp4                              1.510159
#>      CGAGAAGAGGGTTCCCstochastic_pseudotime
#> Pyy                               1.035194
#> Sst                               2.652878
#> Rbp4                              1.522835
#>      GAAATGAAGCTCCTCTstochastic_pseudotime
#> Pyy                               1.064512
#> Sst                               2.816396
#> Rbp4                              1.537458
#>      GGGAGATAGCGAAGGGstochastic_pseudotime
#> Pyy                               1.066091
#> Sst                               2.825176
#> Rbp4                              1.538226
#>      CGGTTAAGTCTAGCCGstochastic_pseudotime
#> Pyy                               1.070407
#> Sst                               2.849154
#> Rbp4                              1.540312
#>      ACGCCAGTCCCAAGTAstochastic_pseudotime
#> Pyy                               1.080947
#> Sst                               2.907596
#> Rbp4                              1.545344
#>      CATGGCGAGAGTCTGGstochastic_pseudotime
#> Pyy                               1.099361
#> Sst                               3.009218
#> Rbp4                              1.553909
#>      CCACGGAGTCAACATCstochastic_pseudotime
#> Pyy                               1.102430
#> Sst                               3.026085
#> Rbp4                              1.555308
#>      GACGCGTGTCGGCTCAstochastic_pseudotime
#> Pyy                               1.127651
#> Sst                               3.163802
#> Rbp4                              1.566495
#>      GTATTCTTCTCTAAGGstochastic_pseudotime
#> Pyy                               1.127716
#> Sst                               3.164155
#> Rbp4                              1.566523
#>      CTTAGGATCAACGGGAstochastic_pseudotime
#> Pyy                               1.162763
#> Sst                               3.352219
#> Rbp4                              1.581119
#>      GGAATAAGTAAGAGAGstochastic_pseudotime
#> Pyy                               1.172718
#> Sst                               3.404777
#> Rbp4                              1.585056
#>      CGGACGTAGACACGACstochastic_pseudotime
#> Pyy                               1.188586
#> Sst                               3.487624
#> Rbp4                              1.591132
#>      GACGTTAAGACGCTTTstochastic_pseudotime
#> Pyy                               1.200080
#> Sst                               3.546862
#> Rbp4                              1.595377
#>      CCGTACTAGGAGTACCstochastic_pseudotime
#> Pyy                               1.279172
#> Sst                               3.931785
#> Rbp4                              1.620680
#>      CCGTACTAGTCCGTATstochastic_pseudotime
#> Pyy                               1.280587
#> Sst                               3.938235
#> Rbp4                              1.621064
#>      ATGTGTGTCAAAGACAstochastic_pseudotime
#> Pyy                               1.289482
#> Sst                               3.978342
#> Rbp4                              1.623421
#>      AACCATGCACAGGAGTstochastic_pseudotime
#> Pyy                               1.308067
#> Sst                               4.059645
#> Rbp4                              1.627995
#>      GAGTCCGAGCGATGACstochastic_pseudotime
#> Pyy                               1.309563
#> Sst                               4.066033
#> Rbp4                              1.628342
#>      GTATCTTGTCGCGGTTstochastic_pseudotime
#> Pyy                               1.312484
#> Sst                               4.078432
#> Rbp4                              1.629009
#>      AGATTGCAGTTCCACAstochastic_pseudotime
#> Pyy                               1.342231
#> Sst                               4.199033
#> Rbp4                              1.635044
#>      CTGAAACTCTGTTTGTstochastic_pseudotime
#> Pyy                               1.344697
#> Sst                               4.208526
#> Rbp4                              1.635478
#>      CCACGGAGTTACGACTstochastic_pseudotime
#> Pyy                               1.350026
#> Sst                               4.228761
#> Rbp4                              1.636380
#>      GTACTTTGTACCGTATstochastic_pseudotime
#> Pyy                               1.367060
#> Sst                               4.290655
#> Rbp4                              1.638907
#>      AGCGTCGCACAGGCCTstochastic_pseudotime
#> Pyy                               1.390984
#> Sst                               4.369608
#> Rbp4                              1.641446
#>      CTGAAACTCGATAGAAstochastic_pseudotime
#> Pyy                               1.394760
#> Sst                               4.381118
#> Rbp4                              1.641727
#>      TTGAACGAGGAGTTTAstochastic_pseudotime
#> Pyy                               1.398352
#> Sst                               4.391808
#> Rbp4                              1.641962
#>      ATCATCTGTGTTTGTGstochastic_pseudotime
#> Pyy                               1.403736
#> Sst                               4.407336
#> Rbp4                              1.642251
#>      GTTCATTTCCGAACGCstochastic_pseudotime
#> Pyy                               1.423189
#> Sst                               4.458075
#> Rbp4                              1.642625
#>      CGTGTAATCATTATCCstochastic_pseudotime
#> Pyy                               1.430658
#> Sst                               4.475093
#> Rbp4                              1.642457
#>      AGGGATGGTTAGGGTGstochastic_pseudotime
#> Pyy                               1.489041
#> Sst                               4.542619
#> Rbp4                              1.632686
#>      TGACTAGGTGTATGGGstochastic_pseudotime
#> Pyy                               1.492600
#> Sst                               4.541551
#> Rbp4                              1.631400
#>      AGCTCTCGTTCAGTACstochastic_pseudotime
#> Pyy                               1.495651
#> Sst                               4.539961
#> Rbp4                              1.630204
#>      TTTCCTCCAGGATCGAstochastic_pseudotime
#> Pyy                               1.502236
#> Sst                               4.534191
#> Rbp4                              1.627299
#>      TCTTTCCGTGCAACTTstochastic_pseudotime
#> Pyy                               1.502686
#> Sst                               4.533670
#> Rbp4                              1.627083
#>      CTGAAACGTATGGTTCstochastic_pseudotime
#> Pyy                               1.514249
#> Sst                               4.513724
#> Rbp4                              1.620595
#>      CGTCAGGTCCGCAAGCstochastic_pseudotime
#> Pyy                               1.528915
#> Sst                               4.463632
#> Rbp4                              1.608751
#>      GTACTCCTCCAACCAAstochastic_pseudotime
#> Pyy                               1.538605
#> Sst                               4.404293
#> Rbp4                              1.596960
#>      AGGCCGTCACGACGAAstochastic_pseudotime
#> Pyy                               1.539062
#> Sst                               4.400642
#> Rbp4                              1.596272
#>      CGAGAAGCACGCATCGstochastic_pseudotime
#> Pyy                               1.540524
#> Sst                               4.388278
#> Rbp4                              1.593964
#>      ATGTGTGCACAAGTAAstochastic_pseudotime
#> Pyy                               1.541808
#> Sst                               4.376429
#> Rbp4                              1.591782
#>      AGTCTTTGTTACGGAGstochastic_pseudotime
#> Pyy                               1.547943
#> Sst                               4.300336
#> Rbp4                              1.578259
#>      TGATTTCGTGAGGGTTstochastic_pseudotime
#> Pyy                               1.553417
#> Sst                               4.129683
#> Rbp4                              1.549485
#>      CCTATTACAGACAAGCstochastic_pseudotime
#> Pyy                               1.553607
#> Sst                               4.088257
#> Rbp4                              1.542643
#>      CCTACCATCTGCCCTAstochastic_pseudotime
#> Pyy                               1.553592
#> Sst                               4.074126
#> Rbp4                              1.540314
#>      AGCGGTCGTCGAATCTstochastic_pseudotime
#> Pyy                               1.552659
#> Sst                               3.984686
#> Rbp4                              1.525614
#>      GGGATGAGTTGTGGAGstochastic_pseudotime
#> Pyy                               1.549727
#> Sst                               3.872839
#> Rbp4                              1.507227
#>      GCTGGGTGTCTAGAGGstochastic_pseudotime
#> Pyy                               1.546370
#> Sst                               3.787769
#> Rbp4                              1.493174
#>      AGAGCTTAGGTTACCTstochastic_pseudotime
#> Pyy                               1.544509
#> Sst                               3.748184
#> Rbp4                              1.486602
#>      GATGAGGAGTGGACGTstochastic_pseudotime
#> Pyy                               1.543688
#> Sst                               3.731838
#> Rbp4                              1.483881
#>      GTCTCGTCATAGTAAGstochastic_pseudotime
#> Pyy                               1.540799
#> Sst                               3.678464
#> Rbp4                              1.474964
#>      CGTCCATGTCTAACGTstochastic_pseudotime
#> Pyy                               1.534028
#> Sst                               3.570475
#> Rbp4                              1.456748
#>      GTGCGGTCAGCGTAAGstochastic_pseudotime
#> Pyy                               1.524673
#> Sst                               3.444700
#> Rbp4                              1.435183
#>      GACCAATGTTGGACCCstochastic_pseudotime
#> Pyy                               1.504389
#> Sst                               3.221079
#> Rbp4                              1.395713
#>      ACATCAGTCCTTGGTCstochastic_pseudotime
#> Pyy                               1.498488
#> Sst                               3.163992
#> Rbp4                              1.385372
#>      AAGTCTGTCCGAGCCAstochastic_pseudotime
#> Pyy                               1.469958
#> Sst                               2.919734
#> Rbp4                              1.339709
#>      TGCCAAACATGACATCstochastic_pseudotime
#> Pyy                               1.417028
#> Sst                               2.553794
#> Rbp4                              1.266219
#>      TGGCTGGGTACCCAATstochastic_pseudotime
#> Pyy                               1.370103
#> Sst                               2.286895
#> Rbp4                              1.207882
#>      TGAGCATGTTAAAGACstochastic_pseudotime
#> Pyy                               1.333032
#> Sst                               2.102148
#> Rbp4                              1.164658
#>      GGCTCGAAGTGACTCTstochastic_pseudotime
#> Pyy                               1.298115
#> Sst                               1.944582
#> Rbp4                              1.125654
#>      GTTTCTAAGAAGATTCstochastic_pseudotime
#> Pyy                               1.268731
#> Sst                               1.822364
#> Rbp4                              1.093866
#>      TTCGGTCAGGGCATGTstochastic_pseudotime
#> Pyy                               1.214684
#> Sst                               1.618384
#> Rbp4                              1.037389
#>      TGACTAGTCCTATTCAstochastic_pseudotime
#> Pyy                               1.198514
#> Sst                               1.561952
#> Rbp4                              1.020918
#>      GTCCTCATCTGGAGCCstochastic_pseudotime
#> Pyy                               1.044239
#> Sst                               1.107002
#> Rbp4                              0.871059
#>      GTATCTTGTGGCCCTAstochastic_pseudotime
#> Pyy                              1.0378733
#> Sst                              1.0909611
#> Rbp4                             0.8651038
#>      AACTCAGGTGCCTGTGstochastic_pseudotime
#> Pyy                              0.7714794
#> Sst                              0.5575509
#> Rbp4                             0.6263670
#>      AAACGGGTCGCATGGCstochastic_pseudotime
#> Pyy                              0.7546940
#> Sst                              0.5313114
#> Rbp4                             0.6118497
#>      AGTGAGGTCTCCAACCstochastic_pseudotime
#> Pyy                              0.7477892
#> Sst                              0.5207261
#> Rbp4                             0.6058916
#>      TGCTGCTAGGATGCGTstochastic_pseudotime
#> Pyy                              0.6981047
#> Sst                              0.4479889
#> Rbp4                             0.5632417
#>      ACTTACTCAATCAGAAstochastic_pseudotime
#> Pyy                              0.6632595
#> Sst                              0.4003835
#> Rbp4                             0.5335449
#>      GAAATGAGTCACCCAGstochastic_pseudotime
#> Pyy                              0.6003933
#> Sst                              0.3210163
#> Rbp4                             0.4803532
#>      AACTTTCCACAGTCGCstochastic_pseudotime
#> Pyy                              0.5504693
#> Sst                              0.2634689
#> Rbp4                             0.4384006
#>      CATGACACATCGATTGstochastic_pseudotime
#> Pyy                             0.36494641
#> Sst                             0.08583078
#> Rbp4                            0.28389919
#>      GGGCATCGTGGTCTCGstochastic_pseudotime
#> Pyy                             0.32070830
#> Sst                             0.05077938
#> Rbp4                            0.24723485
#>      ACTATCTGTCAGGACAstochastic_pseudotime
#> Pyy                             0.22816536
#> Sst                            -0.01477002
#> Rbp4                            0.17057511
#>      TACTCGCCAGATCGGAstochastic_pseudotime
#> Pyy                             0.18456868
#> Sst                            -0.04229669
#> Rbp4                            0.13443736
#>      TATGCCCGTCTAGCCGstochastic_pseudotime
#> Pyy                            -0.05317559
#> Sst                            -0.16042112
#> Rbp4                           -0.06382280
#>      CTCGAAAGTAAGAGGAstochastic_pseudotime
#> Pyy                             -0.2307580
#> Sst                             -0.2198781
#> Rbp4                            -0.2144913
#>      CTAGTGATCCCAGGTGstochastic_pseudotime
#> Pyy                             -0.2930583
#> Sst                             -0.2360700
#> Rbp4                            -0.2681640
#>      CTGTGCTCATTACCTTstochastic_pseudotime
#> Pyy                             -0.3627737
#> Sst                             -0.2517114
#> Rbp4                            -0.3288647
#>      AGCTTGAAGAGTGACCstochastic_pseudotime
#> Pyy                             -0.4062068
#> Sst                             -0.2602427
#> Rbp4                            -0.3670722
#>      GCAGCCACAGCGTTCGstochastic_pseudotime
#> Pyy                             -0.4179805
#> Sst                             -0.2624047
#> Rbp4                            -0.3774860
#>      GGACAGAAGGCCGAATstochastic_pseudotime
#> Pyy                             -0.4207730
#> Sst                             -0.2629083
#> Rbp4                            -0.3799597
#>      TTTACTGCAGCATACTstochastic_pseudotime
#> Pyy                             -0.4291611
#> Sst                             -0.2644002
#> Rbp4                            -0.3873989
#>      ACCAGTATCGTTTAGGstochastic_pseudotime
#> Pyy                             -0.4776734
#> Sst                             -0.2724338
#> Rbp4                            -0.4306910
#>      CTCAGAACACTAGTACstochastic_pseudotime
#> Pyy                             -0.4995540
#> Sst                             -0.2757388
#> Rbp4                            -0.4503762
#>      CAGAATCGTCATTAGCstochastic_pseudotime
#> Pyy                             -0.5012830
#> Sst                             -0.2759918
#> Rbp4                            -0.4519362
#>      TGTGGTATCCAGTAGTstochastic_pseudotime
#> Pyy                             -0.5878483
#> Sst                             -0.2872297
#> Rbp4                            -0.5309518
#>      GGTATTGGTACTCGCGstochastic_pseudotime
#> Pyy                             -0.6030063
#> Sst                             -0.2889252
#> Rbp4                            -0.5449883
#>      TCAGGTATCTGTCCGTstochastic_pseudotime
#> Pyy                             -0.6530435
#> Sst                             -0.2939951
#> Rbp4                            -0.5918052
#>      GGGTCTGCACAGGCCTstochastic_pseudotime
#> Pyy                             -0.6719996
#> Sst                             -0.2957150
#> Rbp4                            -0.6097497
#>      ACTTACTAGAGGGCTTstochastic_pseudotime
#> Pyy                             -0.6835407
#> Sst                             -0.2967107
#> Rbp4                            -0.6207352
#>      AGGGAGTTCCTCAATTstochastic_pseudotime
#> Pyy                             -0.6989962
#> Sst                             -0.2979851
#> Rbp4                            -0.6355215
#>      GATCTAGAGCTAGCCCstochastic_pseudotime
#> Pyy                             -0.7272919
#> Sst                             -0.3001488
#> Rbp4                            -0.6628263
#>      GTATCTTCATCGATGTstochastic_pseudotime
#> Pyy                             -0.7485544
#> Sst                             -0.3016367
#> Rbp4                            -0.6835573
#>      GCGGGTTAGACTTGAAstochastic_pseudotime
#> Pyy                             -0.7954579
#> Sst                             -0.3045264
#> Rbp4                            -0.7300067
#>      TTATGCTCACCCAGTGstochastic_pseudotime
#> Pyy                             -0.8558189
#> Sst                             -0.3075200
#> Rbp4                            -0.7914823
#>      TTCGAAGAGTACGCCCstochastic_pseudotime
#> Pyy                             -0.8892554
#> Sst                             -0.3088618
#> Rbp4                            -0.8265305
#>      GCTGCGAAGCTAGCCCstochastic_pseudotime
#> Pyy                             -0.8951479
#> Sst                             -0.3090766
#> Rbp4                            -0.8327924
#>      GCCTCTAAGTCCGGTCstochastic_pseudotime
#> Pyy                             -0.9185352
#> Sst                             -0.3098683
#> Rbp4                            -0.8579272
#>      CCTCAGTGTTACCAGTstochastic_pseudotime
#> Pyy                             -0.9724446
#> Sst                             -0.3113433
#> Rbp4                            -0.9179135
#>      TCACAAGAGTGACATAstochastic_pseudotime
#> Pyy                             -0.9733598
#> Sst                             -0.3113644
#> Rbp4                            -0.9189617
#>      ATTATCCAGACCCACCstochastic_pseudotime
#> Pyy                             -0.9856169
#> Sst                             -0.3116333
#> Rbp4                            -0.9331160
#>      CTCGAGGTCTGTTGAGstochastic_pseudotime
#> Pyy                             -0.9858904
#> Sst                             -0.3116391
#> Rbp4                            -0.9334344
#>      TTCCCAGAGTTACCCAstochastic_pseudotime
#> Pyy                             -1.0041396
#> Sst                             -0.3119951
#> Rbp4                            -0.9549703
#>      AGCGGTCGTACCGAGAstochastic_pseudotime
#> Pyy                             -1.0046193
#> Sst                             -0.3120037
#> Rbp4                            -0.9555447
#>      TCGGGACAGAGAACAGstochastic_pseudotime
#> Pyy                             -1.0075576
#> Sst                             -0.3120558
#> Rbp4                            -0.9590736
#>      ATTATCCAGAAACCTAstochastic_pseudotime
#> Pyy                             -1.0185317
#> Sst                             -0.3122376
#> Rbp4                            -0.9724233
#>      CGCGTTTCAGTCTTCCstochastic_pseudotime
#> Pyy                             -1.0188231
#> Sst                             -0.3122421
#> Rbp4                            -0.9727818
#>      CCAATCCAGAGTAAGGstochastic_pseudotime
#> Pyy                             -1.0274667
#> Sst                             -0.3123697
#> Rbp4                            -0.9835210
#>      ACTTTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -1.0294854
#> Sst                             -0.3123974
#> Rbp4                            -0.9860614
#>      GCATGATAGCCTTGATstochastic_pseudotime
#> Pyy                             -1.0326627
#> Sst                             -0.3124392
#> Rbp4                            -0.9900875
#>      TATCTCACATTGTGCAstochastic_pseudotime
#> Pyy                             -1.0350218
#> Sst                             -0.3124689
#> Rbp4                            -0.9931001
#>      CAGCTGGGTTCTCATTstochastic_pseudotime
#> Pyy                             -1.0382802
#> Sst                             -0.3125077
#> Rbp4                            -0.9972969
#>      GTCTCGTAGGCTAGACstochastic_pseudotime
#> Pyy                             -1.0404268
#> Sst                             -0.3125318
#> Rbp4                            -1.0000862
#>      GTCCTCATCATTTGGGstochastic_pseudotime
#> Pyy                              -1.056039
#> Sst                              -0.312662
#> Rbp4                             -1.021150
#>      GAAACTCAGAAACCTAstochastic_pseudotime
#> Pyy                             -1.0620946
#> Sst                             -0.3126795
#> Rbp4                            -1.0298752
#>      CGACCTTAGCGCTTATstochastic_pseudotime
#> Pyy                             -1.0649631
#> Sst                             -0.3126762
#> Rbp4                            -1.0341923
#>      CGATCGGCACGCCAGTstochastic_pseudotime
#> Pyy                             -1.0653543
#> Sst                             -0.3126749
#> Rbp4                            -1.0347932
#>      TGATTTCAGTACGCCCstochastic_pseudotime
#> Pyy                              -1.067865
#> Sst                              -0.312661
#> Rbp4                             -1.038737
#>      CGTCAGGTCCTGCCATstochastic_pseudotime
#> Pyy                             -1.0693060
#> Sst                             -0.3126473
#> Rbp4                            -1.0410833
#>      GTGCGGTAGTAGGCCAstochastic_pseudotime
#> Pyy                             -1.0711039
#> Sst                             -0.3126218
#> Rbp4                            -1.0441257
#>      TCAGCTCAGTCGTACTstochastic_pseudotime
#> Pyy                             -1.0730950
#> Sst                             -0.3125776
#> Rbp4                            -1.0477075
#>      AACTTTCCACTGCCAGstochastic_pseudotime
#> Pyy                             -1.0742329
#> Sst                             -0.3125402
#> Rbp4                            -1.0499069
#>      TGCCCATCACTGTGTAstochastic_pseudotime
#> Pyy                             -1.0753655
#> Sst                             -0.3124884
#> Rbp4                            -1.0522747
#>      AGAGTGGGTCACCTAAstochastic_pseudotime
#> Pyy                             -1.0764336
#> Sst                             -0.3124152
#> Rbp4                            -1.0547914
#>      CTTGGCTCATCCCACTstochastic_pseudotime
#> Pyy                             -1.0764645
#> Sst                             -0.3124125
#> Rbp4                            -1.0548708
#>      CGTGTAAGTCTTGATGstochastic_pseudotime
#> Pyy                             -1.0777100
#> Sst                             -0.3122327
#> Rbp4                            -1.0588538
#>      CGGACTGTCTTGAGACstochastic_pseudotime
#> Pyy                             -1.0778816
#> Sst                             -0.3121737
#> Rbp4                            -1.0597697
#>      ACGAGGAAGTTACCCAstochastic_pseudotime
#> Pyy                              -1.078007
#> Sst                              -0.312003
#> Rbp4                             -1.061779
#>      TCTGGAATCGGTGTCGstochastic_pseudotime
#> Pyy                             -1.0780009
#> Sst                             -0.3119921
#> Rbp4                            -1.0618814
#>      TGGACGCAGTGTTGAAstochastic_pseudotime
#> Pyy                             -1.0779803
#> Sst                             -0.3119651
#> Rbp4                            -1.0621270
#>      AAGGCAGTCCGAGCCAstochastic_pseudotime
#> Pyy                             -1.0779787
#> Sst                             -0.3119634
#> Rbp4                            -1.0621419
#>      GGTGTTAAGGATCGCAstochastic_pseudotime
#> Pyy                             -1.0779672
#> Sst                             -0.3119518
#> Rbp4                            -1.0622425
#>      ATTGGACAGTAATCCCstochastic_pseudotime
#> Pyy                             -1.0779566
#> Sst                             -0.3119421
#> Rbp4                            -1.0623246
#>      CGAGCACTCTTCGAGAstochastic_pseudotime
#> Pyy                             -1.0779379
#> Sst                             -0.3119267
#> Rbp4                            -1.0624521
#>      CGACTTCCAGCCTATAstochastic_pseudotime
#> Pyy                             -1.0779238
#> Sst                             -0.3119162
#> Rbp4                            -1.0625363
#>      GGGTCTGCATGAAGTAstochastic_pseudotime
#> Pyy                             -1.0779173
#> Sst                             -0.3119116
#> Rbp4                            -1.0625729
#>      AAACGGGTCAGCTCTCstochastic_pseudotime
#> Pyy                             -1.0779012
#> Sst                             -0.3119007
#> Rbp4                            -1.0626575
#>      TCGTACCCATATACGCstochastic_pseudotime
#> Pyy                             -1.0779008
#> Sst                             -0.3119004
#> Rbp4                            -1.0626596
#>      CCGTTCACACATGACTstochastic_pseudotime
#> Pyy                             -1.0778751
#> Sst                             -0.3118843
#> Rbp4                            -1.0627813
#>      CGAATGTAGTTCCACAstochastic_pseudotime
#> Pyy                             -1.0777836
#> Sst                             -0.3118349
#> Rbp4                            -1.0631290
#>      AGAGCTTAGCCGATTTstochastic_pseudotime
#> Pyy                             -1.0776874
#> Sst                             -0.3117909
#> Rbp4                            -1.0634097
#>      CTTAGGAGTCCGAATTstochastic_pseudotime
#> Pyy                             -1.0776629
#> Sst                             -0.3117806
#> Rbp4                            -1.0634722
#>      CATGACACAACAACCTstochastic_pseudotime
#> Pyy                             -1.0776284
#> Sst                             -0.3117664
#> Rbp4                            -1.0635553
#>      TACGGATGTGTTCGATstochastic_pseudotime
#> Pyy                             -1.0775518
#> Sst                             -0.3117366
#> Rbp4                            -1.0637232
#>      TGAGCCGTCCGAACGCstochastic_pseudotime
#> Pyy                             -1.0774775
#> Sst                             -0.3117093
#> Rbp4                            -1.0638678
#>      ATCCGAAGTGACTCATstochastic_pseudotime
#> Pyy                             -1.0772730
#> Sst                             -0.3116403
#> Rbp4                            -1.0642002
#>      CTCTGGTCAACACCTAstochastic_pseudotime
#> Pyy                             -1.0771701
#> Sst                             -0.3116079
#> Rbp4                            -1.0643406
#>      ACATACGAGTAGGCCAstochastic_pseudotime
#> Pyy                             -1.0770302
#> Sst                             -0.3115659
#> Rbp4                            -1.0645096
#>      CAGATCAAGCTGTTCAstochastic_pseudotime
#> Pyy                             -1.0765611
#> Sst                             -0.3114359
#> Rbp4                            -1.0649486
#>      TGCGGGTTCTCAAACGstochastic_pseudotime
#> Pyy                             -1.0765587
#> Sst                             -0.3114353
#> Rbp4                            -1.0649505
#>      TCTTTCCTCGTCTGCTstochastic_pseudotime
#> Pyy                              -1.076396
#> Sst                              -0.311393
#> Rbp4                             -1.065069
#>      TTCTACAGTACTCTCCstochastic_pseudotime
#> Pyy                             -1.0762285
#> Sst                             -0.3113505
#> Rbp4                            -1.0651783
#>      CACATTTTCGCAGGCTstochastic_pseudotime
#> Pyy                             -1.0761713
#> Sst                             -0.3113362
#> Rbp4                            -1.0652127
#>      GATGAGGTCTGGGCCAstochastic_pseudotime
#> Pyy                             -1.0761592
#> Sst                             -0.3113332
#> Rbp4                            -1.0652198
#>      CTTTGCGAGGAATGGAstochastic_pseudotime
#> Pyy                             -1.0761044
#> Sst                             -0.3113197
#> Rbp4                            -1.0652513
#>      CGCTGGAGTCCATCCTstochastic_pseudotime
#> Pyy                             -1.0760429
#> Sst                             -0.3113045
#> Rbp4                            -1.0652854
#>      CCACTACAGAGTACCGstochastic_pseudotime
#> Pyy                             -1.0760155
#> Sst                             -0.3112978
#> Rbp4                            -1.0653001
#>      GCAGTTATCTTGAGACstochastic_pseudotime
#> Pyy                             -1.0760049
#> Sst                             -0.3112953
#> Rbp4                            -1.0653057
#>      GGCCGATTCACTCTTAstochastic_pseudotime
#> Pyy                             -1.0759943
#> Sst                             -0.3112927
#> Rbp4                            -1.0653113
#>      TGGCCAGCACGAAAGCstochastic_pseudotime
#> Pyy                             -1.0758658
#> Sst                             -0.3112615
#> Rbp4                            -1.0653760
#>      TCAGGTAAGTCAATAGstochastic_pseudotime
#> Pyy                             -1.0756852
#> Sst                             -0.3112184
#> Rbp4                            -1.0654584
#>      AACTTTCGTGACCAAGstochastic_pseudotime
#> Pyy                             -1.0756503
#> Sst                             -0.3112102
#> Rbp4                            -1.0654732
#>      CGGTTAATCGTCTGAAstochastic_pseudotime
#> Pyy                             -1.0756361
#> Sst                             -0.3112068
#> Rbp4                            -1.0654791
#>      GCCTCTATCGCTTAGAstochastic_pseudotime
#> Pyy                             -1.0755631
#> Sst                             -0.3111896
#> Rbp4                            -1.0655088
#>      GTAACGTGTCACTGGCstochastic_pseudotime
#> Pyy                             -1.0754224
#> Sst                             -0.3111568
#> Rbp4                            -1.0655622
#>      TTCTACAAGATCGATAstochastic_pseudotime
#> Pyy                             -1.0753679
#> Sst                             -0.3111442
#> Rbp4                            -1.0655816
#>      TACACGACAACGATCTstochastic_pseudotime
#> Pyy                             -1.0750583
#> Sst                             -0.3110733
#> Rbp4                            -1.0656797
#>      ACTGAGTAGTGACTCTstochastic_pseudotime
#> Pyy                             -1.0750008
#> Sst                             -0.3110603
#> Rbp4                            -1.0656957
#>      GAGGTGAGTGACAAATstochastic_pseudotime
#> Pyy                             -1.0749938
#> Sst                             -0.3110587
#> Rbp4                            -1.0656977
#>      TTGAACGGTGATGCCCstochastic_pseudotime
#> Pyy                             -1.0749246
#> Sst                             -0.3110431
#> Rbp4                            -1.0657161
#>      CAAGGCCAGCTGTTCAstochastic_pseudotime
#> Pyy                             -1.0747715
#> Sst                             -0.3110086
#> Rbp4                            -1.0657539
#>      GTATCTTCAGTGGGATstochastic_pseudotime
#> Pyy                             -1.0747488
#> Sst                             -0.3110036
#> Rbp4                            -1.0657591
#>      AGGCCACAGCTGTCTAstochastic_pseudotime
#> Pyy                             -1.0746380
#> Sst                             -0.3109788
#> Rbp4                            -1.0657836
#>      ACTTACTAGCCGGTAAstochastic_pseudotime
#> Pyy                             -1.0741639
#> Sst                             -0.3108742
#> Rbp4                            -1.0658677
#>      GATGAAACATGGTTGTstochastic_pseudotime
#> Pyy                             -1.0739240
#> Sst                             -0.3108218
#> Rbp4                            -1.0658989
#>      GACGCGTCAGTCGTGCstochastic_pseudotime
#> Pyy                             -1.0737311
#> Sst                             -0.3107799
#> Rbp4                            -1.0659191
#>      TACGGGCAGCGTCAAGstochastic_pseudotime
#> Pyy                             -1.0736231
#> Sst                             -0.3107565
#> Rbp4                            -1.0659286
#>      GACGTTAAGACGCAACstochastic_pseudotime
#> Pyy                             -1.0729863
#> Sst                             -0.3106194
#> Rbp4                            -1.0659615
#>      CCTACACAGATCCGAGstochastic_pseudotime
#> Pyy                             -1.0728520
#> Sst                             -0.3105907
#> Rbp4                            -1.0659638
#>      AGCTCCTGTGTGCCTGstochastic_pseudotime
#> Pyy                             -1.0727784
#> Sst                             -0.3105749
#> Rbp4                            -1.0659644
#>      TGTCCCAAGTTTAGGAstochastic_pseudotime
#> Pyy                             -1.0724028
#> Sst                             -0.3104947
#> Rbp4                            -1.0659612
#>      CGTCACTAGGGTGTGTstochastic_pseudotime
#> Pyy                             -1.0722054
#> Sst                             -0.3104526
#> Rbp4                            -1.0659555
#>      CGTCACTCACGAAATAstochastic_pseudotime
#> Pyy                              -1.072174
#> Sst                              -0.310446
#> Rbp4                             -1.065954
#>      AGCTTGACAATAGAGTstochastic_pseudotime
#> Pyy                             -1.0721689
#> Sst                             -0.3104448
#> Rbp4                            -1.0659542
#>      CTTGGCTAGGACACCAstochastic_pseudotime
#> Pyy                             -1.0719436
#> Sst                             -0.3103968
#> Rbp4                            -1.0659441
#>      CTCAGAATCAGGTTCAstochastic_pseudotime
#> Pyy                             -1.0716892
#> Sst                             -0.3103426
#> Rbp4                            -1.0659290
#>      TATCTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -1.0712619
#> Sst                             -0.3102514
#> Rbp4                            -1.0658955
#>      TTGGAACGTATATGAGstochastic_pseudotime
#> Pyy                             -1.0708010
#> Sst                             -0.3101529
#> Rbp4                            -1.0658492
#>      TTTACTGGTGCTAGCCstochastic_pseudotime
#> Pyy                             -1.0706629
#> Sst                             -0.3101233
#> Rbp4                            -1.0658334
#>      CTCGTACCAATGGATAstochastic_pseudotime
#> Pyy                             -1.0705159
#> Sst                             -0.3100918
#> Rbp4                            -1.0658157
#>      AGTAGTCCACATGACTstochastic_pseudotime
#> Pyy                             -1.0702956
#> Sst                             -0.3100445
#> Rbp4                            -1.0657874
#>      TCAATCTCATCTCCCAstochastic_pseudotime
#> Pyy                             -1.0701249
#> Sst                             -0.3100077
#> Rbp4                            -1.0657642
#>      TTGCCGTTCTCTGTCGstochastic_pseudotime
#> Pyy                             -1.0697519
#> Sst                             -0.3099273
#> Rbp4                            -1.0657097
#>      AGAGTGGGTAGCACGAstochastic_pseudotime
#> Pyy                             -1.0687523
#> Sst                             -0.3097097
#> Rbp4                            -1.0655401
#>      CCACCTATCACCAGGCstochastic_pseudotime
#> Pyy                              -1.068612
#> Sst                              -0.309679
#> Rbp4                             -1.065514
#>      CGTTAGAGTGCAACGAstochastic_pseudotime
#> Pyy                             -1.0685678
#> Sst                             -0.3096692
#> Rbp4                            -1.0655055
#>      GATCTAGGTGCCTTGGstochastic_pseudotime
#> Pyy                             -1.0684184
#> Sst                             -0.3096364
#> Rbp4                            -1.0654767
#>      CGTCAGGCATGGGAACstochastic_pseudotime
#> Pyy                             -1.0682553
#> Sst                             -0.3096004
#> Rbp4                            -1.0654446
#>      TCATTACGTCAGAAGCstochastic_pseudotime
#> Pyy                             -1.0682097
#> Sst                             -0.3095903
#> Rbp4                            -1.0654355
#>      GCGGGTTGTTACTGACstochastic_pseudotime
#> Pyy                             -1.0679698
#> Sst                             -0.3095372
#> Rbp4                            -1.0653867
#>      CCTTACGAGATCACGGstochastic_pseudotime
#> Pyy                             -1.0676861
#> Sst                             -0.3094742
#> Rbp4                            -1.0653271
#>      GCATGATTCAAACGGGstochastic_pseudotime
#> Pyy                             -1.0674816
#> Sst                             -0.3094286
#> Rbp4                            -1.0652829
#>      AGTAGTCGTTCTGTTTstochastic_pseudotime
#> Pyy                             -1.0674488
#> Sst                             -0.3094212
#> Rbp4                            -1.0652757
#>      TGAGCCGTCATGTCCCstochastic_pseudotime
#> Pyy                             -1.0672823
#> Sst                             -0.3093839
#> Rbp4                            -1.0652389
#>      GCAATCATCCAACCAAstochastic_pseudotime
#> Pyy                             -1.0672805
#> Sst                             -0.3093835
#> Rbp4                            -1.0652385
#>      CAGAATCGTGTAAGTAstochastic_pseudotime
#> Pyy                             -1.0663839
#> Sst                             -0.3091809
#> Rbp4                            -1.0650296
#>      GTTCATTAGAAGATTCstochastic_pseudotime
#> Pyy                             -1.0660756
#> Sst                             -0.3091104
#> Rbp4                            -1.0649540
#>      TGCGTGGAGGTAAACTstochastic_pseudotime
#> Pyy                             -1.0653379
#> Sst                             -0.3089402
#> Rbp4                            -1.0647657
#>      TGATTTCAGTGACTCTstochastic_pseudotime
#> Pyy                             -1.0652870
#> Sst                             -0.3089284
#> Rbp4                            -1.0647524
#>      TCAGCAACACCAGGCTstochastic_pseudotime
#> Pyy                             -1.0647653
#> Sst                             -0.3088065
#> Rbp4                            -1.0646130
#>      CGTCAGGAGCCTATGTstochastic_pseudotime
#> Pyy                             -1.0645969
#> Sst                             -0.3087668
#> Rbp4                            -1.0645671
#>      GGGCACTCATTTCAGGstochastic_pseudotime
#> Pyy                             -1.0643440
#> Sst                             -0.3087071
#> Rbp4                            -1.0644973
#>      GAGTCCGAGTTAAGTGstochastic_pseudotime
#> Pyy                             -1.0638875
#> Sst                             -0.3085984
#> Rbp4                            -1.0643687
#>      TATCAGGCACCAGATTstochastic_pseudotime
#> Pyy                             -1.0638503
#> Sst                             -0.3085895
#> Rbp4                            -1.0643581
#>      GTTCTCGCAGACACTTstochastic_pseudotime
#> Pyy                             -1.0634877
#> Sst                             -0.3085025
#> Rbp4                            -1.0642537
#>      GTATTCTCACATTAGCstochastic_pseudotime
#> Pyy                             -1.0633322
#> Sst                             -0.3084649
#> Rbp4                            -1.0642083
#>      AACTCCCTCGAGAGCAstochastic_pseudotime
#> Pyy                              -1.063328
#> Sst                              -0.308464
#> Rbp4                             -1.064207
#>      TGACTTTAGCGTTTACstochastic_pseudotime
#> Pyy                             -1.0632719
#> Sst                             -0.3084504
#> Rbp4                            -1.0641906
#>      CCTTCGACAGCTGTTAstochastic_pseudotime
#> Pyy                             -1.0630917
#> Sst                             -0.3084067
#> Rbp4                            -1.0641375
#>      ACACCGGCATCACGTAstochastic_pseudotime
#> Pyy                             -1.0630341
#> Sst                             -0.3083927
#> Rbp4                            -1.0641204
#>      AAACGGGAGATATGGTstochastic_pseudotime
#> Pyy                             -1.0626388
#> Sst                             -0.3082962
#> Rbp4                            -1.0640021
#>      CTCGTCATCGGTCTAAstochastic_pseudotime
#> Pyy                             -1.0623695
#> Sst                             -0.3082299
#> Rbp4                            -1.0639203
#>      GAAATGAGTCGGATCCstochastic_pseudotime
#> Pyy                              -1.061749
#> Sst                              -0.308076
#> Rbp4                             -1.063729
#>      TCTCTAAGTCGCATCGstochastic_pseudotime
#> Pyy                             -1.0617388
#> Sst                             -0.3080735
#> Rbp4                            -1.0637253
#>      TGGTTCCCACCGCTAGstochastic_pseudotime
#> Pyy                             -1.0617059
#> Sst                             -0.3080652
#> Rbp4                            -1.0637149
#>      GACGGCTAGGTGATATstochastic_pseudotime
#> Pyy                             -1.0615182
#> Sst                             -0.3080182
#> Rbp4                            -1.0636559
#>      GACTACAGTTCGTGATstochastic_pseudotime
#> Pyy                             -1.0613362
#> Sst                             -0.3079725
#> Rbp4                            -1.0635983
#>      TCACGAATCCTATGTTstochastic_pseudotime
#> Pyy                             -1.0612380
#> Sst                             -0.3079477
#> Rbp4                            -1.0635671
#>      TGAAAGAGTTATCACGstochastic_pseudotime
#> Pyy                             -1.0609886
#> Sst                             -0.3078846
#> Rbp4                            -1.0634872
#>      ATAGACCTCTAACTGGstochastic_pseudotime
#> Pyy                             -1.0607761
#> Sst                             -0.3078306
#> Rbp4                            -1.0634187
#>      CAAGGCCAGTGTGGCAstochastic_pseudotime
#> Pyy                             -1.0607018
#> Sst                             -0.3078117
#> Rbp4                            -1.0633946
#>      ACTTTCAAGGATATACstochastic_pseudotime
#> Pyy                              -1.060692
#> Sst                              -0.307809
#> Rbp4                             -1.063391
#>      TACGGATAGTTTAGGAstochastic_pseudotime
#> Pyy                             -1.0606510
#> Sst                             -0.3077987
#> Rbp4                            -1.0633781
#>      CTCTGGTCAGTTCATGstochastic_pseudotime
#> Pyy                             -1.0606199
#> Sst                             -0.3077907
#> Rbp4                            -1.0633680
#>      TGAGGGACAGACTCGCstochastic_pseudotime
#> Pyy                             -1.0605363
#> Sst                             -0.3077693
#> Rbp4                            -1.0633408
#>      GCTGCAGGTTACCGATstochastic_pseudotime
#> Pyy                             -1.0603987
#> Sst                             -0.3077341
#> Rbp4                            -1.0632958
#>      GGGTTGCTCGAACGGAstochastic_pseudotime
#> Pyy                             -1.0603633
#> Sst                             -0.3077249
#> Rbp4                            -1.0632841
#>      CTAAGACAGCCACGCTstochastic_pseudotime
#> Pyy                             -1.0599921
#> Sst                             -0.3076292
#> Rbp4                            -1.0631617
#>      CGTCACTCATAGTAAGstochastic_pseudotime
#> Pyy                             -1.0598886
#> Sst                             -0.3076023
#> Rbp4                            -1.0631273
#>      GGTGTTAGTTCCTCCAstochastic_pseudotime
#> Pyy                             -1.0596814
#> Sst                             -0.3075485
#> Rbp4                            -1.0630581
#>      AAAGATGAGTCTTGCAstochastic_pseudotime
#> Pyy                              -1.059676
#> Sst                              -0.307547
#> Rbp4                             -1.063056
#>      AACGTTGGTGCCTTGGstochastic_pseudotime
#> Pyy                             -1.0593488
#> Sst                             -0.3074614
#> Rbp4                            -1.0629462
#>      CATGACAGTAGCTTGTstochastic_pseudotime
#> Pyy                             -1.0592475
#> Sst                             -0.3074348
#> Rbp4                            -1.0629119
#>      CTGGTCTCAGGAATGCstochastic_pseudotime
#> Pyy                             -1.0589162
#> Sst                             -0.3073474
#> Rbp4                            -1.0627991
#>      CACAGTATCCTCCTAGstochastic_pseudotime
#> Pyy                             -1.0587635
#> Sst                             -0.3073068
#> Rbp4                            -1.0627467
#>      CGTAGGCAGCGTAGTGstochastic_pseudotime
#> Pyy                             -1.0587413
#> Sst                             -0.3073009
#> Rbp4                            -1.0627391
#>      CCTTCCCCAAAGTGCGstochastic_pseudotime
#> Pyy                             -1.0584242
#> Sst                             -0.3072163
#> Rbp4                            -1.0626296
#>      GACACGCCAAAGGAAGstochastic_pseudotime
#> Pyy                             -1.0581594
#> Sst                             -0.3071453
#> Rbp4                            -1.0625376
#>      AGCAGCCTCTTTAGTCstochastic_pseudotime
#> Pyy                             -1.0570692
#> Sst                             -0.3068485
#> Rbp4                            -1.0621522
#>      ACGGCCACACGTGAGAstochastic_pseudotime
#> Pyy                             -1.0567832
#> Sst                             -0.3067695
#> Rbp4                            -1.0620496
#>      TACTTACCAACTGCGCstochastic_pseudotime
#> Pyy                             -1.0567668
#> Sst                             -0.3067649
#> Rbp4                            -1.0620436
#>      TAGTTGGGTACCGTATstochastic_pseudotime
#> Pyy                             -1.0561144
#> Sst                             -0.3065829
#> Rbp4                            -1.0618069
#>      AGACGTTTCTTACCTAstochastic_pseudotime
#> Pyy                             -1.0561026
#> Sst                             -0.3065796
#> Rbp4                            -1.0618026
#>      GAATGAACAATCGAAAstochastic_pseudotime
#> Pyy                             -1.0552744
#> Sst                             -0.3063449
#> Rbp4                            -1.0614974
#>      ACTGAACAGCCTCGTGstochastic_pseudotime
#> Pyy                              -1.055243
#> Sst                              -0.306336
#> Rbp4                             -1.061486
#>      CGAACATCAGACGTAGstochastic_pseudotime
#> Pyy                             -1.0550229
#> Sst                             -0.3062728
#> Rbp4                            -1.0614038
#>      CAAGAAATCGTACCGGstochastic_pseudotime
#> Pyy                             -1.0549860
#> Sst                             -0.3062622
#> Rbp4                            -1.0613900
#>      GTACTCCCACGTCTCTstochastic_pseudotime
#> Pyy                             -1.0547828
#> Sst                             -0.3062036
#> Rbp4                            -1.0613139
#>      CTTAACTCAAGCGTAGstochastic_pseudotime
#> Pyy                             -1.0542903
#> Sst                             -0.3060607
#> Rbp4                            -1.0611284
#>      GGGTTGCCAAGCTGAGstochastic_pseudotime
#> Pyy                             -1.0536086
#> Sst                             -0.3058604
#> Rbp4                            -1.0608689
#>      AAAGCAACACATGGGAstochastic_pseudotime
#> Pyy                             -1.0529774
#> Sst                             -0.3056725
#> Rbp4                            -1.0606260
#>      GGGACCTTCGAATGCTstochastic_pseudotime
#> Pyy                             -1.0528316
#> Sst                             -0.3056287
#> Rbp4                            -1.0605695
#>      TAGACCAAGGAATCGCstochastic_pseudotime
#> Pyy                              -1.052766
#> Sst                              -0.305609
#> Rbp4                             -1.060544
#>      CTAACTTCATTGGGCCstochastic_pseudotime
#> Pyy                             -1.0525829
#> Sst                             -0.3055538
#> Rbp4                            -1.0604728
#>      TGAGCATTCACCCGAGstochastic_pseudotime
#> Pyy                             -1.0523417
#> Sst                             -0.3054808
#> Rbp4                            -1.0603787
#>      GAAGCAGTCAAGGCTTstochastic_pseudotime
#> Pyy                             -1.0520617
#> Sst                             -0.3053955
#> Rbp4                            -1.0602690
#>      AAACCTGAGCCTTGATstochastic_pseudotime
#> Pyy                             -1.0520308
#> Sst                             -0.3053861
#> Rbp4                            -1.0602569
#>      CCTTCCCAGGCAGGTTstochastic_pseudotime
#> Pyy                             -1.0514758
#> Sst                             -0.3052156
#> Rbp4                            -1.0600380
#>      TACAGTGCAATCTACGstochastic_pseudotime
#> Pyy                              -1.051154
#> Sst                              -0.305116
#> Rbp4                             -1.059910
#>      GGACGTCCACGCATCGstochastic_pseudotime
#> Pyy                             -1.0508974
#> Sst                             -0.3050359
#> Rbp4                            -1.0598079
#>      CACAGTAGTTCCTCCAstochastic_pseudotime
#> Pyy                             -1.0495696
#> Sst                             -0.3046155
#> Rbp4                            -1.0592725
#>      TACGGTACATGTTCCCstochastic_pseudotime
#> Pyy                             -1.0493398
#> Sst                             -0.3045416
#> Rbp4                            -1.0591788
#>      CTCGTCACAACTGCTAstochastic_pseudotime
#> Pyy                             -1.0490236
#> Sst                             -0.3044394
#> Rbp4                            -1.0590495
#>      CGGACTGAGGTAGCTGstochastic_pseudotime
#> Pyy                             -1.0486349
#> Sst                             -0.3043129
#> Rbp4                            -1.0588898
#>      AACTGGTTCTTACCGCstochastic_pseudotime
#> Pyy                             -1.0480643
#> Sst                             -0.3041254
#> Rbp4                            -1.0586538
#>      ATTTCTGGTAGTGAATstochastic_pseudotime
#> Pyy                             -1.0478898
#> Sst                             -0.3040677
#> Rbp4                            -1.0585813
#>      GAAATGACACATGACTstochastic_pseudotime
#> Pyy                             -1.0475668
#> Sst                             -0.3039603
#> Rbp4                            -1.0584467
#>      CTCAGAAGTCCTCTTGstochastic_pseudotime
#> Pyy                             -1.0463045
#> Sst                             -0.3035342
#> Rbp4                            -1.0579156
#>      CACCTTGAGAGTAATCstochastic_pseudotime
#> Pyy                             -1.0440679
#> Sst                             -0.3027542
#> Rbp4                            -1.0569552
#>      TACCTTAAGTACATGAstochastic_pseudotime
#> Pyy                             -1.0440638
#> Sst                             -0.3027527
#> Rbp4                            -1.0569534
#>      GAAATGACAACTGGCCstochastic_pseudotime
#> Pyy                             -1.0437292
#> Sst                             -0.3026332
#> Rbp4                            -1.0568076
#>      TACGGATAGGCACATGstochastic_pseudotime
#> Pyy                             -1.0434031
#> Sst                             -0.3025161
#> Rbp4                            -1.0566650
#>      ATTGGTGTCCAGTATGstochastic_pseudotime
#> Pyy                             -1.0425447
#> Sst                             -0.3022045
#> Rbp4                            -1.0562873
#>      CTAGCCTTCAGCGATTstochastic_pseudotime
#> Pyy                             -1.0411091
#> Sst                             -0.3016727
#> Rbp4                            -1.0556480
#>      GGGAATGAGTGTACGGstochastic_pseudotime
#> Pyy                             -1.0407391
#> Sst                             -0.3015335
#> Rbp4                            -1.0554816
#>      CACATTTGTTTCGCTCstochastic_pseudotime
#> Pyy                             -1.0404981
#> Sst                             -0.3014424
#> Rbp4                            -1.0553729
#>      TGCTACCCAGGGCATAstochastic_pseudotime
#> Pyy                             -1.0403312
#> Sst                             -0.3013791
#> Rbp4                            -1.0552975
#>      TAAGTGCAGGGTTTCTstochastic_pseudotime
#> Pyy                             -1.0401748
#> Sst                             -0.3013195
#> Rbp4                            -1.0552267
#>      TTTGCGCAGACAGACCstochastic_pseudotime
#> Pyy                             -1.0391676
#> Sst                             -0.3009325
#> Rbp4                            -1.0547681
#>      TTTGGTTCACCAGATTstochastic_pseudotime
#> Pyy                             -1.0386991
#> Sst                             -0.3007503
#> Rbp4                            -1.0545531
#>      GGACAAGGTGATAAGTstochastic_pseudotime
#> Pyy                             -1.0360077
#> Sst                             -0.2996766
#> Rbp4                            -1.0532987
#>      GGTATTGGTCCCGACAstochastic_pseudotime
#> Pyy                             -1.0357387
#> Sst                             -0.2995668
#> Rbp4                            -1.0531715
#>      CTGAAACAGACTAGGCstochastic_pseudotime
#> Pyy                             -1.0352015
#> Sst                             -0.2993462
#> Rbp4                            -1.0529165
#>      ACACCGGCATGTCCTCstochastic_pseudotime
#> Pyy                              -1.033982
#> Sst                              -0.298839
#> Rbp4                             -1.052333
#>      TTGACTTTCAGCAACTstochastic_pseudotime
#> Pyy                             -1.0332844
#> Sst                             -0.2985445
#> Rbp4                            -1.0519955
#>      AGGGATGTCGCATGATstochastic_pseudotime
#> Pyy                             -1.0331956
#> Sst                             -0.2985068
#> Rbp4                            -1.0519525
#>      CGAGCACTCTACCTGCstochastic_pseudotime
#> Pyy                             -1.0329629
#> Sst                             -0.2984079
#> Rbp4                            -1.0518395
#>      TTCTCCTCAATGGTCTstochastic_pseudotime
#> Pyy                             -1.0318993
#> Sst                             -0.2979514
#> Rbp4                            -1.0513198
#>      GGGCACTGTGGCAAACstochastic_pseudotime
#> Pyy                             -1.0306591
#> Sst                             -0.2974108
#> Rbp4                            -1.0507073
#>      TTGACTTTCCAGTAGTstochastic_pseudotime
#> Pyy                             -1.0305272
#> Sst                             -0.2973528
#> Rbp4                            -1.0506418
#>      CAGCAGCGTTCACGGCstochastic_pseudotime
#> Pyy                             -1.0290584
#> Sst                             -0.2966997
#> Rbp4                            -1.0499066
#>      TCTGGAACACACGCTGstochastic_pseudotime
#> Pyy                             -1.0289414
#> Sst                             -0.2966472
#> Rbp4                            -1.0498476
#>      CTGATCCTCAGTTGACstochastic_pseudotime
#> Pyy                             -1.0256100
#> Sst                             -0.2951187
#> Rbp4                            -1.0481422
#>      AGATTGCCACCAGGCTstochastic_pseudotime
#> Pyy                             -1.0240926
#> Sst                             -0.2944025
#> Rbp4                            -1.0473490
#>      ATCATGGTCCAGAGGAstochastic_pseudotime
#> Pyy                             -1.0239829
#> Sst                             -0.2943503
#> Rbp4                            -1.0472913
#>      AAGGCAGAGCTGAAATstochastic_pseudotime
#> Pyy                              -1.022325
#> Sst                              -0.293553
#> Rbp4                             -1.046412
#>      TTCCCAGGTATTCTCTstochastic_pseudotime
#> Pyy                             -1.0196549
#> Sst                             -0.2922406
#> Rbp4                            -1.0449705
#>      CCTAGCTCAGACTCGCstochastic_pseudotime
#> Pyy                             -1.0182919
#> Sst                             -0.2915575
#> Rbp4                            -1.0442224
#>      CGGGTCAAGGCTCAGAstochastic_pseudotime
#> Pyy                             -1.0180752
#> Sst                             -0.2914482
#> Rbp4                            -1.0441027
#>      TGGCCAGGTCGGCATCstochastic_pseudotime
#> Pyy                              -1.016532
#> Sst                              -0.290663
#> Rbp4                             -1.043244
#>      ACTGATGTCCCTAACCstochastic_pseudotime
#> Pyy                             -1.0164061
#> Sst                             -0.2905986
#> Rbp4                            -1.0431739
#>      CCGGTAGTCGGCTTGGstochastic_pseudotime
#> Pyy                              -1.016027
#> Sst                              -0.290404
#> Rbp4                             -1.042961
#>      TCTCTAAGTTACCGATstochastic_pseudotime
#> Pyy                             -1.0126881
#> Sst                             -0.2886644
#> Rbp4                            -1.0410611
#>      TGACGGCAGGAATGGAstochastic_pseudotime
#> Pyy                              -1.010039
#> Sst                              -0.287254
#> Rbp4                             -1.039519
#>      TAGCCGGTCAACACACstochastic_pseudotime
#> Pyy                             -1.0093448
#> Sst                             -0.2868805
#> Rbp4                            -1.0391095
#>      CGGACACAGTGTCTCAstochastic_pseudotime
#> Pyy                             -1.0082264
#> Sst                             -0.2862754
#> Rbp4                            -1.0384459
#>      TACGGATAGAGCTGGTstochastic_pseudotime
#> Pyy                             -1.0082252
#> Sst                             -0.2862747
#> Rbp4                            -1.0384451
#>      TCAGATGGTACACCGCstochastic_pseudotime
#> Pyy                             -1.0077939
#> Sst                             -0.2860403
#> Rbp4                            -1.0381877
#>      TGGGAAGCATGACATCstochastic_pseudotime
#> Pyy                             -1.0076309
#> Sst                             -0.2859515
#> Rbp4                            -1.0380902
#>      CGATGGCGTTACGGAGstochastic_pseudotime
#> Pyy                             -1.0040747
#> Sst                             -0.2839968
#> Rbp4                            -1.0359347
#>      CTGAAACTCGCCGTGAstochastic_pseudotime
#> Pyy                              -1.001886
#> Sst                              -0.282778
#> Rbp4                             -1.034581
#>      CTACCCAGTGCGCTTGstochastic_pseudotime
#> Pyy                             -0.9980163
#> Sst                             -0.2806000
#> Rbp4                            -1.0321367
#>      CTTAGGACACAGGTTTstochastic_pseudotime
#> Pyy                             -0.9972944
#> Sst                             -0.2801911
#> Rbp4                            -1.0316736
#>      ACGTCAACACTGCCAGstochastic_pseudotime
#> Pyy                             -0.9964937
#> Sst                             -0.2797367
#> Rbp4                            -1.0311573
#>      AAAGATGAGTCATGCTstochastic_pseudotime
#> Pyy                             -0.9959397
#> Sst                             -0.2794219
#> Rbp4                            -1.0307985
#>      TGGTTCCGTCACCTAAstochastic_pseudotime
#> Pyy                             -0.9927033
#> Sst                             -0.2775774
#> Rbp4                            -1.0286761
#>      CGTCTACAGATCTGCTstochastic_pseudotime
#> Pyy                             -0.9908078
#> Sst                             -0.2764943
#> Rbp4                            -1.0274122
#>      CGAGCACAGCTGAAATstochastic_pseudotime
#> Pyy                             -0.9892526
#> Sst                             -0.2756052
#> Rbp4                            -1.0263638
#>      TTCTCAAAGCCGCCTAstochastic_pseudotime
#> Pyy                             -0.9892081
#> Sst                             -0.2755797
#> Rbp4                            -1.0263336
#>      GAACATCGTAATCGTCstochastic_pseudotime
#> Pyy                             -0.9864999
#> Sst                             -0.2740325
#> Rbp4                            -1.0244827
#>      GCGCAACTCACAACGTstochastic_pseudotime
#> Pyy                             -0.9839240
#> Sst                             -0.2725648
#> Rbp4                            -1.0226933
#>      AGCGGTCGTGTATGGGstochastic_pseudotime
#> Pyy                             -0.9780798
#> Sst                             -0.2692647
#> Rbp4                            -1.0185288
#>      AGGGAGTAGTGTCTCAstochastic_pseudotime
#> Pyy                             -0.9775112
#> Sst                             -0.2689467
#> Rbp4                            -1.0181159
#>      CGGAGCTCATTGGGCCstochastic_pseudotime
#> Pyy                             -0.9745236
#> Sst                             -0.2672880
#> Rbp4                            -1.0159239
#>      GCATGATCAAGCGAGTstochastic_pseudotime
#> Pyy                             -0.9734815
#> Sst                             -0.2667146
#> Rbp4                            -1.0151503
#>      GAGTCCGCAATTGCTGstochastic_pseudotime
#> Pyy                             -0.9714827
#> Sst                             -0.2656234
#> Rbp4                            -1.0136539
#>      CTAGCCTGTCTCTTATstochastic_pseudotime
#> Pyy                             -0.9698082
#> Sst                             -0.2647184
#> Rbp4                            -1.0123873
#>      GGAGCAAAGCTTCGCGstochastic_pseudotime
#> Pyy                             -0.9650985
#> Sst                             -0.2622251
#> Rbp4                            -1.0087614
#>      CCCATACCACATCCAAstochastic_pseudotime
#> Pyy                             -0.9604246
#> Sst                             -0.2598376
#> Rbp4                            -1.0050710
#>      CTACACCGTCGAATCTstochastic_pseudotime
#> Pyy                             -0.9581147
#> Sst                             -0.2586937
#> Rbp4                            -1.0032133
#>      CATATGGCATGTTGACstochastic_pseudotime
#> Pyy                             -0.9483967
#> Sst                             -0.2541801
#> Rbp4                            -0.9951539
#>      CACACAACAGCCAGAAstochastic_pseudotime
#> Pyy                             -0.9462284
#> Sst                             -0.2532450
#> Rbp4                            -0.9933021
#>      GGAATAAGTTCGGGCTstochastic_pseudotime
#> Pyy                             -0.9441784
#> Sst                             -0.2523869
#> Rbp4                            -0.9915335
#>      GGGCATCCATGACATCstochastic_pseudotime
#> Pyy                             -0.9433353
#> Sst                             -0.2520414
#> Rbp4                            -0.9908011
#>      GTGTTAGCACCATCCTstochastic_pseudotime
#> Pyy                             -0.9388596
#> Sst                             -0.2502822
#> Rbp4                            -0.9868640
#>      CATCAGACATGCATGTstochastic_pseudotime
#> Pyy                             -0.9370533
#> Sst                             -0.2496087
#> Rbp4                            -0.9852518
#>      CCTCTGAGTATAAACGstochastic_pseudotime
#> Pyy                             -0.9342911
#> Sst                             -0.2486201
#> Rbp4                            -0.9827607
#>      ACAGCTATCAGCATGTstochastic_pseudotime
#> Pyy                             -0.9341881
#> Sst                             -0.2485843
#> Rbp4                            -0.9826672
#>      CAACCAATCGGGAGTAstochastic_pseudotime
#> Pyy                             -0.9305303
#> Sst                             -0.2473562
#> Rbp4                            -0.9793191
#>      GTGCATAAGCTCCTCTstochastic_pseudotime
#> Pyy                             -0.9287763
#> Sst                             -0.2467995
#> Rbp4                            -0.9776944
#>      AGCTTGACAGTAAGATstochastic_pseudotime
#> Pyy                             -0.9191864
#> Sst                             -0.2441283
#> Rbp4                            -0.9685935
#>      CATCGAAGTACGAAATstochastic_pseudotime
#> Pyy                             -0.9175600
#> Sst                             -0.2437380
#> Rbp4                            -0.9670139
#>      ATCCACCGTACCATCAstochastic_pseudotime
#> Pyy                             -0.9126886
#> Sst                             -0.2426769
#> Rbp4                            -0.9622204
#>      CCGTGGATCCAAGCCGstochastic_pseudotime
#> Pyy                             -0.9120270
#> Sst                             -0.2425452
#> Rbp4                            -0.9615622
#>      CTGCCTAGTCTCCATCstochastic_pseudotime
#> Pyy                             -0.9097461
#> Sst                             -0.2421136
#> Rbp4                            -0.9592803
#>      GGAAAGCGTCTTGATGstochastic_pseudotime
#> Pyy                             -0.9094649
#> Sst                             -0.2420628
#> Rbp4                            -0.9589976
#>      AGCTTGAGTGCAGTAGstochastic_pseudotime
#> Pyy                             -0.9069477
#> Sst                             -0.2416314
#> Rbp4                            -0.9564533
#>      CCTACACTCCTTTACAstochastic_pseudotime
#> Pyy                             -0.9008679
#> Sst                             -0.2407592
#> Rbp4                            -0.9502090
#>      CTTCTCTAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.8947241
#> Sst                             -0.2401129
#> Rbp4                            -0.9437590
#>      TCTATTGCACAGGAGTstochastic_pseudotime
#> Pyy                             -0.8883876
#> Sst                             -0.2396800
#> Rbp4                            -0.9369631
#>      CGTTAGAAGTATCTCGstochastic_pseudotime
#> Pyy                             -0.8742600
#> Sst                             -0.2394815
#> Rbp4                            -0.9213087
#>      CACACTCTCATGTGGTstochastic_pseudotime
#> Pyy                             -0.8741537
#> Sst                             -0.2394836
#> Rbp4                            -0.9211885
#>      ATTATCCAGTGAACATstochastic_pseudotime
#> Pyy                             -0.8690564
#> Sst                             -0.2396427
#> Rbp4                            -0.9153765
#>      AGAGCGATCTCTTGATstochastic_pseudotime
#> Pyy                             -0.8579628
#> Sst                             -0.2403344
#> Rbp4                            -0.9024515
#>      AGAATAGCAGTACACTstochastic_pseudotime
#> Pyy                             -0.8541543
#> Sst                             -0.2406656
#> Rbp4                            -0.8979318
#>      CCTAGCTGTGGCCCTAstochastic_pseudotime
#> Pyy                             -0.8519849
#> Sst                             -0.2408728
#> Rbp4                            -0.8953394
#>      ACAGCCGAGTGGTCCCstochastic_pseudotime
#> Pyy                             -0.8518753
#> Sst                             -0.2408836
#> Rbp4                            -0.8952081
#>      GAGCAGAAGCGTGTCCstochastic_pseudotime
#> Pyy                             -0.8468839
#> Sst                             -0.2414074
#> Rbp4                            -0.8891949
#>      CGATGGCTCTGTGCAAstochastic_pseudotime
#> Pyy                             -0.8455443
#> Sst                             -0.2415579
#> Rbp4                            -0.8875701
#>      CAAGGCCTCTTACCTAstochastic_pseudotime
#> Pyy                             -0.8447845
#> Sst                             -0.2416449
#> Rbp4                            -0.8866466
#>      TTGACTTAGACAAGCCstochastic_pseudotime
#> Pyy                             -0.8435498
#> Sst                             -0.2417890
#> Rbp4                            -0.8851429
#>      AGGGAGTAGGGTATCGstochastic_pseudotime
#> Pyy                             -0.8333876
#> Sst                             -0.2430776
#> Rbp4                            -0.8726314
#>      GACTAACTCAGTCCCTstochastic_pseudotime
#> Pyy                             -0.8332765
#> Sst                             -0.2430926
#> Rbp4                            -0.8724933
#>      ACGTCAACAGACAAATstochastic_pseudotime
#> Pyy                             -0.8323775
#> Sst                             -0.2432141
#> Rbp4                            -0.8713754
#>      GTCAAGTGTGGTCCGTstochastic_pseudotime
#> Pyy                             -0.8318530
#> Sst                             -0.2432855
#> Rbp4                            -0.8707223
#>      TGGCGCATCACAACGTstochastic_pseudotime
#> Pyy                             -0.8314111
#> Sst                             -0.2433458
#> Rbp4                            -0.8701717
#>      CAGCGACCACAGGAGTstochastic_pseudotime
#> Pyy                             -0.8268568
#> Sst                             -0.2439802
#> Rbp4                            -0.8644744
#>      AGAGCTTGTTAAAGACstochastic_pseudotime
#> Pyy                             -0.8266846
#> Sst                             -0.2440045
#> Rbp4                            -0.8642582
#>      GCTGCTTGTCTTCAAGstochastic_pseudotime
#> Pyy                             -0.8252777
#> Sst                             -0.2442045
#> Rbp4                            -0.8624898
#>      CACATAGAGCGTGTCCstochastic_pseudotime
#> Pyy                             -0.8241590
#> Sst                             -0.2443647
#> Rbp4                            -0.8610811
#>      TGCTGCTGTCCGTCAGstochastic_pseudotime
#> Pyy                             -0.8235244
#> Sst                             -0.2444559
#> Rbp4                            -0.8602810
#>      CAGCTAAAGTCCATACstochastic_pseudotime
#> Pyy                             -0.8233254
#> Sst                             -0.2444846
#> Rbp4                            -0.8600299
#>      ACTGATGGTGAGGGTTstochastic_pseudotime
#> Pyy                             -0.8218993
#> Sst                             -0.2446907
#> Rbp4                            -0.8582288
#>      CACACTCAGATCCCATstochastic_pseudotime
#> Pyy                             -0.8210218
#> Sst                             -0.2448181
#> Rbp4                            -0.8571187
#>      CCTAGCTCATACTACGstochastic_pseudotime
#> Pyy                             -0.8210214
#> Sst                             -0.2448182
#> Rbp4                            -0.8571183
#>      CTAGTGATCAAAGACAstochastic_pseudotime
#> Pyy                             -0.8197084
#> Sst                             -0.2450096
#> Rbp4                            -0.8554551
#>      TGTATTCCAATGACCTstochastic_pseudotime
#> Pyy                             -0.8193066
#> Sst                             -0.2450683
#> Rbp4                            -0.8549456
#>      AAGACCTCAGGTTTCAstochastic_pseudotime
#> Pyy                             -0.8177976
#> Sst                             -0.2452894
#> Rbp4                            -0.8530297
#>      ACCGTAAGTTATCACGstochastic_pseudotime
#> Pyy                             -0.8161368
#> Sst                             -0.2455336
#> Rbp4                            -0.8509170
#>      ACGGGCTGTGCGATAGstochastic_pseudotime
#> Pyy                             -0.8145088
#> Sst                             -0.2457737
#> Rbp4                            -0.8488420
#>      TCTGAGAGTCAAACTCstochastic_pseudotime
#> Pyy                             -0.8123269
#> Sst                             -0.2460961
#> Rbp4                            -0.8460553
#>      CTAGAGTCACTGTTAGstochastic_pseudotime
#> Pyy                             -0.8110060
#> Sst                             -0.2462914
#> Rbp4                            -0.8443651
#>      CTGCGGATCAGCATGTstochastic_pseudotime
#> Pyy                             -0.8048406
#> Sst                             -0.2472007
#> Rbp4                            -0.8364477
#>      CAGCCGAAGCGATATAstochastic_pseudotime
#> Pyy                             -0.8045266
#> Sst                             -0.2472468
#> Rbp4                            -0.8360434
#>      CAGCAGCAGACAAAGGstochastic_pseudotime
#> Pyy                             -0.8044022
#> Sst                             -0.2472650
#> Rbp4                            -0.8358831
#>      TGTGGTACAAGACGTGstochastic_pseudotime
#> Pyy                             -0.8036366
#> Sst                             -0.2473773
#> Rbp4                            -0.8348967
#>      TCGTACCCACGCATCGstochastic_pseudotime
#> Pyy                             -0.8026449
#> Sst                             -0.2475223
#> Rbp4                            -0.8336179
#>      GATGCTAGTAGAAGGAstochastic_pseudotime
#> Pyy                             -0.8019451
#> Sst                             -0.2476244
#> Rbp4                            -0.8327150
#>      CCTAGCTTCGATGAGGstochastic_pseudotime
#> Pyy                             -0.8013416
#> Sst                             -0.2477123
#> Rbp4                            -0.8319360
#>      GTAACTGTCAGTGCATstochastic_pseudotime
#> Pyy                             -0.8005996
#> Sst                             -0.2478202
#> Rbp4                            -0.8309778
#>      GGACGTCGTCTAGGTTstochastic_pseudotime
#> Pyy                             -0.7973048
#> Sst                             -0.2482956
#> Rbp4                            -0.8267171
#>      ATAAGAGCATCGATTGstochastic_pseudotime
#> Pyy                             -0.7954638
#> Sst                             -0.2485584
#> Rbp4                            -0.8243328
#>      CGGTTAAAGAGCTGCAstochastic_pseudotime
#> Pyy                             -0.7930495
#> Sst                             -0.2489000
#> Rbp4                            -0.8212020
#>      CGCCAAGCACCGGAAAstochastic_pseudotime
#> Pyy                             -0.7915466
#> Sst                             -0.2491109
#> Rbp4                            -0.8192509
#>      AGGCCGTGTACCAGTTstochastic_pseudotime
#> Pyy                             -0.7885615
#> Sst                             -0.2495257
#> Rbp4                            -0.8153706
#>      CGGTTAAGTAAACACAstochastic_pseudotime
#> Pyy                             -0.7880925
#> Sst                             -0.2495905
#> Rbp4                            -0.8147604
#>      CTCTGGTAGTCCCACGstochastic_pseudotime
#> Pyy                             -0.7878810
#> Sst                             -0.2496196
#> Rbp4                            -0.8144851
#>      ACAGCTATCTGGGCCAstochastic_pseudotime
#> Pyy                             -0.7870487
#> Sst                             -0.2497340
#> Rbp4                            -0.8134016
#>      AACTCAGTCTGGTTCCstochastic_pseudotime
#> Pyy                             -0.7844071
#> Sst                             -0.2500945
#> Rbp4                            -0.8099593
#>      TGAGAGGGTAAGGGCTstochastic_pseudotime
#> Pyy                             -0.7821988
#> Sst                             -0.2503929
#> Rbp4                            -0.8070777
#>      CCCAGTTAGCATCATCstochastic_pseudotime
#> Pyy                             -0.7809120
#> Sst                             -0.2505655
#> Rbp4                            -0.8053970
#>      TGCGCAGTCTGCCAGGstochastic_pseudotime
#> Pyy                             -0.7783306
#> Sst                             -0.2509090
#> Rbp4                            -0.8020216
#>      CGTGTAATCCATGCTCstochastic_pseudotime
#> Pyy                             -0.7774086
#> Sst                             -0.2510308
#> Rbp4                            -0.8008149
#>      TGACTAGTCCACGTGGstochastic_pseudotime
#> Pyy                             -0.7764940
#> Sst                             -0.2511511
#> Rbp4                            -0.7996173
#>      TAGAGCTGTATGAAACstochastic_pseudotime
#> Pyy                             -0.7736053
#> Sst                             -0.2515283
#> Rbp4                            -0.7958306
#>      TAGAGCTGTTGTCGCGstochastic_pseudotime
#> Pyy                             -0.7692968
#> Sst                             -0.2520827
#> Rbp4                            -0.7901717
#>      TTTATGCTCCATTCTAstochastic_pseudotime
#> Pyy                             -0.7680181
#> Sst                             -0.2522454
#> Rbp4                            -0.7884896
#>      AAGCCGCAGGCAGGTTstochastic_pseudotime
#> Pyy                             -0.7659111
#> Sst                             -0.2525117
#> Rbp4                            -0.7857155
#>      TTTATGCGTCGAACAGstochastic_pseudotime
#> Pyy                             -0.7649676
#> Sst                             -0.2526302
#> Rbp4                            -0.7844722
#>      CGAACATAGCCAGTAGstochastic_pseudotime
#> Pyy                             -0.7647205
#> Sst                             -0.2526612
#> Rbp4                            -0.7841466
#>      GTCGGGTCAGGCTGAAstochastic_pseudotime
#> Pyy                             -0.7611356
#> Sst                             -0.2531069
#> Rbp4                            -0.7794165
#>      CACATTTTCGAGGTAGstochastic_pseudotime
#> Pyy                             -0.7590886
#> Sst                             -0.2533586
#> Rbp4                            -0.7767115
#>      GGACAGATCCTGCAGGstochastic_pseudotime
#> Pyy                             -0.7574028
#> Sst                             -0.2535644
#> Rbp4                            -0.7744816
#>      GATCGTAGTCCTAGCGstochastic_pseudotime
#> Pyy                             -0.7568716
#> Sst                             -0.2536290
#> Rbp4                            -0.7737786
#>      TGACTTTAGCGTTGCCstochastic_pseudotime
#> Pyy                             -0.7552666
#> Sst                             -0.2538232
#> Rbp4                            -0.7716533
#>      CTCATTACAAGCGCTCstochastic_pseudotime
#> Pyy                             -0.7547187
#> Sst                             -0.2538892
#> Rbp4                            -0.7709273
#>      TACGGGCTCTTCGGTCstochastic_pseudotime
#> Pyy                             -0.7533256
#> Sst                             -0.2540565
#> Rbp4                            -0.7690806
#>      TGCGCAGGTACGCACCstochastic_pseudotime
#> Pyy                             -0.7509395
#> Sst                             -0.2543408
#> Rbp4                            -0.7659144
#>      AACTGGTGTGCCTGTGstochastic_pseudotime
#> Pyy                             -0.7466280
#> Sst                             -0.2548480
#> Rbp4                            -0.7601836
#>      CTACACCGTCGCCATGstochastic_pseudotime
#> Pyy                             -0.7422370
#> Sst                             -0.2553559
#> Rbp4                            -0.7543341
#>      TAGGCATTCCAAAGTCstochastic_pseudotime
#> Pyy                             -0.7422302
#> Sst                             -0.2553567
#> Rbp4                            -0.7543250
#>      CAACCAACAACCGCCAstochastic_pseudotime
#> Pyy                             -0.7362151
#> Sst                             -0.2560388
#> Rbp4                            -0.7462910
#>      CCATTCGGTAGGGACTstochastic_pseudotime
#> Pyy                             -0.7261792
#> Sst                             -0.2571431
#> Rbp4                            -0.7328336
#>      TATGCCCTCTTACCTAstochastic_pseudotime
#> Pyy                             -0.7244128
#> Sst                             -0.2573332
#> Rbp4                            -0.7304582
#>      TTAACTCCAAAGAATCstochastic_pseudotime
#> Pyy                             -0.7239789
#> Sst                             -0.2573797
#> Rbp4                            -0.7298745
#>      GATGAAAAGTTGTAGAstochastic_pseudotime
#> Pyy                             -0.7235611
#> Sst                             -0.2574244
#> Rbp4                            -0.7293122
#>      CCGTACTTCAATAAGGstochastic_pseudotime
#> Pyy                             -0.7233584
#> Sst                             -0.2574461
#> Rbp4                            -0.7290394
#>      ATTACTCTCAGGCGAAstochastic_pseudotime
#> Pyy                             -0.7205002
#> Sst                             -0.2577500
#> Rbp4                            -0.7251896
#>      TGGCCAGAGCTACCTAstochastic_pseudotime
#> Pyy                             -0.7173274
#> Sst                             -0.2580836
#> Rbp4                            -0.7209099
#>      GAATGAACAGCCTTGGstochastic_pseudotime
#> Pyy                             -0.7155707
#> Sst                             -0.2582667
#> Rbp4                            -0.7185375
#>      TTCTACATCGCTTAGAstochastic_pseudotime
#> Pyy                             -0.7149824
#> Sst                             -0.2583277
#> Rbp4                            -0.7177427
#>      TAGTGGTCAACTGCGCstochastic_pseudotime
#> Pyy                             -0.7117786
#> Sst                             -0.2586579
#> Rbp4                            -0.7134099
#>      CTAGCCTGTTCGGGCTstochastic_pseudotime
#> Pyy                             -0.7115896
#> Sst                             -0.2586772
#> Rbp4                            -0.7131541
#>      GGACGTCGTTTAGCTGstochastic_pseudotime
#> Pyy                             -0.7104784
#> Sst                             -0.2587908
#> Rbp4                            -0.7116497
#>      AACTCAGGTTGGAGGTstochastic_pseudotime
#> Pyy                             -0.7102521
#> Sst                             -0.2588139
#> Rbp4                            -0.7113432
#>      TGAGCATTCGGGAGTAstochastic_pseudotime
#> Pyy                             -0.7085501
#> Sst                             -0.2589868
#> Rbp4                            -0.7090373
#>      AGTGAGGGTTACGGAGstochastic_pseudotime
#> Pyy                             -0.7053869
#> Sst                             -0.2593054
#> Rbp4                            -0.7047467
#>      CAGCAGCCACACTGCGstochastic_pseudotime
#> Pyy                             -0.7026262
#> Sst                             -0.2595805
#> Rbp4                            -0.7009971
#>      ATTACTCAGTTGCAGGstochastic_pseudotime
#> Pyy                             -0.7020886
#> Sst                             -0.2596338
#> Rbp4                            -0.7002664
#>      AGTCTTTCACCAGGTCstochastic_pseudotime
#> Pyy                             -0.7002528
#> Sst                             -0.2598149
#> Rbp4                            -0.6977696
#>      AGTGGGACAGGTCGTCstochastic_pseudotime
#> Pyy                             -0.6983939
#> Sst                             -0.2599971
#> Rbp4                            -0.6952395
#>      CAGCCGAAGTAGGTGCstochastic_pseudotime
#> Pyy                             -0.6977001
#> Sst                             -0.2600648
#> Rbp4                            -0.6942946
#>      TATCTCAGTTCGCTAAstochastic_pseudotime
#> Pyy                             -0.6949585
#> Sst                             -0.2603308
#> Rbp4                            -0.6905580
#>      TTCTACATCTGCAGTAstochastic_pseudotime
#> Pyy                             -0.6917225
#> Sst                             -0.2606414
#> Rbp4                            -0.6861416
#>      ATCATCTTCTGTCAAGstochastic_pseudotime
#> Pyy                             -0.6893941
#> Sst                             -0.2608628
#> Rbp4                            -0.6829600
#>      TGGACGCTCAGGTAAAstochastic_pseudotime
#> Pyy                             -0.6888207
#> Sst                             -0.2609171
#> Rbp4                            -0.6821759
#>      TCTCTAAAGCTGGAACstochastic_pseudotime
#> Pyy                             -0.6869043
#> Sst                             -0.2610976
#> Rbp4                            -0.6795543
#>      TGAGAGGAGATCCCATstochastic_pseudotime
#> Pyy                             -0.6815986
#> Sst                             -0.2615913
#> Rbp4                            -0.6722844
#>      GTATTCTGTCCGTCAGstochastic_pseudotime
#> Pyy                             -0.6814586
#> Sst                             -0.2616042
#> Rbp4                            -0.6720923
#>      CCATTCGCATGCCTAAstochastic_pseudotime
#> Pyy                             -0.6724577
#> Sst                             -0.2624212
#> Rbp4                            -0.6597204
#>      AGATTGCAGCTGTCTAstochastic_pseudotime
#> Pyy                             -0.6708611
#> Sst                             -0.2625635
#> Rbp4                            -0.6575209
#>      GTATTCTCACTTGGATstochastic_pseudotime
#> Pyy                             -0.6682525
#> Sst                             -0.2627945
#> Rbp4                            -0.6539240
#>      CGACCTTAGGCTAGACstochastic_pseudotime
#> Pyy                             -0.6618149
#> Sst                             -0.2633561
#> Rbp4                            -0.6450310
#>      CTTGGCTAGACAGACCstochastic_pseudotime
#> Pyy                             -0.6526412
#> Sst                             -0.2641360
#> Rbp4                            -0.6323177
#>      TCAATCTGTCAAAGATstochastic_pseudotime
#> Pyy                             -0.6434849
#> Sst                             -0.2648918
#> Rbp4                            -0.6195814
#>      CTTGGCTGTTAGAACAstochastic_pseudotime
#> Pyy                             -0.6186523
#> Sst                             -0.2668350
#> Rbp4                            -0.5848108
#>      TACTCATCAGATCCATstochastic_pseudotime
#> Pyy                             -0.5129445
#> Sst                             -0.2736898
#> Rbp4                            -0.4333651
#>      TCGCGAGAGTCCGTATstochastic_pseudotime
#> Pyy                             -0.3915534
#> Sst                             -0.2795500
#> Rbp4                            -0.2535790
#>      ACGAGGATCATCACCCstochastic_pseudotime
#> Pyy                             -0.3077774
#> Sst                             -0.2827476
#> Rbp4                            -0.1264929
#>      CCCATACAGACTTTCGstochastic_pseudotime
#> Pyy                           -0.224804162
#> Sst                           -0.285421563
#> Rbp4                           0.001381008
#>      AGGCCGTCACGGTGTCstochastic_pseudotime
#> Pyy                            -0.17846136
#> Sst                            -0.28674322
#> Rbp4                            0.07355036
#>      CCGTGGACATAGAAACstochastic_pseudotime
#> Pyy                             -0.1604011
#> Sst                             -0.2872294
#> Rbp4                             0.1018054
#>      GCATGTAAGAATCTCCstochastic_pseudotime
#> Pyy                             -0.1288978
#> Sst                             -0.2880419
#> Rbp4                             0.1512536
#>      CACTCCATCAACACACstochastic_pseudotime
#> Pyy                            -0.05563876
#> Sst                            -0.28977424
#> Rbp4                            0.26695800
#>      GTGCTTCTCACGACTAstochastic_pseudotime
#> Pyy                             -0.0480042
#> Sst                             -0.2899433
#> Rbp4                             0.2790678
#>      GGCTCGACAGCTTCGGstochastic_pseudotime
#> Pyy                              0.2054249
#> Sst                             -0.2946064
#> Rbp4                             0.6848581
#>      GCTCCTAAGGGCATGTstochastic_pseudotime
#> Pyy                              0.3980365
#> Sst                             -0.2972564
#> Rbp4                             0.9951150
#>      GTCTCGTTCTGCGTAAstochastic_pseudotime
#> Pyy                              0.4054258
#> Sst                             -0.2973469
#> Rbp4                             1.0069941
#>      GTGGGTCCATCGATTGstochastic_pseudotime
#> Pyy                              0.6351619
#> Sst                             -0.2998484
#> Rbp4                             1.3732462
#>      CCACGGAGTCTAGTCAstochastic_pseudotime
#> Pyy                              0.7961802
#> Sst                             -0.3013168
#> Rbp4                             1.6234768
#>      CTAACTTAGATCACGGstochastic_pseudotime
#> Pyy                              0.9402633
#> Sst                             -0.3024902
#> Rbp4                             1.8392243
#>      TGTCCCAGTCCGTGACstochastic_pseudotime
#> Pyy                              1.0671630
#> Sst                             -0.3034454
#> Rbp4                             2.0192069
#>      ATCCACCAGAGCCTAGstochastic_pseudotime
#> Pyy                              1.1403031
#> Sst                             -0.3039739
#> Rbp4                             2.1166767
#>      GTAGGCCGTTACCAGTstochastic_pseudotime
#> Pyy                              1.2395215
#> Sst                             -0.3046805
#> Rbp4                             2.2379978
#>      TCTGGAACACCTATCCstochastic_pseudotime
#> Pyy                              1.2487751
#> Sst                             -0.3047465
#> Rbp4                             2.2484656
#>      ACGTCAAAGGTGTTAAstochastic_pseudotime
#> Pyy                              1.3451852
#> Sst                             -0.3054519
#> Rbp4                             2.3445397
#>      ACTTGTTCATGGTCATstochastic_pseudotime
#> Pyy                               1.398016
#> Sst                              -0.305871
#> Rbp4                              2.381460
#>      GGACAGATCAGGCGAAstochastic_pseudotime
#> Pyy                              1.4362303
#> Sst                             -0.3062132
#> Rbp4                             2.3941317
#>      TACTTGTAGACAAAGGstochastic_pseudotime
#> Pyy                              1.4391259
#> Sst                             -0.3062416
#> Rbp4                             2.3942864
#>      AAATGCCAGACTAAGTstochastic_pseudotime
#> Pyy                              1.4769401
#> Sst                             -0.3066914
#> Rbp4                             2.3728625
#>      AAGTCTGGTTGGGACAstochastic_pseudotime
#> Pyy                              1.4906069
#> Sst                             -0.3070019
#> Rbp4                             2.3238970
#>      AGACGTTCAATGGTCTstochastic_pseudotime
#> Pyy                              1.4918405
#> Sst                             -0.3070879
#> Rbp4                             2.3037895
#>      CGTTCTGTCGATAGAAstochastic_pseudotime
#> Pyy                              1.4893133
#> Sst                             -0.3072824
#> Rbp4                             2.2453549
#>      CTACACCCATCAGTCAstochastic_pseudotime
#> Pyy                              1.4830522
#> Sst                             -0.3074061
#> Rbp4                             2.1972412
#>      GACCAATAGAATGTGTstochastic_pseudotime
#> Pyy                              1.4638998
#> Sst                             -0.3075993
#> Rbp4                             2.1010478
#>      TGACTTTAGCCAGAACstochastic_pseudotime
#> Pyy                              1.4533172
#> Sst                             -0.3076704
#> Rbp4                             2.0580117
#>      TCACAAGCAGACTCGCstochastic_pseudotime
#> Pyy                              1.4108137
#> Sst                             -0.3078676
#> Rbp4                             1.9114212
#>      CTTTGCGTCCACGTGGstochastic_pseudotime
#> Pyy                              1.3370687
#> Sst                             -0.3080747
#> Rbp4                             1.7004612
#>      CTCTGGTCAAGACACGstochastic_pseudotime
#> Pyy                              1.2270487
#> Sst                             -0.3082628
#> Rbp4                             1.4304223
#>      AGCATACGTCCGCTGAstochastic_pseudotime
#> Pyy                              0.9742847
#> Sst                             -0.3084851
#> Rbp4                             0.9074140
#>      GCTGCGATCGGCTACGstochastic_pseudotime
#> Pyy                              0.9580318
#> Sst                             -0.3084938
#> Rbp4                             0.8769259
#>      CCACCTACAGCTCGACstochastic_pseudotime
#> Pyy                              0.8850385
#> Sst                             -0.3085272
#> Rbp4                             0.7437361
#>      GGCAATTGTAGAGCTGstochastic_pseudotime
#> Pyy                              0.7958294
#> Sst                             -0.3085570
#> Rbp4                             0.5886183
#>      ACTTACTCACATGTGTstochastic_pseudotime
#> Pyy                              0.7361342
#> Sst                             -0.3085711
#> Rbp4                             0.4891395
#>      ACGAGCCCATGAGCGAstochastic_pseudotime
#> Pyy                              0.7140539
#> Sst                             -0.3085752
#> Rbp4                             0.4531732
#>      CCCAGTTGTGCGCTTGstochastic_pseudotime
#> Pyy                              0.7051303
#> Sst                             -0.3085767
#> Rbp4                             0.4387614
#>      GACGTTACATCTCGCTstochastic_pseudotime
#> Pyy                              0.4297936
#> Sst                             -0.3085835
#> Rbp4                             0.0266653
#>      AAACGGGGTACAGTTCstochastic_pseudotime
#> Pyy                              0.1209476
#> Sst                             -0.3085153
#> Rbp4                            -0.3676197
#>      GCGCGATGTCCGAATTstochastic_pseudotime
#> Pyy                              0.1119209
#> Sst                             -0.3085123
#> Rbp4                            -0.3781364
#>      TGGCCAGCAAGGACTGstochastic_pseudotime
#> Pyy                             -0.0423344
#> Sst                             -0.3084537
#> Rbp4                            -0.5492933
#>      TAGTTGGTCTTGCAAGstochastic_pseudotime
#> Pyy                            -0.04689558
#> Sst                            -0.30845180
#> Rbp4                           -0.55410979
#>      TGTGGTACACATCCGGstochastic_pseudotime
#> Pyy                            -0.04885095
#> Sst                            -0.30845098
#> Rbp4                           -0.55617036
#> attr(,"scaled:center")
#>       Pyy       Sst      Rbp4 
#> 62.618456  3.767695 19.926528 
#> attr(,"scaled:scale")
#>      Pyy      Sst     Rbp4 
#> 50.72619 11.98590 15.72052 
#> 
#> $cell_order
#> $cell_order$stochastic_pseudotime
#>    [1] "GCATGCGAGGAACTGCstochastic_pseudotime"
#>    [2] "CATCGGGCACGTTGGCstochastic_pseudotime"
#>    [3] "CTCGGAGGTCATTAGCstochastic_pseudotime"
#>    [4] "GGACATTGTCTCCCTAstochastic_pseudotime"
#>    [5] "AGCTCCTTCGCCCTTAstochastic_pseudotime"
#>    [6] "CAGGTGCCATCAGTCAstochastic_pseudotime"
#>    [7] "AGCTCCTCAGTTTACGstochastic_pseudotime"
#>    [8] "GACGGCTTCAACTCTTstochastic_pseudotime"
#>    [9] "GTCTTCGTCTACCAGAstochastic_pseudotime"
#>   [10] "ACGAGGATCGCCGTGAstochastic_pseudotime"
#>   [11] "TCCCGATAGTCGTTTGstochastic_pseudotime"
#>   [12] "GTTCATTAGGCAGTCAstochastic_pseudotime"
#>   [13] "TCATTTGCATTGGCGCstochastic_pseudotime"
#>   [14] "ACATCAGAGAGACTTAstochastic_pseudotime"
#>   [15] "GTACTCCCACAGGTTTstochastic_pseudotime"
#>   [16] "CAGCTGGGTCTAACGTstochastic_pseudotime"
#>   [17] "GCGCCAATCTTATCTGstochastic_pseudotime"
#>   [18] "AATCCAGTCGGTCCGAstochastic_pseudotime"
#>   [19] "AATCCAGTCAGGTTCAstochastic_pseudotime"
#>   [20] "CCTACCAGTTGCTCCTstochastic_pseudotime"
#>   [21] "CAAGTTGTCTAACTGGstochastic_pseudotime"
#>   [22] "GAAGCAGGTCTAGAGGstochastic_pseudotime"
#>   [23] "CTGCTGTAGACTGGGTstochastic_pseudotime"
#>   [24] "CTTAACTTCAGCTTAGstochastic_pseudotime"
#>   [25] "CTAATGGCATGGTTGTstochastic_pseudotime"
#>   [26] "CATCAGATCATCGCTCstochastic_pseudotime"
#>   [27] "AGGTCATGTTTAGGAAstochastic_pseudotime"
#>   [28] "AATCCAGAGGGTTCCCstochastic_pseudotime"
#>   [29] "TCTCATAAGAATGTGTstochastic_pseudotime"
#>   [30] "GTCACGGCATCCGTGGstochastic_pseudotime"
#>   [31] "AGATCTGTCCGTCAAAstochastic_pseudotime"
#>   [32] "TAGTTGGCAGTAACGGstochastic_pseudotime"
#>   [33] "CGACTTCCAGCTGCACstochastic_pseudotime"
#>   [34] "TTAACTCGTAATCGTCstochastic_pseudotime"
#>   [35] "GTAACTGGTCACTGGCstochastic_pseudotime"
#>   [36] "AAAGCAAAGCTACCTAstochastic_pseudotime"
#>   [37] "CTTGGCTCAAAGGAAGstochastic_pseudotime"
#>   [38] "AGGGTGAGTGCGAAACstochastic_pseudotime"
#>   [39] "CATTATCGTCATATCGstochastic_pseudotime"
#>   [40] "TACTCATAGGTGTTAAstochastic_pseudotime"
#>   [41] "TTAGGCAGTTTGTTGGstochastic_pseudotime"
#>   [42] "TCAGCAATCGCCGTGAstochastic_pseudotime"
#>   [43] "CGTCACTTCACTTATCstochastic_pseudotime"
#>   [44] "CGTTGGGCAGGCGATAstochastic_pseudotime"
#>   [45] "AACACGTGTCGACTGCstochastic_pseudotime"
#>   [46] "CCTTCCCTCGTCGTTCstochastic_pseudotime"
#>   [47] "CTACATTGTTACTGACstochastic_pseudotime"
#>   [48] "GAAGCAGGTCGCGAAAstochastic_pseudotime"
#>   [49] "GCAGTTACATTTCACTstochastic_pseudotime"
#>   [50] "GATCGCGCACCAGGCTstochastic_pseudotime"
#>   [51] "TATGCCCAGAACAATCstochastic_pseudotime"
#>   [52] "CAGCTGGCAAGTTGTCstochastic_pseudotime"
#>   [53] "CAACCAAAGCGTTCCGstochastic_pseudotime"
#>   [54] "CTACATTCACCGTTGGstochastic_pseudotime"
#>   [55] "TTGACTTTCACCGGGTstochastic_pseudotime"
#>   [56] "CCTTACGCACAACGCCstochastic_pseudotime"
#>   [57] "GTATCTTTCGTGGACCstochastic_pseudotime"
#>   [58] "CGAGCACAGGCAGGTTstochastic_pseudotime"
#>   [59] "CGGTTAATCACCGGGTstochastic_pseudotime"
#>   [60] "GGGCACTGTCTCCATCstochastic_pseudotime"
#>   [61] "ACATACGGTCTCTCTGstochastic_pseudotime"
#>   [62] "CCCATACTCGTGGACCstochastic_pseudotime"
#>   [63] "TTCCCAGTCATAGCACstochastic_pseudotime"
#>   [64] "TGGCCAGCAATGAAACstochastic_pseudotime"
#>   [65] "CTCGGAGGTGAAGGCTstochastic_pseudotime"
#>   [66] "GACGTTAGTGTGAATAstochastic_pseudotime"
#>   [67] "CGTAGCGTCGGAAATAstochastic_pseudotime"
#>   [68] "GGAGCAAGTATCACCAstochastic_pseudotime"
#>   [69] "CCTCAGTAGTGCGATGstochastic_pseudotime"
#>   [70] "TAGTGGTTCCTCAACCstochastic_pseudotime"
#>   [71] "TGGGAAGCACGTAAGGstochastic_pseudotime"
#>   [72] "AGAGTGGCACCTCGGAstochastic_pseudotime"
#>   [73] "TAGTGGTCATTCACTTstochastic_pseudotime"
#>   [74] "CAAGGCCTCTCAACTTstochastic_pseudotime"
#>   [75] "TTTCCTCGTGGCGAATstochastic_pseudotime"
#>   [76] "CTTTGCGAGACCACGAstochastic_pseudotime"
#>   [77] "ACTGATGCATGTCTCCstochastic_pseudotime"
#>   [78] "CCACTACAGAGCTGGTstochastic_pseudotime"
#>   [79] "CATATGGTCACAGTACstochastic_pseudotime"
#>   [80] "CTGGTCTAGGATTCGGstochastic_pseudotime"
#>   [81] "ATCCGAAGTACACCGCstochastic_pseudotime"
#>   [82] "GTTACAGTCCGATATGstochastic_pseudotime"
#>   [83] "GAGGTGAGTCCGCTGAstochastic_pseudotime"
#>   [84] "AAGGTTCGTCAACTGTstochastic_pseudotime"
#>   [85] "CAGCTAACATGCCTAAstochastic_pseudotime"
#>   [86] "TGCCCTAAGTCCGGTCstochastic_pseudotime"
#>   [87] "GGACGTCGTCAGAAGCstochastic_pseudotime"
#>   [88] "TTCGGTCTCGACCAGCstochastic_pseudotime"
#>   [89] "ACACCCTTCCCAGGTGstochastic_pseudotime"
#>   [90] "CACCACTGTGTGTGCCstochastic_pseudotime"
#>   [91] "GCTGCGATCCCGGATGstochastic_pseudotime"
#>   [92] "ATGAGGGTCTTGACGAstochastic_pseudotime"
#>   [93] "AACCATGGTTGAACTCstochastic_pseudotime"
#>   [94] "AGTGGGAGTGGCAAACstochastic_pseudotime"
#>   [95] "CCAGCGACAGACGTAGstochastic_pseudotime"
#>   [96] "CCACGGAAGTCTCCTCstochastic_pseudotime"
#>   [97] "ACTGCTCAGGTGTTAAstochastic_pseudotime"
#>   [98] "CACACCTGTTATCCGAstochastic_pseudotime"
#>   [99] "GGCGTGTCACCGATATstochastic_pseudotime"
#>  [100] "TTCTACACATAAAGGTstochastic_pseudotime"
#>  [101] "AGTTGGTGTTGTCGCGstochastic_pseudotime"
#>  [102] "TAAGAGATCTGCGACGstochastic_pseudotime"
#>  [103] "GTGTTAGGTGAACCTTstochastic_pseudotime"
#>  [104] "CTGCGGACATTAGCCAstochastic_pseudotime"
#>  [105] "GTCACGGTCTTGCCGTstochastic_pseudotime"
#>  [106] "TGGACGCTCTTGTTTGstochastic_pseudotime"
#>  [107] "CTCCTAGCACAACGTTstochastic_pseudotime"
#>  [108] "TCAATCTCATCCGCGAstochastic_pseudotime"
#>  [109] "TCTCTAATCCCAACGGstochastic_pseudotime"
#>  [110] "GCTCTGTGTCCATGATstochastic_pseudotime"
#>  [111] "TACAGTGGTAAACGCGstochastic_pseudotime"
#>  [112] "TCGAGGCAGGCATGGTstochastic_pseudotime"
#>  [113] "GGCGTGTTCGAGAGCAstochastic_pseudotime"
#>  [114] "GGGTCTGTCATGTGGTstochastic_pseudotime"
#>  [115] "GATGCTAGTGTTAAGAstochastic_pseudotime"
#>  [116] "CTACGTCAGCTAGGCAstochastic_pseudotime"
#>  [117] "TCAGATGCACTCGACGstochastic_pseudotime"
#>  [118] "GTGCATATCAAGGCTTstochastic_pseudotime"
#>  [119] "CTCTGGTCAAGTTGTCstochastic_pseudotime"
#>  [120] "GAGCAGATCAACTCTTstochastic_pseudotime"
#>  [121] "CACAGGCCAAGGCTCCstochastic_pseudotime"
#>  [122] "TGCGTGGGTCCCTTGTstochastic_pseudotime"
#>  [123] "GGCTCGAAGATCACGGstochastic_pseudotime"
#>  [124] "TGTGTTTGTATGAATGstochastic_pseudotime"
#>  [125] "GATCGTAAGTAGGTGCstochastic_pseudotime"
#>  [126] "GATCGATAGAAAGTGGstochastic_pseudotime"
#>  [127] "GTGCATAAGAACAATCstochastic_pseudotime"
#>  [128] "TGGTTAGCAGTCGTGCstochastic_pseudotime"
#>  [129] "AAGGTTCAGATAGTCAstochastic_pseudotime"
#>  [130] "ACTGTCCCAGTTTACGstochastic_pseudotime"
#>  [131] "CACCAGGCATAAAGGTstochastic_pseudotime"
#>  [132] "CTGAAGTTCTGCGGCAstochastic_pseudotime"
#>  [133] "GCATGATGTCACTGGCstochastic_pseudotime"
#>  [134] "TCAACGATCAACCATGstochastic_pseudotime"
#>  [135] "ACAGCCGAGAGAACAGstochastic_pseudotime"
#>  [136] "GACGTGCGTCACCTAAstochastic_pseudotime"
#>  [137] "TCGGGACGTCTGCGGTstochastic_pseudotime"
#>  [138] "TGAGCCGTCCGCTGTTstochastic_pseudotime"
#>  [139] "TCGCGAGAGTACGCGAstochastic_pseudotime"
#>  [140] "GATTCAGAGCGATGACstochastic_pseudotime"
#>  [141] "TCTATTGTCGTATCAGstochastic_pseudotime"
#>  [142] "AGCTTGAAGACACGACstochastic_pseudotime"
#>  [143] "CGAATGTGTCTGGAGAstochastic_pseudotime"
#>  [144] "CGTTCTGAGGACGAAAstochastic_pseudotime"
#>  [145] "ACTATCTTCAAGATCCstochastic_pseudotime"
#>  [146] "ACACTGAGTCCATGATstochastic_pseudotime"
#>  [147] "GCAATCACAGCTGCTGstochastic_pseudotime"
#>  [148] "AAGCCGCCAATGGACGstochastic_pseudotime"
#>  [149] "GTAGGCCAGCTAACAAstochastic_pseudotime"
#>  [150] "GATCAGTAGGTGATATstochastic_pseudotime"
#>  [151] "TGCGTGGAGGTTCCTAstochastic_pseudotime"
#>  [152] "GGACAGACAATCCAACstochastic_pseudotime"
#>  [153] "CATGACATCATGCATGstochastic_pseudotime"
#>  [154] "CATGCCTCACCCATTCstochastic_pseudotime"
#>  [155] "TGCGGGTTCCACTCCAstochastic_pseudotime"
#>  [156] "TGTCCCACATTATCTCstochastic_pseudotime"
#>  [157] "CAAGTTGGTGTGCCTGstochastic_pseudotime"
#>  [158] "CTCGAAAGTGGTGTAGstochastic_pseudotime"
#>  [159] "TCGCGAGCACTTCTGCstochastic_pseudotime"
#>  [160] "CGAGCACGTAGTAGTAstochastic_pseudotime"
#>  [161] "TTGGAACAGAAGATTCstochastic_pseudotime"
#>  [162] "ACGATGTTCGCAGGCTstochastic_pseudotime"
#>  [163] "GGACGTCAGGGTGTTGstochastic_pseudotime"
#>  [164] "TGCGGGTAGACGCACAstochastic_pseudotime"
#>  [165] "CCATGTCTCCGCATAAstochastic_pseudotime"
#>  [166] "CAAGAAAAGAGTACCGstochastic_pseudotime"
#>  [167] "CCGTGGATCTTAGAGCstochastic_pseudotime"
#>  [168] "GGCGTGTCACCAGGTCstochastic_pseudotime"
#>  [169] "GAACCTAAGTAGATGTstochastic_pseudotime"
#>  [170] "CGAGAAGTCGTTACAGstochastic_pseudotime"
#>  [171] "AAGGAGCAGTACATGAstochastic_pseudotime"
#>  [172] "CGTCTACAGGGAGTAAstochastic_pseudotime"
#>  [173] "TTGCCGTCACCGAATTstochastic_pseudotime"
#>  [174] "CGATTGAAGTCAAGGCstochastic_pseudotime"
#>  [175] "AGCTCTCAGACATAACstochastic_pseudotime"
#>  [176] "AAATGCCCATGCTAGTstochastic_pseudotime"
#>  [177] "GACTGCGAGAATCTCCstochastic_pseudotime"
#>  [178] "ACGCCGAAGGGTGTTGstochastic_pseudotime"
#>  [179] "TTGTAGGTCAAAGTAGstochastic_pseudotime"
#>  [180] "CGATGTAAGTTTAGGAstochastic_pseudotime"
#>  [181] "GGAAAGCCAATGTAAGstochastic_pseudotime"
#>  [182] "AAACCTGGTAAGTGGCstochastic_pseudotime"
#>  [183] "AGGGAGTTCGGAGGTAstochastic_pseudotime"
#>  [184] "CTCGAGGGTAAGTAGTstochastic_pseudotime"
#>  [185] "CTTCTCTAGTCCGTATstochastic_pseudotime"
#>  [186] "GAATGAAAGAATGTGTstochastic_pseudotime"
#>  [187] "ATCATGGAGGAATGGAstochastic_pseudotime"
#>  [188] "GATTCAGTCGCCAGCAstochastic_pseudotime"
#>  [189] "CAGGTGCAGCCAGGATstochastic_pseudotime"
#>  [190] "ACGCAGCCATTCTCATstochastic_pseudotime"
#>  [191] "GCGCAGTGTGCTAGCCstochastic_pseudotime"
#>  [192] "GACCAATCAGACGCCTstochastic_pseudotime"
#>  [193] "CGGGTCAAGGCGATACstochastic_pseudotime"
#>  [194] "CTGCTGTCAAAGCAATstochastic_pseudotime"
#>  [195] "CGTGTCTCAGAAGCACstochastic_pseudotime"
#>  [196] "CCTACACAGTACGTTCstochastic_pseudotime"
#>  [197] "AACTCCCCAGTTCATGstochastic_pseudotime"
#>  [198] "AGCTCCTTCAATACCGstochastic_pseudotime"
#>  [199] "TGCGTGGGTTACTGACstochastic_pseudotime"
#>  [200] "AAGGTTCAGGAGCGAGstochastic_pseudotime"
#>  [201] "TGTCCCATCTCTGCTGstochastic_pseudotime"
#>  [202] "TGCTACCCATAAGACAstochastic_pseudotime"
#>  [203] "ACCAGTAGTTCGTTGAstochastic_pseudotime"
#>  [204] "CGATCGGAGCTGCAAGstochastic_pseudotime"
#>  [205] "CCACTACCAACGCACCstochastic_pseudotime"
#>  [206] "CATCAGAGTCTCTTTAstochastic_pseudotime"
#>  [207] "AGTCTTTCACGGACAAstochastic_pseudotime"
#>  [208] "TTTGTCATCTGTTTGTstochastic_pseudotime"
#>  [209] "AACTTTCTCCAGATCAstochastic_pseudotime"
#>  [210] "AAAGCAAGTGCAGTAGstochastic_pseudotime"
#>  [211] "CGTGTCTGTAGCGTCCstochastic_pseudotime"
#>  [212] "CAGTCCTCAGGGCATAstochastic_pseudotime"
#>  [213] "AGGGATGGTCGACTGCstochastic_pseudotime"
#>  [214] "CCACGGAAGATCCTGTstochastic_pseudotime"
#>  [215] "TTGAACGGTGAGGGAGstochastic_pseudotime"
#>  [216] "GTGCTTCAGAGGGATAstochastic_pseudotime"
#>  [217] "TCACAAGCAACGATGGstochastic_pseudotime"
#>  [218] "TTCGGTCGTATAGTAGstochastic_pseudotime"
#>  [219] "CCTTCGATCACAGGCCstochastic_pseudotime"
#>  [220] "AGGCCACCATGACATCstochastic_pseudotime"
#>  [221] "TTTGTCAAGAATGTGTstochastic_pseudotime"
#>  [222] "CCTCAGTCAGACACTTstochastic_pseudotime"
#>  [223] "CGCTTCATCTCTGTCGstochastic_pseudotime"
#>  [224] "GTGCATAGTTCAACCAstochastic_pseudotime"
#>  [225] "ATAGACCAGCGTGTCCstochastic_pseudotime"
#>  [226] "GGGAATGAGCACACAGstochastic_pseudotime"
#>  [227] "GTACTTTTCTTCGGTCstochastic_pseudotime"
#>  [228] "AGCAGCCGTCGCATCGstochastic_pseudotime"
#>  [229] "GTAGTCAGTCCGACGTstochastic_pseudotime"
#>  [230] "CTTCTCTGTCAGAGGTstochastic_pseudotime"
#>  [231] "GTGGGTCGTCCAGTATstochastic_pseudotime"
#>  [232] "GATGAAACAGGCTGAAstochastic_pseudotime"
#>  [233] "ATCTGCCCATGTTGACstochastic_pseudotime"
#>  [234] "TTAACTCCAGGATTGGstochastic_pseudotime"
#>  [235] "AGAATAGCAATCACACstochastic_pseudotime"
#>  [236] "GGGACCTCACTGTCGGstochastic_pseudotime"
#>  [237] "AGGGATGTCGCATGGCstochastic_pseudotime"
#>  [238] "GGACGTCGTCGAGATGstochastic_pseudotime"
#>  [239] "TCAATCTTCGAACGGAstochastic_pseudotime"
#>  [240] "CTCGTCAAGAACAATCstochastic_pseudotime"
#>  [241] "CACAGTAGTGCACTTAstochastic_pseudotime"
#>  [242] "AGTGTCATCCCTTGTGstochastic_pseudotime"
#>  [243] "GACCTGGGTTCGTGATstochastic_pseudotime"
#>  [244] "GCTGCGAAGAACAACTstochastic_pseudotime"
#>  [245] "CAGCCGATCCTATTCAstochastic_pseudotime"
#>  [246] "CACCACTAGCGATCCCstochastic_pseudotime"
#>  [247] "GCTTCCAAGTTAAGTGstochastic_pseudotime"
#>  [248] "ACACCGGTCATGTCTTstochastic_pseudotime"
#>  [249] "TGACGGCCATATGAGAstochastic_pseudotime"
#>  [250] "GCATGTAGTTGCGCACstochastic_pseudotime"
#>  [251] "CATATTCTCCTGTAGAstochastic_pseudotime"
#>  [252] "CGTTCTGTCACAATGCstochastic_pseudotime"
#>  [253] "AATCGGTTCGTATCAGstochastic_pseudotime"
#>  [254] "GAACGGAAGTCAAGCGstochastic_pseudotime"
#>  [255] "GGCGTGTTCAGCTTAGstochastic_pseudotime"
#>  [256] "CTAACTTAGAGGTACCstochastic_pseudotime"
#>  [257] "CTCTAATAGTTATCGCstochastic_pseudotime"
#>  [258] "TACACGAGTGCAACTTstochastic_pseudotime"
#>  [259] "CTCGTACTCTTTAGTCstochastic_pseudotime"
#>  [260] "CAAGGCCTCAGTGCATstochastic_pseudotime"
#>  [261] "AGGTCATCACTTCTGCstochastic_pseudotime"
#>  [262] "ATCACGAGTGCACCACstochastic_pseudotime"
#>  [263] "ACAGCTAAGTCGTTTGstochastic_pseudotime"
#>  [264] "GTGGGTCTCCGAATGTstochastic_pseudotime"
#>  [265] "CCGGGATAGATCTGAAstochastic_pseudotime"
#>  [266] "GGATGTTAGTAAGTACstochastic_pseudotime"
#>  [267] "TTCTTAGTCAGTGTTGstochastic_pseudotime"
#>  [268] "GCGACCATCACCAGGCstochastic_pseudotime"
#>  [269] "ACCAGTAAGTGATCGGstochastic_pseudotime"
#>  [270] "TAAGTGCTCCAAATGCstochastic_pseudotime"
#>  [271] "CGCGGTAAGAGCTGGTstochastic_pseudotime"
#>  [272] "CCGTTCAAGCTAAACAstochastic_pseudotime"
#>  [273] "ATAGACCAGTGATCGGstochastic_pseudotime"
#>  [274] "GACTAACGTCATCCCTstochastic_pseudotime"
#>  [275] "GGCAATTCACGGCTACstochastic_pseudotime"
#>  [276] "CCTACACTCTTTACACstochastic_pseudotime"
#>  [277] "GCATGTATCGTGGACCstochastic_pseudotime"
#>  [278] "CGCCAAGCAAAGTCAAstochastic_pseudotime"
#>  [279] "CGAATGTTCACCGGGTstochastic_pseudotime"
#>  [280] "ACGATACTCGCGCCAAstochastic_pseudotime"
#>  [281] "GGACAAGTCAAAGACAstochastic_pseudotime"
#>  [282] "AAGACCTAGTAGGCCAstochastic_pseudotime"
#>  [283] "CCACCTATCTAAGCCAstochastic_pseudotime"
#>  [284] "AGCATACCACAGCCCAstochastic_pseudotime"
#>  [285] "CCTACCACATCGATTGstochastic_pseudotime"
#>  [286] "ACGATACAGATGAGAGstochastic_pseudotime"
#>  [287] "AGCAGCCGTGTAACGGstochastic_pseudotime"
#>  [288] "CCGGGATTCGACAGCCstochastic_pseudotime"
#>  [289] "ACTTGTTAGTTAAGTGstochastic_pseudotime"
#>  [290] "TGAGAGGAGGTAAACTstochastic_pseudotime"
#>  [291] "AAGACCTAGACCTTTGstochastic_pseudotime"
#>  [292] "TGCGGGTAGGGTTCCCstochastic_pseudotime"
#>  [293] "GATCGATAGGATGGTCstochastic_pseudotime"
#>  [294] "CCCAATCCAACTGGCCstochastic_pseudotime"
#>  [295] "TGTGTTTAGGAATGGAstochastic_pseudotime"
#>  [296] "GATCGATGTTGTTTGGstochastic_pseudotime"
#>  [297] "GTCCTCACACCACGTGstochastic_pseudotime"
#>  [298] "TATCTCACACGCGAAAstochastic_pseudotime"
#>  [299] "CTGGTCTAGGTGATTAstochastic_pseudotime"
#>  [300] "TTCTCCTAGTCCATACstochastic_pseudotime"
#>  [301] "TCGAGGCAGAACAATCstochastic_pseudotime"
#>  [302] "CCCAGTTGTCGTGGCTstochastic_pseudotime"
#>  [303] "TCAGGTAAGGCGTACAstochastic_pseudotime"
#>  [304] "TACTCGCAGTAAGTACstochastic_pseudotime"
#>  [305] "GACCAATGTTATCCGAstochastic_pseudotime"
#>  [306] "TCTGGAAGTCATCGGCstochastic_pseudotime"
#>  [307] "GTACGTAGTCTCTTTAstochastic_pseudotime"
#>  [308] "CGTCACTTCTAACTCTstochastic_pseudotime"
#>  [309] "GGGCATCCATTTCAGGstochastic_pseudotime"
#>  [310] "GCATGATGTTTAGGAAstochastic_pseudotime"
#>  [311] "GTATCTTGTACTCAACstochastic_pseudotime"
#>  [312] "CTAACTTCATCTGGTAstochastic_pseudotime"
#>  [313] "GTCACGGAGCTATGCTstochastic_pseudotime"
#>  [314] "GTTCTCGAGAATTCCCstochastic_pseudotime"
#>  [315] "TCAGGTACAGTCACTAstochastic_pseudotime"
#>  [316] "GAATGAATCATCGGATstochastic_pseudotime"
#>  [317] "CGATGTAGTACTCTCCstochastic_pseudotime"
#>  [318] "ACACCCTAGGACAGCTstochastic_pseudotime"
#>  [319] "TGAGAGGAGCCAGGATstochastic_pseudotime"
#>  [320] "CGCGGTAAGCTGCAAGstochastic_pseudotime"
#>  [321] "GGCTGGTTCATGCATGstochastic_pseudotime"
#>  [322] "TGCGTGGAGTACACCTstochastic_pseudotime"
#>  [323] "TGAGAGGCAGGAACGTstochastic_pseudotime"
#>  [324] "GCACTCTCAGCGTCCAstochastic_pseudotime"
#>  [325] "GGGATGAAGGGCATGTstochastic_pseudotime"
#>  [326] "CCCAGTTTCAGCACATstochastic_pseudotime"
#>  [327] "TTGACTTGTTCCATGAstochastic_pseudotime"
#>  [328] "TGTGTTTGTAGATTAGstochastic_pseudotime"
#>  [329] "TGGGCGTGTGTGCGTCstochastic_pseudotime"
#>  [330] "CAGCCGATCGTTTATCstochastic_pseudotime"
#>  [331] "GGACGTCCATGGATGGstochastic_pseudotime"
#>  [332] "AGGTCATGTCTCCCTAstochastic_pseudotime"
#>  [333] "AGATTGCAGCGTGTCCstochastic_pseudotime"
#>  [334] "CACTCCACATCGTCGGstochastic_pseudotime"
#>  [335] "TTTGGTTAGACAGGCTstochastic_pseudotime"
#>  [336] "GCGCCAAAGCTCCTTCstochastic_pseudotime"
#>  [337] "CTGAAGTGTTTGGGCCstochastic_pseudotime"
#>  [338] "CGGAGTCAGTTTAGGAstochastic_pseudotime"
#>  [339] "CAGTCCTAGGTGGGTTstochastic_pseudotime"
#>  [340] "AGTTGGTAGAGCCTAGstochastic_pseudotime"
#>  [341] "CAGAATCCATTGGGCCstochastic_pseudotime"
#>  [342] "GGGTCTGTCCATGAACstochastic_pseudotime"
#>  [343] "GTGTGCGGTGTAAGTAstochastic_pseudotime"
#>  [344] "GCATGTAAGCGTGTCCstochastic_pseudotime"
#>  [345] "ACCTTTATCCTGCAGGstochastic_pseudotime"
#>  [346] "GGCTGGTAGTCACGCCstochastic_pseudotime"
#>  [347] "ACAGCTATCAAGGCTTstochastic_pseudotime"
#>  [348] "TCACGAAGTCCGAGTCstochastic_pseudotime"
#>  [349] "ACAGCCGAGAGGTACCstochastic_pseudotime"
#>  [350] "TATTACCTCGTTGCCTstochastic_pseudotime"
#>  [351] "GAACCTATCACTTATCstochastic_pseudotime"
#>  [352] "CGATGGCTCAGTGTTGstochastic_pseudotime"
#>  [353] "GTCATTTGTCGTCTTCstochastic_pseudotime"
#>  [354] "TTTATGCAGGTTACCTstochastic_pseudotime"
#>  [355] "CACACCTGTACATGTCstochastic_pseudotime"
#>  [356] "AAATGCCCAATGGAATstochastic_pseudotime"
#>  [357] "ATTACTCTCTGCCAGGstochastic_pseudotime"
#>  [358] "TGTGGTAGTGCCTGCAstochastic_pseudotime"
#>  [359] "GCGGGTTCAAGTAATGstochastic_pseudotime"
#>  [360] "TCTATTGTCCTATTCAstochastic_pseudotime"
#>  [361] "TTTGGTTTCCTTTCGGstochastic_pseudotime"
#>  [362] "GGGAATGCAAAGTGCGstochastic_pseudotime"
#>  [363] "GCGAGAAGTGGGTATGstochastic_pseudotime"
#>  [364] "TACACGACATGTAGTCstochastic_pseudotime"
#>  [365] "CCGTTCAAGGACCACAstochastic_pseudotime"
#>  [366] "AATCGGTTCCAAACTGstochastic_pseudotime"
#>  [367] "CTGAAACGTAGCGCAAstochastic_pseudotime"
#>  [368] "GTCTTCGTCCGTTGCTstochastic_pseudotime"
#>  [369] "CCACCTAAGACTTGAAstochastic_pseudotime"
#>  [370] "ACTTGTTTCACCCGAGstochastic_pseudotime"
#>  [371] "CGCCAAGTCTAAGCCAstochastic_pseudotime"
#>  [372] "CGTGAGCTCTTCTGGCstochastic_pseudotime"
#>  [373] "AAATGCCCACCAACCGstochastic_pseudotime"
#>  [374] "CGTCAGGGTAGGCATGstochastic_pseudotime"
#>  [375] "CGACTTCTCACGCGGTstochastic_pseudotime"
#>  [376] "CCACGGAAGCAATCTCstochastic_pseudotime"
#>  [377] "GCTGCGAAGTATTGGAstochastic_pseudotime"
#>  [378] "CTCGAGGCAAGTCTGTstochastic_pseudotime"
#>  [379] "AACTCAGTCCTCGCATstochastic_pseudotime"
#>  [380] "GGGATGAGTCTGCGGTstochastic_pseudotime"
#>  [381] "GATTCAGAGTACGTTCstochastic_pseudotime"
#>  [382] "AGTAGTCGTGAGGGAGstochastic_pseudotime"
#>  [383] "TTCTCAACAATAGAGTstochastic_pseudotime"
#>  [384] "GTAGTCACATCGATGTstochastic_pseudotime"
#>  [385] "TACAGTGCACAGGCCTstochastic_pseudotime"
#>  [386] "AGTTGGTTCACATAGCstochastic_pseudotime"
#>  [387] "ACGCCAGGTTACAGAAstochastic_pseudotime"
#>  [388] "TTTGTCATCGAATGCTstochastic_pseudotime"
#>  [389] "GCTTCCAAGACATAACstochastic_pseudotime"
#>  [390] "TCCACACCAGGCGATAstochastic_pseudotime"
#>  [391] "ACTGATGCATGGTCTAstochastic_pseudotime"
#>  [392] "AACTCTTCATTTGCCCstochastic_pseudotime"
#>  [393] "CTCAGAACATAGGATAstochastic_pseudotime"
#>  [394] "GGACAGATCCTTTCGGstochastic_pseudotime"
#>  [395] "TACGGGCAGCCAACAGstochastic_pseudotime"
#>  [396] "GTTAAGCAGGATTCGGstochastic_pseudotime"
#>  [397] "GTTCTCGCAAGAAAGGstochastic_pseudotime"
#>  [398] "CTCGAAATCCTGCCATstochastic_pseudotime"
#>  [399] "TTCGGTCAGCTAGCCCstochastic_pseudotime"
#>  [400] "AAATGCCGTGATAAGTstochastic_pseudotime"
#>  [401] "CGTCCATTCGCGCCAAstochastic_pseudotime"
#>  [402] "CAGCTAAGTACCGAGAstochastic_pseudotime"
#>  [403] "GCATGCGCAGCTCCGAstochastic_pseudotime"
#>  [404] "GCTTGAAGTCAACTGTstochastic_pseudotime"
#>  [405] "GGGAGATAGTGACATAstochastic_pseudotime"
#>  [406] "TAGTTGGCACGTAAGGstochastic_pseudotime"
#>  [407] "ACTGATGCATCAGTACstochastic_pseudotime"
#>  [408] "GTTTCTACAGGACCCTstochastic_pseudotime"
#>  [409] "AGAGTGGAGGCGATACstochastic_pseudotime"
#>  [410] "TTGCGTCCACCGATATstochastic_pseudotime"
#>  [411] "ACGCAGCGTAATCACCstochastic_pseudotime"
#>  [412] "CCATGTCGTCGGCTCAstochastic_pseudotime"
#>  [413] "CTGCGGACAGACTCGCstochastic_pseudotime"
#>  [414] "CACAAACAGTAGCGGTstochastic_pseudotime"
#>  [415] "CTCACACGTCACACGCstochastic_pseudotime"
#>  [416] "AACCATGGTTCGGCACstochastic_pseudotime"
#>  [417] "CAGAGAGTCGGCCGATstochastic_pseudotime"
#>  [418] "CAAGGCCCACATTTCTstochastic_pseudotime"
#>  [419] "GACGTGCTCCGGCACAstochastic_pseudotime"
#>  [420] "CCCATACCATGCAACTstochastic_pseudotime"
#>  [421] "CTAAGACGTATCAGTCstochastic_pseudotime"
#>  [422] "ATCATGGGTATTCGTGstochastic_pseudotime"
#>  [423] "AGAGCGACAGTAAGATstochastic_pseudotime"
#>  [424] "CTTCTCTCAGCATACTstochastic_pseudotime"
#>  [425] "CCACGGAAGGGTGTTGstochastic_pseudotime"
#>  [426] "CCTATTATCATCGCTCstochastic_pseudotime"
#>  [427] "GGCAATTCATCTCCCAstochastic_pseudotime"
#>  [428] "GCTGGGTGTCATACTGstochastic_pseudotime"
#>  [429] "GAGCAGAGTTTGCATGstochastic_pseudotime"
#>  [430] "GCAAACTCAGACAAATstochastic_pseudotime"
#>  [431] "ATCATGGGTCCGCTGAstochastic_pseudotime"
#>  [432] "GGGAGATAGATGTGGCstochastic_pseudotime"
#>  [433] "CGTGAGCTCAGCGATTstochastic_pseudotime"
#>  [434] "GAACGGACAAGGACACstochastic_pseudotime"
#>  [435] "CTAAGACTCTGCGGCAstochastic_pseudotime"
#>  [436] "ACAGCCGCATCGTCGGstochastic_pseudotime"
#>  [437] "GTCGTAAAGTTTAGGAstochastic_pseudotime"
#>  [438] "GGGAATGTCACAAACCstochastic_pseudotime"
#>  [439] "GTCACAACACGTCTCTstochastic_pseudotime"
#>  [440] "TTTGGTTCACGAAGCAstochastic_pseudotime"
#>  [441] "GCGCAACGTTTCGCTCstochastic_pseudotime"
#>  [442] "TACCTTATCCGATATGstochastic_pseudotime"
#>  [443] "TGAGAGGAGGCGATACstochastic_pseudotime"
#>  [444] "TCCCGATCAACGATGGstochastic_pseudotime"
#>  [445] "TTGGCAAAGAGTAAGGstochastic_pseudotime"
#>  [446] "CCTTCGAAGATAGCATstochastic_pseudotime"
#>  [447] "GGCGACTGTAAGTTCCstochastic_pseudotime"
#>  [448] "TGAGGGAAGAGGTACCstochastic_pseudotime"
#>  [449] "ATCTGCCCAATCCGATstochastic_pseudotime"
#>  [450] "CACAGGCCACGTCAGCstochastic_pseudotime"
#>  [451] "TACTTACCATGGGACAstochastic_pseudotime"
#>  [452] "AATCGGTCATATGGTCstochastic_pseudotime"
#>  [453] "TTGAACGAGAGCAATTstochastic_pseudotime"
#>  [454] "CTTGGCTCAGCCACCAstochastic_pseudotime"
#>  [455] "CGAATGTTCTTCATGTstochastic_pseudotime"
#>  [456] "TGAGAGGAGCCGTCGTstochastic_pseudotime"
#>  [457] "AGCTTGACATCTATGGstochastic_pseudotime"
#>  [458] "ACTTGTTAGTGAAGTTstochastic_pseudotime"
#>  [459] "TACGGATTCCAAGTACstochastic_pseudotime"
#>  [460] "CTTAGGATCACTCTTAstochastic_pseudotime"
#>  [461] "GGATGTTTCCAGTATGstochastic_pseudotime"
#>  [462] "ACGCCGAGTTTGCATGstochastic_pseudotime"
#>  [463] "TCTCATAGTGGACGATstochastic_pseudotime"
#>  [464] "CGGCTAGGTACTCGCGstochastic_pseudotime"
#>  [465] "CACAAACCACGAAAGCstochastic_pseudotime"
#>  [466] "CTAAGACAGCACGCCTstochastic_pseudotime"
#>  [467] "GTGGGTCCATAAAGGTstochastic_pseudotime"
#>  [468] "ACACCAAGTATAGGTAstochastic_pseudotime"
#>  [469] "CTGGTCTTCGCATGGCstochastic_pseudotime"
#>  [470] "ACTGTCCTCGCTAGCGstochastic_pseudotime"
#>  [471] "TCTCATACAAATTGCCstochastic_pseudotime"
#>  [472] "GATCGATGTAGCGATGstochastic_pseudotime"
#>  [473] "CAACCAACACCTCGTTstochastic_pseudotime"
#>  [474] "CATATGGGTTACTGACstochastic_pseudotime"
#>  [475] "AAACGGGCAAAGAATCstochastic_pseudotime"
#>  [476] "CATCGGGTCTTCGGTCstochastic_pseudotime"
#>  [477] "TCTCTAAGTGGCAAACstochastic_pseudotime"
#>  [478] "CTCTGGTAGCAATCTCstochastic_pseudotime"
#>  [479] "CCGTGGATCTGTCTCGstochastic_pseudotime"
#>  [480] "CGGAGCTAGATAGTCAstochastic_pseudotime"
#>  [481] "GACCAATTCAACGGCCstochastic_pseudotime"
#>  [482] "CATCGAAGTCACACGCstochastic_pseudotime"
#>  [483] "GATTCAGTCGTCCAGGstochastic_pseudotime"
#>  [484] "TACTCGCGTCGGCATCstochastic_pseudotime"
#>  [485] "GTCTCGTTCACTTCATstochastic_pseudotime"
#>  [486] "GCCTCTAAGCTAGTCTstochastic_pseudotime"
#>  [487] "ATCCACCCAGCTCGCAstochastic_pseudotime"
#>  [488] "CGAGAAGAGGGTTCCCstochastic_pseudotime"
#>  [489] "GAAATGAAGCTCCTCTstochastic_pseudotime"
#>  [490] "GGGAGATAGCGAAGGGstochastic_pseudotime"
#>  [491] "CGGTTAAGTCTAGCCGstochastic_pseudotime"
#>  [492] "ACGCCAGTCCCAAGTAstochastic_pseudotime"
#>  [493] "CATGGCGAGAGTCTGGstochastic_pseudotime"
#>  [494] "CCACGGAGTCAACATCstochastic_pseudotime"
#>  [495] "GACGCGTGTCGGCTCAstochastic_pseudotime"
#>  [496] "GTATTCTTCTCTAAGGstochastic_pseudotime"
#>  [497] "CTTAGGATCAACGGGAstochastic_pseudotime"
#>  [498] "GGAATAAGTAAGAGAGstochastic_pseudotime"
#>  [499] "CGGACGTAGACACGACstochastic_pseudotime"
#>  [500] "GACGTTAAGACGCTTTstochastic_pseudotime"
#>  [501] "CCGTACTAGGAGTACCstochastic_pseudotime"
#>  [502] "CCGTACTAGTCCGTATstochastic_pseudotime"
#>  [503] "ATGTGTGTCAAAGACAstochastic_pseudotime"
#>  [504] "AACCATGCACAGGAGTstochastic_pseudotime"
#>  [505] "GAGTCCGAGCGATGACstochastic_pseudotime"
#>  [506] "GTATCTTGTCGCGGTTstochastic_pseudotime"
#>  [507] "AGATTGCAGTTCCACAstochastic_pseudotime"
#>  [508] "CTGAAACTCTGTTTGTstochastic_pseudotime"
#>  [509] "CCACGGAGTTACGACTstochastic_pseudotime"
#>  [510] "GTACTTTGTACCGTATstochastic_pseudotime"
#>  [511] "AGCGTCGCACAGGCCTstochastic_pseudotime"
#>  [512] "CTGAAACTCGATAGAAstochastic_pseudotime"
#>  [513] "TTGAACGAGGAGTTTAstochastic_pseudotime"
#>  [514] "ATCATCTGTGTTTGTGstochastic_pseudotime"
#>  [515] "GTTCATTTCCGAACGCstochastic_pseudotime"
#>  [516] "CGTGTAATCATTATCCstochastic_pseudotime"
#>  [517] "AGGGATGGTTAGGGTGstochastic_pseudotime"
#>  [518] "TGACTAGGTGTATGGGstochastic_pseudotime"
#>  [519] "AGCTCTCGTTCAGTACstochastic_pseudotime"
#>  [520] "TTTCCTCCAGGATCGAstochastic_pseudotime"
#>  [521] "TCTTTCCGTGCAACTTstochastic_pseudotime"
#>  [522] "CTGAAACGTATGGTTCstochastic_pseudotime"
#>  [523] "CGTCAGGTCCGCAAGCstochastic_pseudotime"
#>  [524] "GTACTCCTCCAACCAAstochastic_pseudotime"
#>  [525] "AGGCCGTCACGACGAAstochastic_pseudotime"
#>  [526] "CGAGAAGCACGCATCGstochastic_pseudotime"
#>  [527] "ATGTGTGCACAAGTAAstochastic_pseudotime"
#>  [528] "AGTCTTTGTTACGGAGstochastic_pseudotime"
#>  [529] "TGATTTCGTGAGGGTTstochastic_pseudotime"
#>  [530] "CCTATTACAGACAAGCstochastic_pseudotime"
#>  [531] "CCTACCATCTGCCCTAstochastic_pseudotime"
#>  [532] "AGCGGTCGTCGAATCTstochastic_pseudotime"
#>  [533] "GGGATGAGTTGTGGAGstochastic_pseudotime"
#>  [534] "GCTGGGTGTCTAGAGGstochastic_pseudotime"
#>  [535] "AGAGCTTAGGTTACCTstochastic_pseudotime"
#>  [536] "GATGAGGAGTGGACGTstochastic_pseudotime"
#>  [537] "GTCTCGTCATAGTAAGstochastic_pseudotime"
#>  [538] "CGTCCATGTCTAACGTstochastic_pseudotime"
#>  [539] "GTGCGGTCAGCGTAAGstochastic_pseudotime"
#>  [540] "GACCAATGTTGGACCCstochastic_pseudotime"
#>  [541] "ACATCAGTCCTTGGTCstochastic_pseudotime"
#>  [542] "AAGTCTGTCCGAGCCAstochastic_pseudotime"
#>  [543] "TGCCAAACATGACATCstochastic_pseudotime"
#>  [544] "TGGCTGGGTACCCAATstochastic_pseudotime"
#>  [545] "TGAGCATGTTAAAGACstochastic_pseudotime"
#>  [546] "GGCTCGAAGTGACTCTstochastic_pseudotime"
#>  [547] "GTTTCTAAGAAGATTCstochastic_pseudotime"
#>  [548] "TTCGGTCAGGGCATGTstochastic_pseudotime"
#>  [549] "TGACTAGTCCTATTCAstochastic_pseudotime"
#>  [550] "GTCCTCATCTGGAGCCstochastic_pseudotime"
#>  [551] "GTATCTTGTGGCCCTAstochastic_pseudotime"
#>  [552] "AACTCAGGTGCCTGTGstochastic_pseudotime"
#>  [553] "AAACGGGTCGCATGGCstochastic_pseudotime"
#>  [554] "AGTGAGGTCTCCAACCstochastic_pseudotime"
#>  [555] "TGCTGCTAGGATGCGTstochastic_pseudotime"
#>  [556] "ACTTACTCAATCAGAAstochastic_pseudotime"
#>  [557] "GAAATGAGTCACCCAGstochastic_pseudotime"
#>  [558] "AACTTTCCACAGTCGCstochastic_pseudotime"
#>  [559] "CATGACACATCGATTGstochastic_pseudotime"
#>  [560] "GGGCATCGTGGTCTCGstochastic_pseudotime"
#>  [561] "ACTATCTGTCAGGACAstochastic_pseudotime"
#>  [562] "TACTCGCCAGATCGGAstochastic_pseudotime"
#>  [563] "TATGCCCGTCTAGCCGstochastic_pseudotime"
#>  [564] "CTCGAAAGTAAGAGGAstochastic_pseudotime"
#>  [565] "CTAGTGATCCCAGGTGstochastic_pseudotime"
#>  [566] "CTGTGCTCATTACCTTstochastic_pseudotime"
#>  [567] "AGCTTGAAGAGTGACCstochastic_pseudotime"
#>  [568] "GCAGCCACAGCGTTCGstochastic_pseudotime"
#>  [569] "GGACAGAAGGCCGAATstochastic_pseudotime"
#>  [570] "TTTACTGCAGCATACTstochastic_pseudotime"
#>  [571] "ACCAGTATCGTTTAGGstochastic_pseudotime"
#>  [572] "CTCAGAACACTAGTACstochastic_pseudotime"
#>  [573] "CAGAATCGTCATTAGCstochastic_pseudotime"
#>  [574] "TGTGGTATCCAGTAGTstochastic_pseudotime"
#>  [575] "GGTATTGGTACTCGCGstochastic_pseudotime"
#>  [576] "TCAGGTATCTGTCCGTstochastic_pseudotime"
#>  [577] "GGGTCTGCACAGGCCTstochastic_pseudotime"
#>  [578] "ACTTACTAGAGGGCTTstochastic_pseudotime"
#>  [579] "AGGGAGTTCCTCAATTstochastic_pseudotime"
#>  [580] "GATCTAGAGCTAGCCCstochastic_pseudotime"
#>  [581] "GTATCTTCATCGATGTstochastic_pseudotime"
#>  [582] "GCGGGTTAGACTTGAAstochastic_pseudotime"
#>  [583] "TTATGCTCACCCAGTGstochastic_pseudotime"
#>  [584] "TTCGAAGAGTACGCCCstochastic_pseudotime"
#>  [585] "GCTGCGAAGCTAGCCCstochastic_pseudotime"
#>  [586] "GCCTCTAAGTCCGGTCstochastic_pseudotime"
#>  [587] "CCTCAGTGTTACCAGTstochastic_pseudotime"
#>  [588] "TCACAAGAGTGACATAstochastic_pseudotime"
#>  [589] "ATTATCCAGACCCACCstochastic_pseudotime"
#>  [590] "CTCGAGGTCTGTTGAGstochastic_pseudotime"
#>  [591] "TTCCCAGAGTTACCCAstochastic_pseudotime"
#>  [592] "AGCGGTCGTACCGAGAstochastic_pseudotime"
#>  [593] "TCGGGACAGAGAACAGstochastic_pseudotime"
#>  [594] "ATTATCCAGAAACCTAstochastic_pseudotime"
#>  [595] "CGCGTTTCAGTCTTCCstochastic_pseudotime"
#>  [596] "CCAATCCAGAGTAAGGstochastic_pseudotime"
#>  [597] "ACTTTCACAACTGCTAstochastic_pseudotime"
#>  [598] "GCATGATAGCCTTGATstochastic_pseudotime"
#>  [599] "TATCTCACATTGTGCAstochastic_pseudotime"
#>  [600] "CAGCTGGGTTCTCATTstochastic_pseudotime"
#>  [601] "GTCTCGTAGGCTAGACstochastic_pseudotime"
#>  [602] "GTCCTCATCATTTGGGstochastic_pseudotime"
#>  [603] "GAAACTCAGAAACCTAstochastic_pseudotime"
#>  [604] "CGACCTTAGCGCTTATstochastic_pseudotime"
#>  [605] "CGATCGGCACGCCAGTstochastic_pseudotime"
#>  [606] "TGATTTCAGTACGCCCstochastic_pseudotime"
#>  [607] "CGTCAGGTCCTGCCATstochastic_pseudotime"
#>  [608] "GTGCGGTAGTAGGCCAstochastic_pseudotime"
#>  [609] "TCAGCTCAGTCGTACTstochastic_pseudotime"
#>  [610] "AACTTTCCACTGCCAGstochastic_pseudotime"
#>  [611] "TGCCCATCACTGTGTAstochastic_pseudotime"
#>  [612] "AGAGTGGGTCACCTAAstochastic_pseudotime"
#>  [613] "CTTGGCTCATCCCACTstochastic_pseudotime"
#>  [614] "CGTGTAAGTCTTGATGstochastic_pseudotime"
#>  [615] "CGGACTGTCTTGAGACstochastic_pseudotime"
#>  [616] "ACGAGGAAGTTACCCAstochastic_pseudotime"
#>  [617] "TCTGGAATCGGTGTCGstochastic_pseudotime"
#>  [618] "TGGACGCAGTGTTGAAstochastic_pseudotime"
#>  [619] "AAGGCAGTCCGAGCCAstochastic_pseudotime"
#>  [620] "GGTGTTAAGGATCGCAstochastic_pseudotime"
#>  [621] "ATTGGACAGTAATCCCstochastic_pseudotime"
#>  [622] "CGAGCACTCTTCGAGAstochastic_pseudotime"
#>  [623] "CGACTTCCAGCCTATAstochastic_pseudotime"
#>  [624] "GGGTCTGCATGAAGTAstochastic_pseudotime"
#>  [625] "AAACGGGTCAGCTCTCstochastic_pseudotime"
#>  [626] "TCGTACCCATATACGCstochastic_pseudotime"
#>  [627] "CCGTTCACACATGACTstochastic_pseudotime"
#>  [628] "CGAATGTAGTTCCACAstochastic_pseudotime"
#>  [629] "AGAGCTTAGCCGATTTstochastic_pseudotime"
#>  [630] "CTTAGGAGTCCGAATTstochastic_pseudotime"
#>  [631] "CATGACACAACAACCTstochastic_pseudotime"
#>  [632] "TACGGATGTGTTCGATstochastic_pseudotime"
#>  [633] "TGAGCCGTCCGAACGCstochastic_pseudotime"
#>  [634] "ATCCGAAGTGACTCATstochastic_pseudotime"
#>  [635] "CTCTGGTCAACACCTAstochastic_pseudotime"
#>  [636] "ACATACGAGTAGGCCAstochastic_pseudotime"
#>  [637] "CAGATCAAGCTGTTCAstochastic_pseudotime"
#>  [638] "TGCGGGTTCTCAAACGstochastic_pseudotime"
#>  [639] "TCTTTCCTCGTCTGCTstochastic_pseudotime"
#>  [640] "TTCTACAGTACTCTCCstochastic_pseudotime"
#>  [641] "CACATTTTCGCAGGCTstochastic_pseudotime"
#>  [642] "GATGAGGTCTGGGCCAstochastic_pseudotime"
#>  [643] "CTTTGCGAGGAATGGAstochastic_pseudotime"
#>  [644] "CGCTGGAGTCCATCCTstochastic_pseudotime"
#>  [645] "CCACTACAGAGTACCGstochastic_pseudotime"
#>  [646] "GCAGTTATCTTGAGACstochastic_pseudotime"
#>  [647] "GGCCGATTCACTCTTAstochastic_pseudotime"
#>  [648] "TGGCCAGCACGAAAGCstochastic_pseudotime"
#>  [649] "TCAGGTAAGTCAATAGstochastic_pseudotime"
#>  [650] "AACTTTCGTGACCAAGstochastic_pseudotime"
#>  [651] "CGGTTAATCGTCTGAAstochastic_pseudotime"
#>  [652] "GCCTCTATCGCTTAGAstochastic_pseudotime"
#>  [653] "GTAACGTGTCACTGGCstochastic_pseudotime"
#>  [654] "TTCTACAAGATCGATAstochastic_pseudotime"
#>  [655] "TACACGACAACGATCTstochastic_pseudotime"
#>  [656] "ACTGAGTAGTGACTCTstochastic_pseudotime"
#>  [657] "GAGGTGAGTGACAAATstochastic_pseudotime"
#>  [658] "TTGAACGGTGATGCCCstochastic_pseudotime"
#>  [659] "CAAGGCCAGCTGTTCAstochastic_pseudotime"
#>  [660] "GTATCTTCAGTGGGATstochastic_pseudotime"
#>  [661] "AGGCCACAGCTGTCTAstochastic_pseudotime"
#>  [662] "ACTTACTAGCCGGTAAstochastic_pseudotime"
#>  [663] "GATGAAACATGGTTGTstochastic_pseudotime"
#>  [664] "GACGCGTCAGTCGTGCstochastic_pseudotime"
#>  [665] "TACGGGCAGCGTCAAGstochastic_pseudotime"
#>  [666] "GACGTTAAGACGCAACstochastic_pseudotime"
#>  [667] "CCTACACAGATCCGAGstochastic_pseudotime"
#>  [668] "AGCTCCTGTGTGCCTGstochastic_pseudotime"
#>  [669] "TGTCCCAAGTTTAGGAstochastic_pseudotime"
#>  [670] "CGTCACTAGGGTGTGTstochastic_pseudotime"
#>  [671] "CGTCACTCACGAAATAstochastic_pseudotime"
#>  [672] "AGCTTGACAATAGAGTstochastic_pseudotime"
#>  [673] "CTTGGCTAGGACACCAstochastic_pseudotime"
#>  [674] "CTCAGAATCAGGTTCAstochastic_pseudotime"
#>  [675] "TATCTCACAACTGCTAstochastic_pseudotime"
#>  [676] "TTGGAACGTATATGAGstochastic_pseudotime"
#>  [677] "TTTACTGGTGCTAGCCstochastic_pseudotime"
#>  [678] "CTCGTACCAATGGATAstochastic_pseudotime"
#>  [679] "AGTAGTCCACATGACTstochastic_pseudotime"
#>  [680] "TCAATCTCATCTCCCAstochastic_pseudotime"
#>  [681] "TTGCCGTTCTCTGTCGstochastic_pseudotime"
#>  [682] "AGAGTGGGTAGCACGAstochastic_pseudotime"
#>  [683] "CCACCTATCACCAGGCstochastic_pseudotime"
#>  [684] "CGTTAGAGTGCAACGAstochastic_pseudotime"
#>  [685] "GATCTAGGTGCCTTGGstochastic_pseudotime"
#>  [686] "CGTCAGGCATGGGAACstochastic_pseudotime"
#>  [687] "TCATTACGTCAGAAGCstochastic_pseudotime"
#>  [688] "GCGGGTTGTTACTGACstochastic_pseudotime"
#>  [689] "CCTTACGAGATCACGGstochastic_pseudotime"
#>  [690] "GCATGATTCAAACGGGstochastic_pseudotime"
#>  [691] "AGTAGTCGTTCTGTTTstochastic_pseudotime"
#>  [692] "TGAGCCGTCATGTCCCstochastic_pseudotime"
#>  [693] "GCAATCATCCAACCAAstochastic_pseudotime"
#>  [694] "CAGAATCGTGTAAGTAstochastic_pseudotime"
#>  [695] "GTTCATTAGAAGATTCstochastic_pseudotime"
#>  [696] "TGCGTGGAGGTAAACTstochastic_pseudotime"
#>  [697] "TGATTTCAGTGACTCTstochastic_pseudotime"
#>  [698] "TCAGCAACACCAGGCTstochastic_pseudotime"
#>  [699] "CGTCAGGAGCCTATGTstochastic_pseudotime"
#>  [700] "GGGCACTCATTTCAGGstochastic_pseudotime"
#>  [701] "GAGTCCGAGTTAAGTGstochastic_pseudotime"
#>  [702] "TATCAGGCACCAGATTstochastic_pseudotime"
#>  [703] "GTTCTCGCAGACACTTstochastic_pseudotime"
#>  [704] "GTATTCTCACATTAGCstochastic_pseudotime"
#>  [705] "AACTCCCTCGAGAGCAstochastic_pseudotime"
#>  [706] "TGACTTTAGCGTTTACstochastic_pseudotime"
#>  [707] "CCTTCGACAGCTGTTAstochastic_pseudotime"
#>  [708] "ACACCGGCATCACGTAstochastic_pseudotime"
#>  [709] "AAACGGGAGATATGGTstochastic_pseudotime"
#>  [710] "CTCGTCATCGGTCTAAstochastic_pseudotime"
#>  [711] "GAAATGAGTCGGATCCstochastic_pseudotime"
#>  [712] "TCTCTAAGTCGCATCGstochastic_pseudotime"
#>  [713] "TGGTTCCCACCGCTAGstochastic_pseudotime"
#>  [714] "GACGGCTAGGTGATATstochastic_pseudotime"
#>  [715] "GACTACAGTTCGTGATstochastic_pseudotime"
#>  [716] "TCACGAATCCTATGTTstochastic_pseudotime"
#>  [717] "TGAAAGAGTTATCACGstochastic_pseudotime"
#>  [718] "ATAGACCTCTAACTGGstochastic_pseudotime"
#>  [719] "CAAGGCCAGTGTGGCAstochastic_pseudotime"
#>  [720] "ACTTTCAAGGATATACstochastic_pseudotime"
#>  [721] "TACGGATAGTTTAGGAstochastic_pseudotime"
#>  [722] "CTCTGGTCAGTTCATGstochastic_pseudotime"
#>  [723] "TGAGGGACAGACTCGCstochastic_pseudotime"
#>  [724] "GCTGCAGGTTACCGATstochastic_pseudotime"
#>  [725] "GGGTTGCTCGAACGGAstochastic_pseudotime"
#>  [726] "CTAAGACAGCCACGCTstochastic_pseudotime"
#>  [727] "CGTCACTCATAGTAAGstochastic_pseudotime"
#>  [728] "GGTGTTAGTTCCTCCAstochastic_pseudotime"
#>  [729] "AAAGATGAGTCTTGCAstochastic_pseudotime"
#>  [730] "AACGTTGGTGCCTTGGstochastic_pseudotime"
#>  [731] "CATGACAGTAGCTTGTstochastic_pseudotime"
#>  [732] "CTGGTCTCAGGAATGCstochastic_pseudotime"
#>  [733] "CACAGTATCCTCCTAGstochastic_pseudotime"
#>  [734] "CGTAGGCAGCGTAGTGstochastic_pseudotime"
#>  [735] "CCTTCCCCAAAGTGCGstochastic_pseudotime"
#>  [736] "GACACGCCAAAGGAAGstochastic_pseudotime"
#>  [737] "AGCAGCCTCTTTAGTCstochastic_pseudotime"
#>  [738] "ACGGCCACACGTGAGAstochastic_pseudotime"
#>  [739] "TACTTACCAACTGCGCstochastic_pseudotime"
#>  [740] "TAGTTGGGTACCGTATstochastic_pseudotime"
#>  [741] "AGACGTTTCTTACCTAstochastic_pseudotime"
#>  [742] "GAATGAACAATCGAAAstochastic_pseudotime"
#>  [743] "ACTGAACAGCCTCGTGstochastic_pseudotime"
#>  [744] "CGAACATCAGACGTAGstochastic_pseudotime"
#>  [745] "CAAGAAATCGTACCGGstochastic_pseudotime"
#>  [746] "GTACTCCCACGTCTCTstochastic_pseudotime"
#>  [747] "CTTAACTCAAGCGTAGstochastic_pseudotime"
#>  [748] "GGGTTGCCAAGCTGAGstochastic_pseudotime"
#>  [749] "AAAGCAACACATGGGAstochastic_pseudotime"
#>  [750] "GGGACCTTCGAATGCTstochastic_pseudotime"
#>  [751] "TAGACCAAGGAATCGCstochastic_pseudotime"
#>  [752] "CTAACTTCATTGGGCCstochastic_pseudotime"
#>  [753] "TGAGCATTCACCCGAGstochastic_pseudotime"
#>  [754] "GAAGCAGTCAAGGCTTstochastic_pseudotime"
#>  [755] "AAACCTGAGCCTTGATstochastic_pseudotime"
#>  [756] "CCTTCCCAGGCAGGTTstochastic_pseudotime"
#>  [757] "TACAGTGCAATCTACGstochastic_pseudotime"
#>  [758] "GGACGTCCACGCATCGstochastic_pseudotime"
#>  [759] "CACAGTAGTTCCTCCAstochastic_pseudotime"
#>  [760] "TACGGTACATGTTCCCstochastic_pseudotime"
#>  [761] "CTCGTCACAACTGCTAstochastic_pseudotime"
#>  [762] "CGGACTGAGGTAGCTGstochastic_pseudotime"
#>  [763] "AACTGGTTCTTACCGCstochastic_pseudotime"
#>  [764] "ATTTCTGGTAGTGAATstochastic_pseudotime"
#>  [765] "GAAATGACACATGACTstochastic_pseudotime"
#>  [766] "CTCAGAAGTCCTCTTGstochastic_pseudotime"
#>  [767] "CACCTTGAGAGTAATCstochastic_pseudotime"
#>  [768] "TACCTTAAGTACATGAstochastic_pseudotime"
#>  [769] "GAAATGACAACTGGCCstochastic_pseudotime"
#>  [770] "TACGGATAGGCACATGstochastic_pseudotime"
#>  [771] "ATTGGTGTCCAGTATGstochastic_pseudotime"
#>  [772] "CTAGCCTTCAGCGATTstochastic_pseudotime"
#>  [773] "GGGAATGAGTGTACGGstochastic_pseudotime"
#>  [774] "CACATTTGTTTCGCTCstochastic_pseudotime"
#>  [775] "TGCTACCCAGGGCATAstochastic_pseudotime"
#>  [776] "TAAGTGCAGGGTTTCTstochastic_pseudotime"
#>  [777] "TTTGCGCAGACAGACCstochastic_pseudotime"
#>  [778] "TTTGGTTCACCAGATTstochastic_pseudotime"
#>  [779] "GGACAAGGTGATAAGTstochastic_pseudotime"
#>  [780] "GGTATTGGTCCCGACAstochastic_pseudotime"
#>  [781] "CTGAAACAGACTAGGCstochastic_pseudotime"
#>  [782] "ACACCGGCATGTCCTCstochastic_pseudotime"
#>  [783] "TTGACTTTCAGCAACTstochastic_pseudotime"
#>  [784] "AGGGATGTCGCATGATstochastic_pseudotime"
#>  [785] "CGAGCACTCTACCTGCstochastic_pseudotime"
#>  [786] "TTCTCCTCAATGGTCTstochastic_pseudotime"
#>  [787] "GGGCACTGTGGCAAACstochastic_pseudotime"
#>  [788] "TTGACTTTCCAGTAGTstochastic_pseudotime"
#>  [789] "CAGCAGCGTTCACGGCstochastic_pseudotime"
#>  [790] "TCTGGAACACACGCTGstochastic_pseudotime"
#>  [791] "CTGATCCTCAGTTGACstochastic_pseudotime"
#>  [792] "AGATTGCCACCAGGCTstochastic_pseudotime"
#>  [793] "ATCATGGTCCAGAGGAstochastic_pseudotime"
#>  [794] "AAGGCAGAGCTGAAATstochastic_pseudotime"
#>  [795] "TTCCCAGGTATTCTCTstochastic_pseudotime"
#>  [796] "CCTAGCTCAGACTCGCstochastic_pseudotime"
#>  [797] "CGGGTCAAGGCTCAGAstochastic_pseudotime"
#>  [798] "TGGCCAGGTCGGCATCstochastic_pseudotime"
#>  [799] "ACTGATGTCCCTAACCstochastic_pseudotime"
#>  [800] "CCGGTAGTCGGCTTGGstochastic_pseudotime"
#>  [801] "TCTCTAAGTTACCGATstochastic_pseudotime"
#>  [802] "TGACGGCAGGAATGGAstochastic_pseudotime"
#>  [803] "TAGCCGGTCAACACACstochastic_pseudotime"
#>  [804] "CGGACACAGTGTCTCAstochastic_pseudotime"
#>  [805] "TACGGATAGAGCTGGTstochastic_pseudotime"
#>  [806] "TCAGATGGTACACCGCstochastic_pseudotime"
#>  [807] "TGGGAAGCATGACATCstochastic_pseudotime"
#>  [808] "CGATGGCGTTACGGAGstochastic_pseudotime"
#>  [809] "CTGAAACTCGCCGTGAstochastic_pseudotime"
#>  [810] "CTACCCAGTGCGCTTGstochastic_pseudotime"
#>  [811] "CTTAGGACACAGGTTTstochastic_pseudotime"
#>  [812] "ACGTCAACACTGCCAGstochastic_pseudotime"
#>  [813] "AAAGATGAGTCATGCTstochastic_pseudotime"
#>  [814] "TGGTTCCGTCACCTAAstochastic_pseudotime"
#>  [815] "CGTCTACAGATCTGCTstochastic_pseudotime"
#>  [816] "CGAGCACAGCTGAAATstochastic_pseudotime"
#>  [817] "TTCTCAAAGCCGCCTAstochastic_pseudotime"
#>  [818] "GAACATCGTAATCGTCstochastic_pseudotime"
#>  [819] "GCGCAACTCACAACGTstochastic_pseudotime"
#>  [820] "AGCGGTCGTGTATGGGstochastic_pseudotime"
#>  [821] "AGGGAGTAGTGTCTCAstochastic_pseudotime"
#>  [822] "CGGAGCTCATTGGGCCstochastic_pseudotime"
#>  [823] "GCATGATCAAGCGAGTstochastic_pseudotime"
#>  [824] "GAGTCCGCAATTGCTGstochastic_pseudotime"
#>  [825] "CTAGCCTGTCTCTTATstochastic_pseudotime"
#>  [826] "GGAGCAAAGCTTCGCGstochastic_pseudotime"
#>  [827] "CCCATACCACATCCAAstochastic_pseudotime"
#>  [828] "CTACACCGTCGAATCTstochastic_pseudotime"
#>  [829] "CATATGGCATGTTGACstochastic_pseudotime"
#>  [830] "CACACAACAGCCAGAAstochastic_pseudotime"
#>  [831] "GGAATAAGTTCGGGCTstochastic_pseudotime"
#>  [832] "GGGCATCCATGACATCstochastic_pseudotime"
#>  [833] "GTGTTAGCACCATCCTstochastic_pseudotime"
#>  [834] "CATCAGACATGCATGTstochastic_pseudotime"
#>  [835] "CCTCTGAGTATAAACGstochastic_pseudotime"
#>  [836] "ACAGCTATCAGCATGTstochastic_pseudotime"
#>  [837] "CAACCAATCGGGAGTAstochastic_pseudotime"
#>  [838] "GTGCATAAGCTCCTCTstochastic_pseudotime"
#>  [839] "AGCTTGACAGTAAGATstochastic_pseudotime"
#>  [840] "CATCGAAGTACGAAATstochastic_pseudotime"
#>  [841] "ATCCACCGTACCATCAstochastic_pseudotime"
#>  [842] "CCGTGGATCCAAGCCGstochastic_pseudotime"
#>  [843] "CTGCCTAGTCTCCATCstochastic_pseudotime"
#>  [844] "GGAAAGCGTCTTGATGstochastic_pseudotime"
#>  [845] "AGCTTGAGTGCAGTAGstochastic_pseudotime"
#>  [846] "CCTACACTCCTTTACAstochastic_pseudotime"
#>  [847] "CTTCTCTAGGCAGGTTstochastic_pseudotime"
#>  [848] "TCTATTGCACAGGAGTstochastic_pseudotime"
#>  [849] "CGTTAGAAGTATCTCGstochastic_pseudotime"
#>  [850] "CACACTCTCATGTGGTstochastic_pseudotime"
#>  [851] "ATTATCCAGTGAACATstochastic_pseudotime"
#>  [852] "AGAGCGATCTCTTGATstochastic_pseudotime"
#>  [853] "AGAATAGCAGTACACTstochastic_pseudotime"
#>  [854] "CCTAGCTGTGGCCCTAstochastic_pseudotime"
#>  [855] "ACAGCCGAGTGGTCCCstochastic_pseudotime"
#>  [856] "GAGCAGAAGCGTGTCCstochastic_pseudotime"
#>  [857] "CGATGGCTCTGTGCAAstochastic_pseudotime"
#>  [858] "CAAGGCCTCTTACCTAstochastic_pseudotime"
#>  [859] "TTGACTTAGACAAGCCstochastic_pseudotime"
#>  [860] "AGGGAGTAGGGTATCGstochastic_pseudotime"
#>  [861] "GACTAACTCAGTCCCTstochastic_pseudotime"
#>  [862] "ACGTCAACAGACAAATstochastic_pseudotime"
#>  [863] "GTCAAGTGTGGTCCGTstochastic_pseudotime"
#>  [864] "TGGCGCATCACAACGTstochastic_pseudotime"
#>  [865] "CAGCGACCACAGGAGTstochastic_pseudotime"
#>  [866] "AGAGCTTGTTAAAGACstochastic_pseudotime"
#>  [867] "GCTGCTTGTCTTCAAGstochastic_pseudotime"
#>  [868] "CACATAGAGCGTGTCCstochastic_pseudotime"
#>  [869] "TGCTGCTGTCCGTCAGstochastic_pseudotime"
#>  [870] "CAGCTAAAGTCCATACstochastic_pseudotime"
#>  [871] "ACTGATGGTGAGGGTTstochastic_pseudotime"
#>  [872] "CACACTCAGATCCCATstochastic_pseudotime"
#>  [873] "CCTAGCTCATACTACGstochastic_pseudotime"
#>  [874] "CTAGTGATCAAAGACAstochastic_pseudotime"
#>  [875] "TGTATTCCAATGACCTstochastic_pseudotime"
#>  [876] "AAGACCTCAGGTTTCAstochastic_pseudotime"
#>  [877] "ACCGTAAGTTATCACGstochastic_pseudotime"
#>  [878] "ACGGGCTGTGCGATAGstochastic_pseudotime"
#>  [879] "TCTGAGAGTCAAACTCstochastic_pseudotime"
#>  [880] "CTAGAGTCACTGTTAGstochastic_pseudotime"
#>  [881] "CTGCGGATCAGCATGTstochastic_pseudotime"
#>  [882] "CAGCCGAAGCGATATAstochastic_pseudotime"
#>  [883] "CAGCAGCAGACAAAGGstochastic_pseudotime"
#>  [884] "TGTGGTACAAGACGTGstochastic_pseudotime"
#>  [885] "TCGTACCCACGCATCGstochastic_pseudotime"
#>  [886] "GATGCTAGTAGAAGGAstochastic_pseudotime"
#>  [887] "CCTAGCTTCGATGAGGstochastic_pseudotime"
#>  [888] "GTAACTGTCAGTGCATstochastic_pseudotime"
#>  [889] "GGACGTCGTCTAGGTTstochastic_pseudotime"
#>  [890] "ATAAGAGCATCGATTGstochastic_pseudotime"
#>  [891] "CGGTTAAAGAGCTGCAstochastic_pseudotime"
#>  [892] "CGCCAAGCACCGGAAAstochastic_pseudotime"
#>  [893] "AGGCCGTGTACCAGTTstochastic_pseudotime"
#>  [894] "CGGTTAAGTAAACACAstochastic_pseudotime"
#>  [895] "CTCTGGTAGTCCCACGstochastic_pseudotime"
#>  [896] "ACAGCTATCTGGGCCAstochastic_pseudotime"
#>  [897] "AACTCAGTCTGGTTCCstochastic_pseudotime"
#>  [898] "TGAGAGGGTAAGGGCTstochastic_pseudotime"
#>  [899] "CCCAGTTAGCATCATCstochastic_pseudotime"
#>  [900] "TGCGCAGTCTGCCAGGstochastic_pseudotime"
#>  [901] "CGTGTAATCCATGCTCstochastic_pseudotime"
#>  [902] "TGACTAGTCCACGTGGstochastic_pseudotime"
#>  [903] "TAGAGCTGTATGAAACstochastic_pseudotime"
#>  [904] "TAGAGCTGTTGTCGCGstochastic_pseudotime"
#>  [905] "TTTATGCTCCATTCTAstochastic_pseudotime"
#>  [906] "AAGCCGCAGGCAGGTTstochastic_pseudotime"
#>  [907] "TTTATGCGTCGAACAGstochastic_pseudotime"
#>  [908] "CGAACATAGCCAGTAGstochastic_pseudotime"
#>  [909] "GTCGGGTCAGGCTGAAstochastic_pseudotime"
#>  [910] "CACATTTTCGAGGTAGstochastic_pseudotime"
#>  [911] "GGACAGATCCTGCAGGstochastic_pseudotime"
#>  [912] "GATCGTAGTCCTAGCGstochastic_pseudotime"
#>  [913] "TGACTTTAGCGTTGCCstochastic_pseudotime"
#>  [914] "CTCATTACAAGCGCTCstochastic_pseudotime"
#>  [915] "TACGGGCTCTTCGGTCstochastic_pseudotime"
#>  [916] "TGCGCAGGTACGCACCstochastic_pseudotime"
#>  [917] "AACTGGTGTGCCTGTGstochastic_pseudotime"
#>  [918] "CTACACCGTCGCCATGstochastic_pseudotime"
#>  [919] "TAGGCATTCCAAAGTCstochastic_pseudotime"
#>  [920] "CAACCAACAACCGCCAstochastic_pseudotime"
#>  [921] "CCATTCGGTAGGGACTstochastic_pseudotime"
#>  [922] "TATGCCCTCTTACCTAstochastic_pseudotime"
#>  [923] "TTAACTCCAAAGAATCstochastic_pseudotime"
#>  [924] "GATGAAAAGTTGTAGAstochastic_pseudotime"
#>  [925] "CCGTACTTCAATAAGGstochastic_pseudotime"
#>  [926] "ATTACTCTCAGGCGAAstochastic_pseudotime"
#>  [927] "TGGCCAGAGCTACCTAstochastic_pseudotime"
#>  [928] "GAATGAACAGCCTTGGstochastic_pseudotime"
#>  [929] "TTCTACATCGCTTAGAstochastic_pseudotime"
#>  [930] "TAGTGGTCAACTGCGCstochastic_pseudotime"
#>  [931] "CTAGCCTGTTCGGGCTstochastic_pseudotime"
#>  [932] "GGACGTCGTTTAGCTGstochastic_pseudotime"
#>  [933] "AACTCAGGTTGGAGGTstochastic_pseudotime"
#>  [934] "TGAGCATTCGGGAGTAstochastic_pseudotime"
#>  [935] "AGTGAGGGTTACGGAGstochastic_pseudotime"
#>  [936] "CAGCAGCCACACTGCGstochastic_pseudotime"
#>  [937] "ATTACTCAGTTGCAGGstochastic_pseudotime"
#>  [938] "AGTCTTTCACCAGGTCstochastic_pseudotime"
#>  [939] "AGTGGGACAGGTCGTCstochastic_pseudotime"
#>  [940] "CAGCCGAAGTAGGTGCstochastic_pseudotime"
#>  [941] "TATCTCAGTTCGCTAAstochastic_pseudotime"
#>  [942] "TTCTACATCTGCAGTAstochastic_pseudotime"
#>  [943] "ATCATCTTCTGTCAAGstochastic_pseudotime"
#>  [944] "TGGACGCTCAGGTAAAstochastic_pseudotime"
#>  [945] "TCTCTAAAGCTGGAACstochastic_pseudotime"
#>  [946] "TGAGAGGAGATCCCATstochastic_pseudotime"
#>  [947] "GTATTCTGTCCGTCAGstochastic_pseudotime"
#>  [948] "CCATTCGCATGCCTAAstochastic_pseudotime"
#>  [949] "AGATTGCAGCTGTCTAstochastic_pseudotime"
#>  [950] "GTATTCTCACTTGGATstochastic_pseudotime"
#>  [951] "CGACCTTAGGCTAGACstochastic_pseudotime"
#>  [952] "CTTGGCTAGACAGACCstochastic_pseudotime"
#>  [953] "TCAATCTGTCAAAGATstochastic_pseudotime"
#>  [954] "CTTGGCTGTTAGAACAstochastic_pseudotime"
#>  [955] "TACTCATCAGATCCATstochastic_pseudotime"
#>  [956] "TCGCGAGAGTCCGTATstochastic_pseudotime"
#>  [957] "ACGAGGATCATCACCCstochastic_pseudotime"
#>  [958] "CCCATACAGACTTTCGstochastic_pseudotime"
#>  [959] "AGGCCGTCACGGTGTCstochastic_pseudotime"
#>  [960] "CCGTGGACATAGAAACstochastic_pseudotime"
#>  [961] "GCATGTAAGAATCTCCstochastic_pseudotime"
#>  [962] "CACTCCATCAACACACstochastic_pseudotime"
#>  [963] "GTGCTTCTCACGACTAstochastic_pseudotime"
#>  [964] "GGCTCGACAGCTTCGGstochastic_pseudotime"
#>  [965] "GCTCCTAAGGGCATGTstochastic_pseudotime"
#>  [966] "GTCTCGTTCTGCGTAAstochastic_pseudotime"
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
#>  [981] "CGTTCTGTCGATAGAAstochastic_pseudotime"
#>  [982] "CTACACCCATCAGTCAstochastic_pseudotime"
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
#>  [998] "TGGCCAGCAAGGACTGstochastic_pseudotime"
#>  [999] "TAGTTGGTCTTGCAAGstochastic_pseudotime"
#> [1000] "TGTGGTACACATCCGGstochastic_pseudotime"
#> 
#> 
#> $feature_split
#> NULL
#> 
#> $cell_metadata
#>                             cells Pseudotime_assign stochastic_pseudotime
#> AAACCTGAGCCTTGAT AAACCTGAGCCTTGAT        0.39345182            0.39345182
#> AAACCTGGTAAGTGGC AAACCTGGTAAGTGGC        0.14588694            0.14588694
#> AAACGGGAGATATGGT AAACGGGAGATATGGT        0.39079586            0.39079586
#> AAACGGGCAAAGAATC AAACGGGCAAAGAATC        0.24873873            0.24873873
#> AAACGGGGTACAGTTC AAACGGGGTACAGTTC        0.97053641            0.97053641
#> AAACGGGTCAGCTCTC AAACGGGTCAGCTCTC        0.38303477            0.38303477
#> AAACGGGTCGCATGGC AAACGGGTCGCATGGC        0.31480847            0.31480847
#> AAAGATGAGTCATGCT AAAGATGAGTCATGCT        0.40507686            0.40507686
#> AAAGATGAGTCTTGCA AAAGATGAGTCTTGCA        0.39159208            0.39159208
#> AAAGCAAAGCTACCTA AAAGCAAAGCTACCTA        0.10042365            0.10042365
#> AAAGCAACACATGGGA AAAGCAACACATGGGA        0.39323271            0.39323271
#> AAAGCAAGTGCAGTAG AAAGCAAGTGCAGTAG        0.15416147            0.15416147
#> AAATGCCAGACTAAGT AAATGCCAGACTAAGT        0.70807697            0.70807697
#> AAATGCCCAATGGAAT AAATGCCCAATGGAAT        0.19652438            0.19652438
#> AAATGCCCACCAACCG AAATGCCCACCAACCG        0.20032321            0.20032321
#> AAATGCCCATGCTAGT AAATGCCCATGCTAGT        0.14425109            0.14425109
#> AAATGCCGTGATAAGT AAATGCCGTGATAAGT        0.20673283            0.20673283
#> AACACGTGTCGACTGC AACACGTGTCGACTGC        0.10294021            0.10294021
#> AACCATGCACAGGAGT AACCATGCACAGGAGT        0.27257975            0.27257975
#> AACCATGGTTCGGCAC AACCATGGTTCGGCAC        0.21154504            0.21154504
#> AACCATGGTTGAACTC AACCATGGTTGAACTC        0.11211866            0.11211866
#> AACGTTGGTGCCTTGG AACGTTGGTGCCTTGG        0.39167676            0.39167676
#> AACTCAGGTGCCTGTG AACTCAGGTGCCTGTG        0.31451971            0.31451971
#> AACTCAGGTTGGAGGT AACTCAGGTTGGAGGT        0.45340875            0.45340875
#> AACTCAGTCCTCGCAT AACTCAGTCCTCGCAT        0.20163494            0.20163494
#> AACTCAGTCTGGTTCC AACTCAGTCTGGTTCC        0.44248692            0.44248692
#> AACTCCCCAGTTCATG AACTCCCCAGTTCATG        0.14971548            0.14971548
#> AACTCCCTCGAGAGCA AACTCCCTCGAGAGCA        0.39060195            0.39060195
#> AACTCTTCATTTGCCC AACTCTTCATTTGCCC        0.20478129            0.20478129
#> AACTGGTGTGCCTGTG AACTGGTGTGCCTGTG        0.44818742            0.44818742
#> AACTGGTTCTTACCGC AACTGGTTCTTACCGC        0.39434602            0.39434602
#> AACTTTCCACAGTCGC AACTTTCCACAGTCGC        0.31820825            0.31820825
#> AACTTTCCACTGCCAG AACTTTCCACTGCCAG        0.37780005            0.37780005
#> AACTTTCGTGACCAAG AACTTTCGTGACCAAG        0.38569419            0.38569419
#> AACTTTCTCCAGATCA AACTTTCTCCAGATCA        0.15412362            0.15412362
#> AAGACCTAGACCTTTG AAGACCTAGACCTTTG        0.17576926            0.17576926
#> AAGACCTAGTAGGCCA AAGACCTAGTAGGCCA        0.17296260            0.17296260
#> AAGACCTCAGGTTTCA AAGACCTCAGGTTTCA        0.43717803            0.43717803
#> AAGCCGCAGGCAGGTT AAGCCGCAGGCAGGTT        0.44531585            0.44531585
#> AAGCCGCCAATGGACG AAGCCGCCAATGGACG        0.13665544            0.13665544
#> AAGGAGCAGTACATGA AAGGAGCAGTACATGA        0.14322014            0.14322014
#> AAGGCAGAGCTGAAAT AAGGCAGAGCTGAAAT        0.39976643            0.39976643
#> AAGGCAGTCCGAGCCA AAGGCAGTCCGAGCCA        0.38270817            0.38270817
#> AAGGTTCAGATAGTCA AAGGTTCAGATAGTCA        0.12787275            0.12787275
#> AAGGTTCAGGAGCGAG AAGGTTCAGGAGCGAG        0.15017861            0.15017861
#> AAGGTTCGTCAACTGT AAGGTTCGTCAACTGT        0.10978795            0.10978795
#> AAGTCTGGTTGGGACA AAGTCTGGTTGGGACA        0.72065292            0.72065292
#> AAGTCTGTCCGAGCCA AAGTCTGTCCGAGCCA        0.29757185            0.29757185
#> AATCCAGAGGGTTCCC AATCCAGAGGGTTCCC        0.09858278            0.09858278
#> AATCCAGTCAGGTTCA AATCCAGTCAGGTTCA        0.08717351            0.08717351
#> AATCCAGTCGGTCCGA AATCCAGTCGGTCCGA        0.08589006            0.08589006
#> AATCGGTCATATGGTC AATCGGTCATATGGTC        0.23266474            0.23266474
#> AATCGGTTCCAAACTG AATCGGTTCCAAACTG        0.19908460            0.19908460
#> AATCGGTTCGTATCAG AATCGGTTCGTATCAG        0.16614816            0.16614816
#> ACACCAAGTATAGGTA ACACCAAGTATAGGTA        0.24389748            0.24389748
#> ACACCCTAGGACAGCT ACACCCTAGGACAGCT        0.18529915            0.18529915
#> ACACCCTTCCCAGGTG ACACCCTTCCCAGGTG        0.11084609            0.11084609
#> ACACCGGCATCACGTA ACACCGGCATCACGTA        0.39068513            0.39068513
#> ACACCGGCATGTCCTC ACACCGGCATGTCCTC        0.39735103            0.39735103
#> ACACCGGTCATGTCTT ACACCGGTCATGTCTT        0.16415294            0.16415294
#> ACACTGAGTCCATGAT ACACTGAGTCCATGAT        0.13644643            0.13644643
#> ACAGCCGAGAGAACAG ACAGCCGAGAGAACAG        0.13092033            0.13092033
#> ACAGCCGAGAGGTACC ACAGCCGAGAGGTACC        0.19307210            0.19307210
#> ACAGCCGAGTGGTCCC ACAGCCGAGTGGTCCC        0.43149944            0.43149944
#> ACAGCCGCATCGTCGG ACAGCCGCATCGTCGG        0.22017959            0.22017959
#> ACAGCTAAGTCGTTTG ACAGCTAAGTCGTTTG        0.16875456            0.16875456
#> ACAGCTATCAAGGCTT ACAGCTATCAAGGCTT        0.19275711            0.19275711
#> ACAGCTATCAGCATGT ACAGCTATCAGCATGT        0.41686735            0.41686735
#> ACAGCTATCTGGGCCA ACAGCTATCTGGGCCA        0.44207659            0.44207659
#> ACATACGAGTAGGCCA ACATACGAGTAGGCCA        0.38449519            0.38449519
#> ACATACGGTCTCTCTG ACATACGGTCTCTCTG        0.10658971            0.10658971
#> ACATCAGAGAGACTTA ACATCAGAGAGACTTA        0.07995213            0.07995213
#> ACATCAGTCCTTGGTC ACATCAGTCCTTGGTC        0.29604242            0.29604242
#> ACCAGTAAGTGATCGG ACCAGTAAGTGATCGG        0.17017808            0.17017808
#> ACCAGTAGTTCGTTGA ACCAGTAGTTCGTTGA        0.15163955            0.15163955
#> ACCAGTATCGTTTAGG ACCAGTATCGTTTAGG        0.33759604            0.33759604
#> ACCGTAAGTTATCACG ACCGTAAGTTATCACG        0.43744830            0.43744830
#> ACCTTTATCCTGCAGG ACCTTTATCCTGCAGG        0.19238325            0.19238325
#> ACGAGCCCATGAGCGA ACGAGCCCATGAGCGA        0.88415042            0.88415042
#> ACGAGGAAGTTACCCA ACGAGGAAGTTACCCA        0.38249103            0.38249103
#> ACGAGGATCATCACCC ACGAGGATCATCACCC        0.50058010            0.50058010
#> ACGAGGATCGCCGTGA ACGAGGATCGCCGTGA        0.06777252            0.06777252
#> ACGATACAGATGAGAG ACGATACAGATGAGAG        0.17460040            0.17460040
#> ACGATACTCGCGCCAA ACGATACTCGCGCCAA        0.17223070            0.17223070
#> ACGATGTTCGCAGGCT ACGATGTTCGCAGGCT        0.14044275            0.14044275
#> ACGCAGCCATTCTCAT ACGCAGCCATTCTCAT        0.14797536            0.14797536
#> ACGCAGCGTAATCACC ACGCAGCGTAATCACC        0.21083698            0.21083698
#> ACGCCAGGTTACAGAA ACGCCAGGTTACAGAA        0.20327151            0.20327151
#> ACGCCAGTCCCAAGTA ACGCCAGTCCCAAGTA        0.26412116            0.26412116
#> ACGCCGAAGGGTGTTG ACGCCGAAGGGTGTTG        0.14495407            0.14495407
#> ACGCCGAGTTTGCATG ACGCCGAGTTTGCATG        0.23709715            0.23709715
#> ACGGCCACACGTGAGA ACGGCCACACGTGAGA        0.39232283            0.39232283
#> ACGGGCTGTGCGATAG ACGGGCTGTGCGATAG        0.43771263            0.43771263
#> ACGTCAAAGGTGTTAA ACGTCAAAGGTGTTAA        0.67030736            0.67030736
#> ACGTCAACACTGCCAG ACGTCAACACTGCCAG        0.40496733            0.40496733
#> ACGTCAACAGACAAAT ACGTCAACAGACAAAT        0.43477887            0.43477887
#> ACTATCTGTCAGGACA ACTATCTGTCAGGACA        0.32347529            0.32347529
#> ACTATCTTCAAGATCC ACTATCTTCAAGATCC        0.13640802            0.13640802
#> ACTGAACAGCCTCGTG ACTGAACAGCCTCGTG        0.39269714            0.39269714
#> ACTGAGTAGTGACTCT ACTGAGTAGTGACTCT        0.38612550            0.38612550
#> ACTGATGCATCAGTAC ACTGATGCATCAGTAC        0.20853102            0.20853102
#> ACTGATGCATGGTCTA ACTGATGCATGGTCTA        0.20438578            0.20438578
#> ACTGATGCATGTCTCC ACTGATGCATGTCTCC        0.10877167            0.10877167
#> ACTGATGGTGAGGGTT ACTGATGGTGAGGGTT        0.43650787            0.43650787
#> ACTGATGTCCCTAACC ACTGATGTCCCTAACC        0.40097538            0.40097538
#> ACTGCTCAGGTGTTAA ACTGCTCAGGTGTTAA        0.11356938            0.11356938
#> ACTGTCCCAGTTTACG ACTGTCCCAGTTTACG        0.12809540            0.12809540
#> ACTGTCCTCGCTAGCG ACTGTCCTCGCTAGCG        0.24502888            0.24502888
#> ACTTACTAGAGGGCTT ACTTACTAGAGGGCTT        0.34386273            0.34386273
#> ACTTACTAGCCGGTAA ACTTACTAGCCGGTAA        0.38661403            0.38661403
#> ACTTACTCAATCAGAA ACTTACTCAATCAGAA        0.31635270            0.31635270
#> ACTTACTCACATGTGT ACTTACTCACATGTGT        0.88112383            0.88112383
#> ACTTGTTAGTGAAGTT ACTTGTTAGTGAAGTT        0.23511193            0.23511193
#> ACTTGTTAGTTAAGTG ACTTGTTAGTTAAGTG        0.17543828            0.17543828
#> ACTTGTTCATGGTCAT ACTTGTTCATGGTCAT        0.68140359            0.68140359
#> ACTTGTTTCACCCGAG ACTTGTTTCACCCGAG        0.19935642            0.19935642
#> ACTTTCAAGGATATAC ACTTTCAAGGATATAC        0.39132529            0.39132529
#> ACTTTCACAACTGCTA ACTTTCACAACTGCTA        0.36647111            0.36647111
#> AGAATAGCAATCACAC AGAATAGCAATCACAC        0.16064834            0.16064834
#> AGAATAGCAGTACACT AGAATAGCAGTACACT        0.43111105            0.43111105
#> AGACGTTCAATGGTCT AGACGTTCAATGGTCT        0.72447304            0.72447304
#> AGACGTTTCTTACCTA AGACGTTTCTTACCTA        0.39248940            0.39248940
#> AGAGCGACAGTAAGAT AGAGCGACAGTAAGAT        0.21388426            0.21388426
#> AGAGCGATCTCTTGAT AGAGCGATCTCTTGAT        0.43045973            0.43045973
#> AGAGCTTAGCCGATTT AGAGCTTAGCCGATTT        0.38356057            0.38356057
#> AGAGCTTAGGTTACCT AGAGCTTAGGTTACCT        0.29206395            0.29206395
#> AGAGCTTGTTAAAGAC AGAGCTTGTTAAAGAC        0.43572125            0.43572125
#> AGAGTGGAGGCGATAC AGAGTGGAGGCGATAC        0.20897696            0.20897696
#> AGAGTGGCACCTCGGA AGAGTGGCACCTCGGA        0.10814014            0.10814014
#> AGAGTGGGTAGCACGA AGAGTGGGTAGCACGA        0.38890929            0.38890929
#> AGAGTGGGTCACCTAA AGAGTGGGTCACCTAA        0.37938672            0.37938672
#> AGATCTGTCCGTCAAA AGATCTGTCCGTCAAA        0.09990159            0.09990159
#> AGATTGCAGCGTGTCC AGATTGCAGCGTGTCC        0.18973908            0.18973908
#> AGATTGCAGCTGTCTA AGATTGCAGCTGTCTA        0.45880732            0.45880732
#> AGATTGCAGTTCCACA AGATTGCAGTTCCACA        0.27397637            0.27397637
#> AGATTGCCACCAGGCT AGATTGCCACCAGGCT        0.39940311            0.39940311
#> AGCAGCCGTCGCATCG AGCAGCCGTCGCATCG        0.15802756            0.15802756
#> AGCAGCCGTGTAACGG AGCAGCCGTGTAACGG        0.17493695            0.17493695
#> AGCAGCCTCTTTAGTC AGCAGCCTCTTTAGTC        0.39225229            0.39225229
#> AGCATACCACAGCCCA AGCATACCACAGCCCA        0.17364001            0.17364001
#> AGCATACGTCCGCTGA AGCATACGTCCGCTGA        0.84795229            0.84795229
#> AGCGGTCGTACCGAGA AGCGGTCGTACCGAGA        0.36311571            0.36311571
#> AGCGGTCGTCGAATCT AGCGGTCGTCGAATCT        0.29017655            0.29017655
#> AGCGGTCGTGTATGGG AGCGGTCGTGTATGGG        0.40856859            0.40856859
#> AGCGTCGCACAGGCCT AGCGTCGCACAGGCCT        0.27609829            0.27609829
#> AGCTCCTCAGTTTACG AGCTCCTCAGTTTACG        0.04387964            0.04387964
#> AGCTCCTGTGTGCCTG AGCTCCTGTGTGCCTG        0.38731102            0.38731102
#> AGCTCCTTCAATACCG AGCTCCTTCAATACCG        0.14980052            0.14980052
#> AGCTCCTTCGCCCTTA AGCTCCTTCGCCCTTA        0.03727371            0.03727371
#> AGCTCTCAGACATAAC AGCTCTCAGACATAAC        0.14412494            0.14412494
#> AGCTCTCGTTCAGTAC AGCTCTCGTTCAGTAC        0.28174612            0.28174612
#> AGCTTGAAGACACGAC AGCTTGAAGACACGAC        0.13414127            0.13414127
#> AGCTTGAAGAGTGACC AGCTTGAAGAGTGACC        0.33579448            0.33579448
#> AGCTTGACAATAGAGT AGCTTGACAATAGAGT        0.38758578            0.38758578
#> AGCTTGACAGTAAGAT AGCTTGACAGTAAGAT        0.41962288            0.41962288
#> AGCTTGACATCTATGG AGCTTGACATCTATGG        0.23482479            0.23482479
#> AGCTTGAGTGCAGTAG AGCTTGAGTGCAGTAG        0.42184207            0.42184207
#> AGGCCACAGCTGTCTA AGGCCACAGCTGTCTA        0.38634513            0.38634513
#> AGGCCACCATGACATC AGGCCACCATGACATC        0.15693944            0.15693944
#> AGGCCGTCACGACGAA AGGCCGTCACGACGAA        0.28559593            0.28559593
#> AGGCCGTCACGGTGTC AGGCCGTCACGGTGTC        0.51336715            0.51336715
#> AGGCCGTGTACCAGTT AGGCCGTGTACCAGTT        0.44184087            0.44184087
#> AGGGAGTAGGGTATCG AGGGAGTAGGGTATCG        0.43461094            0.43461094
#> AGGGAGTAGTGTCTCA AGGGAGTAGTGTCTCA        0.40867855            0.40867855
#> AGGGAGTTCCTCAATT AGGGAGTTCCTCAATT        0.34442562            0.34442562
#> AGGGAGTTCGGAGGTA AGGGAGTTCGGAGGTA        0.14601727            0.14601727
#> AGGGATGGTCGACTGC AGGGATGGTCGACTGC        0.15489290            0.15489290
#> AGGGATGGTTAGGGTG AGGGATGGTTAGGGTG        0.28130316            0.28130316
#> AGGGATGTCGCATGAT AGGGATGTCGCATGAT        0.39751556            0.39751556
#> AGGGATGTCGCATGGC AGGGATGTCGCATGGC        0.16133233            0.16133233
#> AGGGTGAGTGCGAAAC AGGGTGAGTGCGAAAC        0.10079224            0.10079224
#> AGGTCATCACTTCTGC AGGTCATCACTTCTGC        0.16799782            0.16799782
#> AGGTCATGTCTCCCTA AGGTCATGTCTCCCTA        0.18896674            0.18896674
#> AGGTCATGTTTAGGAA AGGTCATGTTTAGGAA        0.09797720            0.09797720
#> AGTAGTCCACATGACT AGTAGTCCACATGACT        0.38834807            0.38834807
#> AGTAGTCGTGAGGGAG AGTAGTCGTGAGGGAG        0.20193413            0.20193413
#> AGTAGTCGTTCTGTTT AGTAGTCGTTCTGTTT        0.38934971            0.38934971
#> AGTCTTTCACCAGGTC AGTCTTTCACCAGGTC        0.45480288            0.45480288
#> AGTCTTTCACGGACAA AGTCTTTCACGGACAA        0.15389607            0.15389607
#> AGTCTTTGTTACGGAG AGTCTTTGTTACGGAG        0.28699499            0.28699499
#> AGTGAGGGTTACGGAG AGTGAGGGTTACGGAG        0.45408917            0.45408917
#> AGTGAGGTCTCCAACC AGTGAGGTCTCCAACC        0.31492673            0.31492673
#> AGTGGGACAGGTCGTC AGTGGGACAGGTCGTC        0.45506021            0.45506021
#> AGTGGGAGTGGCAAAC AGTGGGAGTGGCAAAC        0.11212241            0.11212241
#> AGTGTCATCCCTTGTG AGTGTCATCCCTTGTG        0.16340596            0.16340596
#> AGTTGGTAGAGCCTAG AGTTGGTAGAGCCTAG        0.19087481            0.19087481
#> AGTTGGTGTTGTCGCG AGTTGGTGTTGTCGCG        0.11472631            0.11472631
#> AGTTGGTTCACATAGC AGTTGGTTCACATAGC        0.20317528            0.20317528
#> ATAAGAGCATCGATTG ATAAGAGCATCGATTG        0.44075848            0.44075848
#> ATAGACCAGCGTGTCC ATAGACCAGCGTGTCC        0.15765591            0.15765591
#> ATAGACCAGTGATCGG ATAGACCAGTGATCGG        0.17150397            0.17150397
#> ATAGACCTCTAACTGG ATAGACCTCTAACTGG        0.39130284            0.39130284
#> ATCACGAGTGCACCAC ATCACGAGTGCACCAC        0.16799828            0.16799828
#> ATCATCTGTGTTTGTG ATCATCTGTGTTTGTG        0.27668729            0.27668729
#> ATCATCTTCTGTCAAG ATCATCTTCTGTCAAG        0.45629815            0.45629815
#> ATCATGGAGGAATGGA ATCATGGAGGAATGGA        0.14716813            0.14716813
#> ATCATGGGTATTCGTG ATCATGGGTATTCGTG        0.21361069            0.21361069
#> ATCATGGGTCCGCTGA ATCATGGGTCCGCTGA        0.21754489            0.21754489
#> ATCATGGTCCAGAGGA ATCATGGTCCAGAGGA        0.39942569            0.39942569
#> ATCCACCAGAGCCTAG ATCCACCAGAGCCTAG        0.63947551            0.63947551
#> ATCCACCCAGCTCGCA ATCCACCCAGCTCGCA        0.26158907            0.26158907
#> ATCCACCGTACCATCA ATCCACCGTACCATCA        0.42080429            0.42080429
#> ATCCGAAGTACACCGC ATCCGAAGTACACCGC        0.10940613            0.10940613
#> ATCCGAAGTGACTCAT ATCCGAAGTGACTCAT        0.38420459            0.38420459
#> ATCTGCCCAATCCGAT ATCTGCCCAATCCGAT        0.23035978            0.23035978
#> ATCTGCCCATGTTGAC ATCTGCCCATGTTGAC        0.15933504            0.15933504
#> ATGAGGGTCTTGACGA ATGAGGGTCTTGACGA        0.11202036            0.11202036
#> ATGTGTGCACAAGTAA ATGTGTGCACAAGTAA        0.28596812            0.28596812
#> ATGTGTGTCAAAGACA ATGTGTGTCAAAGACA        0.27184281            0.27184281
#> ATTACTCAGTTGCAGG ATTACTCAGTTGCAGG        0.45454817            0.45454817
#> ATTACTCTCAGGCGAA ATTACTCTCAGGCGAA        0.45196211            0.45196211
#> ATTACTCTCTGCCAGG ATTACTCTCTGCCAGG        0.19662352            0.19662352
#> ATTATCCAGAAACCTA ATTATCCAGAAACCTA        0.36488692            0.36488692
#> ATTATCCAGACCCACC ATTATCCAGACCCACC        0.36101405            0.36101405
#> ATTATCCAGTGAACAT ATTATCCAGTGAACAT        0.42854654            0.42854654
#> ATTGGACAGTAATCCC ATTGGACAGTAATCCC        0.38282126            0.38282126
#> ATTGGTGTCCAGTATG ATTGGTGTCCAGTATG        0.39554389            0.39554389
#> ATTTCTGGTAGTGAAT ATTTCTGGTAGTGAAT        0.39438459            0.39438459
#> CAACCAAAGCGTTCCG CAACCAAAGCGTTCCG        0.10499728            0.10499728
#> CAACCAACAACCGCCA CAACCAACAACCGCCA        0.44970718            0.44970718
#> CAACCAACACCTCGTT CAACCAACACCTCGTT        0.24728329            0.24728329
#> CAACCAATCGGGAGTA CAACCAATCGGGAGTA        0.41754283            0.41754283
#> CAAGAAAAGAGTACCG CAAGAAAAGAGTACCG        0.14149803            0.14149803
#> CAAGAAATCGTACCGG CAAGAAATCGTACCGG        0.39275880            0.39275880
#> CAAGGCCAGCTGTTCA CAAGGCCAGCTGTTCA        0.38626588            0.38626588
#> CAAGGCCAGTGTGGCA CAAGGCCAGTGTGGCA        0.39132258            0.39132258
#> CAAGGCCCACATTTCT CAAGGCCCACATTTCT        0.21257353            0.21257353
#> CAAGGCCTCAGTGCAT CAAGGCCTCAGTGCAT        0.16781306            0.16781306
#> CAAGGCCTCTCAACTT CAAGGCCTCTCAACTT        0.10834019            0.10834019
#> CAAGGCCTCTTACCTA CAAGGCCTCTTACCTA        0.43270116            0.43270116
#> CAAGTTGGTGTGCCTG CAAGTTGGTGTGCCTG        0.13927461            0.13927461
#> CAAGTTGTCTAACTGG CAAGTTGTCTAACTGG        0.08966655            0.08966655
#> CACAAACAGTAGCGGT CACAAACAGTAGCGGT        0.21117911            0.21117911
#> CACAAACCACGAAAGC CACAAACCACGAAAGC        0.24080431            0.24080431
#> CACACAACAGCCAGAA CACACAACAGCCAGAA        0.41462709            0.41462709
#> CACACCTGTACATGTC CACACCTGTACATGTC        0.19609395            0.19609395
#> CACACCTGTTATCCGA CACACCTGTTATCCGA        0.11357530            0.11357530
#> CACACTCAGATCCCAT CACACTCAGATCCCAT        0.43665157            0.43665157
#> CACACTCTCATGTGGT CACACTCTCATGTGGT        0.42765967            0.42765967
#> CACAGGCCAAGGCTCC CACAGGCCAAGGCTCC        0.12362352            0.12362352
#> CACAGGCCACGTCAGC CACAGGCCACGTCAGC        0.23176579            0.23176579
#> CACAGTAGTGCACTTA CACAGTAGTGCACTTA        0.16319923            0.16319923
#> CACAGTAGTTCCTCCA CACAGTAGTTCCTCCA        0.39401077            0.39401077
#> CACAGTATCCTCCTAG CACAGTATCCTCCTAG        0.39182690            0.39182690
#> CACATAGAGCGTGTCC CACATAGAGCGTGTCC        0.43613704            0.43613704
#> CACATTTGTTTCGCTC CACATTTGTTTCGCTC        0.39597932            0.39597932
#> CACATTTTCGAGGTAG CACATTTTCGAGGTAG        0.44634059            0.44634059
#> CACATTTTCGCAGGCT CACATTTTCGCAGGCT        0.38530069            0.38530069
#> CACCACTAGCGATCCC CACCACTAGCGATCCC        0.16382067            0.16382067
#> CACCACTGTGTGTGCC CACCACTGTGTGTGCC        0.11101797            0.11101797
#> CACCAGGCATAAAGGT CACCAGGCATAAAGGT        0.12870404            0.12870404
#> CACCTTGAGAGTAATC CACCTTGAGAGTAATC        0.39521731            0.39521731
#> CACTCCACATCGTCGG CACTCCACATCGTCGG        0.19030421            0.19030421
#> CACTCCATCAACACAC CACTCCATCAACACAC        0.52495605            0.52495605
#> CAGAATCCATTGGGCC CAGAATCCATTGGGCC        0.19097175            0.19097175
#> CAGAATCGTCATTAGC CAGAATCGTCATTAGC        0.33822488            0.33822488
#> CAGAATCGTGTAAGTA CAGAATCGTGTAAGTA        0.38969133            0.38969133
#> CAGAGAGTCGGCCGAT CAGAGAGTCGGCCGAT        0.21249825            0.21249825
#> CAGATCAAGCTGTTCA CAGATCAAGCTGTTCA        0.38496653            0.38496653
#> CAGCAGCAGACAAAGG CAGCAGCAGACAAAGG        0.43933972            0.43933972
#> CAGCAGCCACACTGCG CAGCAGCCACACTGCG        0.45447348            0.45447348
#> CAGCAGCGTTCACGGC CAGCAGCGTTCACGGC        0.39837701            0.39837701
#> CAGCCGAAGCGATATA CAGCCGAAGCGATATA        0.43931984            0.43931984
#> CAGCCGAAGTAGGTGC CAGCCGAAGTAGGTGC        0.45515611            0.45515611
#> CAGCCGATCCTATTCA CAGCCGATCCTATTCA        0.16381236            0.16381236
#> CAGCCGATCGTTTATC CAGCCGATCGTTTATC        0.18847844            0.18847844
#> CAGCGACCACAGGAGT CAGCGACCACAGGAGT        0.43569285            0.43569285
#> CAGCTAAAGTCCATAC CAGCTAAAGTCCATAC        0.43627397            0.43627397
#> CAGCTAACATGCCTAA CAGCTAACATGCCTAA        0.11013866            0.11013866
#> CAGCTAAGTACCGAGA CAGCTAAGTACCGAGA        0.20748258            0.20748258
#> CAGCTGGCAAGTTGTC CAGCTGGCAAGTTGTC        0.10483382            0.10483382
#> CAGCTGGGTCTAACGT CAGCTGGGTCTAACGT        0.08342274            0.08342274
#> CAGCTGGGTTCTCATT CAGCTGGGTTCTCATT        0.36790590            0.36790590
#> CAGGTGCAGCCAGGAT CAGGTGCAGCCAGGAT        0.14775203            0.14775203
#> CAGGTGCCATCAGTCA CAGGTGCCATCAGTCA        0.04375749            0.04375749
#> CAGTCCTAGGTGGGTT CAGTCCTAGGTGGGTT        0.19082729            0.19082729
#> CAGTCCTCAGGGCATA CAGTCCTCAGGGCATA        0.15456697            0.15456697
#> CATATGGCATGTTGAC CATATGGCATGTTGAC        0.41422087            0.41422087
#> CATATGGGTTACTGAC CATATGGGTTACTGAC        0.24786291            0.24786291
#> CATATGGTCACAGTAC CATATGGTCACAGTAC        0.10915658            0.10915658
#> CATATTCTCCTGTAGA CATATTCTCCTGTAGA        0.16582947            0.16582947
#> CATCAGACATGCATGT CATCAGACATGCATGT        0.41633659            0.41633659
#> CATCAGAGTCTCTTTA CATCAGAGTCTCTTTA        0.15260797            0.15260797
#> CATCAGATCATCGCTC CATCAGATCATCGCTC        0.09743510            0.09743510
#> CATCGAAGTACGAAAT CATCGAAGTACGAAAT        0.41991926            0.41991926
#> CATCGAAGTCACACGC CATCGAAGTCACACGC        0.25914204            0.25914204
#> CATCGGGCACGTTGGC CATCGGGCACGTTGGC        0.02535511            0.02535511
#> CATCGGGTCTTCGGTC CATCGGGTCTTCGGTC        0.25064301            0.25064301
#> CATGACACAACAACCT CATGACACAACAACCT        0.38367081            0.38367081
#> CATGACACATCGATTG CATGACACATCGATTG        0.32122268            0.32122268
#> CATGACAGTAGCTTGT CATGACAGTAGCTTGT        0.39170288            0.39170288
#> CATGACATCATGCATG CATGACATCATGCATG        0.13738801            0.13738801
#> CATGCCTCACCCATTC CATGCCTCACCCATTC        0.13846934            0.13846934
#> CATGGCGAGAGTCTGG CATGGCGAGAGTCTGG        0.26478476            0.26478476
#> CATTATCGTCATATCG CATTATCGTCATATCG        0.10082155            0.10082155
#> CCAATCCAGAGTAAGG CCAATCCAGAGTAAGG        0.36616402            0.36616402
#> CCACCTAAGACTTGAA CCACCTAAGACTTGAA        0.19931273            0.19931273
#> CCACCTACAGCTCGAC CCACCTACAGCTCGAC        0.86057123            0.86057123
#> CCACCTATCACCAGGC CCACCTATCACCAGGC        0.38895799            0.38895799
#> CCACCTATCTAAGCCA CCACCTATCTAAGCCA        0.17333462            0.17333462
#> CCACGGAAGATCCTGT CCACGGAAGATCCTGT        0.15521576            0.15521576
#> CCACGGAAGCAATCTC CCACGGAAGCAATCTC        0.20097867            0.20097867
#> CCACGGAAGGGTGTTG CCACGGAAGGGTGTTG        0.21434041            0.21434041
#> CCACGGAAGTCTCCTC CCACGGAAGTCTCCTC        0.11325276            0.11325276
#> CCACGGAGTCAACATC CCACGGAGTCAACATC        0.26489541            0.26489541
#> CCACGGAGTCTAGTCA CCACGGAGTCTAGTCA        0.60209526            0.60209526
#> CCACGGAGTTACGACT CCACGGAGTTACGACT        0.27430415            0.27430415
#> CCACTACAGAGCTGGT CCACTACAGAGCTGGT        0.10879556            0.10879556
#> CCACTACAGAGTACCG CCACTACAGAGTACCG        0.38542384            0.38542384
#> CCACTACCAACGCACC CCACTACCAACGCACC        0.15234229            0.15234229
#> CCAGCGACAGACGTAG CCAGCGACAGACGTAG        0.11290047            0.11290047
#> CCATGTCGTCGGCTCA CCATGTCGTCGGCTCA        0.21091389            0.21091389
#> CCATGTCTCCGCATAA CCATGTCTCCGCATAA        0.14130539            0.14130539
#> CCATTCGCATGCCTAA CCATTCGCATGCCTAA        0.45859322            0.45859322
#> CCATTCGGTAGGGACT CCATTCGGTAGGGACT        0.45115246            0.45115246
#> CCCAATCCAACTGGCC CCCAATCCAACTGGCC        0.17606903            0.17606903
#> CCCAGTTAGCATCATC CCCAGTTAGCATCATC        0.44302736            0.44302736
#> CCCAGTTGTCGTGGCT CCCAGTTGTCGTGGCT        0.17879644            0.17879644
#> CCCAGTTGTGCGCTTG CCCAGTTGTGCGCTTG        0.88537346            0.88537346
#> CCCAGTTTCAGCACAT CCCAGTTTCAGCACAT        0.18699051            0.18699051
#> CCCATACAGACTTTCG CCCATACAGACTTTCG        0.50886497            0.50886497
#> CCCATACCACATCCAA CCCATACCACATCCAA        0.41195150            0.41195150
#> CCCATACCATGCAACT CCCATACCATGCAACT        0.21303925            0.21303925
#> CCCATACTCGTGGACC CCCATACTCGTGGACC        0.10659256            0.10659256
#> CCGGGATAGATCTGAA CCGGGATAGATCTGAA        0.16910516            0.16910516
#> CCGGGATTCGACAGCC CCGGGATTCGACAGCC        0.17498718            0.17498718
#> CCGGTAGTCGGCTTGG CCGGTAGTCGGCTTGG        0.40105247            0.40105247
#> CCGTACTAGGAGTACC CCGTACTAGGAGTACC        0.27143976            0.27143976
#> CCGTACTAGTCCGTAT CCGTACTAGTCCGTAT        0.27149488            0.27149488
#> CCGTACTTCAATAAGG CCGTACTTCAATAAGG        0.45155534            0.45155534
#> CCGTGGACATAGAAAC CCGTGGACATAGAAAC        0.51510054            0.51510054
#> CCGTGGATCCAAGCCG CCGTGGATCCAAGCCG        0.42092419            0.42092419
#> CCGTGGATCTGTCTCG CCGTGGATCTGTCTCG        0.25490796            0.25490796
#> CCGTGGATCTTAGAGC CCGTGGATCTTAGAGC        0.14171407            0.14171407
#> CCGTTCAAGCTAAACA CCGTTCAAGCTAAACA        0.17148368            0.17148368
#> CCGTTCAAGGACCACA CCGTTCAAGGACCACA        0.19901805            0.19901805
#> CCGTTCACACATGACT CCGTTCACACATGACT        0.38311688            0.38311688
#> CCTACACAGATCCGAG CCTACACAGATCCGAG        0.38727673            0.38727673
#> CCTACACAGTACGTTC CCTACACAGTACGTTC        0.14961278            0.14961278
#> CCTACACTCCTTTACA CCTACACTCCTTTACA        0.42293489            0.42293489
#> CCTACACTCTTTACAC CCTACACTCTTTACAC        0.17175405            0.17175405
#> CCTACCACATCGATTG CCTACCACATCGATTG        0.17443291            0.17443291
#> CCTACCAGTTGCTCCT CCTACCAGTTGCTCCT        0.08865846            0.08865846
#> CCTACCATCTGCCCTA CCTACCATCTGCCCTA        0.28938104            0.28938104
#> CCTAGCTCAGACTCGC CCTAGCTCAGACTCGC        0.40059138            0.40059138
#> CCTAGCTCATACTACG CCTAGCTCATACTACG        0.43665162            0.43665162
#> CCTAGCTGTGGCCCTA CCTAGCTGTGGCCCTA        0.43148078            0.43148078
#> CCTAGCTTCGATGAGG CCTAGCTTCGATGAGG        0.43982770            0.43982770
#> CCTATTACAGACAAGC CCTATTACAGACAAGC        0.28924970            0.28924970
#> CCTATTATCATCGCTC CCTATTATCATCGCTC        0.21451750            0.21451750
#> CCTCAGTAGTGCGATG CCTCAGTAGTGCGATG        0.10762250            0.10762250
#> CCTCAGTCAGACACTT CCTCAGTCAGACACTT        0.15729527            0.15729527
#> CCTCAGTGTTACCAGT CCTCAGTGTTACCAGT        0.35971968            0.35971968
#> CCTCTGAGTATAAACG CCTCTGAGTATAAACG        0.41684830            0.41684830
#> CCTTACGAGATCACGG CCTTACGAGATCACGG        0.38927147            0.38927147
#> CCTTACGCACAACGCC CCTTACGCACAACGCC        0.10613478            0.10613478
#> CCTTCCCAGGCAGGTT CCTTCCCAGGCAGGTT        0.39357916            0.39357916
#> CCTTCCCCAAAGTGCG CCTTCCCCAAAGTGCG        0.39191318            0.39191318
#> CCTTCCCTCGTCGTTC CCTTCCCTCGTCGTTC        0.10294110            0.10294110
#> CCTTCGAAGATAGCAT CCTTCGAAGATAGCAT        0.22756508            0.22756508
#> CCTTCGACAGCTGTTA CCTTCGACAGCTGTTA        0.39066890            0.39066890
#> CCTTCGATCACAGGCC CCTTCGATCACAGGCC        0.15678939            0.15678939
#> CGAACATAGCCAGTAG CGAACATAGCCAGTAG        0.44549539            0.44549539
#> CGAACATCAGACGTAG CGAACATCAGACGTAG        0.39274996            0.39274996
#> CGAATGTAGTTCCACA CGAATGTAGTTCCACA        0.38335640            0.38335640
#> CGAATGTGTCTGGAGA CGAATGTGTCTGGAGA        0.13423160            0.13423160
#> CGAATGTTCACCGGGT CGAATGTTCACCGGGT        0.17211515            0.17211515
#> CGAATGTTCTTCATGT CGAATGTTCTTCATGT        0.23378519            0.23378519
#> CGACCTTAGCGCTTAT CGACCTTAGCGCTTAT        0.37395056            0.37395056
#> CGACCTTAGGCTAGAC CGACCTTAGGCTAGAC        0.46001335            0.46001335
#> CGACTTCCAGCCTATA CGACTTCCAGCCTATA        0.38295587            0.38295587
#> CGACTTCCAGCTGCAC CGACTTCCAGCTGCAC        0.10000748            0.10000748
#> CGACTTCTCACGCGGT CGACTTCTCACGCGGT        0.20093691            0.20093691
#> CGAGAAGAGGGTTCCC CGAGAAGAGGGTTCCC        0.26247090            0.26247090
#> CGAGAAGCACGCATCG CGAGAAGCACGCATCG        0.28578947            0.28578947
#> CGAGAAGTCGTTACAG CGAGAAGTCGTTACAG        0.14315433            0.14315433
#> CGAGCACAGCTGAAAT CGAGCACAGCTGAAAT        0.40639303            0.40639303
#> CGAGCACAGGCAGGTT CGAGCACAGGCAGGTT        0.10624165            0.10624165
#> CGAGCACGTAGTAGTA CGAGCACGTAGTAGTA        0.14000445            0.14000445
#> CGAGCACTCTACCTGC CGAGCACTCTACCTGC        0.39756417            0.39756417
#> CGAGCACTCTTCGAGA CGAGCACTCTTCGAGA        0.38290180            0.38290180
#> CGATCGGAGCTGCAAG CGATCGGAGCTGCAAG        0.15210360            0.15210360
#> CGATCGGCACGCCAGT CGATCGGCACGCCAGT        0.37407581            0.37407581
#> CGATGGCGTTACGGAG CGATGGCGTTACGGAG        0.40346043            0.40346043
#> CGATGGCTCAGTGTTG CGATGGCTCAGTGTTG        0.19560447            0.19560447
#> CGATGGCTCTGTGCAA CGATGGCTCTGTGCAA        0.43257288            0.43257288
#> CGATGTAAGTTTAGGA CGATGTAAGTTTAGGA        0.14565898            0.14565898
#> CGATGTAGTACTCTCC CGATGTAGTACTCTCC        0.18465529            0.18465529
#> CGATTGAAGTCAAGGC CGATTGAAGTCAAGGC        0.14365630            0.14365630
#> CGCCAAGCAAAGTCAA CGCCAAGCAAAGTCAA        0.17206416            0.17206416
#> CGCCAAGCACCGGAAA CGCCAAGCACCGGAAA        0.44137416            0.44137416
#> CGCCAAGTCTAAGCCA CGCCAAGTCTAAGCCA        0.19942505            0.19942505
#> CGCGGTAAGAGCTGGT CGCGGTAAGAGCTGGT        0.17131652            0.17131652
#> CGCGGTAAGCTGCAAG CGCGGTAAGCTGCAAG        0.18535857            0.18535857
#> CGCGTTTCAGTCTTCC CGCGTTTCAGTCTTCC        0.36492665            0.36492665
#> CGCTGGAGTCCATCCT CGCTGGAGTCCATCCT        0.38540256            0.38540256
#> CGCTTCATCTCTGTCG CGCTTCATCTCTGTCG        0.15733562            0.15733562
#> CGGACACAGTGTCTCA CGGACACAGTGTCTCA        0.40262863            0.40262863
#> CGGACGTAGACACGAC CGGACGTAGACACGAC        0.26802536            0.26802536
#> CGGACTGAGGTAGCTG CGGACTGAGGTAGCTG        0.39421945            0.39421945
#> CGGACTGTCTTGAGAC CGGACTGTCTTGAGAC        0.38143055            0.38143055
#> CGGAGCTAGATAGTCA CGGAGCTAGATAGTCA        0.25535231            0.25535231
#> CGGAGCTCATTGGGCC CGGAGCTCATTGGGCC        0.40925520            0.40925520
#> CGGAGTCAGTTTAGGA CGGAGTCAGTTTAGGA        0.19079149            0.19079149
#> CGGCTAGGTACTCGCG CGGCTAGGTACTCGCG        0.24057477            0.24057477
#> CGGGTCAAGGCGATAC CGGGTCAAGGCGATAC        0.14870864            0.14870864
#> CGGGTCAAGGCTCAGA CGGGTCAAGGCTCAGA        0.40063555            0.40063555
#> CGGTTAAAGAGCTGCA CGGTTAAAGAGCTGCA        0.44113838            0.44113838
#> CGGTTAAGTAAACACA CGGTTAAGTAAACACA        0.44191400            0.44191400
#> CGGTTAAGTCTAGCCG CGGTTAAGTCTAGCCG        0.26374139            0.26374139
#> CGGTTAATCACCGGGT CGGTTAATCACCGGGT        0.10640141            0.10640141
#> CGGTTAATCGTCTGAA CGGTTAATCGTCTGAA        0.38570422            0.38570422
#> CGTAGCGTCGGAAATA CGTAGCGTCGGAAATA        0.10733498            0.10733498
#> CGTAGGCAGCGTAGTG CGTAGGCAGCGTAGTG        0.39183257            0.39183257
#> CGTCACTAGGGTGTGT CGTCACTAGGGTGTGT        0.38756978            0.38756978
#> CGTCACTCACGAAATA CGTCACTCACGAAATA        0.38758347            0.38758347
#> CGTCACTCATAGTAAG CGTCACTCATAGTAAG        0.39153669            0.39153669
#> CGTCACTTCACTTATC CGTCACTTCACTTATC        0.10262502            0.10262502
#> CGTCACTTCTAACTCT CGTCACTTCTAACTCT        0.18155724            0.18155724
#> CGTCAGGAGCCTATGT CGTCAGGAGCCTATGT        0.39023483            0.39023483
#> CGTCAGGCATGGGAAC CGTCAGGCATGGGAAC        0.38908040            0.38908040
#> CGTCAGGGTAGGCATG CGTCAGGGTAGGCATG        0.20093297            0.20093297
#> CGTCAGGTCCGCAAGC CGTCAGGTCCGCAAGC        0.28444865            0.28444865
#> CGTCAGGTCCTGCCAT CGTCAGGTCCTGCCAT        0.37547378            0.37547378
#> CGTCCATGTCTAACGT CGTCCATGTCTAACGT        0.29334747            0.29334747
#> CGTCCATTCGCGCCAA CGTCCATTCGCGCCAA        0.20718311            0.20718311
#> CGTCTACAGATCTGCT CGTCTACAGATCTGCT        0.40608791            0.40608791
#> CGTCTACAGGGAGTAA CGTCTACAGGGAGTAA        0.14362472            0.14362472
#> CGTGAGCTCAGCGATT CGTGAGCTCAGCGATT        0.21816735            0.21816735
#> CGTGAGCTCTTCTGGC CGTGAGCTCTTCTGGC        0.19985774            0.19985774
#> CGTGTAAGTCTTGATG CGTGTAAGTCTTGATG        0.38100505            0.38100505
#> CGTGTAATCATTATCC CGTGTAATCATTATCC        0.27799411            0.27799411
#> CGTGTAATCCATGCTC CGTGTAATCCATGCTC        0.44356629            0.44356629
#> CGTGTCTCAGAAGCAC CGTGTCTCAGAAGCAC        0.14924889            0.14924889
#> CGTGTCTGTAGCGTCC CGTGTCTGTAGCGTCC        0.15448899            0.15448899
#> CGTTAGAAGTATCTCG CGTTAGAAGTATCTCG        0.42764114            0.42764114
#> CGTTAGAGTGCAACGA CGTTAGAGTGCAACGA        0.38897330            0.38897330
#> CGTTCTGAGGACGAAA CGTTCTGAGGACGAAA        0.13458628            0.13458628
#> CGTTCTGTCACAATGC CGTTCTGTCACAATGC        0.16610155            0.16610155
#> CGTTCTGTCGATAGAA CGTTCTGTCGATAGAA        0.73377004            0.73377004
#> CGTTGGGCAGGCGATA CGTTGGGCAGGCGATA        0.10283518            0.10283518
#> CTAACTTAGAGGTACC CTAACTTAGAGGTACC        0.16731071            0.16731071
#> CTAACTTAGATCACGG CTAACTTAGATCACGG        0.61662943            0.61662943
#> CTAACTTCATCTGGTA CTAACTTCATCTGGTA        0.18255059            0.18255059
#> CTAACTTCATTGGGCC CTAACTTCATTGGGCC        0.39332433            0.39332433
#> CTAAGACAGCACGCCT CTAAGACAGCACGCCT        0.24289725            0.24289725
#> CTAAGACAGCCACGCT CTAAGACAGCCACGCT        0.39150966            0.39150966
#> CTAAGACGTATCAGTC CTAAGACGTATCAGTC        0.21328646            0.21328646
#> CTAAGACTCTGCGGCA CTAAGACTCTGCGGCA        0.21903756            0.21903756
#> CTAATGGCATGGTTGT CTAATGGCATGGTTGT        0.09732986            0.09732986
#> CTACACCCATCAGTCA CTACACCCATCAGTCA        0.74023835            0.74023835
#> CTACACCGTCGAATCT CTACACCGTCGAATCT        0.41238947            0.41238947
#> CTACACCGTCGCCATG CTACACCGTCGCCATG        0.44883084            0.44883084
#> CTACATTCACCGTTGG CTACATTCACCGTTGG        0.10582470            0.10582470
#> CTACATTGTTACTGAC CTACATTGTTACTGAC        0.10342900            0.10342900
#> CTACCCAGTGCGCTTG CTACCCAGTGCGCTTG        0.40466587            0.40466587
#> CTACGTCAGCTAGGCA CTACGTCAGCTAGGCA        0.12083497            0.12083497
#> CTAGAGTCACTGTTAG CTAGAGTCACTGTTAG        0.43827926            0.43827926
#> CTAGCCTGTCTCTTAT CTAGCCTGTCTCTTAT        0.41016150            0.41016150
#> CTAGCCTGTTCGGGCT CTAGCCTGTTCGGGCT        0.45322098            0.45322098
#> CTAGCCTTCAGCGATT CTAGCCTTCAGCGATT        0.39584966            0.39584966
#> CTAGTGATCAAAGACA CTAGTGATCAAAGACA        0.43686630            0.43686630
#> CTAGTGATCCCAGGTG CTAGTGATCCCAGGTG        0.33319706            0.33319706
#> CTCACACGTCACACGC CTCACACGTCACACGC        0.21149168            0.21149168
#> CTCAGAACACTAGTAC CTCAGAACACTAGTAC        0.33817820            0.33817820
#> CTCAGAACATAGGATA CTCAGAACATAGGATA        0.20532201            0.20532201
#> CTCAGAAGTCCTCTTG CTCAGAAGTCCTCTTG        0.39473268            0.39473268
#> CTCAGAATCAGGTTCA CTCAGAATCAGGTTCA        0.38779159            0.38779159
#> CTCATTACAAGCGCTC CTCATTACAAGCGCTC        0.44699186            0.44699186
#> CTCCTAGCACAACGTT CTCCTAGCACAACGTT        0.11721684            0.11721684
#> CTCGAAAGTAAGAGGA CTCGAAAGTAAGAGGA        0.33187219            0.33187219
#> CTCGAAAGTGGTGTAG CTCGAAAGTGGTGTAG        0.13974008            0.13974008
#> CTCGAAATCCTGCCAT CTCGAAATCCTGCCAT        0.20598717            0.20598717
#> CTCGAGGCAAGTCTGT CTCGAGGCAAGTCTGT        0.20120653            0.20120653
#> CTCGAGGGTAAGTAGT CTCGAGGGTAAGTAGT        0.14643462            0.14643462
#> CTCGAGGTCTGTTGAG CTCGAGGTCTGTTGAG        0.36104219            0.36104219
#> CTCGGAGGTCATTAGC CTCGGAGGTCATTAGC        0.02950687            0.02950687
#> CTCGGAGGTGAAGGCT CTCGGAGGTGAAGGCT        0.10693042            0.10693042
#> CTCGTACCAATGGATA CTCGTACCAATGGATA        0.38826366            0.38826366
#> CTCGTACTCTTTAGTC CTCGTACTCTTTAGTC        0.16763118            0.16763118
#> CTCGTCAAGAACAATC CTCGTCAAGAACAATC        0.16260073            0.16260073
#> CTCGTCACAACTGCTA CTCGTCACAACTGCTA        0.39413289            0.39413289
#> CTCGTCATCGGTCTAA CTCGTCATCGGTCTAA        0.39087063            0.39087063
#> CTCTAATAGTTATCGC CTCTAATAGTTATCGC        0.16754348            0.16754348
#> CTCTGGTAGCAATCTC CTCTGGTAGCAATCTC        0.25216333            0.25216333
#> CTCTGGTAGTCCCACG CTCTGGTAGTCCCACG        0.44194697            0.44194697
#> CTCTGGTCAACACCTA CTCTGGTCAACACCTA        0.38433302            0.38433302
#> CTCTGGTCAAGACACG CTCTGGTCAAGACACG        0.80871234            0.80871234
#> CTCTGGTCAAGTTGTC CTCTGGTCAAGTTGTC        0.12298779            0.12298779
#> CTCTGGTCAGTTCATG CTCTGGTCAGTTCATG        0.39134432            0.39134432
#> CTGAAACAGACTAGGC CTGAAACAGACTAGGC        0.39709562            0.39709562
#> CTGAAACGTAGCGCAA CTGAAACGTAGCGCAA        0.19909417            0.19909417
#> CTGAAACGTATGGTTC CTGAAACGTATGGTTC        0.28312950            0.28312950
#> CTGAAACTCGATAGAA CTGAAACTCGATAGAA        0.27627094            0.27627094
#> CTGAAACTCGCCGTGA CTGAAACTCGCCGTGA        0.40389703            0.40389703
#> CTGAAACTCTGTTTGT CTGAAACTCTGTTTGT        0.27407964            0.27407964
#> CTGAAGTGTTTGGGCC CTGAAGTGTTTGGGCC        0.19075388            0.19075388
#> CTGAAGTTCTGCGGCA CTGAAGTTCTGCGGCA        0.13008644            0.13008644
#> CTGATCCTCAGTTGAC CTGATCCTCAGTTGAC        0.39909047            0.39909047
#> CTGCCTAGTCTCCATC CTGCCTAGTCTCCATC        0.42133692            0.42133692
#> CTGCGGACAGACTCGC CTGCGGACAGACTCGC        0.21097328            0.21097328
#> CTGCGGACATTAGCCA CTGCGGACATTAGCCA        0.11598431            0.11598431
#> CTGCGGATCAGCATGT CTGCGGATCAGCATGT        0.43926965            0.43926965
#> CTGCTGTAGACTGGGT CTGCTGTAGACTGGGT        0.09388167            0.09388167
#> CTGCTGTCAAAGCAAT CTGCTGTCAAAGCAAT        0.14894398            0.14894398
#> CTGGTCTAGGATTCGG CTGGTCTAGGATTCGG        0.10931741            0.10931741
#> CTGGTCTAGGTGATTA CTGGTCTAGGTGATTA        0.17788128            0.17788128
#> CTGGTCTCAGGAATGC CTGGTCTCAGGAATGC        0.39178790            0.39178790
#> CTGGTCTTCGCATGGC CTGGTCTTCGCATGGC        0.24440935            0.24440935
#> CTGTGCTCATTACCTT CTGTGCTCATTACCTT        0.33476455            0.33476455
#> CTTAACTCAAGCGTAG CTTAACTCAAGCGTAG        0.39292442            0.39292442
#> CTTAACTTCAGCTTAG CTTAACTTCAGCTTAG        0.09667512            0.09667512
#> CTTAGGACACAGGTTT CTTAGGACACAGGTTT        0.40480888            0.40480888
#> CTTAGGAGTCCGAATT CTTAGGAGTCCGAATT        0.38360750            0.38360750
#> CTTAGGATCAACGGGA CTTAGGATCAACGGGA        0.26708023            0.26708023
#> CTTAGGATCACTCTTA CTTAGGATCACTCTTA        0.23645547            0.23645547
#> CTTCTCTAGGCAGGTT CTTCTCTAGGCAGGTT        0.42403270            0.42403270
#> CTTCTCTAGTCCGTAT CTTCTCTAGTCCGTAT        0.14688626            0.14688626
#> CTTCTCTCAGCATACT CTTCTCTCAGCATACT        0.21407507            0.21407507
#> CTTCTCTGTCAGAGGT CTTCTCTGTCAGAGGT        0.15847705            0.15847705
#> CTTGGCTAGACAGACC CTTGGCTAGACAGACC        0.46122434            0.46122434
#> CTTGGCTAGGACACCA CTTGGCTAGGACACCA        0.38768350            0.38768350
#> CTTGGCTCAAAGGAAG CTTGGCTCAAAGGAAG        0.10076573            0.10076573
#> CTTGGCTCAGCCACCA CTTGGCTCAGCCACCA        0.23336847            0.23336847
#> CTTGGCTCATCCCACT CTTGGCTCATCCCACT        0.37941518            0.37941518
#> CTTGGCTGTTAGAACA CTTGGCTGTTAGAACA        0.46561164            0.46561164
#> CTTTGCGAGACCACGA CTTTGCGAGACCACGA        0.10861156            0.10861156
#> CTTTGCGAGGAATGGA CTTTGCGAGGAATGGA        0.38535422            0.38535422
#> CTTTGCGTCCACGTGG CTTTGCGTCCACGTGG        0.78772205            0.78772205
#> GAAACTCAGAAACCTA GAAACTCAGAAACCTA        0.37308643            0.37308643
#> GAAATGAAGCTCCTCT GAAATGAAGCTCCTCT        0.26352895            0.26352895
#> GAAATGACAACTGGCC GAAATGACAACTGGCC        0.39529015            0.39529015
#> GAAATGACACATGACT GAAATGACACATGACT        0.39445586            0.39445586
#> GAAATGAGTCACCCAG GAAATGAGTCACCCAG        0.31739173            0.31739173
#> GAAATGAGTCGGATCC GAAATGAGTCGGATCC        0.39104090            0.39104090
#> GAACATCGTAATCGTC GAACATCGTAATCGTC        0.40693170            0.40693170
#> GAACCTAAGTAGATGT GAACCTAAGTAGATGT        0.14285272            0.14285272
#> GAACCTATCACTTATC GAACCTATCACTTATC        0.19476563            0.19476563
#> GAACGGAAGTCAAGCG GAACGGAAGTCAAGCG        0.16641208            0.16641208
#> GAACGGACAAGGACAC GAACGGACAAGGACAC        0.21850137            0.21850137
#> GAAGCAGGTCGCGAAA GAAGCAGGTCGCGAAA        0.10441518            0.10441518
#> GAAGCAGGTCTAGAGG GAAGCAGGTCTAGAGG        0.08991317            0.08991317
#> GAAGCAGTCAAGGCTT GAAGCAGTCAAGGCTT        0.39344472            0.39344472
#> GAATGAAAGAATGTGT GAATGAAAGAATGTGT        0.14696401            0.14696401
#> GAATGAACAATCGAAA GAATGAACAATCGAAA        0.39268966            0.39268966
#> GAATGAACAGCCTTGG GAATGAACAGCCTTGG        0.45266026            0.45266026
#> GAATGAATCATCGGAT GAATGAATCATCGGAT        0.18446878            0.18446878
#> GACACGCCAAAGGAAG GACACGCCAAAGGAAG        0.39198009            0.39198009
#> GACCAATAGAATGTGT GACCAATAGAATGTGT        0.75141631            0.75141631
#> GACCAATCAGACGCCT GACCAATCAGACGCCT        0.14830944            0.14830944
#> GACCAATGTTATCCGA GACCAATGTTATCCGA        0.18024428            0.18024428
#> GACCAATGTTGGACCC GACCAATGTTGGACCC        0.29567775            0.29567775
#> GACCAATTCAACGGCC GACCAATTCAACGGCC        0.25913059            0.25913059
#> GACCTGGGTTCGTGAT GACCTGGGTTCGTGAT        0.16347853            0.16347853
#> GACGCGTCAGTCGTGC GACGCGTCAGTCGTGC        0.38684445            0.38684445
#> GACGCGTGTCGGCTCA GACGCGTGTCGGCTCA        0.26580596            0.26580596
#> GACGGCTAGGTGATAT GACGGCTAGGTGATAT        0.39110364            0.39110364
#> GACGGCTTCAACTCTT GACGGCTTCAACTCTT        0.05536723            0.05536723
#> GACGTGCGTCACCTAA GACGTGCGTCACCTAA        0.13107786            0.13107786
#> GACGTGCTCCGGCACA GACGTGCTCCGGCACA        0.21286188            0.21286188
#> GACGTTAAGACGCAAC GACGTTAAGACGCAAC        0.38721340            0.38721340
#> GACGTTAAGACGCTTT GACGTTAAGACGCTTT        0.26844904            0.26844904
#> GACGTTACATCTCGCT GACGTTACATCTCGCT        0.92366047            0.92366047
#> GACGTTAGTGTGAATA GACGTTAGTGTGAATA        0.10725862            0.10725862
#> GACTAACGTCATCCCT GACTAACGTCATCCCT        0.17150824            0.17150824
#> GACTAACTCAGTCCCT GACTAACTCAGTCCCT        0.43462943            0.43462943
#> GACTACAGTTCGTGAT GACTACAGTTCGTGAT        0.39115282            0.39115282
#> GACTGCGAGAATCTCC GACTGCGAGAATCTCC        0.14457882            0.14457882
#> GAGCAGAAGCGTGTCC GAGCAGAAGCGTGTCC        0.43234641            0.43234641
#> GAGCAGAGTTTGCATG GAGCAGAGTTTGCATG        0.21642591            0.21642591
#> GAGCAGATCAACTCTT GAGCAGATCAACTCTT        0.12332929            0.12332929
#> GAGGTGAGTCCGCTGA GAGGTGAGTCCGCTGA        0.10962166            0.10962166
#> GAGGTGAGTGACAAAT GAGGTGAGTGACAAAT        0.38612988            0.38612988
#> GAGTCCGAGCGATGAC GAGTCCGAGCGATGAC        0.27263970            0.27263970
#> GAGTCCGAGTTAAGTG GAGTCCGAGTTAAGTG        0.39044182            0.39044182
#> GAGTCCGCAATTGCTG GAGTCCGCAATTGCTG        0.40984019            0.40984019
#> GATCAGTAGGTGATAT GATCAGTAGGTGATAT        0.13703864            0.13703864
#> GATCGATAGAAAGTGG GATCGATAGAAAGTGG        0.12490817            0.12490817
#> GATCGATAGGATGGTC GATCGATAGGATGGTC        0.17599447            0.17599447
#> GATCGATGTAGCGATG GATCGATGTAGCGATG        0.24702510            0.24702510
#> GATCGATGTTGTTTGG GATCGATGTTGTTTGG        0.17636884            0.17636884
#> GATCGCGCACCAGGCT GATCGCGCACCAGGCT        0.10475098            0.10475098
#> GATCGTAAGTAGGTGC GATCGTAAGTAGGTGC        0.12418242            0.12418242
#> GATCGTAGTCCTAGCG GATCGTAGTCCTAGCG        0.44667149            0.44667149
#> GATCTAGAGCTAGCCC GATCTAGAGCTAGCCC        0.34550296            0.34550296
#> GATCTAGGTGCCTTGG GATCTAGGTGCCTTGG        0.38902471            0.38902471
#> GATGAAAAGTTGTAGA GATGAAAAGTTGTAGA        0.45152645            0.45152645
#> GATGAAACAGGCTGAA GATGAAACAGGCTGAA        0.15918134            0.15918134
#> GATGAAACATGGTTGT GATGAAACATGGTTGT        0.38674339            0.38674339
#> GATGAGGAGTGGACGT GATGAGGAGTGGACGT        0.29218590            0.29218590
#> GATGAGGTCTGGGCCA GATGAGGTCTGGGCCA        0.38531049            0.38531049
#> GATGCTAGTAGAAGGA GATGCTAGTAGAAGGA        0.43973166            0.43973166
#> GATGCTAGTGTTAAGA GATGCTAGTGTTAAGA        0.12059156            0.12059156
#> GATTCAGAGCGATGAC GATTCAGAGCGATGAC        0.13310637            0.13310637
#> GATTCAGAGTACGTTC GATTCAGAGTACGTTC        0.20187020            0.20187020
#> GATTCAGTCGCCAGCA GATTCAGTCGCCAGCA        0.14721423            0.14721423
#> GATTCAGTCGTCCAGG GATTCAGTCGTCCAGG        0.26044843            0.26044843
#> GCAAACTCAGACAAAT GCAAACTCAGACAAAT        0.21679054            0.21679054
#> GCAATCACAGCTGCTG GCAATCACAGCTGCTG        0.13654784            0.13654784
#> GCAATCATCCAACCAA GCAATCATCCAACCAA        0.38940470            0.38940470
#> GCACTCTCAGCGTCCA GCACTCTCAGCGTCCA        0.18574978            0.18574978
#> GCAGCCACAGCGTTCG GCAGCCACAGCGTTCG        0.33608162            0.33608162
#> GCAGTTACATTTCACT GCAGTTACATTTCACT        0.10443508            0.10443508
#> GCAGTTATCTTGAGAC GCAGTTATCTTGAGAC        0.38543202            0.38543202
#> GCATGATAGCCTTGAT GCATGATAGCCTTGAT        0.36697029            0.36697029
#> GCATGATCAAGCGAGT GCATGATCAAGCGAGT        0.40945591            0.40945591
#> GCATGATGTCACTGGC GCATGATGTCACTGGC        0.13057757            0.13057757
#> GCATGATGTTTAGGAA GCATGATGTTTAGGAA        0.18187454            0.18187454
#> GCATGATTCAAACGGG GCATGATTCAAACGGG        0.38933894            0.38933894
#> GCATGCGAGGAACTGC GCATGCGAGGAACTGC        0.00000000            0.00000000
#> GCATGCGCAGCTCCGA GCATGCGCAGCTCCGA        0.20753117            0.20753117
#> GCATGTAAGAATCTCC GCATGTAAGAATCTCC        0.51809823            0.51809823
#> GCATGTAAGCGTGTCC GCATGTAAGCGTGTCC        0.19237312            0.19237312
#> GCATGTAGTTGCGCAC GCATGTAGTTGCGCAC        0.16485105            0.16485105
#> GCATGTATCGTGGACC GCATGTATCGTGGACC        0.17182146            0.17182146
#> GCCTCTAAGCTAGTCT GCCTCTAAGCTAGTCT        0.26113891            0.26113891
#> GCCTCTAAGTCCGGTC GCCTCTAAGTCCGGTC        0.35535096            0.35535096
#> GCCTCTATCGCTTAGA GCCTCTATCGCTTAGA        0.38575548            0.38575548
#> GCGACCATCACCAGGC GCGACCATCACCAGGC        0.17001781            0.17001781
#> GCGAGAAGTGGGTATG GCGAGAAGTGGGTATG        0.19838409            0.19838409
#> GCGCAACGTTTCGCTC GCGCAACGTTTCGCTC        0.22273586            0.22273586
#> GCGCAACTCACAACGT GCGCAACTCACAACGT        0.40743415            0.40743415
#> GCGCAGTGTGCTAGCC GCGCAGTGTGCTAGCC        0.14814288            0.14814288
#> GCGCCAAAGCTCCTTC GCGCCAAAGCTCCTTC        0.19070409            0.19070409
#> GCGCCAATCTTATCTG GCGCCAATCTTATCTG        0.08424236            0.08424236
#> GCGCGATGTCCGAATT GCGCGATGTCCGAATT        0.97201795            0.97201795
#> GCGGGTTAGACTTGAA GCGGGTTAGACTTGAA        0.34839977            0.34839977
#> GCGGGTTCAAGTAATG GCGGGTTCAAGTAATG        0.19721160            0.19721160
#> GCGGGTTGTTACTGAC GCGGGTTGTTACTGAC        0.38917689            0.38917689
#> GCTCCTAAGGGCATGT GCTCCTAAGGGCATGT        0.56556551            0.56556551
#> GCTCTGTGTCCATGAT GCTCTGTGTCCATGAT        0.11997197            0.11997197
#> GCTGCAGGTTACCGAT GCTGCAGGTTACCGAT        0.39140283            0.39140283
#> GCTGCGAAGAACAACT GCTGCGAAGAACAACT        0.16356166            0.16356166
#> GCTGCGAAGCTAGCCC GCTGCGAAGCTAGCCC        0.35378184            0.35378184
#> GCTGCGAAGTATTGGA GCTGCGAAGTATTGGA        0.20099111            0.20099111
#> GCTGCGATCCCGGATG GCTGCGATCCCGGATG        0.11162380            0.11162380
#> GCTGCGATCGGCTACG GCTGCGATCGGCTACG        0.85027740            0.85027740
#> GCTGCTTGTCTTCAAG GCTGCTTGTCTTCAAG        0.43595305            0.43595305
#> GCTGGGTGTCATACTG GCTGGGTGTCATACTG        0.21546991            0.21546991
#> GCTGGGTGTCTAGAGG GCTGGGTGTCTAGAGG        0.29176478            0.29176478
#> GCTTCCAAGACATAAC GCTTCCAAGACATAAC        0.20356093            0.20356093
#> GCTTCCAAGTTAAGTG GCTTCCAAGTTAAGTG        0.16384216            0.16384216
#> GCTTGAAGTCAACTGT GCTTGAAGTCAACTGT        0.20755836            0.20755836
#> GGAAAGCCAATGTAAG GGAAAGCCAATGTAAG        0.14573757            0.14573757
#> GGAAAGCGTCTTGATG GGAAAGCGTCTTGATG        0.42138773            0.42138773
#> GGAATAAGTAAGAGAG GGAATAAGTAAGAGAG        0.26744362            0.26744362
#> GGAATAAGTTCGGGCT GGAATAAGTTCGGGCT        0.41501037            0.41501037
#> GGACAAGGTGATAAGT GGACAAGGTGATAAGT        0.39692642            0.39692642
#> GGACAAGTCAAAGACA GGACAAGTCAAAGACA        0.17260807            0.17260807
#> GGACAGAAGGCCGAAT GGACAGAAGGCCGAAT        0.33615025            0.33615025
#> GGACAGACAATCCAAC GGACAGACAATCCAAC        0.13726474            0.13726474
#> GGACAGATCAGGCGAA GGACAGATCAGGCGAA        0.69159319            0.69159319
#> GGACAGATCCTGCAGG GGACAGATCCTGCAGG        0.44659231            0.44659231
#> GGACAGATCCTTTCGG GGACAGATCCTTTCGG        0.20539649            0.20539649
#> GGACATTGTCTCCCTA GGACATTGTCTCCCTA        0.03703934            0.03703934
#> GGACGTCAGGGTGTTG GGACGTCAGGGTGTTG        0.14056809            0.14056809
#> GGACGTCCACGCATCG GGACGTCCACGCATCG        0.39371104            0.39371104
#> GGACGTCCATGGATGG GGACGTCCATGGATGG        0.18849251            0.18849251
#> GGACGTCGTCAGAAGC GGACGTCGTCAGAAGC        0.11046225            0.11046225
#> GGACGTCGTCGAGATG GGACGTCGTCGAGATG        0.16231080            0.16231080
#> GGACGTCGTCTAGGTT GGACGTCGTCTAGGTT        0.44046786            0.44046786
#> GGACGTCGTTTAGCTG GGACGTCGTTTAGCTG        0.45337700            0.45337700
#> GGAGCAAAGCTTCGCG GGAGCAAAGCTTCGCG        0.41106212            0.41106212
#> GGAGCAAGTATCACCA GGAGCAAGTATCACCA        0.10745723            0.10745723
#> GGATGTTAGTAAGTAC GGATGTTAGTAAGTAC        0.16983348            0.16983348
#> GGATGTTTCCAGTATG GGATGTTTCCAGTATG        0.23709567            0.23709567
#> GGCAATTCACGGCTAC GGCAATTCACGGCTAC        0.17165429            0.17165429
#> GGCAATTCATCTCCCA GGCAATTCATCTCCCA        0.21535518            0.21535518
#> GGCAATTGTAGAGCTG GGCAATTGTAGAGCTG        0.87292687            0.87292687
#> GGCCGATTCACTCTTA GGCCGATTCACTCTTA        0.38544021            0.38544021
#> GGCGACTGTAAGTTCC GGCGACTGTAAGTTCC        0.22988875            0.22988875
#> GGCGTGTCACCAGGTC GGCGTGTCACCAGGTC        0.14187285            0.14187285
#> GGCGTGTCACCGATAT GGCGTGTCACCGATAT        0.11389201            0.11389201
#> GGCGTGTTCAGCTTAG GGCGTGTTCAGCTTAG        0.16682853            0.16682853
#> GGCGTGTTCGAGAGCA GGCGTGTTCGAGAGCA        0.12036731            0.12036731
#> GGCTCGAAGATCACGG GGCTCGAAGATCACGG        0.12386101            0.12386101
#> GGCTCGAAGTGACTCT GGCTCGAAGTGACTCT        0.30359542            0.30359542
#> GGCTCGACAGCTTCGG GGCTCGACAGCTTCGG        0.54852550            0.54852550
#> GGCTGGTAGTCACGCC GGCTGGTAGTCACGCC        0.19268553            0.19268553
#> GGCTGGTTCATGCATG GGCTGGTTCATGCATG        0.18545531            0.18545531
#> GGGAATGAGCACACAG GGGAATGAGCACACAG        0.15773709            0.15773709
#> GGGAATGAGTGTACGG GGGAATGAGTGTACGG        0.39592821            0.39592821
#> GGGAATGCAAAGTGCG GGGAATGCAAAGTGCG        0.19805310            0.19805310
#> GGGAATGTCACAAACC GGGAATGTCACAAACC        0.22027971            0.22027971
#> GGGACCTCACTGTCGG GGGACCTCACTGTCGG        0.16131936            0.16131936
#> GGGACCTTCGAATGCT GGGACCTTCGAATGCT        0.39326663            0.39326663
#> GGGAGATAGATGTGGC GGGAGATAGATGTGGC        0.21781342            0.21781342
#> GGGAGATAGCGAAGGG GGGAGATAGCGAAGGG        0.26358587            0.26358587
#> GGGAGATAGTGACATA GGGAGATAGTGACATA        0.20762840            0.20762840
#> GGGATGAAGGGCATGT GGGATGAAGGGCATGT        0.18617771            0.18617771
#> GGGATGAGTCTGCGGT GGGATGAGTCTGCGGT        0.20165154            0.20165154
#> GGGATGAGTTGTGGAG GGGATGAGTTGTGGAG        0.29110129            0.29110129
#> GGGCACTCATTTCAGG GGGCACTCATTTCAGG        0.39030915            0.39030915
#> GGGCACTGTCTCCATC GGGCACTGTCTCCATC        0.10652323            0.10652323
#> GGGCACTGTGGCAAAC GGGCACTGTGGCAAAC        0.39804445            0.39804445
#> GGGCATCCATGACATC GGGCATCCATGACATC        0.41516777            0.41516777
#> GGGCATCCATTTCAGG GGGCATCCATTTCAGG        0.18166048            0.18166048
#> GGGCATCGTGGTCTCG GGGCATCGTGGTCTCG        0.32194573            0.32194573
#> GGGTCTGCACAGGCCT GGGTCTGCACAGGCCT        0.34345321            0.34345321
#> GGGTCTGCATGAAGTA GGGTCTGCATGAAGTA        0.38297955            0.38297955
#> GGGTCTGTCATGTGGT GGGTCTGTCATGTGGT        0.12041901            0.12041901
#> GGGTCTGTCCATGAAC GGGTCTGTCCATGAAC        0.19228655            0.19228655
#> GGGTTGCCAAGCTGAG GGGTTGCCAAGCTGAG        0.39308517            0.39308517
#> GGGTTGCTCGAACGGA GGGTTGCTCGAACGGA        0.39141219            0.39141219
#> GGTATTGGTACTCGCG GGTATTGGTACTCGCG        0.34117291            0.34117291
#> GGTATTGGTCCCGACA GGTATTGGTCCCGACA        0.39698291            0.39698291
#> GGTGTTAAGGATCGCA GGTGTTAAGGATCGCA        0.38277008            0.38277008
#> GGTGTTAGTTCCTCCA GGTGTTAGTTCCTCCA        0.39159064            0.39159064
#> GTAACGTGTCACTGGC GTAACGTGTCACTGGC        0.38585190            0.38585190
#> GTAACTGGTCACTGGC GTAACTGGTCACTGGC        0.10007944            0.10007944
#> GTAACTGTCAGTGCAT GTAACTGTCAGTGCAT        0.43994566            0.43994566
#> GTACGTAGTCTCTTTA GTACGTAGTCTCTTTA        0.18129975            0.18129975
#> GTACTCCCACAGGTTT GTACTCCCACAGGTTT        0.08042772            0.08042772
#> GTACTCCCACGTCTCT GTACTCCCACGTCTCT        0.39280734            0.39280734
#> GTACTCCTCCAACCAA GTACTCCTCCAACCAA        0.28553724            0.28553724
#> GTACTTTGTACCGTAT GTACTTTGTACCGTAT        0.27503444            0.27503444
#> GTACTTTTCTTCGGTC GTACTTTTCTTCGGTC        0.15789618            0.15789618
#> GTAGGCCAGCTAACAA GTAGGCCAGCTAACAA        0.13690831            0.13690831
#> GTAGGCCGTTACCAGT GTAGGCCGTTACCAGT        0.65289548            0.65289548
#> GTAGTCACATCGATGT GTAGTCACATCGATGT        0.20243069            0.20243069
#> GTAGTCAGTCCGACGT GTAGTCAGTCCGACGT        0.15806390            0.15806390
#> GTATCTTCAGTGGGAT GTATCTTCAGTGGGAT        0.38627944            0.38627944
#> GTATCTTCATCGATGT GTATCTTCATCGATGT        0.34635683            0.34635683
#> GTATCTTGTACTCAAC GTATCTTGTACTCAAC        0.18240146            0.18240146
#> GTATCTTGTCGCGGTT GTATCTTGTCGCGGTT        0.27275702            0.27275702
#> GTATCTTGTGGCCCTA GTATCTTGTGGCCCTA        0.30960267            0.30960267
#> GTATCTTTCGTGGACC GTATCTTTCGTGGACC        0.10623687            0.10623687
#> GTATTCTCACATTAGC GTATTCTCACATTAGC        0.39060082            0.39060082
#> GTATTCTCACTTGGAT GTATTCTCACTTGGAT        0.45915634            0.45915634
#> GTATTCTGTCCGTCAG GTATTCTGTCCGTCAG        0.45737900            0.45737900
#> GTATTCTTCTCTAAGG GTATTCTTCTCTAAGG        0.26580831            0.26580831
#> GTCAAGTGTGGTCCGT GTCAAGTGTGGTCCGT        0.43486600            0.43486600
#> GTCACAACACGTCTCT GTCACAACACGTCTCT        0.22150469            0.22150469
#> GTCACGGAGCTATGCT GTCACGGAGCTATGCT        0.18288962            0.18288962
#> GTCACGGCATCCGTGG GTCACGGCATCCGTGG        0.09968496            0.09968496
#> GTCACGGTCTTGCCGT GTCACGGTCTTGCCGT        0.11614353            0.11614353
#> GTCATTTGTCGTCTTC GTCATTTGTCGTCTTC        0.19590647            0.19590647
#> GTCCTCACACCACGTG GTCCTCACACCACGTG        0.17647205            0.17647205
#> GTCCTCATCATTTGGG GTCCTCATCATTTGGG        0.37149983            0.37149983
#> GTCCTCATCTGGAGCC GTCCTCATCTGGAGCC        0.30947466            0.30947466
#> GTCGGGTCAGGCTGAA GTCGGGTCAGGCTGAA        0.44603416            0.44603416
#> GTCGTAAAGTTTAGGA GTCGTAAAGTTTAGGA        0.22018410            0.22018410
#> GTCTCGTAGGCTAGAC GTCTCGTAGGCTAGAC        0.36828369            0.36828369
#> GTCTCGTCATAGTAAG GTCTCGTCATAGTAAG        0.29257818            0.29257818
#> GTCTCGTTCACTTCAT GTCTCGTTCACTTCAT        0.26083130            0.26083130
#> GTCTCGTTCTGCGTAA GTCTCGTTCTGCGTAA        0.56622053            0.56622053
#> GTCTTCGTCCGTTGCT GTCTTCGTCCGTTGCT        0.19917424            0.19917424
#> GTCTTCGTCTACCAGA GTCTTCGTCTACCAGA        0.06474840            0.06474840
#> GTGCATAAGAACAATC GTGCATAAGAACAATC        0.12623075            0.12623075
#> GTGCATAAGCTCCTCT GTGCATAAGCTCCTCT        0.41786591            0.41786591
#> GTGCATAGTTCAACCA GTGCATAGTTCAACCA        0.15743837            0.15743837
#> GTGCATATCAAGGCTT GTGCATATCAAGGCTT        0.12223569            0.12223569
#> GTGCGGTAGTAGGCCA GTGCGGTAGTAGGCCA        0.37621799            0.37621799
#> GTGCGGTCAGCGTAAG GTGCGGTCAGCGTAAG        0.29420985            0.29420985
#> GTGCTTCAGAGGGATA GTGCTTCAGAGGGATA        0.15559466            0.15559466
#> GTGCTTCTCACGACTA GTGCTTCTCACGACTA        0.52566261            0.52566261
#> GTGGGTCCATAAAGGT GTGGGTCCATAAAGGT        0.24304576            0.24304576
#> GTGGGTCCATCGATTG GTGGGTCCATCGATTG        0.58689757            0.58689757
#> GTGGGTCGTCCAGTAT GTGGGTCGTCCAGTAT        0.15867665            0.15867665
#> GTGGGTCTCCGAATGT GTGGGTCTCCGAATGT        0.16896055            0.16896055
#> GTGTGCGGTGTAAGTA GTGTGCGGTGTAAGTA        0.19234889            0.19234889
#> GTGTTAGCACCATCCT GTGTTAGCACCATCCT        0.41600124            0.41600124
#> GTGTTAGGTGAACCTT GTGTTAGGTGAACCTT        0.11552059            0.11552059
#> GTTAAGCAGGATTCGG GTTAAGCAGGATTCGG        0.20567902            0.20567902
#> GTTACAGTCCGATATG GTTACAGTCCGATATG        0.10943257            0.10943257
#> GTTCATTAGAAGATTC GTTCATTAGAAGATTC        0.38978756            0.38978756
#> GTTCATTAGGCAGTCA GTTCATTAGGCAGTCA        0.07399304            0.07399304
#> GTTCATTTCCGAACGC GTTCATTTCCGAACGC        0.27762177            0.27762177
#> GTTCTCGAGAATTCCC GTTCTCGAGAATTCCC        0.18411619            0.18411619
#> GTTCTCGCAAGAAAGG GTTCTCGCAAGAAAGG        0.20586615            0.20586615
#> GTTCTCGCAGACACTT GTTCTCGCAGACACTT        0.39055655            0.39055655
#> GTTTCTAAGAAGATTC GTTTCTAAGAAGATTC        0.30438440            0.30438440
#> GTTTCTACAGGACCCT GTTTCTACAGGACCCT        0.20863331            0.20863331
#> TAAGAGATCTGCGACG TAAGAGATCTGCGACG        0.11491443            0.11491443
#> TAAGTGCAGGGTTTCT TAAGTGCAGGGTTTCT        0.39604784            0.39604784
#> TAAGTGCTCCAAATGC TAAGTGCTCCAAATGC        0.17088580            0.17088580
#> TACACGACAACGATCT TACACGACAACGATCT        0.38608946            0.38608946
#> TACACGACATGTAGTC TACACGACATGTAGTC        0.19895942            0.19895942
#> TACACGAGTGCAACTT TACACGAGTGCAACTT        0.16755987            0.16755987
#> TACAGTGCAATCTACG TACAGTGCAATCTACG        0.39365255            0.39365255
#> TACAGTGCACAGGCCT TACAGTGCACAGGCCT        0.20307221            0.20307221
#> TACAGTGGTAAACGCG TACAGTGGTAAACGCG        0.11998791            0.11998791
#> TACCTTAAGTACATGA TACCTTAAGTACATGA        0.39521819            0.39521819
#> TACCTTATCCGATATG TACCTTATCCGATATG        0.22315621            0.22315621
#> TACGGATAGAGCTGGT TACGGATAGAGCTGGT        0.40262888            0.40262888
#> TACGGATAGGCACATG TACGGATAGGCACATG        0.39536016            0.39536016
#> TACGGATAGTTTAGGA TACGGATAGTTTAGGA        0.39133609            0.39133609
#> TACGGATGTGTTCGAT TACGGATGTGTTCGAT        0.38380219            0.38380219
#> TACGGATTCCAAGTAC TACGGATTCCAAGTAC        0.23512068            0.23512068
#> TACGGGCAGCCAACAG TACGGGCAGCCAACAG        0.20566851            0.20566851
#> TACGGGCAGCGTCAAG TACGGGCAGCGTCAAG        0.38690002            0.38690002
#> TACGGGCTCTTCGGTC TACGGGCTCTTCGGTC        0.44719865            0.44719865
#> TACGGTACATGTTCCC TACGGTACATGTTCCC        0.39406223            0.39406223
#> TACTCATAGGTGTTAA TACTCATAGGTGTTAA        0.10126233            0.10126233
#> TACTCATCAGATCCAT TACTCATCAGATCCAT        0.47840192            0.47840192
#> TACTCGCAGTAAGTAC TACTCGCAGTAAGTAC        0.17983335            0.17983335
#> TACTCGCCAGATCGGA TACTCGCCAGATCGGA        0.32420686            0.32420686
#> TACTCGCGTCGGCATC TACTCGCGTCGGCATC        0.26056453            0.26056453
#> TACTTACCAACTGCGC TACTTACCAACTGCGC        0.39232688            0.39232688
#> TACTTACCATGGGACA TACTTACCATGGGACA        0.23207934            0.23207934
#> TACTTGTAGACAAAGG TACTTGTAGACAAAGG        0.69249290            0.69249290
#> TAGACCAAGGAATCGC TAGACCAAGGAATCGC        0.39328187            0.39328187
#> TAGAGCTGTATGAAAC TAGAGCTGTATGAAAC        0.44414824            0.44414824
#> TAGAGCTGTTGTCGCG TAGAGCTGTTGTCGCG        0.44480364            0.44480364
#> TAGCCGGTCAACACAC TAGCCGGTCAACACAC        0.40240374            0.40240374
#> TAGGCATTCCAAAGTC TAGGCATTCCAAAGTC        0.44883183            0.44883183
#> TAGTGGTCAACTGCGC TAGTGGTCAACTGCGC        0.45319442            0.45319442
#> TAGTGGTCATTCACTT TAGTGGTCATTCACTT        0.10831365            0.10831365
#> TAGTGGTTCCTCAACC TAGTGGTTCCTCAACC        0.10767331            0.10767331
#> TAGTTGGCACGTAAGG TAGTTGGCACGTAAGG        0.20833557            0.20833557
#> TAGTTGGCAGTAACGG TAGTTGGCAGTAACGG        0.09990868            0.09990868
#> TAGTTGGGTACCGTAT TAGTTGGGTACCGTAT        0.39248651            0.39248651
#> TAGTTGGTCTTGCAAG TAGTTGGTCTTGCAAG        0.99963864            0.99963864
#> TATCAGGCACCAGATT TATCAGGCACCAGATT        0.39045258            0.39045258
#> TATCTCACAACTGCTA TATCTCACAACTGCTA        0.38796826            0.38796826
#> TATCTCACACGCGAAA TATCTCACACGCGAAA        0.17774621            0.17774621
#> TATCTCACATTGTGCA TATCTCACATTGTGCA        0.36735443            0.36735443
#> TATCTCAGTTCGCTAA TATCTCAGTTCGCTAA        0.45553430            0.45553430
#> TATGCCCAGAACAATC TATGCCCAGAACAATC        0.10481246            0.10481246
#> TATGCCCGTCTAGCCG TATGCCCGTCTAGCCG        0.32839621            0.32839621
#> TATGCCCTCTTACCTA TATGCCCTCTTACCTA        0.45140492            0.45140492
#> TATTACCTCGTTGCCT TATTACCTCGTTGCCT        0.19319546            0.19319546
#> TCAACGATCAACCATG TCAACGATCAACCATG        0.13088033            0.13088033
#> TCAATCTCATCCGCGA TCAATCTCATCCGCGA        0.11741761            0.11741761
#> TCAATCTCATCTCCCA TCAATCTCATCTCCCA        0.38841266            0.38841266
#> TCAATCTGTCAAAGAT TCAATCTGTCAAAGAT        0.46242134            0.46242134
#> TCAATCTTCGAACGGA TCAATCTTCGAACGGA        0.16238219            0.16238219
#> TCACAAGAGTGACATA TCACAAGAGTGACATA        0.35980599            0.35980599
#> TCACAAGCAACGATGG TCACAAGCAACGATGG        0.15567535            0.15567535
#> TCACAAGCAGACTCGC TCACAAGCAGACTCGC        0.76987840            0.76987840
#> TCACGAAGTCCGAGTC TCACGAAGTCCGAGTC        0.19295452            0.19295452
#> TCACGAATCCTATGTT TCACGAATCCTATGTT        0.39117926            0.39117926
#> TCAGATGCACTCGACG TCAGATGCACTCGACG        0.12206194            0.12206194
#> TCAGATGGTACACCGC TCAGATGGTACACCGC        0.40271551            0.40271551
#> TCAGCAACACCAGGCT TCAGCAACACCAGGCT        0.39018501            0.39018501
#> TCAGCAATCGCCGTGA TCAGCAATCGCCGTGA        0.10233049            0.10233049
#> TCAGCTCAGTCGTACT TCAGCTCAGTCGTACT        0.37716789            0.37716789
#> TCAGGTAAGGCGTACA TCAGGTAAGGCGTACA        0.17933566            0.17933566
#> TCAGGTAAGTCAATAG TCAGGTAAGTCAATAG        0.38566936            0.38566936
#> TCAGGTACAGTCACTA TCAGGTACAGTCACTA        0.18427906            0.18427906
#> TCAGGTATCTGTCCGT TCAGGTATCTGTCCGT        0.34279928            0.34279928
#> TCATTACGTCAGAAGC TCATTACGTCAGAAGC        0.38909591            0.38909591
#> TCATTTGCATTGGCGC TCATTTGCATTGGCGC        0.07811804            0.07811804
#> TCCACACCAGGCGATA TCCACACCAGGCGATA        0.20356548            0.20356548
#> TCCCGATAGTCGTTTG TCCCGATAGTCGTTTG        0.07372476            0.07372476
#> TCCCGATCAACGATGG TCCCGATCAACGATGG        0.22677029            0.22677029
#> TCGAGGCAGAACAATC TCGAGGCAGAACAATC        0.17853003            0.17853003
#> TCGAGGCAGGCATGGT TCGAGGCAGGCATGGT        0.12013944            0.12013944
#> TCGCGAGAGTACGCGA TCGCGAGAGTACGCGA        0.13285679            0.13285679
#> TCGCGAGAGTCCGTAT TCGCGAGAGTCCGTAT        0.49186056            0.49186056
#> TCGCGAGCACTTCTGC TCGCGAGCACTTCTGC        0.13989472            0.13989472
#> TCGGGACAGAGAACAG TCGGGACAGAGAACAG        0.36347096            0.36347096
#> TCGGGACGTCTGCGGT TCGGGACGTCTGCGGT        0.13165609            0.13165609
#> TCGTACCCACGCATCG TCGTACCCACGCATCG        0.43962019            0.43962019
#> TCGTACCCATATACGC TCGTACCCATATACGC        0.38303613            0.38303613
#> TCTATTGCACAGGAGT TCTATTGCACAGGAGT        0.42515801            0.42515801
#> TCTATTGTCCTATTCA TCTATTGTCCTATTCA        0.19751003            0.19751003
#> TCTATTGTCGTATCAG TCTATTGTCGTATCAG        0.13410368            0.13410368
#> TCTCATAAGAATGTGT TCTCATAAGAATGTGT        0.09941768            0.09941768
#> TCTCATACAAATTGCC TCTCATACAAATTGCC        0.24546300            0.24546300
#> TCTCATAGTGGACGAT TCTCATAGTGGACGAT        0.23924247            0.23924247
#> TCTCTAAAGCTGGAAC TCTCTAAAGCTGGAAC        0.45663833            0.45663833
#> TCTCTAAGTCGCATCG TCTCTAAGTCGCATCG        0.39104372            0.39104372
#> TCTCTAAGTGGCAAAC TCTCTAAGTGGCAAAC        0.25123793            0.25123793
#> TCTCTAAGTTACCGAT TCTCTAAGTTACCGAT        0.40172929            0.40172929
#> TCTCTAATCCCAACGG TCTCTAATCCCAACGG        0.11967899            0.11967899
#> TCTGAGAGTCAAACTC TCTGAGAGTCAAACTC        0.43806591            0.43806591
#> TCTGGAACACACGCTG TCTGGAACACACGCTG        0.39840128            0.39840128
#> TCTGGAACACCTATCC TCTGGAACACCTATCC        0.65426080            0.65426080
#> TCTGGAAGTCATCGGC TCTGGAAGTCATCGGC        0.18126525            0.18126525
#> TCTGGAATCGGTGTCG TCTGGAATCGGTGTCG        0.38255141            0.38255141
#> TCTTTCCGTGCAACTT TCTTTCCGTGCAACTT        0.28224240            0.28224240
#> TCTTTCCTCGTCTGCT TCTTTCCTCGTCTGCT        0.38511321            0.38511321
#> TGAAAGAGTTATCACG TGAAAGAGTTATCACG        0.39124616            0.39124616
#> TGACGGCAGGAATGGA TGACGGCAGGAATGGA        0.40226399            0.40226399
#> TGACGGCCATATGAGA TGACGGCCATATGAGA        0.16450645            0.16450645
#> TGACTAGGTGTATGGG TGACTAGGTGTATGGG        0.28153908            0.28153908
#> TGACTAGTCCACGTGG TGACTAGTCCACGTGG        0.44370652            0.44370652
#> TGACTAGTCCTATTCA TGACTAGTCCTATTCA        0.30612641            0.30612641
#> TGACTTTAGCCAGAAC TGACTTTAGCCAGAAC        0.75592340            0.75592340
#> TGACTTTAGCGTTGCC TGACTTTAGCGTTGCC        0.44691041            0.44691041
#> TGACTTTAGCGTTTAC TGACTTTAGCGTTTAC        0.39061794            0.39061794
#> TGAGAGGAGATCCCAT TGAGAGGAGATCCCAT        0.45736001            0.45736001
#> TGAGAGGAGCCAGGAT TGAGAGGAGCCAGGAT        0.18533543            0.18533543
#> TGAGAGGAGCCGTCGT TGAGAGGAGCCGTCGT        0.23480848            0.23480848
#> TGAGAGGAGGCGATAC TGAGAGGAGGCGATAC        0.22496278            0.22496278
#> TGAGAGGAGGTAAACT TGAGAGGAGGTAAACT        0.17559506            0.17559506
#> TGAGAGGCAGGAACGT TGAGAGGCAGGAACGT        0.18571191            0.18571191
#> TGAGAGGGTAAGGGCT TGAGAGGGTAAGGGCT        0.44282872            0.44282872
#> TGAGCATGTTAAAGAC TGAGCATGTTAAAGAC        0.30259819            0.30259819
#> TGAGCATTCACCCGAG TGAGCATTCACCCGAG        0.39338014            0.39338014
#> TGAGCATTCGGGAGTA TGAGCATTCGGGAGTA        0.45364723            0.45364723
#> TGAGCCGTCATGTCCC TGAGCCGTCATGTCCC        0.38940412            0.38940412
#> TGAGCCGTCCGAACGC TGAGCCGTCCGAACGC        0.38391924            0.38391924
#> TGAGCCGTCCGCTGTT TGAGCCGTCCGCTGTT        0.13179411            0.13179411
#> TGAGGGAAGAGGTACC TGAGGGAAGAGGTACC        0.22991504            0.22991504
#> TGAGGGACAGACTCGC TGAGGGACAGACTCGC        0.39136648            0.39136648
#> TGATTTCAGTACGCCC TGATTTCAGTACGCCC        0.37493239            0.37493239
#> TGATTTCAGTGACTCT TGATTTCAGTGACTCT        0.39002891            0.39002891
#> TGATTTCGTGAGGGTT TGATTTCGTGAGGGTT        0.28885414            0.28885414
#> TGCCAAACATGACATC TGCCAAACATGACATC        0.29981422            0.29981422
#> TGCCCATCACTGTGTA TGCCCATCACTGTGTA        0.37853276            0.37853276
#> TGCCCTAAGTCCGGTC TGCCCTAAGTCCGGTC        0.11045788            0.11045788
#> TGCGCAGGTACGCACC TGCGCAGGTACGCACC        0.44755194            0.44755194
#> TGCGCAGTCTGCCAGG TGCGCAGTCTGCCAGG        0.44342473            0.44342473
#> TGCGGGTAGACGCACA TGCGGGTAGACGCACA        0.14129548            0.14129548
#> TGCGGGTAGGGTTCCC TGCGGGTAGGGTTCCC        0.17596371            0.17596371
#> TGCGGGTTCCACTCCA TGCGGGTTCCACTCCA        0.13855331            0.13855331
#> TGCGGGTTCTCAAACG TGCGGGTTCTCAAACG        0.38496871            0.38496871
#> TGCGTGGAGGTAAACT TGCGTGGAGGTAAACT        0.39001354            0.39001354
#> TGCGTGGAGGTTCCTA TGCGTGGAGGTTCCTA        0.13706876            0.13706876
#> TGCGTGGAGTACACCT TGCGTGGAGTACACCT        0.18565147            0.18565147
#> TGCGTGGGTCCCTTGT TGCGTGGGTCCCTTGT        0.12364395            0.12364395
#> TGCGTGGGTTACTGAC TGCGTGGGTTACTGAC        0.14988764            0.14988764
#> TGCTACCCAGGGCATA TGCTACCCAGGGCATA        0.39601469            0.39601469
#> TGCTACCCATAAGACA TGCTACCCATAAGACA        0.15092498            0.15092498
#> TGCTGCTAGGATGCGT TGCTGCTAGGATGCGT        0.31576947            0.31576947
#> TGCTGCTGTCCGTCAG TGCTGCTGTCCGTCAG        0.43624129            0.43624129
#> TGGACGCAGTGTTGAA TGGACGCAGTGTTGAA        0.38269902            0.38269902
#> TGGACGCTCAGGTAAA TGGACGCTCAGGTAAA        0.45637658            0.45637658
#> TGGACGCTCTTGTTTG TGGACGCTCTTGTTTG        0.11648748            0.11648748
#> TGGCCAGAGCTACCTA TGGCCAGAGCTACCTA        0.45241195            0.45241195
#> TGGCCAGCAAGGACTG TGGCCAGCAAGGACTG        0.99879803            0.99879803
#> TGGCCAGCAATGAAAC TGGCCAGCAATGAAAC        0.10669177            0.10669177
#> TGGCCAGCACGAAAGC TGGCCAGCACGAAAGC        0.38553754            0.38553754
#> TGGCCAGGTCGGCATC TGGCCAGGTCGGCATC        0.40094984            0.40094984
#> TGGCGCATCACAACGT TGGCGCATCACAACGT        0.43493934            0.43493934
#> TGGCTGGGTACCCAAT TGGCTGGGTACCCAAT        0.30145002            0.30145002
#> TGGGAAGCACGTAAGG TGGGAAGCACGTAAGG        0.10787200            0.10787200
#> TGGGAAGCATGACATC TGGGAAGCATGACATC        0.40274824            0.40274824
#> TGGGCGTGTGTGCGTC TGGGCGTGTGTGCGTC        0.18754150            0.18754150
#> TGGTTAGCAGTCGTGC TGGTTAGCAGTCGTGC        0.12763465            0.12763465
#> TGGTTCCCACCGCTAG TGGTTCCCACCGCTAG        0.39105268            0.39105268
#> TGGTTCCGTCACCTAA TGGTTCCGTCACCTAA        0.40571521            0.40571521
#> TGTATTCCAATGACCT TGTATTCCAATGACCT        0.43693192            0.43693192
#> TGTCCCAAGTTTAGGA TGTCCCAAGTTTAGGA        0.38748223            0.38748223
#> TGTCCCACATTATCTC TGTCCCACATTATCTC        0.13885077            0.13885077
#> TGTCCCAGTCCGTGAC TGTCCCAGTCCGTGAC        0.63062828            0.63062828
#> TGTCCCATCTCTGCTG TGTCCCATCTCTGCTG        0.15073939            0.15073939
#> TGTGGTACAAGACGTG TGTGGTACAAGACGTG        0.43946200            0.43946200
#> TGTGGTACACATCCGG TGTGGTACACATCCGG        1.00000000            1.00000000
#> TGTGGTAGTGCCTGCA TGTGGTAGTGCCTGCA        0.19709420            0.19709420
#> TGTGGTATCCAGTAGT TGTGGTATCCAGTAGT        0.34070568            0.34070568
#> TGTGTTTAGGAATGGA TGTGTTTAGGAATGGA        0.17615413            0.17615413
#> TGTGTTTGTAGATTAG TGTGTTTGTAGATTAG        0.18724131            0.18724131
#> TGTGTTTGTATGAATG TGTGTTTGTATGAATG        0.12409521            0.12409521
#> TTAACTCCAAAGAATC TTAACTCCAAAGAATC        0.45146685            0.45146685
#> TTAACTCCAGGATTGG TTAACTCCAGGATTGG        0.16048138            0.16048138
#> TTAACTCGTAATCGTC TTAACTCGTAATCGTC        0.10007480            0.10007480
#> TTAGGCAGTTTGTTGG TTAGGCAGTTTGTTGG        0.10183060            0.10183060
#> TTATGCTCACCCAGTG TTATGCTCACCCAGTG        0.35144609            0.35144609
#> TTCCCAGAGTTACCCA TTCCCAGAGTTACCCA        0.36305858            0.36305858
#> TTCCCAGGTATTCTCT TTCCCAGGTATTCTCT        0.40031313            0.40031313
#> TTCCCAGTCATAGCAC TTCCCAGTCATAGCAC        0.10666207            0.10666207
#> TTCGAAGAGTACGCCC TTCGAAGAGTACGCCC        0.35340972            0.35340972
#> TTCGGTCAGCTAGCCC TTCGGTCAGCTAGCCC        0.20626969            0.20626969
#> TTCGGTCAGGGCATGT TTCGGTCAGGGCATGT        0.30574064            0.30574064
#> TTCGGTCGTATAGTAG TTCGGTCGTATAGTAG        0.15668028            0.15668028
#> TTCGGTCTCGACCAGC TTCGGTCTCGACCAGC        0.11072727            0.11072727
#> TTCTACAAGATCGATA TTCTACAAGATCGATA        0.38588853            0.38588853
#> TTCTACACATAAAGGT TTCTACACATAAAGGT        0.11427302            0.11427302
#> TTCTACAGTACTCTCC TTCTACAGTACTCTCC        0.38525417            0.38525417
#> TTCTACATCGCTTAGA TTCTACATCGCTTAGA        0.45274329            0.45274329
#> TTCTACATCTGCAGTA TTCTACATCTGCAGTA        0.45597912            0.45597912
#> TTCTCAAAGCCGCCTA TTCTCAAAGCCGCCTA        0.40640175            0.40640175
#> TTCTCAACAATAGAGT TTCTCAACAATAGAGT        0.20215902            0.20215902
#> TTCTCCTAGTCCATAC TTCTCCTAGTCCATAC        0.17820049            0.17820049
#> TTCTCCTCAATGGTCT TTCTCCTCAATGGTCT        0.39778615            0.39778615
#> TTCTTAGTCAGTGTTG TTCTTAGTCAGTGTTG        0.16999740            0.16999740
#> TTGAACGAGAGCAATT TTGAACGAGAGCAATT        0.23279176            0.23279176
#> TTGAACGAGGAGTTTA TTGAACGAGGAGTTTA        0.27643654            0.27643654
#> TTGAACGGTGAGGGAG TTGAACGGTGAGGGAG        0.15550984            0.15550984
#> TTGAACGGTGATGCCC TTGAACGGTGATGCCC        0.38617273            0.38617273
#> TTGACTTAGACAAGCC TTGACTTAGACAAGCC        0.43290936            0.43290936
#> TTGACTTGTTCCATGA TTGACTTGTTCCATGA        0.18719247            0.18719247
#> TTGACTTTCACCGGGT TTGACTTTCACCGGGT        0.10610177            0.10610177
#> TTGACTTTCAGCAACT TTGACTTTCAGCAACT        0.39749700            0.39749700
#> TTGACTTTCCAGTAGT TTGACTTTCCAGTAGT        0.39807188            0.39807188
#> TTGCCGTCACCGAATT TTGCCGTCACCGAATT        0.14363331            0.14363331
#> TTGCCGTTCTCTGTCG TTGCCGTTCTCTGTCG        0.38855146            0.38855146
#> TTGCGTCCACCGATAT TTGCGTCCACCGATAT        0.21063840            0.21063840
#> TTGGAACAGAAGATTC TTGGAACAGAAGATTC        0.14015561            0.14015561
#> TTGGAACGTATATGAG TTGGAACGTATATGAG        0.38815256            0.38815256
#> TTGGCAAAGAGTAAGG TTGGCAAAGAGTAAGG        0.22727657            0.22727657
#> TTGTAGGTCAAAGTAG TTGTAGGTCAAAGTAG        0.14531083            0.14531083
#> TTTACTGCAGCATACT TTTACTGCAGCATACT        0.33635765            0.33635765
#> TTTACTGGTGCTAGCC TTTACTGGTGCTAGCC        0.38820663            0.38820663
#> TTTATGCAGGTTACCT TTTATGCAGGTTACCT        0.19591461            0.19591461
#> TTTATGCGTCGAACAG TTTATGCGTCGAACAG        0.44545815            0.44545815
#> TTTATGCTCCATTCTA TTTATGCTCCATTCTA        0.44499739            0.44499739
#> TTTCCTCCAGGATCGA TTTCCTCCAGGATCGA        0.28220979            0.28220979
#> TTTCCTCGTGGCGAAT TTTCCTCGTGGCGAAT        0.10858573            0.10858573
#> TTTGCGCAGACAGACC TTTGCGCAGACAGACC        0.39626086            0.39626086
#> TTTGGTTAGACAGGCT TTTGGTTAGACAGGCT        0.19051740            0.19051740
#> TTTGGTTCACCAGATT TTTGGTTCACCAGATT        0.39635978            0.39635978
#> TTTGGTTCACGAAGCA TTTGGTTCACGAAGCA        0.22168078            0.22168078
#> TTTGGTTTCCTTTCGG TTTGGTTTCCTTTCGG        0.19771378            0.19771378
#> TTTGTCAAGAATGTGT TTTGTCAAGAATGTGT        0.15706795            0.15706795
#> TTTGTCATCGAATGCT TTTGTCATCGAATGCT        0.20337776            0.20337776
#> TTTGTCATCTGTTTGT TTTGTCATCTGTTTGT        0.15401696            0.15401696
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
#> Pyy       Pyy                     0.1700076                             665
#> Sst       Sst                     0.2816426                              64
#> Rbp4     Rbp4                     0.6920430                             647
#>      stochastic_pseudotimer.sq stochastic_pseudotimedev.expl
#> Pyy                 0.13243812                     0.1458385
#> Sst                 0.07437209                     0.4873345
#> Rbp4                0.16038102                     0.1794935
#>      stochastic_pseudotimepadjust  peaktime feature_split index
#> Pyy                             0 0.1700076            NA     1
#> Sst                             0 0.2816426            NA     2
#> Rbp4                            0 0.6920430            NA     3
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
