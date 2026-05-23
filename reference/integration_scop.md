# The integration workflow

Integrate single-cell RNA-seq data using various integration methods.
For `ChromatinAssay`, the current workflow uses `TFIDF + SVD/LSI`
preprocessing. In this setting, `Uncorrected` is supported directly and
`Harmony5` will be automatically redirected to the legacy `Harmony`
workflow. `Seurat` and `RPCA` are currently not supported for
`ChromatinAssay`.

## Usage

``` r
integration_scop(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_method = c("Uncorrected", "Seurat", "CCA", "RPCA", "scVI", "PeakVI",
    "PoissonVI", "WNN", "MultiMAP", "GLUE", "scVI5", "MNN", "fastMNN", "fastMNN5",
    "Harmony", "Harmony5", "Scanorama", "BBKNN", "CSS", "Coralysis", "LIGER", "Conos",
    "ComBat"),
  compute_lisi = FALSE,
  lisi_label_colnames = NULL,
  lisi_reduction = NULL,
  lisi_dims = NULL,
  lisi_prefix = NULL,
  lisi_tool_name = NULL,
  lisi_perplexity = 30,
  lisi_nn_method = c("auto", "exact"),
  lisi_tol = 1e-05,
  lisi_max_iter = 50,
  compute_metrics = FALSE,
  metrics_batch_col = NULL,
  metrics_celltype_col = NULL,
  metrics_reduction = NULL,
  metrics_cluster_col = NULL,
  metrics_tool_name = NULL,
  metrics_k_graph = 15,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt_merge:

  A merged \`Seurat\` object that includes the batch information.

- batch:

  A character string specifying the batch variable name.

- append:

  Whether the integrated data will be appended to the original Seurat
  object (`srt_merge`). Default is `TRUE`.

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- integration_method:

  A character vector specifying the integration method to use. Supported
  methods are: `"Uncorrected"`, `"Seurat"`, `"CCA"`, `"RPCA"`, `"scVI"`,
  `"PeakVI"`, `"PoissonVI"`, `"WNN"`, `"MultiMAP"`, `"GLUE"`, `"scVI5"`,
  `"MNN"`, `"fastMNN"`, `"fastMNN5"`, `"Harmony"`, `"Harmony5"`,
  `"Scanorama"`, `"BBKNN"`, `"CSS"`, `"Coralysis"`, `"LIGER"`,
  `"Conos"`, `"ComBat"`. Default is `"Uncorrected"`. For
  `ChromatinAssay`, prefer `"Uncorrected"` or `"Harmony5"`; the latter
  is automatically switched to `"Harmony"`.

- compute_lisi:

  Whether to compute LISI scores on the integrated result. Default is
  `FALSE`.

- lisi_label_colnames:

  Character vector of metadata columns used to compute LISI. If `NULL`
  and `compute_lisi = TRUE`, `batch` will be used when it is a single
  metadata column name.

- lisi_reduction:

  Dimensional reduction used for LISI computation. Default is `NULL`,
  which uses
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  from the integrated object.

- lisi_dims:

  Dimensions used from `lisi_reduction`. Default is `NULL`, which uses
  all available dimensions.

- lisi_prefix:

  Prefix used when storing LISI metadata columns. Default is `NULL`,
  which uses `lisi_reduction`.

- lisi_tool_name:

  Name of the tool entry used to store LISI results. Default is `NULL`,
  which uses `paste0(lisi_prefix, "_LISI")`.

- lisi_perplexity:

  Effective neighborhood size used by LISI. Default is `30`.

- lisi_nn_method:

  Nearest-neighbor backend passed to
  [`RunLISI()`](https://mengxu98.github.io/scop/reference/RunLISI.md).
  One of `"auto"` or `"exact"`. Default is `"auto"`.

- lisi_tol:

  Tolerance used in the LISI binary search. Default is `1e-5`.

- lisi_max_iter:

  Maximum iterations used in the LISI binary search. Default is `50`.

- compute_metrics:

  Whether to compute integration summary metrics on the selected
  reduction. Default is `FALSE`.

- metrics_batch_col:

  Metadata column used for batch-mixing metrics. Default is `NULL`,
  which uses `batch` when it is a single metadata column name.

- metrics_celltype_col:

  Metadata column used for biological conservation metrics. Default is
  `NULL`.

- metrics_reduction:

  Reduction used for integration metric computation. Default is `NULL`,
  which uses
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  from the integrated object.

- metrics_cluster_col:

  Metadata column used as cluster labels for `celltype_NMI`,
  `celltype_ARI`, and `celltype_purity`. Default is `NULL`, which
  resolves the integrated cluster column automatically when possible.

- metrics_tool_name:

  Name of the tool entry used to store integration metrics. Default is
  `NULL`, which uses `paste0(integration_method, "_metrics")`.

- metrics_k_graph:

  Number of neighbors used for graph-connectivity computation. Default
  is `15`.

- do_normalization:

  Whether data normalization should be performed. Default is `TRUE`.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, and `"scran"`. Default is
  `"LogNormalize"`.

- do_HVF_finding:

  Whether to perform high variable feature finding. If `TRUE`, the
  function will force to find the highly variable features (HVF) using
  the specified HVF method.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`. Default is `"separate"`.

- HVF_method:

  The method to use for finding highly variable features. Options are
  `"vst"`, `"mvp"`, `"disp"`, or `"scran"`. Default is `"vst"`.

- nHVF:

  The number of highly variable features to select. If NULL, all highly
  variable features will be used. Default is `2000`.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable. Default is `1`.

- HVF:

  A vector of feature names to use as highly variable features. If NULL,
  the function will use the highly variable features identified by the
  HVF method.

- do_scaling:

  Whether to perform scaling. If `TRUE`, the function will force to
  scale the data using the
  [Seurat::ScaleData](https://satijalab.org/seurat/reference/ScaleData.html)
  function.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables. Default is `NULL`.

- regression_model:

  The regression model to use for scaling. Options are `"linear"`,
  `"poisson"`, or `"negativebinomial"`. Default is `"linear"`.

- scale_within_batch:

  Whether to scale data within each batch. Only valid when the
  `integration_method` is one of `"Uncorrected"`, `"Seurat"`, `"MNN"`,
  `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.

- linear_reduction:

  The linear dimensionality reduction method to use. Options are
  `"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, or `"glmpca"`. Default is
  `"pca"`.

- linear_reduction_dims:

  The number of dimensions to keep after linear dimensionality
  reduction. Default is `50`.

- linear_reduction_dims_use:

  The dimensions to use for downstream analysis. If `NULL`, estimated
  dimensions stored in the linear reduction will be used when available;
  otherwise, the first up to `50` dimensions will be used as a fallback.

- linear_reduction_params:

  A list of parameters to pass to the linear dimensionality reduction
  method.

- force_linear_reduction:

  Whether to force linear dimensionality reduction even if the specified
  reduction is already present in the Seurat object.

- nonlinear_reduction:

  The nonlinear dimensionality reduction method to use. Options are
  `"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`,
  `"trimap"`, `"largevis"`, or `"fr"`. Default is `"umap"`.

- nonlinear_reduction_dims:

  The number of dimensions to keep after nonlinear dimensionality
  reduction. If a vector is provided, different numbers of dimensions
  can be specified for each method. Default is `c(2, 3)`.

- nonlinear_reduction_params:

  A list of parameters to pass to the nonlinear dimensionality reduction
  method.

- force_nonlinear_reduction:

  Whether to force nonlinear dimensionality reduction even if the
  specified reduction is already present in the Seurat object. Default
  is `TRUE`.

- neighbor_metric:

  The distance metric to use for finding neighbors. Options are
  `"euclidean"`, `"cosine"`, `"manhattan"`, or `"hamming"`. Default is
  `"euclidean"`.

- neighbor_k:

  The number of nearest neighbors to use for finding neighbors. Default
  is `20`.

- cluster_algorithm:

  The clustering algorithm to use. Options are `"louvain"`, `"slm"`, or
  `"leiden"`. Default is `"louvain"`.

- cluster_resolution:

  The resolution parameter to use for clustering. Larger values result
  in fewer clusters. Default is `0.6`.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to the integration method functions.

## Value

A `Seurat` object.

For `ChromatinAssay`, integrated outputs are additionally normalized to
the ATAC naming convention used in `scop`, including `*lsi`, `*UMAP2D`,
cluster aliases, and `ATAC_default_*` metadata when available.

## See also

[Seurat_integrate](https://mengxu98.github.io/scop/reference/Seurat_integrate.md),
[scVI_integrate](https://mengxu98.github.io/scop/reference/scVI_integrate.md),
[MultiMAP_integrate](https://mengxu98.github.io/scop/reference/MultiMAP_integrate.md),
[GLUE_integrate](https://mengxu98.github.io/scop/reference/GLUE_integrate.md),
[MNN_integrate](https://mengxu98.github.io/scop/reference/MNN_integrate.md),
[fastMNN_integrate](https://mengxu98.github.io/scop/reference/fastMNN_integrate.md),
[Harmony_integrate](https://mengxu98.github.io/scop/reference/Harmony_integrate.md),
[Scanorama_integrate](https://mengxu98.github.io/scop/reference/Scanorama_integrate.md),
[BBKNN_integrate](https://mengxu98.github.io/scop/reference/BBKNN_integrate.md),
[CSS_integrate](https://mengxu98.github.io/scop/reference/CSS_integrate.md),
[Coralysis_integrate](https://mengxu98.github.io/scop/reference/Coralysis_integrate.md),
[LIGER_integrate](https://mengxu98.github.io/scop/reference/LIGER_integrate.md),
[Conos_integrate](https://mengxu98.github.io/scop/reference/Conos_integrate.md),
[ComBat_integrate](https://mengxu98.github.io/scop/reference/ComBat_integrate.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected",
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10,
  nonlinear_reduction_dims = 2,
  compute_lisi = TRUE,
  lisi_label_colnames = "tech",
  lisi_perplexity = 10
)
#> ◌ [2026-05-23 09:23:52] Run integration workflow...
#> ℹ [2026-05-23 09:23:53] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-23 09:23:53] Checking a list of <Seurat>...
#> ! [2026-05-23 09:23:53] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:23:53] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-05-23 09:23:55] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-05-23 09:23:56] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:23:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-05-23 09:23:58] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-05-23 09:23:58] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:23:58] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:00] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-05-23 09:24:00] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:24:00] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:02] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-05-23 09:24:03] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:24:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:05] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:05] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:24:05] Number of available HVF: 500
#> ℹ [2026-05-23 09:24:05] Finished check
#> ℹ [2026-05-23 09:24:07] Perform Uncorrected integration
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-23 09:24:07] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-23 09:24:08] Perform "pca" linear dimension reduction
#> ℹ [2026-05-23 09:24:08] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-23 09:24:09] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 09:24:09] Reorder clusters...
#> ℹ [2026-05-23 09:24:09] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:24:09] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:10)
#> ℹ [2026-05-23 09:24:16] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:10)
#> ℹ [2026-05-23 09:24:23] Compute LISI scores from reduction "UncorrectedpcaUMAP2D"
#> ℹ [2026-05-23 09:24:23] Compute LISI scores from reduction "UncorrectedUMAP2D"
#> ✔ [2026-05-23 09:24:23] Stored LISI scores in metadata: "UncorrectedpcaUMAP2D_tech_LISI" and "UncorrectedUMAP2D_tech_LISI"
#> ✔ [2026-05-23 09:24:24] Uncorrected integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "UncorrectedUMAP2D"
)

LISIPlot(
  panc8_sub,
  features = c("UncorrectedpcaUMAP2D_tech_LISI", "UncorrectedUMAP2D_tech_LISI")
)


panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "LIGER"
)
#> ◌ [2026-05-23 09:24:25] Run integration workflow...
#> ℹ [2026-05-23 09:24:25] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-23 09:24:26] Checking a list of <Seurat>...
#> ℹ [2026-05-23 09:24:27] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:24:27] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:27] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:24:27] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:28] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:24:28] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:28] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:24:28] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:29] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:24:29] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:24:29] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:24:30] Number of available HVF: 2000
#> ℹ [2026-05-23 09:24:30] Finished check
#> Warning: Layer ‘ligerScaleData’ is empty
#> ℹ [2026-05-23 09:24:36] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-05-23 09:24:36] Perform LIGER integration
#> ℹ [2026-05-23 09:24:42] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-23 09:24:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 09:24:43] Reorder clusters...
#> ℹ [2026-05-23 09:24:44] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:24:44] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-05-23 09:24:51] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-05-23 09:24:57] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:7)
#> ✔ [2026-05-23 09:25:06] LIGER integration completed
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Harmony",
  compute_lisi = TRUE,
  lisi_label_colnames = "tech"
)
#> ◌ [2026-05-23 09:25:06] Run integration workflow...
#> ℹ [2026-05-23 09:25:06] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-23 09:25:08] Checking a list of <Seurat>...
#> ℹ [2026-05-23 09:25:08] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:25:08] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-05-23 09:25:08] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:25:08] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-05-23 09:25:09] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:25:09] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-05-23 09:25:10] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:25:10] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-05-23 09:25:10] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:25:10] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:25:11] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:25:11] Number of available HVF: 2000
#> ℹ [2026-05-23 09:25:11] Finished check
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-23 09:25:36] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-23 09:25:37] Perform linear dimension reduction("pca")
#> ℹ [2026-05-23 09:25:37] Perform Harmony integration
#> ℹ [2026-05-23 09:25:37] Using "Harmonypca" (1:20) as input
#> ℹ [2026-05-23 09:25:40] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-23 09:25:41] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 09:25:41] Reorder clusters...
#> ℹ [2026-05-23 09:25:42] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:25:42] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-05-23 09:25:48] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-05-23 09:25:55] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:7)
#> ℹ [2026-05-23 09:26:01] Compute LISI scores from reduction "UncorrectedpcaUMAP2D"
#> ℹ [2026-05-23 09:26:01] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-05-23 09:26:02] Stored LISI scores in metadata: "UncorrectedpcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-05-23 09:26:08] Harmony integration completed
LISIPlot(
  panc8_sub,
  features = c("HarmonypcaUMAP2D_tech_LISI", "HarmonyUMAP2D_tech_LISI")
)
#> Error in benchmark_feature_plot(srt = srt, features = features, tool_name = tool_name,     reduction = reduction, plot_type = plot_type, plot_boxplot = plot_boxplot,     boxplot_jitter = boxplot_jitter, combine = combine, nrow = nrow,     ncol = ncol, byrow = byrow, pt.size = pt.size, pt.alpha = pt.alpha,     palette = palette, palcolor = palcolor, theme_use = theme_use,     theme_args = theme_args, verbose = verbose, ...): The following benchmark columns are missing:
#> "HarmonypcaUMAP2D_tech_LISI"

if (requireNamespace("Signac", quietly = TRUE)) {
  data("pbmcmultiome_sub", package = "scop")
  pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
  pbmcmultiome_sub <- integration_scop(
    pbmcmultiome_sub,
    batch = "batch",
    assay = "peaks",
    integration_method = "Harmony5",
    normalization_method = "TFIDF"
  )

integration_methods <- c(
  "Uncorrected", "Seurat", "CCA", "RPCA", "scVI", "scVI5",
  "MNN", "fastMNN", "fastMNN5", "Harmony", "Harmony5",
  "Scanorama", "BBKNN", "CSS", "Coralysis", "LIGER", "Conos", "ComBat"
)
p_list <- list()
for (method in integration_methods) {
  panc8_sub <- integration_scop(
    panc8_sub,
    batch = "tech",
    integration_method = method,
    linear_reduction_dims_use = 1:50,
    nonlinear_reduction = "umap"
  )
  p_list[[method]] <- CellDimPlot(
    panc8_sub,
    group.by = c("tech", "celltype"),
    reduction = paste0(method, "UMAP2D"),
    xlab = "", ylab = "",
    title = method,
    legend.position = "none",
    theme_use = "theme_blank"
  )
}

nonlinear_reductions <- c(
  "umap", "tsne", "dm", "phate",
  "pacmap", "trimap", "largevis", "fr"
)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat",
  linear_reduction_dims_use = 1:50,
  nonlinear_reduction = nonlinear_reductions
)
for (nr in nonlinear_reductions) {
  print(
    CellDimPlot(
      panc8_sub,
      group.by = c("tech", "celltype"),
      reduction = paste0("Seurat", nr, "2D"),
      xlab = "", ylab = "", title = nr,
      legend.position = "none", theme_use = "theme_blank"
    )
  )
}
}
#> ◌ [2026-05-23 09:26:09] Run integration workflow...
#> ! [2026-05-23 09:26:09] `integration_method = 'Harmony5'` is not compatible with <ChromatinAssay> in current Seurat v5 workflow. Automatically switch to "Harmony"
#> ℹ [2026-05-23 09:26:09] Split `srt_merge` into `srt_list` by "batch"
#> ℹ [2026-05-23 09:26:09] Checking a list of <Seurat>...
#> ! [2026-05-23 09:26:09] Data 1/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-05-23 09:26:09] Perform `RunTFIDF()` on 1/2 of `srt_list`...
#> ℹ [2026-05-23 09:26:09] Perform `FindTopFeatures()` on 1/2 of `srt_list`...
#> ! [2026-05-23 09:26:10] Data 2/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-05-23 09:26:10] Perform `RunTFIDF()` on 2/2 of `srt_list`...
#> ℹ [2026-05-23 09:26:10] Perform `FindTopFeatures()` on 2/2 of `srt_list`...
#> ℹ [2026-05-23 09:26:10] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:26:10] Number of available HVF: 2000
#> ℹ [2026-05-23 09:26:10] Finished check
#> ℹ [2026-05-23 09:26:10] `normalization_method` is "TFIDF". Use lsi workflow...
#> ℹ [2026-05-23 09:26:10] Perform linear dimension reduction("svd")
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-05-23 09:26:11] Perform Harmony integration
#> ℹ [2026-05-23 09:26:11] Using "Harmonysvd" (2:10) as input
#> ℹ [2026-05-23 09:26:11] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 09:26:11] Reorder clusters...
#> ℹ [2026-05-23 09:26:11] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:26:11] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ℹ [2026-05-23 09:26:17] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ✔ [2026-05-23 09:26:23] Harmony integration completed
#> ◌ [2026-05-23 09:26:23] Run integration workflow...
#> ℹ [2026-05-23 09:26:23] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-23 09:26:24] Checking a list of <Seurat>...
#> ℹ [2026-05-23 09:26:25] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:26:25] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-05-23 09:26:25] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:26:25] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-05-23 09:26:26] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:26:26] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-05-23 09:26:27] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:26:27] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-05-23 09:26:27] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:26:27] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:26:28] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:26:28] Number of available HVF: 2000
#> ℹ [2026-05-23 09:26:29] Finished check
#> ℹ [2026-05-23 09:28:19] Perform Uncorrected integration
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-23 09:28:20] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-23 09:28:20] Perform "pca" linear dimension reduction
#> Warning: Number of dimensions changing from 20 to 50
#> ℹ [2026-05-23 09:28:21] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-23 09:28:21] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-23 09:28:21] Reorder clusters...
#> ℹ [2026-05-23 09:28:22] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:28:22] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-05-23 09:28:29] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-05-23 09:28:35] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> Warning: Number of dimensions changing from 20 to 50
#> ✔ [2026-05-23 09:29:08] Uncorrected integration completed
#> ◌ [2026-05-23 09:29:08] Run integration workflow...
#> ℹ [2026-05-23 09:29:09] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-23 09:29:13] Checking a list of <Seurat>...
#> ℹ [2026-05-23 09:29:13] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:29:13] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-05-23 09:29:14] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:29:14] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-05-23 09:29:15] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:29:15] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-05-23 09:29:15] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:29:15] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-05-23 09:29:16] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-05-23 09:29:16] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-23 09:29:17] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:29:17] Number of available HVF: 2000
#> ℹ [2026-05-23 09:29:18] Finished check
#> ℹ [2026-05-23 09:39:01] Perform FindIntegrationAnchors
#> Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 2 globals exported for future expression (‘FUN()’) is 2.04 GiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are two globals: ‘FUN’ (2.04 GiB of class ‘function’) and ‘anchor.features’ (25.71 KiB of class ‘character’)
```
