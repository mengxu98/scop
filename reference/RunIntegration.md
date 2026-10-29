# Integration workflow

Integrate single-cell data with one or more methods. For
`ChromatinAssay`, the workflow uses TFIDF + SVD/LSI; `Uncorrected` is
supported directly, `Harmony5` is redirected to `Harmony`, and
`Seurat`/`RPCA` are not supported.

## Usage

``` r
RunIntegration(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_methods = c("Uncorrected", "Seurat", "CCA", "RPCA", "scVI", "PeakVI",
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
  lisi_tol = 1e-05,
  lisi_max_iter = 50,
  lisi_knn_algorithm = c("auto", "brute_force", "clustered"),
  lisi_cores = NULL,
  lisi_max_dense_bytes = Inf,
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
  integration_method = NULL,
  ...
)
```

## Arguments

- srt_merge:

  A merged \`Seurat\` object that includes the batch information.

- batch:

  Batch variable name.

- append:

  Append integrated results to `srt_merge`.

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- assay:

  Assay to use. `NULL` uses the default assay.

- integration_methods:

  Method(s) to run. Multiple methods require `append = TRUE` and are
  applied sequentially. For `ChromatinAssay`, prefer `"Uncorrected"` or
  `"Harmony5"`.

- compute_lisi, lisi_label_colnames, lisi_reduction, lisi_dims,
  lisi_prefix, lisi_tool_name, lisi_perplexity, lisi_tol, lisi_max_iter,
  lisi_knn_algorithm, lisi_cores, lisi_max_dense_bytes:

  LISI scores. `lisi_label_colnames = NULL` uses `batch` when it is a
  single metadata column; `lisi_reduction = NULL` uses
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- compute_metrics, metrics_batch_col, metrics_celltype_col,
  metrics_reduction, metrics_cluster_col, metrics_tool_name,
  metrics_k_graph:

  Integration summary metrics on the selected reduction.

- do_normalization:

  Whether data normalization should be performed.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, and `"scran"`.

- do_HVF_finding, HVF_method, nHVF, HVF:

  Highly variable features. `HVF_method` is `"vst"`, `"mvp"`, `"disp"`,
  or `"scran"`.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable.

- do_scaling:

  Force scaling via ScaleData.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables.

- regression_model:

  `"linear"`, `"poisson"`, or `"negativebinomial"`.

- scale_within_batch:

  Scale within each batch. Only used by `"Uncorrected"`, `"Seurat"`,
  `"MNN"`, `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.

- linear_reduction, linear_reduction_dims, linear_reduction_dims_use,
  linear_reduction_params, force_linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

- nonlinear_reduction_dims:

  Number(s) of nonlinear dimensions to compute.

- nonlinear_reduction_params:

  Shared parameters, optionally overridden by lists named after
  individual nonlinear reduction methods.

- force_nonlinear_reduction:

  Recompute nonlinear reductions even when stored results exist.

- neighbor_metric, neighbor_k:

  Neighbor graph (`"euclidean"`, `"cosine"`, `"manhattan"`,
  `"hamming"`).

- cluster_algorithm, cluster_resolution:

  Clustering (`"louvain"`, `"slm"`, `"leiden"`). Larger
  `cluster_resolution` yields fewer clusters.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

- integration_method:

  Deprecated alias of `integration_methods`.

- ...:

  Passed to the integration method functions.

## Value

A `Seurat` object. For `ChromatinAssay`, names follow the ATAC
convention (`*lsi`, `*UMAP2D`, cluster aliases, `ATAC_default_*`).

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
[ComBat_integrate](https://mengxu98.github.io/scop/reference/ComBat_integrate.md),
[RunIntegrationBenchmark](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md),
[IntegrationBenchmarkPlot](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10,
  nonlinear_reduction_dims = 2,
  compute_lisi = TRUE,
  lisi_label_colnames = "tech",
  lisi_perplexity = 10
)
#> ◌ [2026-10-29 12:26:25] Run integration workflow...
#> ℹ [2026-10-29 12:26:26] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:26:26] Checking a list of <Seurat>...
#> ! [2026-10-29 12:26:26] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-29 12:26:26] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:26] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-10-29 12:26:27] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-29 12:26:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:27] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-10-29 12:26:27] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-29 12:26:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:27] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-10-29 12:26:27] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-29 12:26:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:27] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-10-29 12:26:27] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-29 12:26:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:27] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:27] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:26:27] Number of available HVF: 500
#> ℹ [2026-10-29 12:26:27] Finished check
#> ℹ [2026-10-29 12:26:30] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-29 12:26:30] Perform linear dimension reduction("pca")
#> ℹ [2026-10-29 12:26:30] Perform Harmony integration
#> ℹ [2026-10-29 12:26:30] Using "Harmonypca" (1:10) as input
#> ℹ [2026-10-29 12:26:31] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:26:31] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:26:31] Reorder clusters...
#> ℹ [2026-10-29 12:26:31] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:26:31] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ℹ [2026-10-29 12:26:38] Perform umap nonlinear dimension reduction using Harmonypca (1:10)
#> ℹ [2026-10-29 12:26:45] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-29 12:26:45] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-29 12:26:45] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-29 12:26:45] Harmony integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "HarmonyUMAP2D"
)


IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")


panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "LIGER"
)
#> ◌ [2026-10-29 12:26:46] Run integration workflow...
#> ℹ [2026-10-29 12:26:46] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:26:47] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:26:48] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:26:48] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:48] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:26:48] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:48] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:26:48] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:48] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:26:48] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:48] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:26:48] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:26:49] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:26:49] Number of available HVF: 2000
#> ℹ [2026-10-29 12:26:49] Finished check
#> ℹ [2026-10-29 12:26:49] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-29 12:26:49] Perform LIGER integration
#> ℹ [2026-10-29 12:26:54] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:26:54] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:26:55] Reorder clusters...
#> ℹ [2026-10-29 12:26:55] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:26:55] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-29 12:27:01] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-29 12:27:08] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ✔ [2026-10-29 12:27:15] LIGER integration completed
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  compute_lisi = TRUE,
  lisi_label_colnames = "tech"
)
#> ◌ [2026-10-29 12:27:15] Run integration workflow...
#> ℹ [2026-10-29 12:27:16] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:27:17] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:27:17] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:17] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:18] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:18] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:18] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:18] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:18] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:18] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:18] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:18] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:18] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:27:18] Number of available HVF: 2000
#> ℹ [2026-10-29 12:27:18] Finished check
#> ℹ [2026-10-29 12:27:18] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-29 12:27:19] Perform linear dimension reduction("pca")
#> ℹ [2026-10-29 12:27:19] Perform Harmony integration
#> ℹ [2026-10-29 12:27:19] Using "Harmonypca" (1:20) as input
#> ℹ [2026-10-29 12:27:20] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:27:20] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:27:20] Reorder clusters...
#> ℹ [2026-10-29 12:27:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:27:20] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-29 12:27:27] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-29 12:27:34] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ℹ [2026-10-29 12:27:40] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-29 12:27:40] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-29 12:27:40] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-29 12:27:41] Harmony integration completed
IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")


data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
pbmcmultiome_sub <- RunIntegration(
  pbmcmultiome_sub,
  batch = "batch",
  assay = "peaks",
  integration_methods = "Harmony5",
  normalization_method = "TFIDF"
)
#> ◌ [2026-10-29 12:27:41] Run integration workflow...
#> ! [2026-10-29 12:27:41] `integration_method = 'Harmony5'` is not compatible with <ChromatinAssay> in current Seurat v5 workflow. Automatically switch to "Harmony"
#> ℹ [2026-10-29 12:27:41] Split `srt_merge` into `srt_list` by "batch"
#> ℹ [2026-10-29 12:27:41] Checking a list of <Seurat>...
#> ! [2026-10-29 12:27:42] Data 1/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-29 12:27:42] Perform `RunTFIDF()` on 1/2 of `srt_list`...
#> ℹ [2026-10-29 12:27:42] Perform `FindTopFeatures()` on 1/2 of `srt_list`...
#> ! [2026-10-29 12:27:42] Data 2/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-29 12:27:42] Perform `RunTFIDF()` on 2/2 of `srt_list`...
#> ℹ [2026-10-29 12:27:42] Perform `FindTopFeatures()` on 2/2 of `srt_list`...
#> ℹ [2026-10-29 12:27:42] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:27:42] Number of available HVF: 2000
#> ℹ [2026-10-29 12:27:42] Finished check
#> ℹ [2026-10-29 12:27:42] `normalization_method` is "TFIDF". Use lsi workflow...
#> ℹ [2026-10-29 12:27:42] Perform linear dimension reduction("svd")
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-10-29 12:27:43] Perform Harmony integration
#> ℹ [2026-10-29 12:27:43] Using "Harmonysvd" (2:5) as input
#> ℹ [2026-10-29 12:27:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:27:43] Reorder clusters...
#> ℹ [2026-10-29 12:27:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:27:43] Perform umap nonlinear dimension reduction using Harmony (1:4)
#> ℹ [2026-10-29 12:27:49] Perform umap nonlinear dimension reduction using Harmony (1:4)
#> ✔ [2026-10-29 12:27:55] Harmony integration completed

integration_methods <- c(
  "Uncorrected", "Seurat", "CCA", "RPCA",
  "MNN", "fastMNN", "Harmony", "Harmony5",
  "Coralysis", "LIGER", "Conos", "ComBat"
)
p_list <- list()
for (method in integration_methods) {
  panc8_sub <- RunIntegration(
    panc8_sub,
    batch = "tech",
    integration_methods = method,
    linear_reduction_dims_use = 1:50,
    nonlinear_reduction = "umap",
    RunParallelDivisiveICP_params = list(L = 10L)
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
#> ◌ [2026-10-29 12:27:55] Run integration workflow...
#> ℹ [2026-10-29 12:27:56] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:27:56] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:27:57] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:57] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:57] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:57] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:57] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:57] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:57] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:57] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:58] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:27:58] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:27:58] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:27:58] Number of available HVF: 2000
#> ℹ [2026-10-29 12:27:58] Finished check
#> ℹ [2026-10-29 12:27:58] Perform Uncorrected integration
#> ℹ [2026-10-29 12:27:58] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-29 12:27:58] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:27:59] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:27:59] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:27:59] Reorder clusters...
#> ℹ [2026-10-29 12:27:59] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:27:59] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-29 12:28:06] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-29 12:28:12] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-10-29 12:28:19] Uncorrected integration completed
#> ◌ [2026-10-29 12:28:19] Run integration workflow...
#> ℹ [2026-10-29 12:28:20] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:28:21] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:28:21] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:28:21] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:28:21] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:28:21] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:28:22] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:28:22] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:28:22] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:28:22] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:28:22] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:28:22] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:28:22] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:28:22] Number of available HVF: 2000
#> ℹ [2026-10-29 12:28:22] Finished check
#> ℹ [2026-10-29 12:28:22] Perform FindIntegrationAnchors
#> ℹ [2026-10-29 12:28:45] Perform Seurat integration
#> ℹ [2026-10-29 12:29:03] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-29 12:29:03] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:29:04] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:29:04] Reorder clusters...
#> ℹ [2026-10-29 12:29:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:29:04] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:29:11] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:29:18] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-29 12:29:27] Seurat integration completed
#> ◌ [2026-10-29 12:29:28] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-29 12:29:28] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-29 12:29:28] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-29 12:29:29] Number of available HVF: 2000
#> ℹ [2026-10-29 12:29:30] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-29 12:29:30] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-29 12:29:31] Perform Seurat v5 integration with `CCAIntegration()`
#> ℹ [2026-10-29 12:29:37] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:29:38] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:29:38] Reorder clusters...
#> ℹ [2026-10-29 12:29:38] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:29:38] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-29 12:29:45] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-29 12:29:51] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-29 12:30:00] CCA integration completed
#> ◌ [2026-10-29 12:30:00] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-29 12:30:01] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-29 12:30:01] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-29 12:30:01] Number of available HVF: 2000
#> ℹ [2026-10-29 12:30:02] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-29 12:30:02] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-29 12:30:03] Perform Seurat v5 integration with `RPCAIntegration()`
#> ℹ [2026-10-29 12:30:14] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:30:15] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:30:15] Reorder clusters...
#> ℹ [2026-10-29 12:30:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:30:15] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-29 12:30:22] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-29 12:30:28] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-29 12:30:37] RPCA integration completed
#> ◌ [2026-10-29 12:30:37] Run integration workflow...
#> ℹ [2026-10-29 12:30:38] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:30:38] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:30:38] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:30:39] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:30:39] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:30:39] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:30:39] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:30:39] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:30:39] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:30:39] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:30:39] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:30:40] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:30:40] Number of available HVF: 2000
#> ℹ [2026-10-29 12:30:40] Finished check
#> ℹ [2026-10-29 12:30:40] Perform MNN integration
#> ℹ [2026-10-29 12:31:23] Perform ScaleData
#> ℹ [2026-10-29 12:31:23] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:31:24] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:31:24] Reorder clusters...
#> ℹ [2026-10-29 12:31:24] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:31:24] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-29 12:31:31] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-29 12:31:38] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ✔ [2026-10-29 12:31:47] MNN integration completed
#> ◌ [2026-10-29 12:31:47] Run integration workflow...
#> ℹ [2026-10-29 12:31:49] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:31:49] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:31:49] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:31:49] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:31:49] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:31:50] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:31:50] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:31:50] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:31:50] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:31:50] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:31:50] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:31:50] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:31:50] Number of available HVF: 2000
#> ℹ [2026-10-29 12:31:50] Finished check
#> ℹ [2026-10-29 12:31:51] Perform fastMNN integration
#> ℹ [2026-10-29 12:31:55] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:31:55] Reorder clusters...
#> ℹ [2026-10-29 12:31:55] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:31:55] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-29 12:32:02] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-29 12:32:09] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ✔ [2026-10-29 12:32:17] fastMNN integration completed
#> ◌ [2026-10-29 12:32:17] Run integration workflow...
#> ℹ [2026-10-29 12:32:19] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:32:20] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:32:20] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:32:20] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:32:20] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:32:20] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:32:21] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:32:21] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:32:21] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:32:21] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:32:21] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:32:21] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:32:21] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:32:21] Number of available HVF: 2000
#> ℹ [2026-10-29 12:32:21] Finished check
#> ℹ [2026-10-29 12:32:21] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-29 12:32:22] Perform linear dimension reduction("pca")
#> ℹ [2026-10-29 12:32:22] Perform Harmony integration
#> ℹ [2026-10-29 12:32:22] Using "Harmonypca" (1:50) as input
#> ℹ [2026-10-29 12:32:23] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:32:23] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:32:23] Reorder clusters...
#> ℹ [2026-10-29 12:32:23] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:32:23] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-29 12:32:30] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-29 12:32:37] Perform umap nonlinear dimension reduction using Harmonypca (1:50)
#> ✔ [2026-10-29 12:32:45] Harmony integration completed
#> ◌ [2026-10-29 12:32:46] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-29 12:32:47] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-29 12:32:47] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-29 12:32:48] Number of available HVF: 2000
#> ℹ [2026-10-29 12:32:49] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-29 12:32:49] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-29 12:32:50] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ℹ [2026-10-29 12:32:51] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:32:51] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:32:51] Reorder clusters...
#> ℹ [2026-10-29 12:32:51] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:32:51] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-29 12:32:58] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-29 12:33:05] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-29 12:33:13] Harmony5 integration completed
#> ◌ [2026-10-29 12:33:14] Run integration workflow...
#> ℹ [2026-10-29 12:33:15] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:33:16] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:33:16] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:33:16] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:33:17] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:33:17] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:33:17] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:33:17] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:33:17] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:33:17] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:33:17] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:33:17] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:33:17] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:33:17] Number of available HVF: 2000
#> ℹ [2026-10-29 12:33:18] Finished check
#> ℹ [2026-10-29 12:33:18] Perform Coralysis integration
#> Data in `logcounts` slot already of `dgCMatrix` class...
#> 2000/2000 features remain after filtering features with only zero values.
#> 
#> Building training set...
#> No. of cells is lower or equal than the no. of 'nclusters': 200 <= 500
#> All the cells will be included and the clustering step will be skipped!
#> Found more than one class "kcca" in cache; using the first, from namespace 'kernlab'
#> Also defined by ‘flexclust’
#> Found more than one class "kcca" in cache; using the first, from namespace 'kernlab'
#> Also defined by ‘flexclust’
#> No. of cells is lower or equal than the no. of 'nclusters': 200 <= 500
#> All the cells will be included and the clustering step will be skipped!
#> No. of cells is lower or equal than the no. of 'nclusters': 200 <= 500
#> All the cells will be included and the clustering step will be skipped!
#> No. of cells is lower or equal than the no. of 'nclusters': 200 <= 500
#> All the cells will be included and the clustering step will be skipped!
#> Training set successfully built.
#> 
#> Computing cluster seed.
#> 
#> Initializing divisive ICP clustering...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%
#> 
#> 
#> Divisive ICP clustering completed successfully.
#> 
#> Predicting cell cluster probabilities using ICP models...
#> Prediction of cell cluster probabilities completed successfully.
#> 
#> Multi-level integration completed successfully.
#> Divisive ICP: selecting ICP tables multiple of 4
#> ℹ [2026-10-29 12:34:46] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:34:46] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:34:47] Reorder clusters...
#> ℹ [2026-10-29 12:34:47] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:34:47] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-29 12:34:54] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-29 12:35:01] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ✔ [2026-10-29 12:35:09] Coralysis integration completed
#> ◌ [2026-10-29 12:35:09] Run integration workflow...
#> ℹ [2026-10-29 12:35:11] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:35:12] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:35:12] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:12] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:13] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:13] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:13] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:13] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:13] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:13] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:13] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:13] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:14] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:35:14] Number of available HVF: 2000
#> ℹ [2026-10-29 12:35:14] Finished check
#> ℹ [2026-10-29 12:35:14] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-29 12:35:14] Perform LIGER integration
#> ℹ [2026-10-29 12:35:19] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:35:19] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:35:20] Reorder clusters...
#> ℹ [2026-10-29 12:35:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:35:20] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-29 12:35:27] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-29 12:35:34] Perform umap nonlinear dimension reduction using LIGER (1:50)
#> ✔ [2026-10-29 12:35:43] LIGER integration completed
#> ◌ [2026-10-29 12:35:43] Run integration workflow...
#> ℹ [2026-10-29 12:35:45] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:35:46] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:35:46] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:46] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:46] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:46] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:46] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:46] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:46] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:46] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:47] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:35:47] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:35:47] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:35:47] Number of available HVF: 2000
#> ℹ [2026-10-29 12:35:47] Finished check
#> ℹ [2026-10-29 12:35:47] Perform ScaleData on the data 1 ...
#> ℹ [2026-10-29 12:35:47] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:35:47] Perform ScaleData on the data 2 ...
#> ℹ [2026-10-29 12:35:47] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:35:47] Perform ScaleData on the data 3 ...
#> ℹ [2026-10-29 12:35:48] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:35:48] Perform ScaleData on the data 4 ...
#> ℹ [2026-10-29 12:35:48] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:35:48] Perform ScaleData on the data 5 ...
#> ℹ [2026-10-29 12:35:48] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:35:49]  Perform Conos integration
#> ℹ [2026-10-29 12:35:49] Conos integration using pca (1:50) as input
#> ℹ [2026-10-29 12:35:53] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-29 12:35:53] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:35:54] Reorder clusters...
#> ℹ [2026-10-29 12:35:54] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:35:54] Perform umap nonlinear dimension reduction using Conos
#> ℹ [2026-10-29 12:36:01] Perform umap nonlinear dimension reduction using Conos
#> ℹ [2026-10-29 12:36:08] Perform umap nonlinear dimension reduction using ConosUMAP2D (1:50)
#> ✔ [2026-10-29 12:36:17] Conos integration completed
#> ◌ [2026-10-29 12:36:17] Run integration workflow...
#> ℹ [2026-10-29 12:36:19] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:36:21] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:36:21] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:21] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:21] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:21] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:21] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:21] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:22] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:22] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:22] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:22] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:22] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:36:22] Number of available HVF: 2000
#> ℹ [2026-10-29 12:36:22] Finished check
#> ℹ [2026-10-29 12:36:22] Perform ComBat integration
#> Found 138 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
#> ℹ [2026-10-29 12:36:23] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-29 12:36:23] Perform linear dimension reduction ("pca")
#> ℹ [2026-10-29 12:36:24] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:36:24] Reorder clusters...
#> ℹ [2026-10-29 12:36:24] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:36:24] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ℹ [2026-10-29 12:36:31] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ℹ [2026-10-29 12:36:38] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ✔ [2026-10-29 12:36:47] ComBat integration completed

# Python-backed methods prepare a scVI/scvi-tools environment and run model
# training, so keep them separate from ordinary example checks.
if (FALSE) { # \dontrun{
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "scVI",
  train_params = list(max_epochs = 2L),
  nonlinear_reduction = "umap"
)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "scVI5",
  IntegrateLayers_params = list(max_epochs = 2L),
  nonlinear_reduction = "umap"
)
} # }

nonlinear_reductions <- c(
  "umap", "tsne", "dm", "phate",
  "pacmap", "trimap", "largevis", "fr"
)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Seurat",
  linear_reduction_dims_use = 1:50,
  nonlinear_reduction = nonlinear_reductions,
  nonlinear_reduction_params = list(largevis = list(n_epochs = 1000L))
)
#> ◌ [2026-10-29 12:36:47] Run integration workflow...
#> ℹ [2026-10-29 12:36:50] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-29 12:36:51] Checking a list of <Seurat>...
#> ℹ [2026-10-29 12:36:51] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:51] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:51] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:51] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:51] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:51] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:52] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:52] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:52] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-29 12:36:52] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-29 12:36:52] Use the separate HVF from `srt_list`
#> ℹ [2026-10-29 12:36:52] Number of available HVF: 2000
#> ℹ [2026-10-29 12:36:52] Finished check
#> ℹ [2026-10-29 12:36:52] Perform FindIntegrationAnchors
#> ℹ [2026-10-29 12:37:27] Perform Seurat integration
#> ℹ [2026-10-29 12:38:20] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-29 12:38:20] Perform "pca" linear dimension reduction
#> ℹ [2026-10-29 12:38:21] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-29 12:38:21] Reorder clusters...
#> ℹ [2026-10-29 12:38:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-29 12:38:21] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:38:28] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:38:35] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:38:38] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:38:43] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:12] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:13] Perform phate nonlinear dimension reduction using Seuratpca (1:50)
#> ◌ [2026-10-29 12:40:16] Preparing python environment...
#> ℹ [2026-10-29 12:40:16] Environment name: 'scprep_env' and python version: 3.10-1
#> ℹ [2026-10-29 12:40:16] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-10-29 12:40:16] Using micromamba executable: '/home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba'
#> ℹ [2026-10-29 12:40:17] Creating Python environment with python "3.10-1" using micromamba...
#> ℹ [2026-10-29 12:40:17] Skipping conda Terms of Service acceptance for micromamba.
#> + /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba 'create' '--yes' '--name' 'scprep_env' 'python=3.10-1' 'pip' 'setuptools<81' 'wheel' '--quiet' '-c' 'conda-forge'
#> ✔ [2026-10-29 12:40:20] Environment created successfully: '/home/runner/work/_temp/pkgdown-python/envs/scprep_env'
#> ℹ [2026-10-29 12:40:20] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-10-29 12:40:20] Checking 3 conda requirements
#> ℹ [2026-10-29 12:40:20] Installing 3 packages into environment: 'scprep_env' using micromamba
#> ℹ [2026-10-29 12:40:20] Installing packages via micromamba...
#> ℹ [2026-10-29 12:40:20] Using channel: "conda-forge"
#> ℹ [2026-10-29 12:40:20] Installing 3 packages...
#> + /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba 'install' '--yes' '--name' 'scprep_env' '-c' 'conda-forge' 'leidenalg=0.10.2' 'tbb=2022.2.0' 'python-igraph=0.11.9'
#> Using Cached Shard Index for conda-forge/linux-64                                                   ✔ Done
#> Using Cached Shard Index for conda-forge/noarch                                                     ✔ Done
#> Fetching and Parsing Packages' Shards                                                           ⧖ Starting
#> Fetching and Parsing Packages' Shards                                                     ✔ Done (0.2 sec)
#> 
#> Pinned packages:
#> 
#>   - python=3.10
#> 
#> Resolving Environment                                                                           ⧖ Starting
#> Resolving Environment                                                                     ✔ Done (0.3 sec)
#> 
#> Transaction
#> 
#>   Prefix: /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> 
#>   Updating specs:
#> 
#>    - leidenalg=0.10.2
#>    - tbb=2022.2.0
#>    - python-igraph=0.11.9
#> 
#> 
#>   Package           Version  Build                  Channel           Size
#> ────────────────────────────────────────────────────────────────────────────
#>   Install:
#> ────────────────────────────────────────────────────────────────────────────
#> 
#>   + arpack            3.9.1  nompi_hf03ea27_102     conda-forge     Cached
#>   + glpk                5.0  h445213a_0             conda-forge     Cached
#>   + gmp               6.3.0  hfd2156b_3             conda-forge     Cached
#>   + igraph          0.10.16  hfe3e89f_1             conda-forge     Cached
#>   + leidenalg        0.10.2  py310hea6c23e_1        conda-forge     Cached
#>   + libblas          3.11.0  11_h4a7cf45_openblas   conda-forge     Cached
#>   + libgfortran      16.2.0  h69a702a_7             conda-forge     Cached
#>   + libgfortran5     16.2.0  h6b99dfc_7             conda-forge     Cached
#>   + libhwloc         2.12.1  default_hafda6a7_1003  conda-forge     Cached
#>   + libiconv           1.18  h0cb94f2_3             conda-forge     Cached
#>   + liblapack        3.11.0  11_h47877c9_openblas   conda-forge     Cached
#>   + libleidenalg     0.11.1  h00ab1b0_0             conda-forge     Cached
#>   + libopenblas      0.3.34  pthreads_hf13c14d_2    conda-forge     Cached
#>   + libstdcxx        16.2.0  h934c35e_7             conda-forge     Cached
#>   + libstdcxx-ng     16.2.0  hdf11a46_7             conda-forge     Cached
#>   + libxml2          2.15.4  ha7a5771_1             conda-forge     Cached
#>   + libxml2-16       2.15.4  hbdfff7e_1             conda-forge     Cached
#>   + libxml2-devel    2.15.4  ha7a5771_1             conda-forge     Cached
#>   + python-igraph    0.11.9  py310hea6c23e_3        conda-forge     Cached
#>   + python_abi         3.10  9_cp310                conda-forge     Cached
#>   + tbb            2022.2.0  hb60516a_1             conda-forge     Cached
#>   + texttable         1.7.0  pyhd8ed1ab_1           conda-forge     Cached
#> 
#>   Summary:
#> 
#>   Install: 22 packages
#> 
#>   Total download: 0 B
#> 
#> ─────────────────────────────────────────────────warning  libmamba Security Warning: This transaction includes executing package scripts (pre/post-link/unlink) if present. These scripts can contain arbitrary code. Please ensure you trust the package sources.
#> ───────────────────────────
#> 
#> 
#> 
#> Transaction starting
#> Linking libstdcxx-16.2.0-h934c35e_7
#> Linking libiconv-1.18-h0cb94f2_3
#> Linking libgfortran5-16.2.0-h6b99dfc_7
#> Linking gmp-6.3.0-hfd2156b_3
#> Linking libstdcxx-ng-16.2.0-hdf11a46_7
#> Linking libxml2-16-2.15.4-hbdfff7e_1
#> Linking libgfortran-16.2.0-h69a702a_7
#> Linking glpk-5.0-h445213a_0
#> Linking libxml2-2.15.4-ha7a5771_1
#> Linking libopenblas-0.3.34-pthreads_hf13c14d_2
#> Linking libhwloc-2.12.1-default_hafda6a7_1003
#> Linking libxml2-devel-2.15.4-ha7a5771_1
#> Linking libblas-3.11.0-11_h4a7cf45_openblas
#> Linking tbb-2022.2.0-hb60516a_1
#> Linking liblapack-3.11.0-11_h47877c9_openblas
#> Linking arpack-3.9.1-nompi_hf03ea27_102
#> Linking igraph-0.10.16-hfe3e89f_1
#> Linking libleidenalg-0.11.1-h00ab1b0_0
#> Linking python_abi-3.10-9_cp310
#> Linking texttable-1.7.0-pyhd8ed1ab_1
#> Linking python-igraph-0.11.9-py310hea6c23e_3
#> Linking leidenalg-0.10.2-py310hea6c23e_1
#> 
#> Transaction finished
#> 
#> ✔ [2026-10-29 12:40:21] micromamba installation completed successfully
#> ℹ [2026-10-29 12:40:22] Checking 12 requirements using uv
#> ℹ [2026-10-29 12:40:22] Installing 12 packages into environment: 'scprep_env' using micromamba
#> ℹ [2026-10-29 12:40:22] Installing packages via uv...
#> + /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/uv 'pip' 'install' '--python' '/home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10' '--force-reinstall' 'setuptools<81' 'matplotlib==3.10.8' 'numba==0.66.0' 'llvmlite==0.48.0' 'numpy==1.26.4' 'packaging>=24.0' 'pandas==2.0.3' 'scikit-learn==1.7.0' 'scipy==1.15.3' 'phate==1.0.11' 'leidenalg==0.10.2' 'igraph==0.11.9'
#> Using Python 3.10.1 environment at: /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Resolved 35 packages in 94ms
#> Downloading pandas (11.7MiB)
#> Downloading pygsp (1.8MiB)
#>  Downloaded pygsp
#>  Downloaded pandas
#> Prepared 35 packages in 187ms
#> Uninstalled 5 packages in 16ms
#> Installed 35 packages in 62ms
#>  + cloudpickle==3.1.2
#>  + contourpy==1.3.2
#>  + cycler==0.12.1
#>  + decorator==5.3.1
#>  + deprecated==1.3.1
#>  + fonttools==4.65.0
#>  + future==1.0.0
#>  + graphtools==2.1.0
#>  - igraph==0.11.9 (from file:///home/conda/feedstock_root/build_artifacts/python-igraph_1761728289959/work)
#>  + igraph==0.11.9
#>  + joblib==1.6.0
#>  + kiwisolver==1.5.1
#>  - leidenalg==0.10.2 (from file:///home/conda/feedstock_root/build_artifacts/leidenalg_1756310397572/work)
#>  + leidenalg==0.10.2
#>  + llvmlite==0.48.0
#>  + matplotlib==3.10.8
#>  + numba==0.66.0
#>  + numpy==1.26.4
#>  - packaging==26.3 (from file:///home/conda/feedstock_root/build_artifacts/bld/rattler-build_packaging_1785888127/work)
#>  + packaging==26.3
#>  + pandas==2.0.3
#>  + phate==1.0.11
#>  + pillow==12.3.0
#>  + pygsp==0.6.1
#>  + pyparsing==3.3.3
#>  + python-dateutil==3.0.0
#>  + pytz==2026.6
#>  + s-gd2==1.8.1
#>  + scikit-learn==1.7.0
#>  + scipy==1.15.3
#>  + scprep==1.2.3
#>  ~ setuptools==80.10.2
#>  + six==1.17.0
#>  + tasklogger==1.2.0
#>  - texttable==1.7.0 (from file:///home/conda/feedstock_root/build_artifacts/texttable_1734145585512/work)
#>  + texttable==1.7.2
#>  + threadpoolctl==3.7.0
#>  + tzdata==2026.6
#>  + wrapt==2.5.1
#> ✔ [2026-10-29 12:40:22] uv installation completed
#> + /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10 -m pip check
#> No broken requirements found.
#> ✔ [2026-10-29 12:40:25] Python environment ready
#> ✔                       Until next loading, the environment will be cached
#> Environment manager config:
#>  Manager:       micromamba
#>  Executable:    /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#>  Environment:   /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Python config:
#>  python:        /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10
#>  version:       Python 3.10.1
#> 
#> ℹ [2026-10-29 12:40:37] Perform phate nonlinear dimension reduction using Seuratpca (1:50)
#> ◌ [2026-10-29 12:40:40] Preparing python environment...
#> ℹ [2026-10-29 12:40:40] Environment name: 'scprep_env' and python version: 3.10-1
#> ℹ [2026-10-29 12:40:40] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-10-29 12:40:40] Using micromamba executable: '/home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba'
#> ℹ [2026-10-29 12:40:40] Using existing environment: '/home/runner/work/_temp/pkgdown-python/envs/scprep_env'
#> ℹ [2026-10-29 12:40:40] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-10-29 12:40:40] Checking 3 conda requirements
#> ℹ [2026-10-29 12:40:41] Checking 12 requirements using uv
#> + /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10 -m pip check
#> No broken requirements found.
#> ✔ [2026-10-29 12:40:43] Python environment ready
#> ✔                       Until next loading, the environment will be cached
#> Environment manager config:
#>  Manager:       micromamba
#>  Executable:    /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#>  Environment:   /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Python config:
#>  python:        /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10
#>  version:       Python 3.10.1
#> 
#> ℹ [2026-10-29 12:40:55] Perform pacmap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:56] Perform pacmap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:56] Perform trimap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:57] Perform trimap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:40:57] Perform largevis nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:41:07] Perform largevis nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-29 12:41:16] Perform fr nonlinear dimension reduction using Seurat_SNN
#> ℹ [2026-10-29 12:41:16] Perform fr nonlinear dimension reduction using Seurat_SNN
#> ℹ [2026-10-29 12:41:20] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-29 12:41:39] Seurat integration completed
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







```
