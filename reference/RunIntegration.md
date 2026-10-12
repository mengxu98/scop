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

- nonlinear_reduction, nonlinear_reduction_dims,
  nonlinear_reduction_params, force_nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

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
#> ◌ [2026-10-11 23:21:07] Run integration workflow...
#> ℹ [2026-10-11 23:21:08] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:21:08] Checking a list of <Seurat>...
#> ! [2026-10-11 23:21:08] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-11 23:21:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:08] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-10-11 23:21:08] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-11 23:21:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:08] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-10-11 23:21:09] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-11 23:21:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:09] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-10-11 23:21:09] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-11 23:21:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:09] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-10-11 23:21:09] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-11 23:21:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:09] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:09] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:21:09] Number of available HVF: 500
#> ℹ [2026-10-11 23:21:09] Finished check
#> ℹ [2026-10-11 23:21:13] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-11 23:21:13] Perform linear dimension reduction("pca")
#> ℹ [2026-10-11 23:21:13] Perform Harmony integration
#> ℹ [2026-10-11 23:21:13] Using "Harmonypca" (1:10) as input
#> ℹ [2026-10-11 23:21:14] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:21:14] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:21:14] Reorder clusters...
#> ℹ [2026-10-11 23:21:14] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:21:14] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ℹ [2026-10-11 23:21:20] Perform umap nonlinear dimension reduction using Harmonypca (1:10)
#> ℹ [2026-10-11 23:21:26] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-11 23:21:26] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-11 23:21:26] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-11 23:21:26] Harmony integration completed
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
#> ◌ [2026-10-11 23:21:27] Run integration workflow...
#> ℹ [2026-10-11 23:21:27] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:21:28] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:21:28] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:28] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:29] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:29] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:29] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:29] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:29] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:29] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:29] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:29] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:29] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:21:30] Number of available HVF: 2000
#> ℹ [2026-10-11 23:21:30] Finished check
#> ℹ [2026-10-11 23:21:30] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-11 23:21:30] Perform LIGER integration
#> ℹ [2026-10-11 23:21:35] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:21:36] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:21:36] Reorder clusters...
#> ℹ [2026-10-11 23:21:36] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:21:36] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-11 23:21:42] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-11 23:21:47] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ✔ [2026-10-11 23:21:54] LIGER integration completed
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  compute_lisi = TRUE,
  lisi_label_colnames = "tech"
)
#> ◌ [2026-10-11 23:21:54] Run integration workflow...
#> ℹ [2026-10-11 23:21:54] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:21:55] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:21:55] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:55] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:56] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:56] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:56] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:56] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:56] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:56] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:56] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:21:56] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:21:56] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:21:56] Number of available HVF: 2000
#> ℹ [2026-10-11 23:21:57] Finished check
#> ℹ [2026-10-11 23:21:57] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-11 23:21:57] Perform linear dimension reduction("pca")
#> ℹ [2026-10-11 23:21:58] Perform Harmony integration
#> ℹ [2026-10-11 23:21:58] Using "Harmonypca" (1:20) as input
#> ℹ [2026-10-11 23:21:58] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:21:58] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:21:58] Reorder clusters...
#> ℹ [2026-10-11 23:21:58] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:21:58] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-11 23:22:04] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-11 23:22:10] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ℹ [2026-10-11 23:22:16] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-11 23:22:16] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-11 23:22:16] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-11 23:22:16] Harmony integration completed
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
#> ◌ [2026-10-11 23:22:17] Run integration workflow...
#> ! [2026-10-11 23:22:17] `integration_method = 'Harmony5'` is not compatible with <ChromatinAssay> in current Seurat v5 workflow. Automatically switch to "Harmony"
#> ℹ [2026-10-11 23:22:17] Split `srt_merge` into `srt_list` by "batch"
#> ℹ [2026-10-11 23:22:17] Checking a list of <Seurat>...
#> ! [2026-10-11 23:22:17] Data 1/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-11 23:22:17] Perform `RunTFIDF()` on 1/2 of `srt_list`...
#> ℹ [2026-10-11 23:22:17] Perform `FindTopFeatures()` on 1/2 of `srt_list`...
#> ! [2026-10-11 23:22:17] Data 2/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-11 23:22:17] Perform `RunTFIDF()` on 2/2 of `srt_list`...
#> ℹ [2026-10-11 23:22:17] Perform `FindTopFeatures()` on 2/2 of `srt_list`...
#> ℹ [2026-10-11 23:22:17] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:22:17] Number of available HVF: 2000
#> ℹ [2026-10-11 23:22:18] Finished check
#> ℹ [2026-10-11 23:22:18] `normalization_method` is "TFIDF". Use lsi workflow...
#> ℹ [2026-10-11 23:22:18] Perform linear dimension reduction("svd")
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-10-11 23:22:18] Perform Harmony integration
#> ℹ [2026-10-11 23:22:18] Using "Harmonysvd" (2:10) as input
#> ℹ [2026-10-11 23:22:18] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:22:18] Reorder clusters...
#> ℹ [2026-10-11 23:22:18] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:22:18] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ℹ [2026-10-11 23:22:23] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ✔ [2026-10-11 23:22:28] Harmony integration completed

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
#> ◌ [2026-10-11 23:22:28] Run integration workflow...
#> ℹ [2026-10-11 23:22:29] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:22:29] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:22:30] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:30] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:30] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:30] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:30] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:30] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:30] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:30] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:31] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:31] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:31] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:22:31] Number of available HVF: 2000
#> ℹ [2026-10-11 23:22:31] Finished check
#> ℹ [2026-10-11 23:22:31] Perform Uncorrected integration
#> ℹ [2026-10-11 23:22:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-11 23:22:31] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:22:32] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:22:32] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:22:32] Reorder clusters...
#> ℹ [2026-10-11 23:22:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:22:32] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-11 23:22:38] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-11 23:22:43] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-10-11 23:22:49] Uncorrected integration completed
#> ◌ [2026-10-11 23:22:49] Run integration workflow...
#> ℹ [2026-10-11 23:22:50] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:22:51] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:22:51] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:51] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:51] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:51] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:51] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:51] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:52] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:52] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:52] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:22:52] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:22:52] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:22:52] Number of available HVF: 2000
#> ℹ [2026-10-11 23:22:52] Finished check
#> ℹ [2026-10-11 23:22:52] Perform FindIntegrationAnchors
#> ℹ [2026-10-11 23:23:14] Perform Seurat integration
#> ℹ [2026-10-11 23:23:32] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-11 23:23:32] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:23:33] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:23:33] Reorder clusters...
#> ℹ [2026-10-11 23:23:33] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:23:34] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:23:39] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:23:45] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-11 23:23:53] Seurat integration completed
#> ◌ [2026-10-11 23:23:53] Run integration workflow...
#> ℹ [2026-10-11 23:23:54] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-11 23:23:54] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-11 23:23:54] Number of available HVF: 2000
#> ℹ [2026-10-11 23:23:56] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-11 23:23:56] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-11 23:23:56] Perform Seurat v5 integration with `CCAIntegration()`
#> ℹ [2026-10-11 23:24:03] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:24:03] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:24:03] Reorder clusters...
#> ℹ [2026-10-11 23:24:03] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:24:03] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-11 23:24:09] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-11 23:24:15] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-11 23:24:22] CCA integration completed
#> ◌ [2026-10-11 23:24:22] Run integration workflow...
#> ℹ [2026-10-11 23:24:23] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-11 23:24:23] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-11 23:24:23] Number of available HVF: 2000
#> ℹ [2026-10-11 23:24:24] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-11 23:24:24] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-11 23:24:25] Perform Seurat v5 integration with `RPCAIntegration()`
#> ℹ [2026-10-11 23:24:36] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:24:36] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:24:37] Reorder clusters...
#> ℹ [2026-10-11 23:24:37] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:24:37] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-11 23:24:42] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-11 23:24:48] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-11 23:24:55] RPCA integration completed
#> ◌ [2026-10-11 23:24:55] Run integration workflow...
#> ℹ [2026-10-11 23:24:56] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:24:57] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:24:57] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:24:57] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:24:57] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:24:57] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:24:57] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:24:57] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:24:57] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:24:58] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:24:58] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:24:58] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:24:58] Number of available HVF: 2000
#> ℹ [2026-10-11 23:24:58] Finished check
#> ℹ [2026-10-11 23:24:58] Perform MNN integration
#> ℹ [2026-10-11 23:25:40] Perform ScaleData
#> ℹ [2026-10-11 23:25:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:25:41] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:25:41] Reorder clusters...
#> ℹ [2026-10-11 23:25:41] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:25:41] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-11 23:25:47] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-11 23:25:53] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ✔ [2026-10-11 23:26:01] MNN integration completed
#> ◌ [2026-10-11 23:26:01] Run integration workflow...
#> ℹ [2026-10-11 23:26:02] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:26:03] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:03] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:03] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:03] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:03] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:03] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:03] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:03] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:04] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:04] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:04] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:26:04] Number of available HVF: 2000
#> ℹ [2026-10-11 23:26:04] Finished check
#> ℹ [2026-10-11 23:26:04] Perform fastMNN integration
#> ℹ [2026-10-11 23:26:08] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:26:08] Reorder clusters...
#> ℹ [2026-10-11 23:26:08] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:26:08] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-11 23:26:14] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-11 23:26:19] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ✔ [2026-10-11 23:26:27] fastMNN integration completed
#> ◌ [2026-10-11 23:26:27] Run integration workflow...
#> ℹ [2026-10-11 23:26:28] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:26:29] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:26:29] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:29] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:30] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:30] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:30] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:30] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:30] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:30] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:30] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:26:30] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:26:30] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:26:30] Number of available HVF: 2000
#> ℹ [2026-10-11 23:26:30] Finished check
#> ℹ [2026-10-11 23:26:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-11 23:26:31] Perform linear dimension reduction("pca")
#> ℹ [2026-10-11 23:26:31] Perform Harmony integration
#> ℹ [2026-10-11 23:26:31] Using "Harmonypca" (1:50) as input
#> ℹ [2026-10-11 23:26:32] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:26:32] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:26:32] Reorder clusters...
#> ℹ [2026-10-11 23:26:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:26:32] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-11 23:26:38] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-11 23:26:44] Perform umap nonlinear dimension reduction using Harmonypca (1:50)
#> ✔ [2026-10-11 23:26:51] Harmony integration completed
#> ◌ [2026-10-11 23:26:51] Run integration workflow...
#> ℹ [2026-10-11 23:26:52] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-11 23:26:52] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-11 23:26:53] Number of available HVF: 2000
#> ℹ [2026-10-11 23:26:54] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-11 23:26:54] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-11 23:26:55] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ℹ [2026-10-11 23:26:56] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:26:56] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:26:56] Reorder clusters...
#> ℹ [2026-10-11 23:26:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:26:56] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-11 23:27:02] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-11 23:27:08] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-11 23:27:15] Harmony5 integration completed
#> ◌ [2026-10-11 23:27:15] Run integration workflow...
#> ℹ [2026-10-11 23:27:17] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:27:18] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:27:18] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:27:18] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:27:18] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:27:18] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:27:18] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:27:18] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:27:18] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:27:18] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:27:19] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:27:19] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:27:19] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:27:19] Number of available HVF: 2000
#> ℹ [2026-10-11 23:27:19] Finished check
#> ℹ [2026-10-11 23:27:19] Perform Coralysis integration
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
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  76%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
#> 
#> 
#> Divisive ICP clustering completed successfully.
#> 
#> Predicting cell cluster probabilities using ICP models...
#> Prediction of cell cluster probabilities completed successfully.
#> 
#> Multi-level integration completed successfully.
#> Divisive ICP: selecting ICP tables multiple of 4
#> ℹ [2026-10-11 23:33:39] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:33:40] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:33:40] Reorder clusters...
#> ℹ [2026-10-11 23:33:40] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:33:40] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-11 23:33:45] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-11 23:33:51] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ✔ [2026-10-11 23:33:58] Coralysis integration completed
#> ◌ [2026-10-11 23:33:59] Run integration workflow...
#> ℹ [2026-10-11 23:34:00] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:34:01] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:34:02] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:02] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:02] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:02] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:02] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:02] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:02] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:02] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:02] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:02] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:03] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:34:03] Number of available HVF: 2000
#> ℹ [2026-10-11 23:34:03] Finished check
#> ℹ [2026-10-11 23:34:03] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-11 23:34:03] Perform LIGER integration
#> ℹ [2026-10-11 23:34:08] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-11 23:34:09] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:34:09] Reorder clusters...
#> ℹ [2026-10-11 23:34:09] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:34:09] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-11 23:34:15] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-11 23:34:21] Perform umap nonlinear dimension reduction using LIGER (1:50)
#> ✔ [2026-10-11 23:34:28] LIGER integration completed
#> ◌ [2026-10-11 23:34:28] Run integration workflow...
#> ℹ [2026-10-11 23:34:37] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:34:38] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:34:38] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:38] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:38] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:38] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:38] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:38] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:39] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:39] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:39] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:39] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:39] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:34:39] Number of available HVF: 2000
#> ℹ [2026-10-11 23:34:39] Finished check
#> ℹ [2026-10-11 23:34:39] Perform ScaleData on the data 1 ...
#> ℹ [2026-10-11 23:34:39] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:34:39] Perform ScaleData on the data 2 ...
#> ℹ [2026-10-11 23:34:39] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:34:40] Perform ScaleData on the data 3 ...
#> ℹ [2026-10-11 23:34:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:34:40] Perform ScaleData on the data 4 ...
#> ℹ [2026-10-11 23:34:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:34:40] Perform ScaleData on the data 5 ...
#> ℹ [2026-10-11 23:34:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:34:41]  Perform Conos integration
#> ℹ [2026-10-11 23:34:41] Conos integration using pca (1:50) as input
#> Error in sccore::plapply(which(is.na(mi)), function(i) {    if (space == "CPCA") {        xcp <- quickCPCA(self$samples[sn.pairs[, i]], data.type = data.type,             ncomps = ncomps, n.odgenes = n.odgenes, verbose = FALSE,             var.scale = var.scale, score.component.variance = score.component.variance)    }    else if (space == "JNMF") {        xcp <- quickJNMF(self$samples[sn.pairs[, i]], data.type = data.type,             n.comps = ncomps, n.odgenes = n.odgenes, var.scale = var.scale,             verbose = FALSE, max.iter = 3000)    }    else if (space == "genes") {        xcp <- quickNULL(p2.objs = self$samples[sn.pairs[, i]],             data.type = data.type, n.odgenes = n.odgenes, var.scale = var.scale,             verbose = FALSE)    }    else if (space == "PCA") {        xcp <- quickPlainPCA(self$samples[sn.pairs[, i]], data.type = data.type,             ncomps = ncomps, n.odgenes = n.odgenes, verbose = FALSE,             var.scale = var.scale, score.component.variance = score.component.variance)    }    else if (space == "CCA" || space == "PMA") {        xcp <- quickCCA(self$samples[sn.pairs[, i]], data.type = data.type,             ncomps = ncomps, n.odgenes = n.odgenes, verbose = FALSE,             var.scale = var.scale, score.component.variance = score.component.variance,             PMA = (space == "PMA"))    }    if (verbose)         cat(".")    xcp}, n.cores = self$n.cores, mc.preschedule = (space == "PCA"),     progress = FALSE, fail.on.error = TRUE): Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.
#> Errors in plapply: Error in (function (cond)  : 
#>   error in evaluating the argument 'x' in selecting a method for function 't': The `slot` argument of `GetAssayData()` was deprecated in SeuratObject
#> 5.0.0 and is now defunct.
#> ℹ Please use the `layer` argument instead.

# Python-backed methods prepare a scVI/scvi-tools environment and run model
# training, so keep them separate from ordinary example checks.
if (FALSE) { # \dontrun{
if (reticulate::py_module_available("scvi")) {
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
}
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
  nonlinear_reduction = nonlinear_reductions
)
#> ◌ [2026-10-11 23:34:42] Run integration workflow...
#> ℹ [2026-10-11 23:34:44] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-11 23:34:45] Checking a list of <Seurat>...
#> ℹ [2026-10-11 23:34:45] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:45] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:45] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:45] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:45] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:45] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:46] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:46] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:46] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-11 23:34:46] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-11 23:34:46] Use the separate HVF from `srt_list`
#> ℹ [2026-10-11 23:34:46] Number of available HVF: 2000
#> ℹ [2026-10-11 23:34:46] Finished check
#> ℹ [2026-10-11 23:34:46] Perform FindIntegrationAnchors
#> ℹ [2026-10-11 23:35:21] Perform Seurat integration
#> ℹ [2026-10-11 23:36:06] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-11 23:36:06] Perform "pca" linear dimension reduction
#> ℹ [2026-10-11 23:36:07] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-11 23:36:07] Reorder clusters...
#> ℹ [2026-10-11 23:36:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-11 23:36:07] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:36:13] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:36:19] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:36:22] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-11 23:36:28] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> Registered S3 methods overwritten by 'RcppEigen':
#>   method               from         
#>   predict.fastLm       RcppArmadillo
#>   print.fastLm         RcppArmadillo
#>   summary.fastLm       RcppArmadillo
#>   print.summary.fastLm RcppArmadillo
#> ! [2026-10-11 23:37:46] (converted from warning) 'as(<dsCMatrix>, "dgTMatrix")' is deprecated.
#> !                       Use 'as(as(., "generalMatrix"), "TsparseMatrix")' instead.
#> !                       See help("Deprecated") and help("Matrix-deprecated").
#> ! [2026-10-11 23:37:46] Error when performing nonlinear dimension reduction. Skip this step
#> ℹ [2026-10-11 23:37:46] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-11 23:38:01] Seurat integration completed
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
