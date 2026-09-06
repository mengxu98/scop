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
#> ◌ [2026-09-06 22:10:38] Run integration workflow...
#> ℹ [2026-09-06 22:10:38] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:10:39] Checking a list of <Seurat>...
#> ! [2026-09-06 22:10:39] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:10:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:39] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-09-06 22:10:39] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:10:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:39] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-09-06 22:10:39] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:10:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:40] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-09-06 22:10:40] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:10:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:40] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-09-06 22:10:40] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:10:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:40] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:40] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:10:40] Number of available HVF: 500
#> ℹ [2026-09-06 22:10:40] Finished check
#> ℹ [2026-09-06 22:10:42] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 22:10:43] Perform linear dimension reduction("pca")
#> ℹ [2026-09-06 22:10:43] Perform Harmony integration
#> ℹ [2026-09-06 22:10:43] Using "Harmonypca" (1:10) as input
#> ℹ [2026-09-06 22:10:43] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:10:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:10:44] Reorder clusters...
#> ℹ [2026-09-06 22:10:44] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:10:44] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ℹ [2026-09-06 22:10:49] Perform umap nonlinear dimension reduction using Harmonypca (1:10)
#> ℹ [2026-09-06 22:10:54] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-09-06 22:10:54] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-09-06 22:10:54] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-09-06 22:10:55] Harmony integration completed
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
#> ◌ [2026-09-06 22:10:55] Run integration workflow...
#> ℹ [2026-09-06 22:10:56] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:10:57] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:10:57] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:10:57] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:57] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:10:57] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:57] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:10:57] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:58] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:10:58] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:58] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:10:58] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:10:58] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:10:58] Number of available HVF: 2000
#> ℹ [2026-09-06 22:10:58] Finished check
#> ℹ [2026-09-06 22:10:58] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-09-06 22:10:58] Perform LIGER integration
#> ℹ [2026-09-06 22:11:03] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:11:04] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:11:04] Reorder clusters...
#> ℹ [2026-09-06 22:11:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:11:04] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-09-06 22:11:09] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-09-06 22:11:15] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ✔ [2026-09-06 22:11:20] LIGER integration completed
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  compute_lisi = TRUE,
  lisi_label_colnames = "tech"
)
#> ◌ [2026-09-06 22:11:20] Run integration workflow...
#> ℹ [2026-09-06 22:11:21] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:11:22] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:11:22] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:22] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:22] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:22] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:22] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:22] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:23] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:23] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:23] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:23] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:23] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:11:23] Number of available HVF: 2000
#> ℹ [2026-09-06 22:11:23] Finished check
#> ℹ [2026-09-06 22:11:23] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 22:11:23] Perform linear dimension reduction("pca")
#> ℹ [2026-09-06 22:11:24] Perform Harmony integration
#> ℹ [2026-09-06 22:11:24] Using "Harmonypca" (1:20) as input
#> ℹ [2026-09-06 22:11:24] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:11:24] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:11:25] Reorder clusters...
#> ℹ [2026-09-06 22:11:25] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:11:25] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-06 22:11:30] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-06 22:11:35] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ℹ [2026-09-06 22:11:40] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-09-06 22:11:40] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-09-06 22:11:40] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-09-06 22:11:41] Harmony integration completed
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
#> ◌ [2026-09-06 22:11:42] Run integration workflow...
#> ! [2026-09-06 22:11:42] `integration_method = 'Harmony5'` is not compatible with <ChromatinAssay> in current Seurat v5 workflow. Automatically switch to "Harmony"
#> ℹ [2026-09-06 22:11:42] Split `srt_merge` into `srt_list` by "batch"
#> ℹ [2026-09-06 22:11:42] Checking a list of <Seurat>...
#> ! [2026-09-06 22:11:42] Data 1/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-09-06 22:11:42] Perform `RunTFIDF()` on 1/2 of `srt_list`...
#> ℹ [2026-09-06 22:11:42] Perform `FindTopFeatures()` on 1/2 of `srt_list`...
#> ! [2026-09-06 22:11:42] Data 2/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-09-06 22:11:42] Perform `RunTFIDF()` on 2/2 of `srt_list`...
#> ℹ [2026-09-06 22:11:42] Perform `FindTopFeatures()` on 2/2 of `srt_list`...
#> ℹ [2026-09-06 22:11:42] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:11:42] Number of available HVF: 2000
#> ℹ [2026-09-06 22:11:42] Finished check
#> ℹ [2026-09-06 22:11:43] `normalization_method` is "TFIDF". Use lsi workflow...
#> ℹ [2026-09-06 22:11:43] Perform linear dimension reduction("svd")
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-09-06 22:11:43] Perform Harmony integration
#> ℹ [2026-09-06 22:11:43] Using "Harmonysvd" (2:10) as input
#> ℹ [2026-09-06 22:11:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:11:43] Reorder clusters...
#> ℹ [2026-09-06 22:11:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:11:43] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ℹ [2026-09-06 22:11:48] Perform umap nonlinear dimension reduction using Harmony (1:9)
#> ✔ [2026-09-06 22:11:52] Harmony integration completed

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
#> ◌ [2026-09-06 22:11:52] Run integration workflow...
#> ℹ [2026-09-06 22:11:53] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:11:54] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:11:54] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:54] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:54] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:54] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:54] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:54] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:55] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:55] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:55] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:11:55] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:11:55] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:11:55] Number of available HVF: 2000
#> ℹ [2026-09-06 22:11:55] Finished check
#> ℹ [2026-09-06 22:11:55] Perform Uncorrected integration
#> ℹ [2026-09-06 22:11:55] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 22:11:55] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:11:56] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:11:56] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:11:56] Reorder clusters...
#> ℹ [2026-09-06 22:11:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:11:56] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-09-06 22:12:01] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-09-06 22:12:07] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-09-06 22:12:12] Uncorrected integration completed
#> ◌ [2026-09-06 22:12:12] Run integration workflow...
#> ℹ [2026-09-06 22:12:13] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:12:14] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:12:14] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:12:14] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:12:14] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:12:14] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:12:14] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:12:14] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:12:15] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:12:15] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:12:15] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:12:15] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:12:15] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:12:15] Number of available HVF: 2000
#> ℹ [2026-09-06 22:12:15] Finished check
#> ℹ [2026-09-06 22:12:15] Perform FindIntegrationAnchors
#> ℹ [2026-09-06 22:12:43] Perform Seurat integration
#> ℹ [2026-09-06 22:13:00] Perform ScaleData on `srt_integrated`
#> ℹ [2026-09-06 22:13:00] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:13:01] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:13:01] Reorder clusters...
#> ℹ [2026-09-06 22:13:01] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:13:01] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:13:06] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:13:12] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-09-06 22:13:19] Seurat integration completed
#> ◌ [2026-09-06 22:13:19] Run integration workflow...
#> ℹ [2026-09-06 22:13:20] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-09-06 22:13:20] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-09-06 22:13:21] Number of available HVF: 2000
#> ℹ [2026-09-06 22:13:22] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-09-06 22:13:22] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-09-06 22:13:22] Perform Seurat v5 integration with `CCAIntegration()`
#> ℹ [2026-09-06 22:13:30] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:13:31] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:13:31] Reorder clusters...
#> ℹ [2026-09-06 22:13:31] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:13:31] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-09-06 22:13:37] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-09-06 22:13:42] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-09-06 22:13:48] CCA integration completed
#> ◌ [2026-09-06 22:13:48] Run integration workflow...
#> ℹ [2026-09-06 22:13:49] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-09-06 22:13:49] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-09-06 22:13:50] Number of available HVF: 2000
#> ℹ [2026-09-06 22:13:51] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-09-06 22:13:51] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-09-06 22:13:51] Perform Seurat v5 integration with `RPCAIntegration()`
#> ℹ [2026-09-06 22:14:03] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:14:04] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:14:04] Reorder clusters...
#> ℹ [2026-09-06 22:14:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:14:04] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-09-06 22:14:09] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-09-06 22:14:14] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-09-06 22:14:21] RPCA integration completed
#> ◌ [2026-09-06 22:14:21] Run integration workflow...
#> ℹ [2026-09-06 22:14:22] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:14:23] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:14:23] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:14:23] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:14:23] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:14:23] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:14:23] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:14:23] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:14:23] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:14:24] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:14:24] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:14:24] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:14:24] Number of available HVF: 2000
#> ℹ [2026-09-06 22:14:24] Finished check
#> ℹ [2026-09-06 22:14:24] Perform MNN integration
#> ℹ [2026-09-06 22:15:06] Perform ScaleData
#> ℹ [2026-09-06 22:15:06] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:15:07] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:15:07] Reorder clusters...
#> ℹ [2026-09-06 22:15:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:15:07] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-09-06 22:15:13] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-09-06 22:15:19] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ✔ [2026-09-06 22:15:26] MNN integration completed
#> ◌ [2026-09-06 22:15:26] Run integration workflow...
#> ℹ [2026-09-06 22:15:27] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:15:27] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:27] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:28] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:28] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:28] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:28] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:28] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:28] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:28] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:28] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:28] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:15:29] Number of available HVF: 2000
#> ℹ [2026-09-06 22:15:29] Finished check
#> ℹ [2026-09-06 22:15:29] Perform fastMNN integration
#> ℹ [2026-09-06 22:15:31] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:15:31] Reorder clusters...
#> ℹ [2026-09-06 22:15:31] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:15:31] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-09-06 22:15:37] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-09-06 22:15:42] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ✔ [2026-09-06 22:15:49] fastMNN integration completed
#> ◌ [2026-09-06 22:15:49] Run integration workflow...
#> ℹ [2026-09-06 22:15:51] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:15:52] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:15:52] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:52] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:52] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:52] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:52] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:52] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:52] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:52] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:53] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:15:53] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:15:53] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:15:53] Number of available HVF: 2000
#> ℹ [2026-09-06 22:15:53] Finished check
#> ℹ [2026-09-06 22:15:53] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 22:15:53] Perform linear dimension reduction("pca")
#> ℹ [2026-09-06 22:15:54] Perform Harmony integration
#> ℹ [2026-09-06 22:15:54] Using "Harmonypca" (1:50) as input
#> ℹ [2026-09-06 22:15:54] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:15:54] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:15:54] Reorder clusters...
#> ℹ [2026-09-06 22:15:55] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:15:55] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-09-06 22:16:00] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-09-06 22:16:05] Perform umap nonlinear dimension reduction using Harmonypca (1:50)
#> ✔ [2026-09-06 22:16:12] Harmony integration completed
#> ◌ [2026-09-06 22:16:12] Run integration workflow...
#> ℹ [2026-09-06 22:16:14] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-09-06 22:16:14] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-09-06 22:16:14] Number of available HVF: 2000
#> ℹ [2026-09-06 22:16:15] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-09-06 22:16:15] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-09-06 22:16:16] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ℹ [2026-09-06 22:16:17] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:16:17] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:16:17] Reorder clusters...
#> ℹ [2026-09-06 22:16:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:16:17] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-09-06 22:16:23] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-09-06 22:16:28] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-09-06 22:16:35] Harmony5 integration completed
#> ◌ [2026-09-06 22:16:35] Run integration workflow...
#> ℹ [2026-09-06 22:16:37] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:16:38] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:16:38] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:16:38] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:16:38] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:16:38] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:16:38] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:16:38] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:16:38] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:16:38] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:16:39] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:16:39] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:16:39] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:16:39] Number of available HVF: 2000
#> ℹ [2026-09-06 22:16:39] Finished check
#> ℹ [2026-09-06 22:16:39] Perform Coralysis integration
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
#> ℹ [2026-09-06 22:22:43] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:22:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:22:43] Reorder clusters...
#> ℹ [2026-09-06 22:22:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:22:43] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-09-06 22:22:49] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-09-06 22:22:54] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ✔ [2026-09-06 22:23:01] Coralysis integration completed
#> ◌ [2026-09-06 22:23:01] Run integration workflow...
#> ℹ [2026-09-06 22:23:03] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:23:04] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:23:04] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:04] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:04] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:04] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:04] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:04] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:05] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:05] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:05] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:05] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:05] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:23:05] Number of available HVF: 2000
#> ℹ [2026-09-06 22:23:05] Finished check
#> ℹ [2026-09-06 22:23:05] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-09-06 22:23:05] Perform LIGER integration
#> ℹ [2026-09-06 22:23:11] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:23:11] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:23:11] Reorder clusters...
#> ℹ [2026-09-06 22:23:11] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:23:11] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-09-06 22:23:17] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-09-06 22:23:22] Perform umap nonlinear dimension reduction using LIGER (1:50)
#> ✔ [2026-09-06 22:23:29] LIGER integration completed
#> ◌ [2026-09-06 22:23:30] Run integration workflow...
#> ℹ [2026-09-06 22:23:38] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:23:39] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:23:39] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:39] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:40] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:40] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:40] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:40] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:40] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:40] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:40] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:40] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:40] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:23:40] Number of available HVF: 2000
#> ℹ [2026-09-06 22:23:40] Finished check
#> ℹ [2026-09-06 22:23:40] Perform ScaleData on the data 1 ...
#> ℹ [2026-09-06 22:23:41] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:23:41] Perform ScaleData on the data 2 ...
#> ℹ [2026-09-06 22:23:41] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:23:41] Perform ScaleData on the data 3 ...
#> ℹ [2026-09-06 22:23:41] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:23:41] Perform ScaleData on the data 4 ...
#> ℹ [2026-09-06 22:23:41] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:23:41] Perform ScaleData on the data 5 ...
#> ℹ [2026-09-06 22:23:41] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:23:41]  Perform Conos integration
#> ℹ [2026-09-06 22:23:41] Conos integration using pca (1:50) as input
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
#> ◌ [2026-09-06 22:23:42] Run integration workflow...
#> ℹ [2026-09-06 22:23:44] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:23:45] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:23:45] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:45] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:46] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:46] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:46] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:46] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:46] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:46] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:46] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:23:46] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 22:23:46] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:23:47] Number of available HVF: 2000
#> ℹ [2026-09-06 22:23:47] Finished check
#> ℹ [2026-09-06 22:23:47] Perform FindIntegrationAnchors
#> ℹ [2026-09-06 22:24:26] Perform Seurat integration
#> ℹ [2026-09-06 22:25:14] Perform ScaleData on `srt_integrated`
#> ℹ [2026-09-06 22:25:14] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 22:25:15] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:25:15] Reorder clusters...
#> ℹ [2026-09-06 22:25:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:25:15] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:25:21] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:25:26] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:25:29] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-09-06 22:25:34] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> Registered S3 methods overwritten by 'RcppEigen':
#>   method               from         
#>   predict.fastLm       RcppArmadillo
#>   print.fastLm         RcppArmadillo
#>   summary.fastLm       RcppArmadillo
#>   print.summary.fastLm RcppArmadillo
#> ! [2026-09-06 22:26:57] (converted from warning) 'as(<dsCMatrix>, "dgTMatrix")' is deprecated.
#> !                       Use 'as(as(., "generalMatrix"), "TsparseMatrix")' instead.
#> !                       See help("Deprecated") and help("Matrix-deprecated").
#> ! [2026-09-06 22:26:57] Error when performing nonlinear dimension reduction. Skip this step
#> ℹ [2026-09-06 22:26:57] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-09-06 22:27:11] Seurat integration completed
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
