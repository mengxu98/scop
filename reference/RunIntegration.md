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
#> ◌ [2026-10-25 21:13:47] Run integration workflow...
#> ℹ [2026-10-25 21:13:48] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:13:48] Checking a list of <Seurat>...
#> ! [2026-10-25 21:13:48] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:13:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:48] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-10-25 21:13:49] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:13:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:49] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-10-25 21:13:49] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:13:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:49] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-10-25 21:13:49] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:13:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:49] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-10-25 21:13:49] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:13:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:49] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:13:49] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:13:49] Number of available HVF: 500
#> ℹ [2026-10-25 21:13:49] Finished check
#> ℹ [2026-10-25 21:13:51] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:13:52] Perform linear dimension reduction("pca")
#> ℹ [2026-10-25 21:13:52] Perform Harmony integration
#> ℹ [2026-10-25 21:13:52] Using "Harmonypca" (1:10) as input
#> ℹ [2026-10-25 21:13:52] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:13:53] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:13:53] Reorder clusters...
#> ℹ [2026-10-25 21:13:53] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:13:53] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ℹ [2026-10-25 21:13:59] Perform umap nonlinear dimension reduction using Harmonypca (1:10)
#> ℹ [2026-10-25 21:14:05] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-25 21:14:05] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-25 21:14:05] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-25 21:14:06] Harmony integration completed
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
#> ◌ [2026-10-25 21:14:06] Run integration workflow...
#> ℹ [2026-10-25 21:14:07] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:14:08] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:14:08] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:08] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:08] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:08] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:08] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:08] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:09] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:09] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:09] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:09] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:09] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:14:09] Number of available HVF: 2000
#> ℹ [2026-10-25 21:14:09] Finished check
#> ℹ [2026-10-25 21:14:09] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-25 21:14:09] Perform LIGER integration
#> ℹ [2026-10-25 21:14:14] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:14:15] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:14:15] Reorder clusters...
#> ℹ [2026-10-25 21:14:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:14:15] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-25 21:14:21] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-25 21:14:27] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ✔ [2026-10-25 21:14:34] LIGER integration completed
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  compute_lisi = TRUE,
  lisi_label_colnames = "tech"
)
#> ◌ [2026-10-25 21:14:34] Run integration workflow...
#> ℹ [2026-10-25 21:14:35] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:14:36] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:14:36] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:36] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:36] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:36] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:36] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:36] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:37] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:37] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:37] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:14:37] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:14:37] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:14:37] Number of available HVF: 2000
#> ℹ [2026-10-25 21:14:37] Finished check
#> ℹ [2026-10-25 21:14:37] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:14:37] Perform linear dimension reduction("pca")
#> ℹ [2026-10-25 21:14:38] Perform Harmony integration
#> ℹ [2026-10-25 21:14:38] Using "Harmonypca" (1:20) as input
#> ℹ [2026-10-25 21:14:38] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:14:39] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:14:39] Reorder clusters...
#> ℹ [2026-10-25 21:14:39] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:14:39] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-25 21:14:45] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-25 21:14:51] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ℹ [2026-10-25 21:14:57] Compute LISI scores from reduction "HarmonypcaUMAP2D"
#> ℹ [2026-10-25 21:14:57] Compute LISI scores from reduction "HarmonyUMAP2D"
#> ✔ [2026-10-25 21:14:57] Stored LISI scores in metadata: "HarmonypcaUMAP2D_tech_LISI" and "HarmonyUMAP2D_tech_LISI"
#> ✔ [2026-10-25 21:14:58] Harmony integration completed
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
#> ◌ [2026-10-25 21:14:58] Run integration workflow...
#> ! [2026-10-25 21:14:58] `integration_method = 'Harmony5'` is not compatible with <ChromatinAssay> in current Seurat v5 workflow. Automatically switch to "Harmony"
#> ℹ [2026-10-25 21:14:58] Split `srt_merge` into `srt_list` by "batch"
#> ℹ [2026-10-25 21:14:58] Checking a list of <Seurat>...
#> ! [2026-10-25 21:14:58] Data 1/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-25 21:14:58] Perform `RunTFIDF()` on 1/2 of `srt_list`...
#> ℹ [2026-10-25 21:14:58] Perform `FindTopFeatures()` on 1/2 of `srt_list`...
#> ! [2026-10-25 21:14:58] Data 2/2 of the `srt_list` is "raw_counts"
#> ℹ [2026-10-25 21:14:58] Perform `RunTFIDF()` on 2/2 of `srt_list`...
#> ℹ [2026-10-25 21:14:59] Perform `FindTopFeatures()` on 2/2 of `srt_list`...
#> ℹ [2026-10-25 21:14:59] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:14:59] Number of available HVF: 2000
#> ℹ [2026-10-25 21:14:59] Finished check
#> ℹ [2026-10-25 21:14:59] `normalization_method` is "TFIDF". Use lsi workflow...
#> ℹ [2026-10-25 21:14:59] Perform linear dimension reduction("svd")
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-10-25 21:14:59] Perform Harmony integration
#> ℹ [2026-10-25 21:14:59] Using "Harmonysvd" (2:5) as input
#> ℹ [2026-10-25 21:15:00] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:15:00] Reorder clusters...
#> ℹ [2026-10-25 21:15:00] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:15:00] Perform umap nonlinear dimension reduction using Harmony (1:4)
#> ℹ [2026-10-25 21:15:05] Perform umap nonlinear dimension reduction using Harmony (1:4)
#> ✔ [2026-10-25 21:15:10] Harmony integration completed

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
#> ◌ [2026-10-25 21:15:10] Run integration workflow...
#> ℹ [2026-10-25 21:15:11] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:15:12] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:15:12] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:12] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:12] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:12] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:12] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:12] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:13] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:13] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:13] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:13] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:13] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:15:13] Number of available HVF: 2000
#> ℹ [2026-10-25 21:15:13] Finished check
#> ℹ [2026-10-25 21:15:13] Perform Uncorrected integration
#> ℹ [2026-10-25 21:15:13] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:15:13] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:15:14] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:15:14] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:15:14] Reorder clusters...
#> ℹ [2026-10-25 21:15:14] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:15:14] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-25 21:15:20] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-10-25 21:15:26] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-10-25 21:15:33] Uncorrected integration completed
#> ◌ [2026-10-25 21:15:33] Run integration workflow...
#> ℹ [2026-10-25 21:15:34] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:15:35] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:15:35] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:35] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:35] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:35] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:35] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:35] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:35] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:35] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:36] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:15:36] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:15:36] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:15:36] Number of available HVF: 2000
#> ℹ [2026-10-25 21:15:36] Finished check
#> ℹ [2026-10-25 21:15:36] Perform FindIntegrationAnchors
#> ℹ [2026-10-25 21:16:00] Perform Seurat integration
#> ℹ [2026-10-25 21:16:14] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-25 21:16:14] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:16:15] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:16:15] Reorder clusters...
#> ℹ [2026-10-25 21:16:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:16:15] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:16:22] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:16:28] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-25 21:16:37] Seurat integration completed
#> ◌ [2026-10-25 21:16:37] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-25 21:16:38] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-25 21:16:38] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-25 21:16:38] Number of available HVF: 2000
#> ℹ [2026-10-25 21:16:39] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-25 21:16:39] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-25 21:16:40] Perform Seurat v5 integration with `CCAIntegration()`
#> ℹ [2026-10-25 21:16:47] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:16:47] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:16:47] Reorder clusters...
#> ℹ [2026-10-25 21:16:47] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:16:47] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-25 21:16:53] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-10-25 21:16:59] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-25 21:17:07] CCA integration completed
#> ◌ [2026-10-25 21:17:07] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-25 21:17:08] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-25 21:17:08] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-25 21:17:08] Number of available HVF: 2000
#> ℹ [2026-10-25 21:17:10] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-25 21:17:10] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-25 21:17:11] Perform Seurat v5 integration with `RPCAIntegration()`
#> ℹ [2026-10-25 21:17:22] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:17:22] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:17:22] Reorder clusters...
#> ℹ [2026-10-25 21:17:22] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:17:22] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-25 21:17:29] Perform umap nonlinear dimension reduction using RPCA (1:50)
#> ℹ [2026-10-25 21:17:35] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-25 21:17:43] RPCA integration completed
#> ◌ [2026-10-25 21:17:43] Run integration workflow...
#> ℹ [2026-10-25 21:17:44] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:17:44] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:17:44] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:17:45] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:17:45] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:17:45] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:17:45] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:17:45] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:17:45] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:17:45] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:17:45] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:17:45] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:17:45] Number of available HVF: 2000
#> ℹ [2026-10-25 21:17:46] Finished check
#> ℹ [2026-10-25 21:17:46] Perform MNN integration
#> ℹ [2026-10-25 21:18:31] Perform ScaleData
#> ℹ [2026-10-25 21:18:31] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:18:32] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:18:32] Reorder clusters...
#> ℹ [2026-10-25 21:18:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:18:32] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-25 21:18:38] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ℹ [2026-10-25 21:18:45] Perform umap nonlinear dimension reduction using MNNpca (1:50)
#> ✔ [2026-10-25 21:18:53] MNN integration completed
#> ◌ [2026-10-25 21:18:53] Run integration workflow...
#> ℹ [2026-10-25 21:18:55] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:18:55] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:18:55] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:18:55] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:18:55] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:18:56] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:18:56] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:18:56] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:18:56] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:18:56] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:18:56] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:18:56] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:18:56] Number of available HVF: 2000
#> ℹ [2026-10-25 21:18:56] Finished check
#> ℹ [2026-10-25 21:18:57] Perform fastMNN integration
#> ℹ [2026-10-25 21:19:00] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:19:01] Reorder clusters...
#> ℹ [2026-10-25 21:19:01] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:19:01] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-25 21:19:07] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ℹ [2026-10-25 21:19:13] Perform umap nonlinear dimension reduction using fastMNN (1:50)
#> ✔ [2026-10-25 21:19:21] fastMNN integration completed
#> ◌ [2026-10-25 21:19:21] Run integration workflow...
#> ℹ [2026-10-25 21:19:23] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:19:24] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:19:24] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:19:24] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:19:25] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:19:25] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:19:25] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:19:25] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:19:25] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:19:25] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:19:25] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:19:25] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:19:25] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:19:25] Number of available HVF: 2000
#> ℹ [2026-10-25 21:19:26] Finished check
#> ℹ [2026-10-25 21:19:26] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:19:26] Perform linear dimension reduction("pca")
#> ℹ [2026-10-25 21:19:26] Perform Harmony integration
#> ℹ [2026-10-25 21:19:26] Using "Harmonypca" (1:50) as input
#> ℹ [2026-10-25 21:19:27] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:19:27] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:19:27] Reorder clusters...
#> ℹ [2026-10-25 21:19:27] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:19:27] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-25 21:19:34] Perform umap nonlinear dimension reduction using Harmony (1:50)
#> ℹ [2026-10-25 21:19:40] Perform umap nonlinear dimension reduction using Harmonypca (1:50)
#> ✔ [2026-10-25 21:19:48] Harmony integration completed
#> ◌ [2026-10-25 21:19:48] Run integration workflow...
#> Splitting ‘counts’, ‘data’ layers. Not splitting . If you would like to split other layers, set in `layers` argument.
#> ℹ [2026-10-25 21:19:49] Data appears already log-normalized. Skipping `Seurat::NormalizeData()`
#> ℹ [2026-10-25 21:19:50] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-10-25 21:19:50] Number of available HVF: 2000
#> ℹ [2026-10-25 21:19:51] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-10-25 21:19:51] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-10-25 21:19:52] Perform Seurat v5 integration with `HarmonyIntegration()`
#> The `features` argument is ignored by `HarmonyIntegration`.
#> This message is displayed once per session.
#> ℹ [2026-10-25 21:19:53] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:19:53] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:19:53] Reorder clusters...
#> ℹ [2026-10-25 21:19:53] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:19:53] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-25 21:20:00] Perform umap nonlinear dimension reduction using Harmony5 (1:50)
#> ℹ [2026-10-25 21:20:06] Perform umap nonlinear dimension reduction using pca (1:50)
#> ✔ [2026-10-25 21:20:14] Harmony5 integration completed
#> ◌ [2026-10-25 21:20:14] Run integration workflow...
#> ℹ [2026-10-25 21:20:15] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:20:16] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:20:17] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:20:17] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:20:17] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:20:17] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:20:17] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:20:17] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:20:17] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:20:17] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:20:18] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:20:18] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:20:18] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:20:18] Number of available HVF: 2000
#> ℹ [2026-10-25 21:20:18] Finished check
#> ℹ [2026-10-25 21:20:18] Perform Coralysis integration
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
#> ℹ [2026-10-25 21:21:43] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:21:43] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:21:43] Reorder clusters...
#> ℹ [2026-10-25 21:21:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:21:43] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-25 21:21:49] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ℹ [2026-10-25 21:21:55] Perform umap nonlinear dimension reduction using Coralysis (1:50)
#> ✔ [2026-10-25 21:22:03] Coralysis integration completed
#> ◌ [2026-10-25 21:22:03] Run integration workflow...
#> ℹ [2026-10-25 21:22:05] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:22:06] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:22:06] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:06] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:07] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:07] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:07] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:07] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:07] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:07] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:07] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:07] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:07] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:22:07] Number of available HVF: 2000
#> ℹ [2026-10-25 21:22:08] Finished check
#> ℹ [2026-10-25 21:22:08] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-10-25 21:22:08] Perform LIGER integration
#> ℹ [2026-10-25 21:22:13] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:22:13] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:22:13] Reorder clusters...
#> ℹ [2026-10-25 21:22:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:22:14] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-25 21:22:20] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-10-25 21:22:26] Perform umap nonlinear dimension reduction using LIGER (1:50)
#> ✔ [2026-10-25 21:22:34] LIGER integration completed
#> ◌ [2026-10-25 21:22:35] Run integration workflow...
#> ℹ [2026-10-25 21:22:36] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:22:37] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:22:38] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:38] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:38] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:38] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:38] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:38] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:38] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:38] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:39] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:22:39] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:22:39] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:22:39] Number of available HVF: 2000
#> ℹ [2026-10-25 21:22:39] Finished check
#> ℹ [2026-10-25 21:22:39] Perform ScaleData on the data 1 ...
#> ℹ [2026-10-25 21:22:39] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:22:39] Perform ScaleData on the data 2 ...
#> ℹ [2026-10-25 21:22:39] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:22:39] Perform ScaleData on the data 3 ...
#> ℹ [2026-10-25 21:22:39] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:22:40] Perform ScaleData on the data 4 ...
#> ℹ [2026-10-25 21:22:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:22:40] Perform ScaleData on the data 5 ...
#> ℹ [2026-10-25 21:22:40] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:22:40]  Perform Conos integration
#> ℹ [2026-10-25 21:22:40] Conos integration using pca (1:50) as input
#> ℹ [2026-10-25 21:22:44] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:22:44] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:22:44] Reorder clusters...
#> ℹ [2026-10-25 21:22:44] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:22:44] Perform umap nonlinear dimension reduction using Conos
#> ℹ [2026-10-25 21:22:51] Perform umap nonlinear dimension reduction using Conos
#> ℹ [2026-10-25 21:22:57] Perform umap nonlinear dimension reduction using ConosUMAP2D (1:50)
#> ✔ [2026-10-25 21:23:05] Conos integration completed
#> ◌ [2026-10-25 21:23:05] Run integration workflow...
#> ℹ [2026-10-25 21:23:07] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:23:08] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:23:08] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:08] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:08] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:08] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:09] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:09] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:09] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:09] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:09] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:09] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:09] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:23:09] Number of available HVF: 2000
#> ℹ [2026-10-25 21:23:09] Finished check
#> ℹ [2026-10-25 21:23:09] Perform ComBat integration
#> Found 138 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
#> ℹ [2026-10-25 21:23:10] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:23:10] Perform linear dimension reduction ("pca")
#> ℹ [2026-10-25 21:23:11] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:23:11] Reorder clusters...
#> ℹ [2026-10-25 21:23:11] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:23:11] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ℹ [2026-10-25 21:23:18] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ℹ [2026-10-25 21:23:24] Perform umap nonlinear dimension reduction using ComBatpca (1:50)
#> ✔ [2026-10-25 21:23:32] ComBat integration completed

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
#> ◌ [2026-10-25 21:23:33] Run integration workflow...
#> ℹ [2026-10-25 21:23:35] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:23:36] Checking a list of <Seurat>...
#> ℹ [2026-10-25 21:23:36] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:36] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:36] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:36] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:36] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:36] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:37] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:37] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:37] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-10-25 21:23:37] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:23:37] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:23:37] Number of available HVF: 2000
#> ℹ [2026-10-25 21:23:37] Finished check
#> ℹ [2026-10-25 21:23:37] Perform FindIntegrationAnchors
#> ℹ [2026-10-25 21:24:12] Perform Seurat integration
#> ℹ [2026-10-25 21:25:01] Perform ScaleData on `srt_integrated`
#> ℹ [2026-10-25 21:25:01] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:25:04] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:25:04] Reorder clusters...
#> ℹ [2026-10-25 21:25:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:25:04] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:25:11] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:25:17] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:25:20] Perform tsne nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:25:26] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:00] Perform dm nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:01] Perform phate nonlinear dimension reduction using Seuratpca (1:50)
#> ◌ [2026-10-25 21:27:04] Preparing python environment...
#> ℹ [2026-10-25 21:27:04] Environment name: 'scprep_env' and python version: 3.10-1
#> ℹ [2026-10-25 21:27:04] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-10-25 21:27:04] Using micromamba executable: '/home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba'
#> ℹ [2026-10-25 21:27:04] Creating Python environment with python "3.10-1" using micromamba...
#> ℹ [2026-10-25 21:27:04] Skipping conda Terms of Service acceptance for micromamba.
#> + /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba 'create' '--yes' '--name' 'scprep_env' 'python=3.10-1' 'pip' 'setuptools<81' 'wheel' '--quiet' '-c' 'conda-forge'
#> ✔ [2026-10-25 21:27:12] Environment created successfully: '/home/runner/work/_temp/pkgdown-python/envs/scprep_env'
#> ℹ [2026-10-25 21:27:12] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-10-25 21:27:12] Checking 3 conda requirements
#> ℹ [2026-10-25 21:27:12] Installing 3 packages into environment: 'scprep_env' using micromamba
#> ℹ [2026-10-25 21:27:12] Installing packages via micromamba...
#> ℹ [2026-10-25 21:27:12] Using channel: "conda-forge"
#> ℹ [2026-10-25 21:27:12] Installing 3 packages...
#> + /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba 'install' '--yes' '--name' 'scprep_env' '-c' 'conda-forge' 'leidenalg=0.10.2' 'tbb=2022.2.0' 'python-igraph=0.11.9'
#> Using Cached Shard Index for conda-forge/linux-64                                                   ✔ Done
#> Using Cached Shard Index for conda-forge/noarch                                                     ✔ Done
#> Fetching and Parsing Packages' Shards                                                           ⧖ Starting
#> Fetching and Parsing Packages' Shards                                                     ✔ Done (1.0 sec)
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
#>   Package           Version  Build                  Channel          Size
#> ───────────────────────────────────────────────────────────────────────────
#>   Install:
#> ───────────────────────────────────────────────────────────────────────────
#> 
#>   + arpack            3.9.1  nompi_hf03ea27_102     conda-forge     130kB
#>   + glpk                5.0  h445213a_0             conda-forge       1MB
#>   + gmp               6.3.0  hfd2156b_3             conda-forge     493kB
#>   + igraph          0.10.16  hfe3e89f_1             conda-forge       1MB
#>   + leidenalg        0.10.2  py310hea6c23e_1        conda-forge      81kB
#>   + libblas          3.11.0  11_h4a7cf45_openblas   conda-forge      18kB
#>   + libgfortran      16.2.0  h69a702a_7             conda-forge      29kB
#>   + libgfortran5     16.2.0  h6b99dfc_7             conda-forge       3MB
#>   + libhwloc         2.12.1  default_hafda6a7_1003  conda-forge       2MB
#>   + libiconv           1.18  h0cb94f2_3             conda-forge     789kB
#>   + liblapack        3.11.0  11_h47877c9_openblas   conda-forge      18kB
#>   + libleidenalg     0.11.1  h00ab1b0_0             conda-forge      83kB
#>   + libopenblas      0.3.34  pthreads_hf13c14d_2    conda-forge       7MB
#>   + libstdcxx        16.2.0  h934c35e_7             conda-forge       7MB
#>   + libstdcxx-ng     16.2.0  hdf11a46_7             conda-forge      29kB
#>   + libxml2          2.15.4  ha7a5771_1             conda-forge      47kB
#>   + libxml2-16       2.15.4  hbdfff7e_1             conda-forge     566kB
#>   + libxml2-devel    2.15.4  ha7a5771_1             conda-forge      80kB
#>   + python-igraph    0.11.9  py310hea6c23e_3        conda-forge     573kB
#>   + python_abi         3.10  9_cp310                conda-forge       7kB
#>   + tbb            2022.2.0  hb60516a_1             conda-forge     183kB
#>   + texttable         1.7.0  pyhd8ed1ab_1           conda-forge      20kB
#> 
#>   Summary:
#> 
#>   Install: 22 packages
#> 
#>   Total download: 24MB
#> 
#> ─────────────────────────────────────────────────────────────────────────warning  libmamba Security Warning: This transaction includes executing package scripts (pre/post-link/unlink) if present. These scripts can contain arbitrary code. Please ensure you trust the package sources.
#> ──
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
#> ✔ [2026-10-25 21:27:16] micromamba installation completed successfully
#> ℹ [2026-10-25 21:27:17] Checking 12 requirements using uv
#> ℹ [2026-10-25 21:27:17] Installing 12 packages into environment: 'scprep_env' using micromamba
#> ℹ [2026-10-25 21:27:17] Installing packages via uv...
#> + /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/uv 'pip' 'install' '--python' '/home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10' '--force-reinstall' 'setuptools<81' 'matplotlib==3.10.8' 'numba==0.66.0' 'llvmlite==0.48.0' 'numpy==1.26.4' 'packaging>=24.0' 'pandas==2.0.3' 'scikit-learn==1.7.0' 'scipy==1.15.3' 'phate==1.0.11' 'leidenalg==0.10.2' 'igraph==0.11.9'
#> Using Python 3.10.1 environment at: /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Resolved 35 packages in 413ms
#> Downloading scipy (35.9MiB)
#> Downloading fonttools (5.0MiB)
#> Downloading setuptools (1.0MiB)
#> Downloading igraph (4.2MiB)
#> Downloading numpy (17.4MiB)
#> Downloading leidenalg (1.9MiB)
#> Downloading kiwisolver (1.6MiB)
#> Downloading pandas (11.7MiB)
#> Downloading pygsp (1.8MiB)
#> Downloading pillow (6.6MiB)
#> Downloading llvmlite (57.1MiB)
#> Downloading scikit-learn (12.3MiB)
#> Downloading numba (3.6MiB)
#> Downloading matplotlib (8.3MiB)
#>  Downloaded leidenalg
#>  Downloaded kiwisolver
#>  Downloaded pygsp
#>  Downloaded setuptools
#>  Downloaded igraph
#>  Downloaded fonttools
#>  Downloaded pillow
#>  Downloaded numba
#>  Downloaded matplotlib
#>  Downloaded scikit-learn
#>  Downloaded numpy
#>  Downloaded pandas
#>  Downloaded scipy
#>  Downloaded llvmlite
#> Prepared 35 packages in 1.80s
#> Uninstalled 5 packages in 38ms
#> Installed 35 packages in 70ms
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
#>  + pytz==2026.5
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
#>  + tzdata==2026.5
#>  + wrapt==2.5.1
#> ✔ [2026-10-25 21:27:20] uv installation completed
#> + /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10 -m pip check
#> No broken requirements found.
#> ✔ [2026-10-25 21:27:22] Python environment ready
#> ✔                       Until next loading, the environment will be cached
#> Environment manager config:
#>  Manager:       micromamba
#>  Executable:    /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#>  Environment:   /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Python config:
#>  python:        /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10
#>  version:       Python 3.10.1
#> 
#> ℹ [2026-10-25 21:27:34] Perform phate nonlinear dimension reduction using Seuratpca (1:50)
#> ◌ [2026-10-25 21:27:37] Preparing python environment...
#> ℹ [2026-10-25 21:27:37] Environment name: 'scprep_env' and python version: 3.10-1
#> ℹ [2026-10-25 21:27:37] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-10-25 21:27:37] Using micromamba executable: '/home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba'
#> ℹ [2026-10-25 21:27:37] Using existing environment: '/home/runner/work/_temp/pkgdown-python/envs/scprep_env'
#> ℹ [2026-10-25 21:27:37] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-10-25 21:27:37] Checking 3 conda requirements
#> ℹ [2026-10-25 21:27:38] Checking 12 requirements using uv
#> + /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10 -m pip check
#> No broken requirements found.
#> ✔ [2026-10-25 21:27:40] Python environment ready
#> ✔                       Until next loading, the environment will be cached
#> Environment manager config:
#>  Manager:       micromamba
#>  Executable:    /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#>  Environment:   /home/runner/work/_temp/pkgdown-python/envs/scprep_env
#> Python config:
#>  python:        /home/runner/work/_temp/pkgdown-python/envs/scprep_env/bin/python3.10
#>  version:       Python 3.10.1
#> 
#> ℹ [2026-10-25 21:27:52] Perform pacmap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:52] Perform pacmap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:53] Perform trimap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:53] Perform trimap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:27:54] Perform largevis nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:28:02] Perform largevis nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-10-25 21:28:10] Perform fr nonlinear dimension reduction using Seurat_SNN
#> ℹ [2026-10-25 21:28:11] Perform fr nonlinear dimension reduction using Seurat_SNN
#> ℹ [2026-10-25 21:28:14] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-10-25 21:28:32] Seurat integration completed
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
