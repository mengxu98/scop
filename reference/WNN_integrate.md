# WNN integration function

WNN integration function

## Usage

``` r
WNN_integrate(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
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
  verbose = TRUE,
  seed = 11
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

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

## Examples

``` r
data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
pbmcmultiome_sub <- WNN_integrate(
  srt_merge = pbmcmultiome_sub,
  batch = "batch",
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10
)
#> ℹ [2026-09-06 22:44:20] Start standard processing workflow...
#> ℹ [2026-09-06 22:44:20] Auto preprocess assays: "RNA" and "peaks"
#> ℹ [2026-09-06 22:44:20] Start standard processing workflow...
#> ℹ [2026-09-06 22:44:20] Checking a list of <Seurat>...
#> ! [2026-09-06 22:44:20] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:44:20] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:44:20] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:44:20] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:44:20] Number of available HVF: 2000
#> ℹ [2026-09-06 22:44:20] Finished check
#> ℹ [2026-09-06 22:44:20] Perform `ScaleData()`
#> ℹ [2026-09-06 22:44:20] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:44:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:44:21] Reorder clusters...
#> ℹ [2026-09-06 22:44:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:44:21] Perform umap nonlinear dimension reduction
#> ℹ [2026-09-06 22:44:28] Perform umap nonlinear dimension reduction using RNApca (1:10)
#> ✔ [2026-09-06 22:44:34] Standard processing workflow completed
#> ℹ [2026-09-06 22:44:34] Start standard processing workflow...
#> ℹ [2026-09-06 22:44:34] Checking a list of <Seurat>...
#> ! [2026-09-06 22:44:34] Data 1/1 of the `srt_list` is "raw_counts"
#> ℹ [2026-09-06 22:44:34] Perform `RunTFIDF()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:44:34] Perform `FindTopFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:44:34] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:44:34] Number of available HVF: 11413
#> ℹ [2026-09-06 22:44:34] Finished check
#> ℹ [2026-09-06 22:44:34] `normalization_method` is TFIDF. Use lsi workflow
#> ℹ [2026-09-06 22:44:34] Perform svd linear dimension reduction
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-09-06 22:44:35] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:44:35] Reorder clusters...
#> ℹ [2026-09-06 22:44:35] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:44:36] Perform umap nonlinear dimension reduction
#> ℹ [2026-09-06 22:44:43] Perform umap nonlinear dimension reduction using ATACsvd (1:10)
#> ✔ [2026-09-06 22:44:48] Standard processing workflow completed
#> ℹ [2026-09-06 22:44:48] Adjust neighbor k from 20 to 20 for small-sample WNN graph construction
#> ℹ [2026-09-06 22:44:48] Adjust WNN knn.range to 80 for small-sample graph construction
#> ℹ [2026-09-06 22:44:48] Perform WNN integration using RNApca and ATAClsi
#> Calculating cell-specific modality weights
#> Finding 20 nearest neighbors for each modality.
#> Calculating kernel bandwidths
#> Finding multimodal neighbors
#> Constructing multimodal KNN graph
#> Constructing multimodal SNN graph
#> ℹ [2026-09-06 22:44:50] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:44:50] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:44:50] Reorder clusters...
#> ℹ [2026-09-06 22:44:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:44:50] Perform umap nonlinear dimension reduction using WNN
#> ℹ [2026-09-06 22:44:55] Perform umap nonlinear dimension reduction using WNN
```
