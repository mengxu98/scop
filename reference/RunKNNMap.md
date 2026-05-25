# Single-cell reference mapping with KNN method

This function performs single-cell reference mapping using the K-nearest
neighbor (KNN) method. It takes two single-cell datasets as input:
srt_query and srt_ref. The function maps cells from the srt_query
dataset to the srt_ref dataset based on their similarity or distance.

## Usage

``` r
RunKNNMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = NULL,
  ref_umap = NULL,
  ref_group = NULL,
  features = NULL,
  nfeatures = 2000,
  query_reduction = NULL,
  ref_reduction = NULL,
  query_dims = 1:30,
  ref_dims = 1:30,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean"
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  A Seurat object or count matrix representing the reference object. If
  provided, the similarities will be calculated between cells from the
  query and reference objects. If not provided, the similarities will be
  calculated within the query object.

- query_assay:

  The assay to use for the query object. If not provided, the default
  assay of the query object will be used.

- ref_assay:

  The assay to use for the reference object. If not provided, the
  default assay of the reference object will be used.

- ref_umap:

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  The grouping variable in the reference object. This variable will be
  used to group cells in the heatmap columns. If not provided, all cells
  will be treated as one group.

- features:

  A vector of feature names to include in the heatmap. If not provided,
  highly variable features (HVF) will be used.

- nfeatures:

  The maximum number of features to include in the heatmap. Default is
  `2000`.

- query_reduction:

  A character string specifying the name of a dimensionality reduction
  in the query object to use for calculating the distance metric.

- ref_reduction:

  A character string specifying the name of a dimensionality reduction
  in the reference object to use for calculating the distance metric.

- query_dims:

  A numeric vector specifying the dimension indices from the query
  reduction to be used for calculating the distance metric.

- ref_dims:

  A numeric vector specifying the dimension indices from the reference
  reduction to be used for calculating the distance metric.

- projection_method:

  A character string specifying the projection method to use. Options
  are "model" and "knn". If "model" is selected, the function will try
  to use a pre-trained UMAP model in the reference object for
  projection. If "knn" is selected, the function will directly find the
  nearest neighbors using the distance metric.

- nn_method:

  A character string specifying the nearest neighbor search method to
  use. Options are "raw", "annoy", "rann", and "cpp". If "raw" is
  selected, the function will use the brute-force method to find the
  nearest neighbors. If "annoy" is selected, the function will use the
  Annoy library for approximate nearest neighbor search. If "rann" is
  selected, the function will use the RANN library for approximate
  nearest neighbor search. If "cpp" is selected, the function will use
  the compiled exact top-k search. If not provided, the function will
  use "cpp" for Euclidean or cosine distance, otherwise it will choose
  the search method based on the size of the query and reference
  datasets.

- k:

  A number of nearest neighbors to find for each cell in the query
  object.

- distance_metric:

  The distance metric to use for calculating the pairwise distances
  between cells. Options include: "pearson", "spearman", "cosine",
  "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman",
  "simple matching", and "faith". Additional distance metrics can also
  be used, such as "euclidean", "manhattan", "hamming", etc.

- vote_fun:

  A character string specifying the function to be used for aggregating
  the nearest neighbors in the reference object. Options are "mean",
  "median", "sum", "min", "max", "sd", "var", etc. If not provided, the
  default is "mean".

## Value

A Seurat object with the projection results stored in the
"ref.embeddings" reduction. If `ref_group` is provided, the function
will also add a new metadata column called "predicted_ref_group" to the
query object.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunSingleR](https://mengxu98.github.io/scop/reference/RunSingleR.md),
[CellCorHeatmap](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-05-25 05:28:52] Start standard processing workflow...
#> ℹ [2026-05-25 05:28:52] Checking a list of <Seurat>...
#> ! [2026-05-25 05:28:52] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:28:52] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:28:54] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:28:55] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:28:55] Number of available HVF: 2000
#> ℹ [2026-05-25 05:28:55] Finished check
#> ℹ [2026-05-25 05:28:55] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:28:55] Perform pca linear dimension reduction
#> ℹ [2026-05-25 05:28:56] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-05-25 05:28:57] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 05:28:57] Reorder clusters...
#> ℹ [2026-05-25 05:28:57] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:28:57] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 05:28:57] Perform umap nonlinear dimension reduction using Standardpca (1:27)
#> ℹ [2026-05-25 05:29:02] Perform umap nonlinear dimension reduction using Standardpca (1:27)
#> ✔ [2026-05-25 05:29:08] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-05-25 05:29:09] Run integration workflow...
#> ℹ [2026-05-25 05:29:09] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-25 05:29:10] Checking a list of <Seurat>...
#> ℹ [2026-05-25 05:29:10] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-05-25 05:29:10] Perform `Seurat::FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-05-25 05:29:10] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-05-25 05:29:10] Perform `Seurat::FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-05-25 05:29:11] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-05-25 05:29:11] Perform `Seurat::FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-05-25 05:29:12] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-05-25 05:29:12] Perform `Seurat::FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-05-25 05:29:12] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:29:12] Number of available HVF: 2000
#> ℹ [2026-05-25 05:29:13] Finished check
#> ℹ [2026-05-25 05:29:14] Perform Uncorrected integration
#> ℹ [2026-05-25 05:29:15] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:29:16] Perform "pca" linear dimension reduction
#> ℹ [2026-05-25 05:29:16] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-25 05:29:17] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-25 05:29:17] Reorder clusters...
#> ℹ [2026-05-25 05:29:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:29:17] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:19)
#> ℹ [2026-05-25 05:29:23] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:19)
#> ℹ [2026-05-25 05:29:29] Perform umap nonlinear dimension reduction using Standardpca (1:27)
#> Warning: Key ‘StandardpcaUMAP2D_’ taken, using ‘standardpcaumap2d_’ instead
#> ✔ [2026-05-25 05:29:36] Uncorrected integration completed
CellDimPlot(
  srt_ref,
  group.by = c("celltype", "tech")
)


# Set the number of threads for RcppParallel
# details see: ?RcppParallel::setThreadOptions
# if (requireNamespace("RcppParallel", quietly = TRUE)) {
#   RcppParallel::setThreadOptions()
# }
# Projection
srt_query <- RunKNNMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_umap = "UncorrectedUMAP2D"
)
#> ℹ [2026-05-25 05:29:36] Use the features to calculate distance metric
#> ℹ [2026-05-25 05:29:36] Data type is log-normalized
#> ℹ [2026-05-25 05:29:37] Data type is log-normalized
#> ℹ [2026-05-25 05:29:37] Use 636 features to calculate distance
#> ℹ [2026-05-25 05:29:37] Use cpp method to find neighbors
#> ℹ [2026-05-25 05:29:37] Running UMAP projection
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
```
