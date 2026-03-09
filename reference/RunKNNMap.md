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

  A Seurat object storing the query cells.

- srt_ref:

  A Seurat object storing the reference cells.

- query_assay:

  A character string specifying the assay name for the query cells. If
  not provided, the default assay for the query object will be used.

- ref_assay:

  A character string specifying the assay name for the reference cells.
  If not provided, the default assay for the reference object will be
  used.

- ref_umap:

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  A character string specifying a metadata column name in the reference
  object to use for grouping.

- features:

  A vector of feature names to include in the heatmap. If not provided,
  highly variable features (HVF) will be used.

- nfeatures:

  A integer specifying the number of highly variable features to be
  calculated if `features` is not provided.

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
  use. Options are "raw", "annoy", and "rann". If "raw" is selected, the
  function will use the brute-force method to find the nearest
  neighbors. If "annoy" is selected, the function will use the Annoy
  library for approximate nearest neighbor search. If "rann" is
  selected, the function will use the RANN library for approximate
  nearest neighbor search. If not provided, the function will choose the
  search method based on the size of the query and reference datasets.

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
#> ℹ [2026-03-09 08:34:46] Start standard scop workflow...
#> ℹ [2026-03-09 08:34:46] Checking a list of <Seurat>...
#> ! [2026-03-09 08:34:46] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-09 08:34:46] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:34:49] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:34:49] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:34:49] Number of available HVF: 2000
#> ℹ [2026-03-09 08:34:50] Finished check
#> ℹ [2026-03-09 08:34:50] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:34:50] Perform pca linear dimension reduction
#> ℹ [2026-03-09 08:34:52] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-09 08:34:52] Reorder clusters...
#> ℹ [2026-03-09 08:34:52] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-09 08:34:52] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-09 08:34:57] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-09 08:35:02] Run scop standard workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-03-09 08:35:03] Run Uncorrected integration...
#> ℹ [2026-03-09 08:35:03] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-03-09 08:35:04] Checking a list of <Seurat>...
#> ℹ [2026-03-09 08:35:04] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:35:04] Perform `Seurat::FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-03-09 08:35:04] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:35:04] Perform `Seurat::FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-03-09 08:35:05] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:35:05] Perform `Seurat::FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-03-09 08:35:06] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:35:06] Perform `Seurat::FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-03-09 08:35:06] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:35:06] Number of available HVF: 2000
#> ℹ [2026-03-09 08:35:07] Finished check
#> ℹ [2026-03-09 08:35:08] Perform Uncorrected integration
#> ℹ [2026-03-09 08:35:10] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:35:11] Perform linear dimension reduction("pca")
#> ℹ [2026-03-09 08:35:12] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-03-09 08:35:12] Reorder clusters...
#> ℹ [2026-03-09 08:35:13] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:10)
#> ℹ [2026-03-09 08:35:17] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:10)
#> ✔ [2026-03-09 08:35:24] Run Uncorrected integration done
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
#> ℹ [2026-03-09 08:35:24] Use the features to calculate distance metric
#> ℹ [2026-03-09 08:35:25] Data type is log-normalized
#> ℹ [2026-03-09 08:35:25] Data type is log-normalized
#> ℹ [2026-03-09 08:35:25] Use 2000 features to calculate distance
#> ℹ [2026-03-09 08:35:26] Use raw method to find neighbors
#> ℹ [2026-03-09 08:35:26] Running UMAP projection
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
