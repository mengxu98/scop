# Single-cell reference mapping with PCA method

Single-cell reference mapping with PCA method

## Usage

``` r
RunPCAMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = srt_ref[[ref_pca]]@assay.used,
  ref_pca = NULL,
  ref_dims = 1:30,
  ref_umap = NULL,
  ref_group = NULL,
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

- ref_pca:

  A character string specifying the name of a PCA reduction in the
  reference object to use for calculating the distance metric. If NULL
  (default), it will be automatically detected as the first PCA
  reduction.

- ref_dims:

  A numeric vector specifying the dimension indices from the reference
  reduction to be used for calculating the distance metric.

- ref_umap:

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  A character string specifying a metadata column name in the reference
  object to use for grouping.

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

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-01-30 17:23:21] Start standard scop workflow...
#> ℹ [2026-01-30 17:23:22] Checking a list of <Seurat>...
#> ! [2026-01-30 17:23:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:23:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:23:24] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:23:25] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:23:25] Number of available HVF: 2000
#> ℹ [2026-01-30 17:23:25] Finished check
#> ℹ [2026-01-30 17:23:26] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:23:26] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:23:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:23:29] Reorder clusters...
#> ℹ [2026-01-30 17:23:29] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:23:29] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:23:35] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:23:40] Run scop standard workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-01-30 17:23:41] Run Uncorrected integration...
#> ℹ [2026-01-30 17:23:41] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-30 17:23:42] Checking a list of <Seurat>...
#> ℹ [2026-01-30 17:23:42] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:23:42] Perform `Seurat::FindVariableFeatures()` on the data 1/4 of the `srt_list`...
#> ℹ [2026-01-30 17:23:43] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:23:43] Perform `Seurat::FindVariableFeatures()` on the data 2/4 of the `srt_list`...
#> ℹ [2026-01-30 17:23:43] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:23:43] Perform `Seurat::FindVariableFeatures()` on the data 3/4 of the `srt_list`...
#> ℹ [2026-01-30 17:23:44] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:23:44] Perform `Seurat::FindVariableFeatures()` on the data 4/4 of the `srt_list`...
#> ℹ [2026-01-30 17:23:44] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:23:45] Number of available HVF: 2000
#> ℹ [2026-01-30 17:23:45] Finished check
#> ℹ [2026-01-30 17:23:46] Perform Uncorrected integration
#> ℹ [2026-01-30 17:23:47] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:23:47] Perform linear dimension reduction("pca")
#> ℹ [2026-01-30 17:23:49] Perform Seurat::FindClusters ("louvain")
#> ℹ [2026-01-30 17:23:49] Reorder clusters...
#> ℹ [2026-01-30 17:23:49] Perform nonlinear dimension reduction ("umap")
#> ℹ [2026-01-30 17:23:49] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ℹ [2026-01-30 17:23:55] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ✔ [2026-01-30 17:24:01] Run Uncorrected integration done
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunPCAMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Uncorrectedpca",
  ref_umap = "UncorrectedUMAP2D"
)
#> ℹ [2026-01-30 17:24:02] Data type is log-normalized
#> ℹ [2026-01-30 17:24:02] Detected srt_query data type: log_normalized_counts
#> ℹ [2026-01-30 17:24:03] Data type is log-normalized
#> ℹ [2026-01-30 17:24:03] Detected srt_ref data type: log_normalized_counts
#> ℹ [2026-01-30 17:24:03] Run PCA projection
#> ℹ [2026-01-30 17:24:03] Use [1] 2000 features to calculate PC.
#> ℹ [2026-01-30 17:24:04] Run UMAP projection
#> ℹ [2026-01-30 17:24:04] Use the reduction to calculate distance metric
#> ℹ [2026-01-30 17:24:04] Use raw method to find neighbors
#> ℹ [2026-01-30 17:24:04] Running UMAP projection
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
