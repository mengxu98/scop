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
  vote_fun = "mean",
  verbose = TRUE
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

- ref_pca:

  Name of a PCA reduction in the reference object to use for calculating
  the distance metric. If NULL (default), it will be automatically
  detected as the first PCA reduction.

- ref_dims:

  Dimension indices from the reference reduction to be used for
  calculating the distance metric.

- ref_umap:

  Name of the UMAP reduction in the reference object. If not provided,
  the first UMAP reduction found in the reference object will be used.

- ref_group:

  The grouping variable in the reference object. This variable will be
  used to group cells in the heatmap columns. If not provided, all cells
  will be treated as one group.

- projection_method:

  Projection method to use. Options are "model" and "knn". If "model" is
  selected, the function will try to use a pre-trained UMAP model in the
  reference object for projection. If "knn" is selected, the function
  will directly find the nearest neighbors using the distance metric.

- nn_method:

  Nearest neighbor search method to use. Options are "raw", "annoy", and
  "rann". If "raw" is selected, the function will use the brute-force
  method to find the nearest neighbors. If "annoy" is selected, the
  function will use the Annoy library for approximate nearest neighbor
  search. If "rann" is selected, the function will use the RANN library
  for approximate nearest neighbor search. If not provided, the function
  will choose the search method based on the size of the query and
  reference datasets. For finite dense inputs using cosine or Euclidean
  distance, the raw path uses a bounded-memory compiled top-k kernel
  unless the full distance matrix is requested. Other metrics and sparse
  inputs retain the existing proxyC path.

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

  Function to be used for aggregating the nearest neighbors in the
  reference object. Options are "mean", "median", "sum", "min", "max",
  "sd", "var", etc. If not provided, the default is "mean".

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- RunStandardWorkflow(panc8_sub)
#> ℹ [2026-09-06 22:29:37] Start standard processing workflow...
#> ℹ [2026-09-06 22:29:38] Checking a list of <Seurat>...
#> ! [2026-09-06 22:29:38] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:29:38] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:29:38] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:29:38] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:29:38] Number of available HVF: 2000
#> ℹ [2026-09-06 22:29:38] Finished check
#> ℹ [2026-09-06 22:29:38] Perform `ScaleData()`
#> ℹ [2026-09-06 22:29:38] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:29:39] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-09-06 22:29:39] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:29:39] Reorder clusters...
#> ℹ [2026-09-06 22:29:39] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:29:39] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:29:47] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- RunIntegration(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-09-06 22:29:48] Run integration workflow...
#> ℹ [2026-09-06 22:29:48] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 22:29:49] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:29:49] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:29:49] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-09-06 22:29:49] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:29:49] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-09-06 22:29:50] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:29:50] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-09-06 22:29:50] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:29:50] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-09-06 22:29:50] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:29:50] Number of available HVF: 2000
#> ℹ [2026-09-06 22:29:50] Finished check
#> ℹ [2026-09-06 22:29:50] Perform Uncorrected integration
#> ℹ [2026-09-06 22:29:50] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 22:29:50] Perform "pca" linear dimension reduction
#> ! [2026-09-06 22:29:50] Some PCA features are absent from scale.data and will be dropped: "NPTX2", "MMP14", "VCAN", "LTB", "ALDH1A3", "NOTCH3", "DEFB1", "COL18A1", "COL5A3", "S100A6", "HABP2", "PLAT", "SERPINF1", "GPX2", "TFPI2", "SRXN1", "CYR61", "NCAM1", …, "AARS2", and "POMC"
#> ℹ [2026-09-06 22:29:51] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 22:29:51] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 22:29:51] Reorder clusters...
#> ℹ [2026-09-06 22:29:51] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:29:51] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:24)
#> ℹ [2026-09-06 22:29:57] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:24)
#> ℹ [2026-09-06 22:30:02] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:24)
#> ✔ [2026-09-06 22:30:08] Uncorrected integration completed
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunPCAMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Uncorrectedpca",
  ref_umap = "UncorrectedUMAP2D"
)
#> ℹ [2026-09-06 22:30:09] Data type is log-normalized
#> ℹ [2026-09-06 22:30:09] Detected srt_query data type: log_normalized_counts
#> ℹ [2026-09-06 22:30:09] Data type is log-normalized
#> ℹ [2026-09-06 22:30:09] Detected srt_ref data type: log_normalized_counts
#> ℹ [2026-09-06 22:30:09] Run PCA projection
#> ℹ [2026-09-06 22:30:09] Use [1] 675 features to calculate PC.
#> ℹ [2026-09-06 22:30:09] Run UMAP projection
#> ℹ [2026-09-06 22:30:09] Use the reduction to calculate distance metric
#> ℹ [2026-09-06 22:30:09] Use raw method to find neighbors
#> ℹ [2026-09-06 22:30:09] Running UMAP projection
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
