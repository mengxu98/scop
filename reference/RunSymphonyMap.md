# Single-cell reference mapping with Symphony method

Single-cell reference mapping with Symphony method

## Usage

``` r
RunSymphonyMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = srt_ref[[ref_pca]]@assay.used,
  ref_pca = NULL,
  ref_harmony = NULL,
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

  The PCA reduction in the reference object to use for calculating the
  distance metric.

- ref_harmony:

  The Harmony reduction in the reference object to use for calculating
  the distance metric.

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
#> ℹ [2026-09-20 22:52:15] Start standard processing workflow...
#> ℹ [2026-09-20 22:52:15] Checking a list of <Seurat>...
#> ! [2026-09-20 22:52:15] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:52:15] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:52:15] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:52:15] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:52:16] Number of available HVF: 2000
#> ℹ [2026-09-20 22:52:16] Finished check
#> ℹ [2026-09-20 22:52:16] Perform `ScaleData()`
#> ℹ [2026-09-20 22:52:16] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:52:16] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-09-20 22:52:17] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:52:17] Reorder clusters...
#> ℹ [2026-09-20 22:52:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:52:17] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:52:26] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- RunIntegration(
  srt_ref,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-09-20 22:52:26] Run integration workflow...
#> ℹ [2026-09-20 22:52:27] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-20 22:52:28] Checking a list of <Seurat>...
#> ℹ [2026-09-20 22:52:28] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-20 22:52:28] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-09-20 22:52:28] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-20 22:52:28] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-09-20 22:52:28] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-20 22:52:28] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-09-20 22:52:28] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-09-20 22:52:28] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-09-20 22:52:29] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:52:29] Number of available HVF: 2000
#> ℹ [2026-09-20 22:52:29] Finished check
#> ℹ [2026-09-20 22:52:29] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-20 22:52:29] Perform linear dimension reduction("pca")
#> ℹ [2026-09-20 22:52:29] Perform Harmony integration
#> ℹ [2026-09-20 22:52:29] Using "Harmonypca" (1:20) as input
#> ℹ [2026-09-20 22:52:29] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-20 22:52:30] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-20 22:52:30] Reorder clusters...
#> ℹ [2026-09-20 22:52:30] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:52:30] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-20 22:52:36] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-20 22:52:42] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ✔ [2026-09-20 22:52:49] Harmony integration completed
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSymphonyMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Harmonypca",
  ref_harmony = "Harmony",
  ref_umap = "HarmonyUMAP2D"
)
#> ℹ [2026-09-20 22:53:20] Data type is log-normalized
#> ℹ [2026-09-20 22:53:20] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-09-20 22:53:20] Data type is log-normalized
#> ℹ [2026-09-20 22:53:20] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-09-20 22:53:20] Build reference
#> ℹ [2026-09-20 22:53:20] Saved embeddings
#> ℹ [2026-09-20 22:53:20] Saved soft cluster assignments
#> ℹ [2026-09-20 22:53:20] Saved variable gene information for 2000 genes
#> ℹ [2026-09-20 22:53:20] Saved PCA loadings
#> ℹ [2026-09-20 22:53:20] Saved metadata
#> ℹ [2026-09-20 22:53:20] Calculate final L2 normalized reference centroids (Y_cos)
#> ! [2026-09-20 22:53:20] Function "cosine_normalize" not found in symphony namespace
#> ℹ [2026-09-20 22:53:20] Calculate reference compression terms (Nr and C)
#> ℹ [2026-09-20 22:53:20] Run mapQuery
#> ℹ [2026-09-20 22:53:20] Scaling and synchronizing query gene expression
#> ℹ [2026-09-20 22:53:20] Found 2000 reference variable genes in query object
#> ℹ [2026-09-20 22:53:20] Project query cells using reference gene loadings
#> ℹ [2026-09-20 22:53:21] Clustering query cells to reference centroids
#> ! [2026-09-20 22:53:21] Function "cosine_normalize" not found in symphony namespace
#> ℹ [2026-09-20 22:53:21] Correcting query batch effects
#> ℹ [2026-09-20 22:53:21] Run UMAP projection
#> ℹ [2026-09-20 22:53:21] Use the reduction to calculate distance metric
#> ℹ [2026-09-20 22:53:21] Use raw method to find neighbors
#> ℹ [2026-09-20 22:53:21] Running UMAP projection
#> ℹ [2026-09-20 22:53:21] Run SymphonyMap finished
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
