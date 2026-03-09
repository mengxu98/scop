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

  The PCA reduction in the reference object to use for calculating the
  distance metric.

- ref_harmony:

  The Harmony reduction in the reference object to use for calculating
  the distance metric.

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
#> ℹ [2026-03-09 08:48:15] Start standard scop workflow...
#> ℹ [2026-03-09 08:48:15] Checking a list of <Seurat>...
#> ! [2026-03-09 08:48:15] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-09 08:48:15] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:48:18] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:48:18] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:48:18] Number of available HVF: 2000
#> ℹ [2026-03-09 08:48:19] Finished check
#> ℹ [2026-03-09 08:48:19] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:48:19] Perform pca linear dimension reduction
#> ℹ [2026-03-09 08:48:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-09 08:48:21] Reorder clusters...
#> ℹ [2026-03-09 08:48:21] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-09 08:48:21] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-09 08:48:26] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-09 08:48:31] Run scop standard workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-03-09 08:48:31] Run Harmony integration...
#> ℹ [2026-03-09 08:48:31] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-03-09 08:48:32] Checking a list of <Seurat>...
#> ℹ [2026-03-09 08:48:32] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:48:32] Perform `Seurat::FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-03-09 08:48:33] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:48:33] Perform `Seurat::FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-03-09 08:48:33] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:48:33] Perform `Seurat::FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-03-09 08:48:34] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-03-09 08:48:34] Perform `Seurat::FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-03-09 08:48:34] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:48:35] Number of available HVF: 2000
#> ℹ [2026-03-09 08:48:35] Finished check
#> ℹ [2026-03-09 08:48:37] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:48:37] Perform linear dimension reduction("pca")
#> ℹ [2026-03-09 08:48:38] Perform Harmony integration
#> ℹ [2026-03-09 08:48:38] Using "CSSpca" (1:10) as input
#> ℹ [2026-03-09 08:48:40] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-03-09 08:48:40] Reorder clusters...
#> ℹ [2026-03-09 08:48:40] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ℹ [2026-03-09 08:48:45] Perform umap nonlinear dimension reduction using Harmony (1:10)
#> ✔ [2026-03-09 08:48:51] Run Harmony integration done
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSymphonyMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Harmonypca",
  ref_harmony = "Harmony",
  ref_umap = "HarmonyUMAP2D"
)
#> ℹ [2026-03-09 08:49:20] Data type is log-normalized
#> ℹ [2026-03-09 08:49:20] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-09 08:49:21] Data type is log-normalized
#> ℹ [2026-03-09 08:49:21] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-09 08:49:21] Build reference
#> ℹ [2026-03-09 08:49:21] Saved embeddings
#> ℹ [2026-03-09 08:49:21] Saved soft cluster assignments
#> ℹ [2026-03-09 08:49:21] Saved variable gene information for 2000 genes
#> ℹ [2026-03-09 08:49:21] Saved PCA loadings
#> ℹ [2026-03-09 08:49:21] Saved metadata
#> ℹ [2026-03-09 08:49:21] Calculate final L2 normalized reference centroids (Y_cos)
#> ℹ [2026-03-09 08:49:21] Calculate reference compression terms (Nr and C)
#> ℹ [2026-03-09 08:49:21] Run mapQuery
#> ℹ [2026-03-09 08:49:21] Scaling and synchronizing query gene expression
#> ℹ [2026-03-09 08:49:21] Found 2000 reference variable genes in query dataset
#> ℹ [2026-03-09 08:49:21] Project query cells using reference gene loadings
#> ℹ [2026-03-09 08:49:21] Clustering query cells to reference centroids
#> ℹ [2026-03-09 08:49:21] Correcting query batch effects
#> ℹ [2026-03-09 08:49:21] Run UMAP projection
#> ℹ [2026-03-09 08:49:21] Use the reduction to calculate distance metric
#> ℹ [2026-03-09 08:49:21] Use raw method to find neighbors
#> ℹ [2026-03-09 08:49:22] Running UMAP projection
#> ℹ [2026-03-09 08:49:22] Run SymphonyMap finished
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
