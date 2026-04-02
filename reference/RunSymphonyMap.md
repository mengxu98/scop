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

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  The grouping variable in the reference object. This variable will be
  used to group cells in the heatmap columns. If not provided, all cells
  will be treated as one group.

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
#> ℹ [2026-04-02 16:55:11] Start standard processing workflow...
#> ℹ [2026-04-02 16:55:12] Checking a list of <Seurat>...
#> ! [2026-04-02 16:55:12] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:55:12] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:55:14] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:55:14] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:55:14] Number of available HVF: 2000
#> ℹ [2026-04-02 16:55:15] Finished check
#> ℹ [2026-04-02 16:55:15] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:55:15] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:55:19] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:55:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:55:20] Reorder clusters...
#> ℹ [2026-04-02 16:55:20] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:55:20] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:55:20] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:55:20] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:55:20] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:55:23] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:55:27] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-04-02 16:55:27] Run integration workflow...
#> ℹ [2026-04-02 16:55:27] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-02 16:55:28] Checking a list of <Seurat>...
#> ℹ [2026-04-02 16:55:28] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-02 16:55:28] Perform `Seurat::FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-04-02 16:55:29] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-02 16:55:29] Perform `Seurat::FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-04-02 16:55:29] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-02 16:55:29] Perform `Seurat::FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-04-02 16:55:30] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-02 16:55:30] Perform `Seurat::FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-04-02 16:55:30] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:55:31] Number of available HVF: 2000
#> ℹ [2026-04-02 16:55:31] Finished check
#> ℹ [2026-04-02 16:55:33] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:55:33] Perform linear dimension reduction("pca")
#> ℹ [2026-04-02 16:55:37] Perform Harmony integration
#> ℹ [2026-04-02 16:55:37] Using "CSSpca" (1:50) as input
#> Error in loadNamespace(x): there is no package called ‘harmony’
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSymphonyMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Harmonypca",
  ref_harmony = "Harmony",
  ref_umap = "HarmonyUMAP2D"
)
#> Error in srt_ref[[ref_pca]]: ‘Harmonypca’ not found in this Seurat object
#>  
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
#> Error in srt_query[[query_reduction]]: ‘ref.embeddings’ not found in this Seurat object
#>  
```
