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

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-06-29 04:40:00] Start standard processing workflow...
#> ℹ [2026-06-29 04:40:01] Checking a list of <Seurat>...
#> ! [2026-06-29 04:40:01] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-29 04:40:01] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-29 04:40:01] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-29 04:40:01] Use the separate HVF from `srt_list`
#> ℹ [2026-06-29 04:40:01] Number of available HVF: 2000
#> ℹ [2026-06-29 04:40:01] Finished check
#> ℹ [2026-06-29 04:40:01] Perform `ScaleData()`
#> ℹ [2026-06-29 04:40:01] Perform pca linear dimension reduction
#> ℹ [2026-06-29 04:40:03] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-06-29 04:40:03] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-29 04:40:03] Reorder clusters...
#> ℹ [2026-06-29 04:40:03] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-29 04:40:03] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-29 04:40:12] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-06-29 04:40:12] Run integration workflow...
#> ℹ [2026-06-29 04:40:13] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-06-29 04:40:13] Checking a list of <Seurat>...
#> ℹ [2026-06-29 04:40:14] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-29 04:40:14] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-06-29 04:40:14] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-29 04:40:14] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-06-29 04:40:14] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-29 04:40:14] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-06-29 04:40:15] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-29 04:40:15] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-06-29 04:40:15] Use the separate HVF from `srt_list`
#> ℹ [2026-06-29 04:40:15] Number of available HVF: 2000
#> ℹ [2026-06-29 04:40:16] Finished check
#> ℹ [2026-06-29 04:40:16] Perform `Seurat::ScaleData()`
#> ℹ [2026-06-29 04:40:16] Perform linear dimension reduction("pca")
#> ! [2026-06-29 04:40:16] Some PCA features are absent from scale.data and will be dropped: "G0S2", "MRC2", "COL18A1", "COL5A3", "NOTCH3", "CRLF1", "PTGR1", "IL11", "LIF", "PDLIM3", "HTRA1", "TFPI2", "NREP", "ENG", "AQP3", "SEMA7A", "NPTX2", "SNAI2", …, "PCDH18", and "SCN11A"
#> ℹ [2026-06-29 04:40:16] Perform Harmony integration
#> ℹ [2026-06-29 04:40:16] Using "Harmonypca" (1:25) as input
#> ℹ [2026-06-29 04:40:16] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-06-29 04:40:17] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-06-29 04:40:17] Reorder clusters...
#> ℹ [2026-06-29 04:40:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-29 04:40:17] Perform umap nonlinear dimension reduction using Harmony (1:25)
#> ℹ [2026-06-29 04:40:23] Perform umap nonlinear dimension reduction using Harmony (1:25)
#> ℹ [2026-06-29 04:40:29] Perform umap nonlinear dimension reduction using Harmonypca (1:25)
#> ✔ [2026-06-29 04:40:35] Harmony integration completed
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSymphonyMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Harmonypca",
  ref_harmony = "Harmony",
  ref_umap = "HarmonyUMAP2D"
)
#> ℹ [2026-06-29 04:41:08] Data type is log-normalized
#> ℹ [2026-06-29 04:41:08] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-06-29 04:41:09] Data type is log-normalized
#> ℹ [2026-06-29 04:41:09] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-06-29 04:41:09] Build reference
#> ℹ [2026-06-29 04:41:09] Saved embeddings
#> ℹ [2026-06-29 04:41:09] Saved soft cluster assignments
#> ℹ [2026-06-29 04:41:09] Saved variable gene information for 636 genes
#> ℹ [2026-06-29 04:41:09] Saved PCA loadings
#> ℹ [2026-06-29 04:41:09] Saved metadata
#> ℹ [2026-06-29 04:41:09] Calculate final L2 normalized reference centroids (Y_cos)
#> ! [2026-06-29 04:41:09] Function "cosine_normalize" not found in symphony namespace
#> ℹ [2026-06-29 04:41:09] Calculate reference compression terms (Nr and C)
#> ℹ [2026-06-29 04:41:09] Run mapQuery
#> ℹ [2026-06-29 04:41:09] Scaling and synchronizing query gene expression
#> ℹ [2026-06-29 04:41:09] Found 636 reference variable genes in query object
#> ℹ [2026-06-29 04:41:09] Project query cells using reference gene loadings
#> ℹ [2026-06-29 04:41:09] Clustering query cells to reference centroids
#> ! [2026-06-29 04:41:09] Function "cosine_normalize" not found in symphony namespace
#> ℹ [2026-06-29 04:41:09] Correcting query batch effects
#> ℹ [2026-06-29 04:41:09] Run UMAP projection
#> ℹ [2026-06-29 04:41:09] Use the reduction to calculate distance metric
#> ℹ [2026-06-29 04:41:09] Use raw method to find neighbors
#> ℹ [2026-06-29 04:41:09] Running UMAP projection
#> ℹ [2026-06-29 04:41:09] Run SymphonyMap finished
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
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```
