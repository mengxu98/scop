# Single-cell reference mapping with Seurat method

Single-cell reference mapping with Seurat method

## Usage

``` r
RunSeuratMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_pca = NULL,
  ref_assay = srt_ref[[ref_pca]]@assay.used,
  ref_dims = 1:30,
  ref_umap = NULL,
  ref_group = NULL,
  normalization.method = "LogNormalize",
  reduction_project_method = "pcaproject",
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  k.weight = 100,
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

- ref_pca:

  A character string specifying the name of the PCA reduction in the
  reference object to use for calculating the distance metric.

- ref_assay:

  A character string specifying the assay name for the reference cells.
  If not provided, the default assay for the reference object will be
  used.

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

- normalization.method:

  The normalization method to use. Default is \`"LogNormalize"\`.

- reduction_project_method:

  Dimensional reduction to perform when finding anchors. Default is
  \`"pcaproject"\`.

- k.anchor:

  How many neighbors (k) to use when finding anchors. Default is \`5\`.

- k.filter:

  How many neighbors (k) to use when filtering anchors. Set to NA to
  turn off filtering. Default is \`200\`.

- k.score:

  How many neighbors (k) to use when scoring anchors. Default is \`30\`.

- k.weight:

  Number of neighbors to consider when weighting anchors. Default is
  \`100\`.

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

## See also

\[RunKNNMap\]

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-01-29 13:33:32] Start standard scop workflow...
#> ℹ [2026-01-29 13:33:32] Checking a list of <Seurat>...
#> ! [2026-01-29 13:33:33] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:33:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-29 13:33:35] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-29 13:33:35] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:33:36] Number of available HVF: 2000
#> ℹ [2026-01-29 13:33:36] Finished check
#> ℹ [2026-01-29 13:33:36] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:33:37] Perform pca linear dimension reduction
#> ℹ [2026-01-29 13:33:38] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-29 13:33:38] Reorder clusters...
#> ℹ [2026-01-29 13:33:38] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-29 13:33:38] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-29 13:33:43] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-29 13:33:48] Run scop standard workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-01-29 13:33:48] Run Uncorrected integration...
#> ℹ [2026-01-29 13:33:48] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-29 13:33:49] Checking a list of <Seurat>...
#> ℹ [2026-01-29 13:33:49] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:33:49] Perform `Seurat::FindVariableFeatures()` on the data 1/4 of the `srt_list`...
#> ℹ [2026-01-29 13:33:50] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:33:50] Perform `Seurat::FindVariableFeatures()` on the data 2/4 of the `srt_list`...
#> ℹ [2026-01-29 13:33:50] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:33:50] Perform `Seurat::FindVariableFeatures()` on the data 3/4 of the `srt_list`...
#> ℹ [2026-01-29 13:33:51] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:33:51] Perform `Seurat::FindVariableFeatures()` on the data 4/4 of the `srt_list`...
#> ℹ [2026-01-29 13:33:51] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:33:52] Number of available HVF: 2000
#> ℹ [2026-01-29 13:33:52] Finished check
#> ℹ [2026-01-29 13:33:53] Perform Uncorrected integration
#> ℹ [2026-01-29 13:33:54] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:33:54] Perform linear dimension reduction("pca")
#> ℹ [2026-01-29 13:33:56] Perform Seurat::FindClusters ("louvain")
#> ℹ [2026-01-29 13:33:56] Reorder clusters...
#> ℹ [2026-01-29 13:33:56] Perform nonlinear dimension reduction ("umap")
#> ℹ [2026-01-29 13:33:56] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ℹ [2026-01-29 13:34:01] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ✔ [2026-01-29 13:34:07] Run Uncorrected integration done
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSeuratMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Uncorrectedpca",
  ref_umap = "UncorrectedUMAP2D",
  k.weight = 50
)
#> ℹ [2026-01-29 13:34:08] Data type is log-normalized
#> ℹ [2026-01-29 13:34:08] Detected srt_query data type: log_normalized_counts
#> ℹ [2026-01-29 13:34:08] Data type is log-normalized
#> ℹ [2026-01-29 13:34:08] Detected srt_ref data type: log_normalized_counts
#> ℹ [2026-01-29 13:34:08] Run FindTransferAnchors
#> Projecting cell embeddings
#> Finding neighborhoods
#> Finding anchors
#>  Found 471 anchors
#> Filtering anchors
#>  Retained 471 anchors
#> Requested to reuse weights matrix, but no weights found. Computing new weights.
#> 
#> Integrating dataset 2 with reference dataset
#> Finding integration vectors
#> Finding integration vector weights
#> Integrating data
#> ℹ [2026-01-29 13:34:13] Run UMAP projection
#> ℹ [2026-01-29 13:34:13] Use the reduction to calculate distance metric
#> ℹ [2026-01-29 13:34:13] Use raw method to find neighbors
#> ℹ [2026-01-29 13:34:14] Running UMAP projection
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
