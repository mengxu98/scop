# Single-cell reference mapping with CSS method

Single-cell reference mapping with CSS method

## Usage

``` r
RunCSSMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = srt_ref[[ref_css]]@assay.used,
  ref_css = NULL,
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

- ref_css:

  CSS reduction in the reference object to use for calculating the
  distance metric.

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

## See also

\[RunKNNMap\]

## Examples

``` r
data(panc8_sub)
panc8_sub <- RunStandardWorkflow(panc8_sub)
#> ℹ [2026-10-04 22:10:53] Start standard processing workflow...
#> ℹ [2026-10-04 22:10:53] Checking a list of <Seurat>...
#> ! [2026-10-04 22:10:53] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 22:10:53] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:10:53] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:10:53] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:10:53] Number of available HVF: 2000
#> ℹ [2026-10-04 22:10:53] Finished check
#> ℹ [2026-10-04 22:10:53] Perform `ScaleData()`
#> ℹ [2026-10-04 22:10:54] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:10:54] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-10-04 22:10:54] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:10:55] Reorder clusters...
#> ℹ [2026-10-04 22:10:55] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:10:55] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 22:11:03] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- RunIntegration(
  srt_ref,
  batch = "tech",
  integration_method = "CSS"
)
#> ◌ [2026-10-04 22:11:03] Run integration workflow...
#> ℹ [2026-10-04 22:11:04] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-04 22:11:04] Checking a list of <Seurat>...
#> ℹ [2026-10-04 22:11:04] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:11:04] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-10-04 22:11:05] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:11:05] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-10-04 22:11:05] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:11:05] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-10-04 22:11:05] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:11:05] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-10-04 22:11:05] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:11:05] Number of available HVF: 2000
#> ℹ [2026-10-04 22:11:05] Finished check
#> ℹ [2026-10-04 22:11:14] Perform ScaleData
#> ℹ [2026-10-04 22:11:15] Perform "pca" linear dimension reduction
#> ℹ [2026-10-04 22:11:15] Perform CSS integration
#> ℹ [2026-10-04 22:11:15] Using "CSSpca" (1:20) as input
#> Loading required package: Seurat
#> Loading required package: SeuratObject
#> Loading required package: sp
#> ‘SeuratObject’ was built under R 4.6.0 but the current version is
#> 4.6.1; it is recomended that you reinstall ‘SeuratObject’ as the ABI
#> for R may have changed
#> 
#> Attaching package: ‘SeuratObject’
#> The following object is masked from ‘package:BiocGenerics’:
#> 
#>     intersect
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t
#> v5.6 adds significant improvements in speed and efficiency, particularly for large datasets.
#> For more information, see https://satijalab.org/seurat/articles/announcements.
#> 
#> Attaching package: ‘Seurat’
#> The following objects are masked from ‘package:scop’:
#> 
#>     AddModuleScore, AggregateExpression, AverageExpression,
#>     CellCycleScoring, FindAllMarkers, FindClusters, FindMarkers,
#>     FindMultiModalNeighbors, FindNeighbors, FindVariableFeatures,
#>     FoldChange, NormalizeData, RunCCA, RunPCA, RunUMAP, SCTransform,
#>     ScaleData, SelectIntegrationFeatures
#> Error in .Call("_Seurat_ComputeSNN", PACKAGE = "Seurat", nn, jaccard_prune): Incorrect number of arguments (2), expecting 3 for '_Seurat_ComputeSNN'
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunCSSMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_css = "CSS",
  ref_umap = "CSSUMAP2D"
)
#> Error in srt_ref[[ref_css]]: ‘CSS’ not found in this Seurat object
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
