# Reorder idents by the gene expression

Reorder idents by the gene expression

## Usage

``` r
srt_reorder(
  object,
  features = NULL,
  reorder_by = NULL,
  layer = "data",
  assay = NULL,
  log = TRUE,
  distance_metric = "euclidean",
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- features:

  Features to plot: a character vector or a named list of assay gene
  names or numeric metadata columns.

- reorder_by:

  Reorder groups instead of idents.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- log:

  Whether [log1p](https://rdrr.io/r/base/Log.html) transformation needs
  to be applied.

- distance_metric:

  Metric to compute distance.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-13 22:57:18] Start standard processing workflow...
#> ℹ [2026-09-13 22:57:18] Checking a list of <Seurat>...
#> ! [2026-09-13 22:57:18] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-13 22:57:18] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:57:18] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:57:19] Use the separate HVF from `srt_list`
#> ℹ [2026-09-13 22:57:19] Number of available HVF: 2000
#> ℹ [2026-09-13 22:57:19] Finished check
#> ℹ [2026-09-13 22:57:19] Perform `ScaleData()`
#> ℹ [2026-09-13 22:57:19] Perform pca linear dimension reduction
#> ℹ [2026-09-13 22:57:19] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-13 22:57:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-13 22:57:20] Reorder clusters...
#> ℹ [2026-09-13 22:57:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-13 22:57:20] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-13 22:57:28] Standard processing workflow completed
pancreas_sub <- srt_reorder(
  pancreas_sub,
  reorder_by = "SubCellType",
  layer = "data"
)
#> ℹ [2026-09-13 22:57:28] Skip `log1p()` because `layer = data` is not "counts"
```
