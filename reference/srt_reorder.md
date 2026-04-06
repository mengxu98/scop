# Reorder idents by the gene expression

Reorder idents by the gene expression

## Usage

``` r
srt_reorder(
  srt,
  features = NULL,
  reorder_by = NULL,
  layer = "data",
  assay = NULL,
  log = TRUE,
  distance_metric = "euclidean",
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A character vector or a named list of features to plot. Features can
  be gene names in Assay or names of numeric columns in meta.data.

- reorder_by:

  Reorder groups instead of idents.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- log:

  Whether [log1p](https://rdrr.io/r/base/Log.html) transformation needs
  to be applied. Default is `TRUE`.

- distance_metric:

  Metric to compute distance. Default is `"euclidean"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-06 05:26:58] Start standard processing workflow...
#> ℹ [2026-04-06 05:26:59] Checking a list of <Seurat>...
#> ! [2026-04-06 05:26:59] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:26:59] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 05:27:02] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 05:27:02] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:27:02] Number of available HVF: 2000
#> ℹ [2026-04-06 05:27:03] Finished check
#> ℹ [2026-04-06 05:27:03] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 05:27:03] Perform pca linear dimension reduction
#> ℹ [2026-04-06 05:27:04] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-06 05:27:04] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-06 05:27:04] Reorder clusters...
#> ℹ [2026-04-06 05:27:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:27:04] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-06 05:27:04] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-06 05:27:09] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-06 05:27:14] Standard processing workflow completed
pancreas_sub <- srt_reorder(
  pancreas_sub,
  reorder_by = "SubCellType",
  layer = "data"
)
#> ℹ [2026-04-06 05:27:14] Skip `log1p()` because `layer = data` is not "counts"
```
