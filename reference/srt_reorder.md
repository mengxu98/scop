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
#> ℹ [2026-01-30 17:40:20] Start standard scop workflow...
#> ℹ [2026-01-30 17:40:21] Checking a list of <Seurat>...
#> ! [2026-01-30 17:40:21] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:40:21] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:40:23] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:40:24] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:40:24] Number of available HVF: 2000
#> ℹ [2026-01-30 17:40:24] Finished check
#> ℹ [2026-01-30 17:40:24] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:40:25] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:40:26] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:40:26] Reorder clusters...
#> ℹ [2026-01-30 17:40:26] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:40:26] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:40:31] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:40:36] Run scop standard workflow completed
pancreas_sub <- srt_reorder(
  pancreas_sub,
  reorder_by = "SubCellType",
  layer = "data"
)
```
