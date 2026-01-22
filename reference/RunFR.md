# Run Force-Directed Layout (Fruchterman-Reingold algorithm)

Run Force-Directed Layout (Fruchterman-Reingold algorithm)

## Usage

``` r
RunFR(object, ...)

# S3 method for class 'Seurat'
RunFR(
  object,
  reduction = NULL,
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  graph = NULL,
  neighbor = NULL,
  k.param = 20,
  ndim = 2,
  niter = 500,
  reduction.name = "FR",
  reduction.key = "FR_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)

# Default S3 method
RunFR(
  object,
  assay = NULL,
  ndim = 2,
  niter = 500,
  reduction.key = "FR_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, a Neighbor object, or a Graph
  object. Default is `NULL`.

- ...:

  Additional arguments to be passed to
  [igraph::layout_with_fr](https://r.igraph.org/reference/layout_with_fr.html).

- reduction:

  Which dimensionality reduction to use. Default is `"pca"`.

- dims:

  The dimensions to be used. Default is `NULL`.

- features:

  A character vector of features to use. Default is `NULL`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `data`.

- graph:

  The name of the Graph object to be used. Default is `NULL`.

- neighbor:

  The name of the Neighbor object to be used. Default is `NULL`.

- k.param:

  The number of nearest neighbors to consider. Default is `20`.

- ndim:

  The number of dimensions for the force-directed layout. Default is
  `2`.

- niter:

  The number of iterations for the force-directed layout. Default is
  `500`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"fr"`.

- reduction.key:

  The prefix for the column names of the force-directed layout
  embeddings. Default is `"FR_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-22 04:04:18] Start standard scop workflow...
#> ℹ [2026-01-22 04:04:19] Checking a list of <Seurat>...
#> ! [2026-01-22 04:04:19] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 04:04:19] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:04:22] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:04:22] Use the separate HVF from srt_list
#> ℹ [2026-01-22 04:04:22] Number of available HVF: 2000
#> ℹ [2026-01-22 04:04:22] Finished check
#> ℹ [2026-01-22 04:04:23] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 04:04:23] Perform pca linear dimension reduction
#> ℹ [2026-01-22 04:04:24] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 04:04:24] Reorder clusters...
#> ℹ [2026-01-22 04:04:24] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 04:04:24] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 04:04:29] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 04:04:33] Run scop standard workflow completed
pancreas_sub <- RunFR(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> ℹ [2026-01-22 04:04:33] Running force-directed layout
#> ℹ [2026-01-22 04:04:33] Computing nearest neighbor graph and SNN
#> ℹ [2026-01-22 04:04:39] Force-directed layout computed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "fr"
)
```
