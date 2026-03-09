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
#> ℹ [2026-03-09 08:33:08] Start standard scop workflow...
#> ℹ [2026-03-09 08:33:09] Checking a list of <Seurat>...
#> ! [2026-03-09 08:33:09] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-09 08:33:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:33:11] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:33:12] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:33:12] Number of available HVF: 2000
#> ℹ [2026-03-09 08:33:12] Finished check
#> ℹ [2026-03-09 08:33:12] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:33:12] Perform pca linear dimension reduction
#> ℹ [2026-03-09 08:33:14] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-09 08:33:14] Reorder clusters...
#> ℹ [2026-03-09 08:33:14] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-09 08:33:14] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-09 08:33:19] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-09 08:33:23] Run scop standard workflow completed
pancreas_sub <- RunFR(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> ℹ [2026-03-09 08:33:23] Running force-directed layout
#> ℹ [2026-03-09 08:33:24] Computing nearest neighbor graph and SNN
#> ℹ [2026-03-09 08:33:29] Force-directed layout computed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "fr"
)
```
