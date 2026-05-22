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
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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
#> ℹ [2026-05-22 17:31:47] Start standard processing workflow...
#> ℹ [2026-05-22 17:31:47] Checking a list of <Seurat>...
#> ! [2026-05-22 17:31:47] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-22 17:31:47] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-22 17:31:49] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-22 17:31:50] Use the separate HVF from `srt_list`
#> ℹ [2026-05-22 17:31:50] Number of available HVF: 2000
#> ℹ [2026-05-22 17:31:50] Finished check
#> ℹ [2026-05-22 17:31:50] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-22 17:31:51] Perform pca linear dimension reduction
#> ℹ [2026-05-22 17:31:51] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-22 17:31:51] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-22 17:31:52] Reorder clusters...
#> ℹ [2026-05-22 17:31:52] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-22 17:31:52] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-22 17:31:52] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-22 17:31:57] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-22 17:32:02] Standard processing workflow completed
pancreas_sub <- RunFR(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> ℹ [2026-05-22 17:32:02] Running force-directed layout
#> ℹ [2026-05-22 17:32:02] Computing nearest neighbor graph and SNN
#> ℹ [2026-05-22 17:32:14] Force-directed layout computed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "fr"
)
```
