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
  object.

- ...:

  Passed to
  [igraph::layout_with_fr](https://r.igraph.org/reference/layout_with_fr.html).

- reduction:

  Linear reduction used as input.

- dims:

  Dimensions to use. Supply only one of `dims`, `features`, `neighbor`,
  or `graph`.

- features:

  Features used instead of a reduction.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- graph:

  Neighbor-graph edges.

- neighbor:

  Neighbor object to be used.

- k.param:

  The number of nearest neighbors to consider.

- ndim:

  The number of dimensions for the force-directed layout.

- niter:

  The number of iterations for the force-directed layout.

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the force-directed layout
  embeddings.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:08:11] Start standard processing workflow...
#> ℹ [2026-09-06 22:08:12] Checking a list of <Seurat>...
#> ! [2026-09-06 22:08:12] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:08:12] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:08:12] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:08:12] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:08:12] Number of available HVF: 2000
#> ℹ [2026-09-06 22:08:12] Finished check
#> ℹ [2026-09-06 22:08:12] Perform `ScaleData()`
#> ℹ [2026-09-06 22:08:12] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:08:13] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:08:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:08:13] Reorder clusters...
#> ℹ [2026-09-06 22:08:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:08:13] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:08:20] Standard processing workflow completed
pancreas_sub <- RunFR(
  object = pancreas_sub,
  graph = "Standardpca_SNN",
  niter = 100
)
#> ℹ [2026-09-06 22:08:20] Running force-directed layout
#> ℹ [2026-09-06 22:08:20] Force-directed layout computed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "fr"
)
```
