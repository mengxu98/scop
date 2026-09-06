# Run UMAP

Run UMAP

## Usage

``` r
RunUMAP2(object, ...)

# S3 method for class 'Seurat'
RunUMAP2(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  neighbor = NULL,
  graph = NULL,
  assay = NULL,
  layer = "data",
  umap.method = "uwot",
  reduction.model = NULL,
  n_threads = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = 200L,
  cores = 1,
  min.dist = 0.3,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  learning.rate = 1,
  repulsion.strength = 1,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunUMAP2(
  object,
  assay = NULL,
  umap.method = "uwot",
  reduction.model = NULL,
  n_threads = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = 200L,
  cores = 1,
  min.dist = 0.3,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  learning.rate = 1,
  repulsion.strength = 1,
  reduction.key = "UMAP_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
```

## Arguments

- object:

  A `Seurat` object, matrix-like object, `Neighbor`, or `Graph`.

- ...:

  Passed to the UMAP implementation.

- reduction:

  Linear reduction used as input.

- dims:

  Dimensions to use. Supply only one of `dims`, `features`, `neighbor`,
  or `graph`.

- features:

  Features used instead of a reduction.

- neighbor, graph:

  Existing `Neighbor` or `Graph` object name.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- umap.method:

  `"uwot"` or `"naive"`.

- reduction.model:

  Pre-trained UMAP `DimReduc` used to embed new data.

- n_threads:

  Number of threads.

- return.model:

  Store the UMAP model.

- n.neighbors, n.components, metric, n.epochs:

  UMAP layout parameters. String `metric` values include `"euclidean"`,
  `"manhattan"`, `"cosine"`, `"pearson"`, and `"pearson2"`.

- cores:

  Number of CPU cores.

- min.dist:

  Minimum embedding distance (how tightly points pack).

- set.op.mix.ratio:

  Mix of fuzzy union (`1`) and intersection (`0`).

- local.connectivity, negative.sample.rate:

  Fuzzy simplicial set and optimization sampling.

- a, b:

  UMAP curve parameters. `NULL` estimates them automatically.

- learning.rate, repulsion.strength:

  Layout optimization rates.

- reduction.name, reduction.key:

  Stored reduction name and embedding prefix.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:38:12] Start standard processing workflow...
#> ℹ [2026-09-06 22:38:13] Checking a list of <Seurat>...
#> ! [2026-09-06 22:38:13] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:38:13] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:38:13] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:38:13] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:38:13] Number of available HVF: 2000
#> ℹ [2026-09-06 22:38:13] Finished check
#> ℹ [2026-09-06 22:38:13] Perform `ScaleData()`
#> ℹ [2026-09-06 22:38:13] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:38:13] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:38:14] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:38:14] Reorder clusters...
#> ℹ [2026-09-06 22:38:14] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:38:14] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:38:22] Standard processing workflow completed
pancreas_sub <- RunUMAP2(pancreas_sub, dims = 1:30)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "umap"
)
```
