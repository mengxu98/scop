# Run MDS (multi-dimensional scaling)

Run MDS (multi-dimensional scaling)

## Usage

``` r
RunMDS(object, ...)

# S3 method for class 'Seurat'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.name = "mds",
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay5'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, an assay object, or a
  matrix-like object.

- ...:

  Passed to [stats::cmdscale](https://rdrr.io/r/stats/cmdscale.html),
  [MASS::isoMDS](https://rdrr.io/pkg/MASS/man/isoMDS.html) or
  [MASS::sammon](https://rdrr.io/pkg/MASS/man/sammon.html).

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- features:

  Features used instead of a reduction.

- nmds:

  The number of dimensions to be computed.

- dist.method:

  The distance metric to be used. Currently supported values are
  `"euclidean"`, `"chisquared"`, `"kullback"`, `"jeffreys"`, `"jensen"`,
  `"manhattan"`, `"maximum"`, `"canberra"`, `"minkowski"`, and
  `"hamming"`.

- mds.method:

  The MDS algorithm to be used. Currently supported values are
  `"cmdscale"`, `"isoMDS"`, and `"sammon"`.

- rev.mds:

  Whether to perform reverse MDS (i.e., transpose the input matrix)
  before running the analysis.

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the basis vectors.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-13 22:40:10] Start standard processing workflow...
#> ℹ [2026-09-13 22:40:10] Checking a list of <Seurat>...
#> ! [2026-09-13 22:40:10] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-13 22:40:10] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:40:10] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:40:11] Use the separate HVF from `srt_list`
#> ℹ [2026-09-13 22:40:11] Number of available HVF: 2000
#> ℹ [2026-09-13 22:40:11] Finished check
#> ℹ [2026-09-13 22:40:11] Perform `ScaleData()`
#> ℹ [2026-09-13 22:40:11] Perform pca linear dimension reduction
#> ℹ [2026-09-13 22:40:11] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-13 22:40:11] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-13 22:40:12] Reorder clusters...
#> ℹ [2026-09-13 22:40:12] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-13 22:40:12] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-13 22:40:19] Standard processing workflow completed
pancreas_sub <- RunMDS(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "mds"
)
```
