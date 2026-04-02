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

  Additional arguments to be passed to
  [stats::cmdscale](https://rdrr.io/r/stats/cmdscale.html),
  [MASS::isoMDS](https://rdrr.io/pkg/MASS/man/isoMDS.html) or
  [MASS::sammon](https://rdrr.io/pkg/MASS/man/sammon.html).

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `data`.

- features:

  A character vector of features to use. Default is `NULL`.

- nmds:

  The number of dimensions to be computed. Default is `50`.

- dist.method:

  The distance metric to be used. Currently supported values are
  `"euclidean"`, `"chisquared"`, `"kullback"`, `"jeffreys"`, `"jensen"`,
  `"manhattan"`, `"maximum"`, `"canberra"`, `"minkowski"`, and
  `"hamming"`. Default is `"euclidean"`.

- mds.method:

  The MDS algorithm to be used. Currently supported values are
  `"cmdscale"`, `"isoMDS"`, and `"sammon"`. Default is `"cmdscale"`.

- rev.mds:

  Whether to perform reverse MDS (i.e., transpose the input matrix)
  before running the analysis. Default is `FALSE`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"mds"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"MDS_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:44:21] Start standard processing workflow...
#> ℹ [2026-04-02 16:44:22] Checking a list of <Seurat>...
#> ! [2026-04-02 16:44:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:44:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:44:24] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:44:25] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:44:25] Number of available HVF: 2000
#> ℹ [2026-04-02 16:44:25] Finished check
#> ℹ [2026-04-02 16:44:25] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:44:25] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:44:29] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:44:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:44:29] Reorder clusters...
#> ℹ [2026-04-02 16:44:30] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:44:30] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:44:30] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:44:30] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:44:30] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:44:33] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:44:36] Standard processing workflow completed
pancreas_sub <- RunMDS(pancreas_sub)
#> Error in loadNamespace(x): there is no package called ‘proxyC’
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "mds"
)
```
