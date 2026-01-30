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
#> ℹ [2026-01-30 17:22:41] Start standard scop workflow...
#> ℹ [2026-01-30 17:22:42] Checking a list of <Seurat>...
#> ! [2026-01-30 17:22:42] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:22:42] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:22:44] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:22:45] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:22:45] Number of available HVF: 2000
#> ℹ [2026-01-30 17:22:45] Finished check
#> ℹ [2026-01-30 17:22:46] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:22:46] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:22:47] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:22:47] Reorder clusters...
#> ℹ [2026-01-30 17:22:47] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:22:47] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:22:52] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:22:57] Run scop standard workflow completed
pancreas_sub <- RunMDS(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "mds"
)
```
