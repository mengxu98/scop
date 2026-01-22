# Run diffusion map (DM)

Run diffusion map (DM)

## Usage

``` r
RunDM(object, ...)

# S3 method for class 'Seurat'
RunDM(
  object,
  reduction = "pca",
  dims = 1:30,
  features = NULL,
  assay = NULL,
  layer = "data",
  ndcs = 2,
  sigma = "local",
  k = 30,
  dist.method = "euclidean",
  npcs = NULL,
  reduction.name = "dm",
  reduction.key = "DM_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunDM(
  object,
  assay = NULL,
  layer = "data",
  ndcs = 2,
  sigma = "local",
  k = 30,
  dist.method = "euclidean",
  npcs = NULL,
  reduction.key = "DM_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object or a matrix-like object.

- ...:

  Additional arguments to be passed to `destiny::DiffusionMap`.

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

- ndcs:

  A number of diffusion components (dimensions) to be computed. Default
  is `2`.

- sigma:

  The diffusion scale parameter of the Gaussian kernel. Currently
  supported values are `"local"` (default) and `"global"`.

- k:

  A number of nearest neighbors to be used for the construction of the
  graph. Default is `30`.

- dist.method:

  The distance metric to be used for the construction of the knn graph.
  Currently supported values are `"euclidean"` and `"cosine"`. Default
  is `"euclidean"`.

- npcs:

  Number of principal components to use for dimensionality reduction
  before computing diffusion map. This can speed up computation when
  using many features. Default is `NULL` (auto-determined based on the
  number of features).

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"dm"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"DM_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-22 03:56:42] Start standard scop workflow...
#> ℹ [2026-01-22 03:56:43] Checking a list of <Seurat>...
#> ! [2026-01-22 03:56:43] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 03:56:43] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 03:56:46] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 03:56:46] Use the separate HVF from srt_list
#> ℹ [2026-01-22 03:56:46] Number of available HVF: 2000
#> ℹ [2026-01-22 03:56:46] Finished check
#> ℹ [2026-01-22 03:56:47] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 03:56:47] Perform pca linear dimension reduction
#> ℹ [2026-01-22 03:56:48] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 03:56:48] Reorder clusters...
#> ℹ [2026-01-22 03:56:48] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 03:56:48] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 03:56:53] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 03:56:57] Run scop standard workflow completed
pancreas_sub <- RunDM(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> ◌ [2026-01-22 03:56:57] Running destiny::DiffusionMap
#> ℹ [2026-01-22 03:57:04] Using 50 principal components to speed up computation (provided 2000 features)
#> Error in loadNamespace(name): there is no package called ‘destiny’

CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "dm"
)
```
