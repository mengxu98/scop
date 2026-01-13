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
pancreas_sub <- RunDM(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> Error in loadNamespace(name): there is no package called ‘destiny’

CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "dm"
)
```
