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

  Passed to `destiny::DiffusionMap`.

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

- ndcs:

  A number of diffusion components (dimensions) to be computed.

- sigma:

  The diffusion scale parameter of the Gaussian kernel. Currently
  supported values are `"local"` (default) and `"global"`.

- k:

  A number of nearest neighbors to be used for the construction of the
  graph.

- dist.method:

  The distance metric to be used for the construction of the knn graph.
  Currently supported values are `"euclidean"` and `"cosine"`.

- npcs:

  Number of principal components to use for dimensionality reduction
  before computing diffusion map. This can speed up computation when
  using many features. Default is `NULL` (auto-determined based on the
  number of features).

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the basis vectors.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.
