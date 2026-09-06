# Run TriMap

Run TriMap

## Usage

``` r
RunTriMap(object, ...)

# S3 method for class 'Seurat'
RunTriMap(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  n_inliers = 12,
  n_outliers = 4,
  n_random = 3,
  distance_method = "euclidean",
  lr = 0.1,
  n_iters = 400,
  apply_pca = TRUE,
  opt_method = "dbd",
  reduction.name = "trimap",
  reduction.key = "TriMap_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
)

# Default S3 method
RunTriMap(
  object,
  assay = NULL,
  n_components = 2,
  n_inliers = 12,
  n_outliers = 4,
  n_random = 3,
  distance_method = "euclidean",
  lr = 0.1,
  n_iters = 400,
  apply_pca = TRUE,
  opt_method = "dbd",
  reduction.key = "TriMap_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
)
```

## Arguments

- object:

  A `Seurat` object, matrix-like object, `Neighbor`, or `Graph`.

- ...:

  Passed to the trimap.TRIMAP function.

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

- n_components:

  A number of TriMap components.

- n_inliers:

  A number of nearest neighbors for forming the nearest neighbor
  triplets.

- n_outliers:

  A number of outliers for forming the nearest neighbor triplets.

- n_random:

  A number of random triplets per point.

- distance_method:

  Distance metric for TriMap. Options are: `"euclidean"`, `"manhattan"`,
  `"angular"`, `"cosine"`, `"hamming"`.

- lr:

  The learning rate for TriMap.

- n_iters:

  A number of iterations for TriMap.

- apply_pca:

  Whether to apply PCA before the nearest-neighbor calculation.

- opt_method:

  Optimization method for TriMap. Options are: `"dbd"`, `"sd"`,
  `"momentum"`.

- reduction.name:

  Name of the reduction to be stored in the Seurat object.

- reduction.key:

  Prefix for the column names of the TriMap embeddings.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

- backend:

  TriMap backend. `"cpp"` uses a compiled triplet sampler and optimizer;
  `"python"` retains the official trimap package.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunTriMap(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "trimap"
)
} # }
```
