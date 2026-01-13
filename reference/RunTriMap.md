# Run TriMap (Large-scale Dimensionality Reduction Using Triplets)

Run TriMap (Large-scale Dimensionality Reduction Using Triplets)

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
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, a matrix-like object, a
  Neighbor object, or a Graph object.

- ...:

  Additional arguments to be passed to the trimap.TRIMAP function.

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

- n_components:

  A number of TriMap components. Default is `2`.

- n_inliers:

  A number of nearest neighbors for forming the nearest neighbor
  triplets. Default is `12`.

- n_outliers:

  A number of outliers for forming the nearest neighbor triplets.
  Default is `4`.

- n_random:

  A number of random triplets per point. Default is `3`.

- distance_method:

  A character string specifying the distance metric for TriMap. Options
  are: `"euclidean"`, `"manhattan"`, `"angular"`, `"cosine"`,
  `"hamming"`. Default is `"euclidean"`.

- lr:

  The learning rate for TriMap. Default is `0.1`.

- n_iters:

  A number of iterations for TriMap. Default is `400`.

- apply_pca:

  Whether to apply PCA before the nearest-neighbor calculation. Default
  is `TRUE`.

- opt_method:

  A character string specifying the optimization method for TriMap.
  Options are: `"dbd"`, `"sd"`, `"momentum"`. Default is `"dbd"`.

- reduction.name:

  A character string specifying the name of the reduction to be stored
  in the Seurat object. Default is `"trimap"`.

- reduction.key:

  A character string specifying the prefix for the column names of the
  TriMap embeddings. Default is `"TriMap_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
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
