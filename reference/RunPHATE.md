# Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)

Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory
Embedding)

## Usage

``` r
RunPHATE(object, ...)

# S3 method for class 'Seurat'
RunPHATE(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  knn = 5,
  decay = 40,
  n_landmark = 2000,
  t = "auto",
  gamma = 1,
  n_pca = 100,
  knn_dist = "euclidean",
  knn_max = NULL,
  t_max = 100,
  do_cluster = FALSE,
  n_clusters = "auto",
  max_clusters = 100,
  reduction.name = "phate",
  reduction.key = "PHATE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunPHATE(
  object,
  assay = NULL,
  n_components = 2,
  knn = 5,
  decay = 40,
  n_landmark = 2000,
  t = "auto",
  gamma = 1,
  n_pca = 100,
  knn_dist = "euclidean",
  knn_max = NULL,
  t_max = 100,
  do_cluster = FALSE,
  n_clusters = "auto",
  max_clusters = 100,
  reduction.key = "PHATE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, a matrix-like object, a
  Neighbor object, or a Graph object.

- ...:

  Additional arguments to be passed to phate.PHATE.

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

  The number of PHATE components. Default is `2`.

- knn:

  A number of nearest neighbors on which to build kernel. Default is
  `5`.

- decay:

  The sets decay rate of kernel tails. Default is `40`.

- n_landmark:

  A number of landmarks to use in fast PHATE. Default is `2000`.

- t:

  The power to which the diffusion operator is powered. This sets the
  level of diffusion. If `"auto"`, `t` is selected according to the knee
  point in the Von Neumann Entropy of the diffusion operator. Default is
  `"auto"`.

- gamma:

  The informational distance constant between `-1` and `1`. `gamma=1`
  gives the PHATE log potential, `gamma=0` gives a square root
  potential. Default is `1`.

- n_pca:

  A number of principal components to use for calculating neighborhoods.
  For extremely large datasets, using `n_pca < 20` allows neighborhoods
  to be calculated in roughly `log(n_samples)` time. Default is `100`.

- knn_dist:

  The distance metric for k-nearest neighbors. Recommended values:
  `"euclidean"`, `"cosine"`, `"precomputed"`. Default is `"euclidean"`.

- knn_max:

  The maximum number of neighbors for which alpha decaying kernel is
  computed for each point. For very large datasets, setting `knn_max` to
  a small multiple of `knn` can speed up computation significantly.
  Default is `NULL`.

- t_max:

  The maximum `t` to test. Default is `100`.

- do_cluster:

  Whether to perform clustering on the PHATE embeddings. Default is
  `FALSE`.

- n_clusters:

  A number of clusters to be identified. Default is `"auto"`.

- max_clusters:

  The maximum number of clusters to test. Default is `100`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"phate"`.

- reduction.key:

  The prefix for the column names of the PHATE embeddings. Default is
  `"PHATE_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunPHATE(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "phate"
)
} # }
```
