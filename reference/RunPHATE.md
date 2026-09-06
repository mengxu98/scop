# Run PHATE

Run PHATE

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
  backend = c("python", "cpp"),
  mds = "metric",
  mds_solver = "sgd",
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
  backend = c("python", "cpp"),
  mds = "metric",
  mds_solver = "sgd",
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

  A `Seurat` object, matrix-like object, `Neighbor`, or `Graph`.

- ...:

  Passed to phate.PHATE.

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

  The number of PHATE components.

- knn:

  A number of nearest neighbors on which to build kernel.

- decay:

  The sets decay rate of kernel tails.

- n_landmark:

  A number of landmarks to use in fast PHATE.

- t:

  The power to which the diffusion operator is powered. This sets the
  level of diffusion. If `"auto"`, `t` is selected according to the knee
  point in the Von Neumann Entropy of the diffusion operator.

- gamma:

  The informational distance constant between `-1` and `1`. `gamma=1`
  gives the PHATE log potential, `gamma=0` gives a square root
  potential.

- n_pca:

  A number of principal components to use for calculating neighborhoods.
  For extremely large datasets, using `n_pca < 20` allows neighborhoods
  to be calculated in roughly `log(n_samples)` time.

- knn_dist:

  The distance metric for k-nearest neighbors. Recommended values:
  `"euclidean"`, `"cosine"`, `"precomputed"`.

- knn_max:

  The maximum number of neighbors for which alpha decaying kernel is
  computed for each point. For very large datasets, setting `knn_max` to
  a small multiple of `knn` can speed up computation significantly.

- t_max:

  The maximum `t` to test.

- do_cluster:

  Whether to perform clustering on the PHATE embeddings.

- backend:

  PHATE backend. `"python"` calls the upstream `phate` Python package
  and `"cpp"` uses the native C++ helper path. The selected backend is
  stable; native arguments never trigger an implicit Python workflow.

- mds:

  MDS algorithm passed to PHATE. The native backend uses its C++
  classical MDS or R's native metric MDS implementation.

- mds_solver:

  Metric MDS solver; recorded by the native backend and passed through
  to Python when that backend is selected.

- n_clusters:

  A number of clusters to be identified.

- max_clusters:

  The maximum number of clusters to test.

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the PHATE embeddings.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
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
