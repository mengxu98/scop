# Run PaCMAP

Run PaCMAP

## Usage

``` r
RunPaCMAP(object, ...)

# S3 method for class 'Seurat'
RunPaCMAP(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  n.neighbors = NULL,
  MN_ratio = 0.5,
  FP_ratio = 2,
  distance_method = "euclidean",
  lr = 1,
  num_iters = 450L,
  apply_pca = TRUE,
  init = "random",
  reduction.name = "pacmap",
  reduction.key = "PaCMAP_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
)

# Default S3 method
RunPaCMAP(
  object,
  assay = NULL,
  n_components = 2,
  n.neighbors = NULL,
  MN_ratio = 0.5,
  FP_ratio = 2,
  distance_method = "euclidean",
  lr = 1,
  num_iters = 450L,
  apply_pca = TRUE,
  init = "random",
  reduction.key = "PaCMAP_",
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

  Passed to pacmap.PaCMAP.

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

  The number of PaCMAP components.

- n.neighbors:

  A number of neighbors considered in the k-Nearest Neighbor graph.
  Default is `10` for dataset whose sample size is smaller than 10000.
  For large dataset whose sample size (n) is larger than 10000, the
  default value is: `10 + 15 * (log10(n) - 4)`.

- MN_ratio:

  The ratio of the ratio of the number of mid-near pairs to the number
  of neighbors.

- FP_ratio:

  The ratio of the ratio of the number of further pairs to the number of
  neighbors.

- distance_method:

  The distance metric to be used.

- lr:

  The learning rate of the Adam optimizer.

- num_iters:

  The number of iterations for PaCMAP optimization.

- apply_pca:

  Whether pacmap should apply PCA to the data before constructing the
  k-Nearest Neighbor graph. Using PCA to preprocess the data can largely
  accelerate the DR process without losing too much accuracy. Notice
  that this option does not affect the initialization of the
  optimization process.

- init:

  The initialization of the lower dimensional embedding. One of `"pca"`
  or `"random"`.

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the PaCMAP embeddings.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

- backend:

  PaCMAP backend. `"cpp"` uses a compiled pair sampler and Adam
  optimizer; `"python"` retains the official pacmap package.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunPaCMAP(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "pacmap"
)
} # }
```
