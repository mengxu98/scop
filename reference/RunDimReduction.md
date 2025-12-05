# Run dimensionality reduction

Run dimensionality reduction

## Usage

``` r
RunDimReduction(
  srt,
  prefix = "",
  features = NULL,
  assay = NULL,
  layer = "data",
  linear_reduction = NULL,
  linear_reduction_dims = 50,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = NULL,
  nonlinear_reduction_dims = 2,
  reduction_use = NULL,
  reduction_dims = NULL,
  graph_use = NULL,
  neighbor_use = NULL,
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- prefix:

  The prefix used to name the result.

- features:

  Use features expression data to run linear or nonlinear dimensionality
  reduction.

- assay:

  Specific assay to get data from.

- layer:

  Specific layer to get data from.

- linear_reduction:

  Method of linear dimensionality reduction. Options are `"pca"`,
  `"ica"`, `"nmf"`, `"mds"`, `"glmpca"`.

- linear_reduction_dims:

  Total number of dimensions to compute and store for
  `linear_reduction`.

- linear_reduction_params:

  Other parameters passed to the `linear_reduction` method.

- force_linear_reduction:

  Whether force to do linear dimensionality reduction.

- nonlinear_reduction:

  Method of nonlinear dimensionality reduction. Options are `"umap"`,
  `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`, `"trimap"`,
  `"largevis"`.

- nonlinear_reduction_dims:

  Total number of dimensions to compute and store for
  `nonlinear_reduction`.

- reduction_use:

  Which dimensional reduction to use as input for `nonlinear_reduction`.

- reduction_dims:

  Which dimensions to use as input for `nonlinear_reduction`, used only
  if `features` is `NULL`.

- graph_use:

  Name of graph to use for the `nonlinear_reduction`.

- neighbor_use:

  Name of neighbor to use for the `nonlinear_reduction`.

- nonlinear_reduction_params:

  Other parameters passed to the `nonlinear_reduction` method.

- force_nonlinear_reduction:

  Whether force to do nonlinear dimensionality reduction.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Set a seed. Default is `11`.

## See also

[DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
