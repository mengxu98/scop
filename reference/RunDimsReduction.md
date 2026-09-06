# Run dimension reduction

Run dimension reduction

## Usage

``` r
RunDimsReduction(
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

  A `Seurat` object.

- prefix:

  The prefix used to name the result.

- features:

  Features to plot.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- linear_reduction_dims:

  Total number of dimensions to compute and store for
  `linear_reduction`.

- linear_reduction_params:

  Other parameters passed to the `linear_reduction` method.

- force_linear_reduction:

  Whether force to do linear dimension reduction.

- nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

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

  Whether force to do nonlinear dimension reduction.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Optional integer seed. When supplied, every input receives a
  deterministic independent L'Ecuyer-CMRG random-number stream, making
  results reproducible across worker counts and scheduling order. The
  caller's random number state is restored when the call finishes.

## See also

[DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
