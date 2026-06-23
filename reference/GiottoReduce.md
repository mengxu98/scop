# Run Giotto dimensional reduction

Run Giotto dimensional reduction

## Usage

``` r
GiottoReduce(
  x,
  reduction = c("pca", "umap"),
  dims = 1:20,
  name = NULL,
  features = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- reduction:

  Dimensional reduction to run.

- dims:

  Dimensions to use.

- name:

  Name for the Giotto reduction.

- features:

  Features used for the reduction.

- params:

  Additional parameters passed to the Giotto reduction function.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
