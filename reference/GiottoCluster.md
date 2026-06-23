# Run Giotto nearest-network clustering

Run Giotto nearest-network clustering

## Usage

``` r
GiottoCluster(
  x,
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  network_name = "scop_NN",
  cluster_name = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- method:

  Giotto clustering method.

- dims:

  Dimensions used to build the nearest-neighbor network.

- k:

  Number of nearest neighbors.

- resolution:

  Clustering resolution.

- network_name:

  Name for the Giotto nearest-neighbor network.

- cluster_name:

  Name for the Giotto cluster result.

- params:

  Additional parameters passed to the Giotto clustering function.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
