# Run Giotto spatial co-expression modules

Run Giotto spatial co-expression modules

## Usage

``` r
GiottoSpatialModules(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  detect_params = list(),
  cluster_params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- features:

  Features to test.

- network_method:

  Spatial network method.

- network_name:

  Name for the Giotto spatial network.

- cor_method:

  Correlation method used by Giotto.

- k:

  Number of spatial co-expression modules.

- detect_params:

  Additional parameters passed to \`Giotto::detectSpatialCorFeats()\`.

- cluster_params:

  Additional parameters passed to \`Giotto::clusterSpatialCorFeats()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
