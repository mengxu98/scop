# Run Giotto spatial gene detection

Run Giotto spatial gene detection

## Usage

``` r
GiottoSpatialGenes(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  top_n = 100,
  params = list(),
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

- bin_method:

  Binarization method passed to \`Giotto::binSpect()\`.

- top_n:

  Number of top spatial genes to store.

- params:

  Additional parameters passed to \`Giotto::binSpect()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
