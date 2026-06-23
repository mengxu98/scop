# Create a Giotto spatial network

Create a Giotto spatial network

## Usage

``` r
GiottoSpatialNetwork(
  x,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  params = list(),
  verbose = TRUE
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- network_method:

  Spatial network method.

- network_name:

  Name for the Giotto spatial network.

- params:

  Additional parameters passed to \`Giotto::createSpatialNetwork()\`.

- verbose:

  Whether to print progress messages.

## Value

A \`giotto2\` workflow object.
