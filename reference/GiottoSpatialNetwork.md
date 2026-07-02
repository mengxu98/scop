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

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
coords <- data.frame(
  cell_ID = colnames(spatial),
  sdimx = spatial$x,
  sdimy = spatial$y,
  row.names = colnames(spatial)
)
edges <- data.frame(
  from = colnames(spatial)[1:79],
  to = colnames(spatial)[2:80]
)
g <- structure(
  list(
    source = list(cells = colnames(spatial), coordinates = coords),
    results = list(spatial_network = list(name = "Delaunay_network", table = edges))
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "network")


if (
  isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE))
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoSpatialNetwork(g, network_method = "Delaunay")
}
#> Error in check_r("giotto-suite/Giotto", verbose = FALSE): could not find function "check_r"
```
