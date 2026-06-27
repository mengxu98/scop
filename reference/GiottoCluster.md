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

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:80],
  features = rownames(visium_human_pancreas_sub)[1:200]
)
#> Warning: Not validating Centroids objects
#> Warning: Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating Seurat objects
coords <- data.frame(
  cell_ID = colnames(spatial),
  sdimx = spatial$x,
  sdimy = spatial$y,
  row.names = colnames(spatial)
)
g <- structure(
  list(
    source = list(cells = colnames(spatial), coordinates = coords),
    results = list(
      cluster = list(
        table = data.frame(
          cell = colnames(spatial),
          cluster = spatial$coda_label,
          row.names = colnames(spatial)
        )
      )
    ),
    parameters = list(k = 8, resolution = 0.4)
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "cluster")


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoPreprocess(g)
  g <- GiottoReduce(g, reduction = "pca", dims = 1:10)
  g <- GiottoCluster(g, dims = 1:10, k = 8, resolution = 0.4)
}
```
