# Run Giotto cell proximity enrichment

Run Giotto cell proximity enrichment

## Usage

``` r
GiottoCellProximity(
  x,
  group.by,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  number_of_simulations = 1000,
  adjust_method = "fdr",
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- group.by:

  Metadata column containing cell or spot groups.

- network_method:

  Spatial network method.

- network_name:

  Name for the Giotto spatial network.

- number_of_simulations:

  Number of label simulations used by Giotto.

- adjust_method:

  Multiple-testing correction method.

- params:

  Additional parameters passed to \`Giotto::cellProximityEnrichment()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.

## Examples

``` r
proximity <- data.frame(
  group_1 = c("Ductal", "Ductal", "Endocrine", "Stromal"),
  group_2 = c("Endocrine", "Stromal", "Stromal", "Ductal"),
  enrichment = c(1.6, 0.8, 1.3, 0.7),
  p.adj = c(0.01, 0.08, 0.03, 0.12)
)
g <- structure(
  list(
    results = list(cell_proximity = list(table = proximity)),
    parameters = list(network_method = "Delaunay", number_of_simulations = 100)
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "cell_proximity")


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
  data(visium_human_pancreas_sub)
  spatial <- subset(
    visium_human_pancreas_sub,
    cells = colnames(visium_human_pancreas_sub)[1:80],
    features = rownames(visium_human_pancreas_sub)[1:200]
  )
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoSpatialNetwork(g)
  g <- GiottoCellProximity(g, group.by = "coda_label", number_of_simulations = 100)
}
```
