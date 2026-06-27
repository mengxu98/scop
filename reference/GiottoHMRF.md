# Run Giotto HMRF spatial domains

Run Giotto HMRF spatial domains

## Usage

``` r
GiottoHMRF(
  x,
  spatial_genes = NULL,
  network_name = "Delaunay_full",
  k = 20,
  betas = c(0, 10, 20),
  hmrf_name = "scop_HMRF",
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- spatial_genes:

  Spatial genes used by Giotto HMRF.

- network_name:

  Name for the Giotto spatial network.

- k:

  Number of HMRF domains.

- betas:

  HMRF beta values.

- hmrf_name:

  Name for the HMRF result.

- params:

  Additional parameters passed to \`Giotto::doHMRF()\`.

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
hmrf_meta <- data.frame(
  cell_ID = colnames(spatial),
  scop_HMRF_k4_b10 = paste0("domain_", as.integer(factor(spatial$coda_label)) %% 4 + 1),
  row.names = colnames(spatial)
)
g <- structure(
  list(
    source = list(cells = colnames(spatial), coordinates = coords),
    results = list(
      hmrf = list(
        table = hmrf_meta["scop_HMRF_k4_b10"],
        metadata = hmrf_meta
      )
    )
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "hmrf")


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoPreprocess(g)
  g <- GiottoSpatialNetwork(g)
  g <- GiottoHMRF(g, spatial_genes = rownames(spatial)[1:30], k = 4)
}
```
