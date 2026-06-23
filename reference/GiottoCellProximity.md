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
