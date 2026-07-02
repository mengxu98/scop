# Run Giotto cell proximity enrichment

Use Giotto as a temporary backend to build a spatial network and test
pairwise enrichment between metadata groups. The complete Giotto object
and result tables are returned as a standalone result; the input
`Seurat` object is not modified.

## Usage

``` r
RunGiottoCellProximity(
  srt,
  group.by,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  number_of_simulations = 1000,
  adjust_method = "fdr",
  tool_name = "GiottoCellProximity",
  store_giotto = TRUE,
  conversion_params = list(),
  network_params = list(),
  enrichment_params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Seurat metadata column containing cell or spot groups.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the expression matrix.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- network_method:

  Spatial network method passed to `Giotto::createSpatialNetwork()`.

- network_name:

  Name for the Giotto spatial network.

- number_of_simulations:

  Number of label simulations used by Giotto.

- adjust_method:

  Multiple-testing correction method.

- tool_name:

  Result name recorded in returned parameters. This function does not
  write to `srt@tools`.

- store_giotto:

  Deprecated compatibility argument. The complete Giotto object is
  always returned in the `giotto` element.

- conversion_params:

  Additional parameters passed to `Giotto::createGiottoObject()`.

- network_params:

  Additional parameters passed to `Giotto::createSpatialNetwork()`.

- enrichment_params:

  Additional parameters passed to `Giotto::cellProximityEnrichment()`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A `giotto2_result` list containing the full Giotto object, enrichment
table, raw Giotto result, parameters, features, and cells.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial$region <- ifelse(
  spatial$x > stats::median(spatial$x),
  "right",
  "left"
)
proximity <- list(
  enrichment = data.frame(
    group_1 = c("left", "left", "right", "right"),
    group_2 = c("left", "right", "left", "right"),
    enrichment = c(1.1, -0.7, -0.5, 1.3),
    type_int = c("enriched", "depleted", "depleted", "enriched")
  ),
  parameters = list(network_method = "Delaunay", number_of_simulations = 100)
)
class(proximity) <- c("giotto2_cell_proximity", "giotto2_result", "list")

head(proximity$enrichment)
#>   group_1 group_2 enrichment type_int
#> 1    left    left        1.1 enriched
#> 2    left   right       -0.7 depleted
#> 3   right    left       -0.5 depleted
#> 4   right   right        1.3 enriched
GiottoPlot(proximity)


if (
  isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE))
) {
spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
proximity <- RunGiottoCellProximity(
  spatial,
  group.by = "region",
  assay = "Spatial",
  layer = "data",
  coord.cols = c("x", "y"),
  network_method = "Delaunay",
  number_of_simulations = 100
)
}
#> Error in check_r("giotto-suite/Giotto", verbose = FALSE): could not find function "check_r"
```
