# Preprocess an internal Giotto workflow object

Preprocess an internal Giotto workflow object

## Usage

``` r
GiottoPreprocess(
  x,
  filter_params = list(),
  norm_params = list(),
  stat_params = list(),
  hvf_params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- filter_params:

  Additional parameters reserved for future filtering.

- norm_params:

  Additional parameters passed to \`Giotto::normalizeGiotto()\`.

- stat_params:

  Additional parameters passed to \`Giotto::addStatistics()\`.

- hvf_params:

  Additional parameters passed to \`Giotto::calculateHVF()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
g <- structure(
  list(
    source = list(
      cells = colnames(spatial),
      features = rownames(spatial)[1:100],
      coordinates = data.frame(cell_ID = colnames(spatial), sdimx = spatial$x, sdimy = spatial$y)
    ),
    results = list(),
    active = NULL
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "spatial")


if (
  isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE))
) {
g <- SeuratToScopGiotto(spatial, assay = "Spatial", coord.cols = c("x", "y"), verbose = FALSE)
g <- GiottoPreprocess(g, verbose = FALSE)
}
#> Error in check_r("giotto-suite/Giotto", verbose = FALSE): could not find function "check_r"
```
