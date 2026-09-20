# Plot stored COMMOT results

Plot stored COMMOT group communication, cluster matrices, or
raw-coordinate direction vectors without rerunning Python.

## Usage

``` r
COMMOTPlot(
  object,
  result.name = NULL,
  plot_type = c("network", "matrix", "direction"),
  key = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
)
```

## Arguments

- object:

  A \`Seurat\` object with COMMOT results.

- result.name:

  Stored COMMOT result name.

- plot_type:

  Network, cluster matrix, or direction-vector view.

- key:

  Stored cluster or direction selection key.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments. Applied to every
  plot type, including the network view.

- ...:

  Arguments passed to the SCOP network plot for \`"network"\`.

## Value

A \`ggplot\` or compatible plot object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- RunCOMMOT(
  visium_human_pancreas_sub,
  group.by = "CellType",
  coord.cols = c("col", "row"),
  cluster = TRUE,
  direction = TRUE,
  backend = "r"
)
COMMOTPlot(spatial, plot_type = "matrix")
COMMOTPlot(spatial, plot_type = "direction")
COMMOTPlot(spatial, plot_type = "network")
} # }
```
