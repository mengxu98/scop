# Plot stored CellOracle results

Display official CellOracle displacements, grid flow, predicted
expression changes or a cluster-specific regulatory network without
rerunning Python.

## Usage

``` r
CellOraclePlot(
  object,
  tool_name = "CellOracle",
  plot_type = c("quiver", "flow", "effect", "network"),
  randomized = FALSE,
  arrow_scale = 1,
  pt.size = 0.5,
  pt.alpha = 0.4,
  top_n = 20L,
  cluster = NULL,
  theme_use = "theme_scop",
  theme_args = list()
)
```

## Arguments

- object:

  A Seurat object returned by
  [`RunCellOracle()`](https://mengxu98.github.io/scop/reference/RunCellOracle.md).

- tool_name:

  Result name in `object@tools`.

- plot_type:

  One of `"quiver"`, `"flow"`, `"effect"`, or `"network"`.

- randomized:

  Show the official randomized displacement reference in quiver/flow
  plots. This is not an independent biological replicate.

- arrow_scale:

  Positive multiplier of stored displacements in embedding units. This
  is a display parameter, not matplotlib's inverse `scale` argument.

- pt.size:

  Background cell point size.

- pt.alpha:

  Background cell opacity.

- top_n:

  Number of genes in effect plots or edges in network plots.

- cluster:

  Cluster name required for network plots when multiple networks are
  stored. Networks are never silently pooled across clusters.

- theme_use:

  Theme function or function name.

- theme_args:

  Named list of arguments passed to the theme function.

## Value

A ggplot object. Flow plots apply the stored upstream `mass_filter`
(`TRUE` means excluded), without smoothing or resampling. Effect bars
display mean absolute `delta_X` across cells, not fold changes or
p-values. Network arrow widths reflect absolute upstream coefficients,
not activation/repression.

## See also

[`RunCellOracle()`](https://mengxu98.github.io/scop/reference/RunCellOracle.md)
