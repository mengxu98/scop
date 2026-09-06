# Plot stored SpatialDM results

Plot stored SpatialDM results

## Usage

``` r
SpatialDMPlot(
  object,
  result.name = NULL,
  plot_type = c("weights", "global", "local"),
  pair = NULL,
  spot = NULL,
  signaling = c("secreted", "contact"),
  highlight = NULL,
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- object:

  A Seurat object with SpatialDM results.

- result.name:

  Stored result name.

- plot_type:

  Plot type (\`"weights"\`, \`"global"\`, or \`"local"\`).

- pair:

  LR interaction name for local plotting.

- spot:

  Spot identifier for weight plotting.

- signaling:

  Which stored weight matrix to use.

- highlight:

  LR interactions to label in the global plot.

- palette, palcolor, theme_use, theme_args:

  SCOP plotting controls.

- ...:

  Additional arguments passed to SCOP spatial plotting functions.

- image.scale:

  Image scale factor matching the selected raster.

## Value

A \`ggplot\` or \`patchwork\` object.
