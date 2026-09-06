# Plot stored SpatialCellChat results

Plot normalized SpatialCellChat results without rerunning the backend.
The network view is a group-aggregated network over group centroids, not
a materialized individual cell-by-cell edge table.

## Usage

``` r
SpatialCellChatPlot(
  object,
  result.name = NULL,
  sample = NULL,
  plot_type = c("spatial_network", "lr_spatial", "pathway", "incoming", "outgoing"),
  signaling = NULL,
  pairLR.use = NULL,
  direction = c("outgoing", "incoming"),
  top_n = 30,
  point.size = 1.5,
  palette = "RdBu",
  palcolor = NULL,
  title = NULL
)
```

## Arguments

- object:

  A \`Seurat\` object with SpatialCellChat results.

- result.name:

  Stored result name.

- sample:

  Spatial sample. Required when multiple samples are stored.

- plot_type:

  Spatial result view.

- signaling:

  Optional pathway name.

- pairLR.use:

  Optional interaction name.

- direction:

  Score direction for spatial score plots.

- top_n:

  Maximum network edges.

- point.size:

  Coordinate point size.

- palette, palcolor:

  Color palette.

- title:

  Optional title.

## Value

A \`ggplot\` object.

## Examples

``` r
data(visium_human_pancreas_results_sub)
SpatialCellChatPlot(
  visium_human_pancreas_results_sub,
  plot_type = "incoming",
  sample = "ALL",
  top_n = 10
)
```
