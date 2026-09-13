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
  pt.size = 1.5,
  point.size = NULL,
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

- pt.size:

  Coordinate point size.

- point.size:

  Deprecated alias(es) for \`pt.size\`; supply exactly one of the two.
  It will be removed in scop 1.0.0.

- palette, palcolor:

  Color palette.

- title:

  Optional title.

## Value

A \`ggplot\` object.

## See also

\[RunSpatialCellChat()\]
