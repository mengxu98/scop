# Spatial neighborhood plot

Visualize results produced by
[`RunSpatialNeighborhood()`](https://mengxu98.github.io/scop/reference/RunSpatialNeighborhood.md)
using scop spatial, statistical, and network plotting conventions.

## Usage

``` r
SpatialNeighborhoodPlot(
  srt,
  method = NULL,
  plot_type = c("heatmap", "network", "stat", "spatial"),
  comparison = NULL,
  condition = NULL,
  value = c("estimate", "fraction", "count"),
  FDR_threshold = 0.05,
  top_n = 30,
  layout = c("fr", "nicely", "kk", "circle", "mds"),
  edge_size = c(0.4, 2),
  cols.enriched = "#d7301f",
  cols.depleted = "#2b8cbe",
  cols.ns = "grey75",
  pair = NULL,
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  split.by = NULL,
  palette = "RdBu",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object containing
  [`RunSpatialNeighborhood()`](https://mengxu98.github.io/scop/reference/RunSpatialNeighborhood.md)
  results.

- method:

  Stored neighborhood method to plot. If `NULL`, the active method is
  used.

- plot_type:

  Plot type. One of `"heatmap"`, `"network"`, `"stat"`, or `"spatial"`.

- comparison, condition:

  Optional filters for stored result tables.

- value:

  Column used as the plotted effect value.

- FDR_threshold:

  FDR cutoff used to mark significant pairs.

- top_n:

  Number of pairs to show for network and statistic plots.

- layout:

  Network layout.

- edge_size:

  Network edge size range.

- cols.enriched, cols.depleted, cols.ns:

  Colors for direction categories.

- pair:

  Pair to visualize for `plot_type = "spatial"`, either `"from|to"` or a
  length-2 character vector.

- image, overlay_image, coord.cols, split.by, palette, palcolor,
  legend.position, legend.direction, legend.title, theme_use,
  theme_args, combine, nrow, ncol, byrow, verbose, ...:

  Arguments forwarded to the corresponding plotting backend.

- seed:

  Random seed used by network layout.

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.
