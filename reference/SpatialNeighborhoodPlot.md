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

  A `Seurat` object.

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

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- overlay_image:

  Whether to draw the spatial image beneath spots.

- coord.cols:

  Metadata coordinate columns used when no Seurat image coordinates are
  available.

- split.by:

  Optional metadata column identifying conditions for differential
  neighborhood statistics.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- legend.title:

  Title for the legend. Default is `NULL`, which uses the group name.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

- seed:

  Random seed used by layouts and jittered plot layers.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the selected backend.

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.
