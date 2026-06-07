# Plot spatial gradient screening results

Visualize normalized results produced by
[`RunSpatialGradientFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialGradientFeatures.md)
without requiring the original SPATA2 object.

## Usage

``` r
SpatialGradientPlot(
  srt,
  result_name = NULL,
  plot_type = c("summary", "surface", "line", "model", "combined"),
  features = NULL,
  nfeatures = 4,
  assay = NULL,
  layer = "data",
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  pt.size = NULL,
  pt.alpha = 0.9,
  stroke = 0.1,
  palette = "Spectral",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  line_size = 1,
  line_alpha = 0.35,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- result_name:

  Stored spatial gradient result name. If `NULL`, the latest stored
  result is used.

- plot_type:

  Plot type: `"summary"`, `"surface"`, `"line"`, `"model"`, or
  `"combined"`.

- features:

  Variables to plot. If `NULL`, top variables from the stored result are
  used.

- nfeatures:

  Number of top variables used when `features = NULL`.

- assay:

  Assay used for `features`. If `NULL`, the default assay is used.

- layer:

  Assay layer used for `features`.

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- overlay_image:

  Whether to draw the spatial image beneath spots.

- image.alpha:

  Transparency of the spatial image.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- flip.y:

  Whether to reverse the y axis for metadata coordinates.

- pt.size:

  Point size.

- pt.alpha:

  Point alpha.

- stroke:

  Point border width.

- palette, palcolor:

  Color palette passed to SCOP plotting helpers.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- line_size:

  Size of fitted gradient lines.

- line_alpha:

  Alpha for raw value points.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

## Value

A `ggplot` or `patchwork` object.
