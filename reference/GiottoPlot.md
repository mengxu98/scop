# Plot Giotto backend results

Plot standalone Giotto backend results with scop plotting conventions.
The input `Seurat` object, when supplied, is copied internally for
plotting and is not modified.

## Usage

``` r
GiottoPlot(x, ...)

# S3 method for class 'giotto2_cluster'
GiottoPlot(
  x,
  srt,
  image = x$parameters$image %||% NULL,
  coord.cols = x$parameters$coord.cols %||% c("x", "y"),
  overlay_image = TRUE,
  crop = TRUE,
  pt.size = NULL,
  pt.alpha = 0.95,
  stroke = 0.08,
  palette = "Chinese",
  feature_palette = "Spectral",
  bg_color = "grey25",
  legend.position = "right",
  theme_use = "theme_blank",
  theme_args = list(),
  title = "Giotto Leiden clusters",
  subtitle = NULL,
  ...
)

# S3 method for class 'giotto2_cell_proximity'
GiottoPlot(
  x,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  title = "Giotto cell proximity enrichment",
  subtitle = NULL,
  ...
)

# S3 method for class 'giotto2_spatial_genes'
GiottoPlot(
  x,
  srt = NULL,
  plot_type = c("ranking", "feature"),
  feature = NULL,
  top_n = 20,
  assay = x$parameters$assay %||% NULL,
  layer = x$parameters$layer %||% "data",
  image = x$parameters$image %||% NULL,
  coord.cols = x$parameters$coord.cols %||% c("x", "y"),
  overlay_image = TRUE,
  crop = TRUE,
  pt.size = NULL,
  pt.alpha = 0.95,
  palette = "Chinese",
  feature_palette = "Spectral",
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  title = NULL,
  subtitle = NULL,
  ...
)

# S3 method for class 'giotto2_spatial_modules'
GiottoPlot(
  x,
  features = NULL,
  top_n = 20,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  title = "Giotto spatial co-expression",
  subtitle = "Spatial correlation among top features",
  ...
)
```

## Arguments

- x:

  A result returned by one of the `RunGiotto*()` functions.

- ...:

  Arguments passed to S3 methods.

- srt:

  Original \`Seurat\` object used to create the Giotto result. Required
  for spatial spot plots.

- image:

  Name of the Seurat spatial image. If \`NULL\`, the first image is used
  when available.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- overlay_image:

  Whether to draw the spatial image beneath spots.

- crop:

  Whether to crop spatial panels to plotted spots.

- pt.size:

  Point size for spatial plots.

- pt.alpha:

  Point alpha for spatial plots.

- stroke:

  Point border width for discrete spatial plots.

- palette:

  Discrete palette used for groups.

- feature_palette:

  Continuous palette used for spatial expression plots.

- bg_color:

  Point border color for discrete spatial plots.

- legend.position:

  Legend position.

- theme_use:

  Theme function name used by scop plots.

- theme_args:

  Additional arguments passed to \`theme_use\`.

- title, subtitle:

  Plot title and subtitle. If \`NULL\`, sensible defaults are used.

- heatmap_palette:

  Continuous palette used for heatmaps.

- heatmap_palcolor:

  Optional custom colors used to create \`heatmap_palette\`.

- plot_type:

  Plot type for spatial gene results. \`"ranking"\` plots the
  feature-level table; \`"feature"\` plots expression of one feature on
  spatial coordinates.

- feature:

  Feature to draw for \`plot_type = "feature"\`. If \`NULL\`, the top
  Giotto feature is used.

- top_n:

  Number of rows shown in ranking plots.

- assay:

  Assay used for spatial feature expression plots.

- layer:

  Assay layer used for spatial feature expression plots.

- features:

  Features used for spatial co-expression heatmaps. If \`NULL\`, top
  features from the Giotto result are used.

## Value

A `ggplot` object.
