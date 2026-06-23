# Plot spatial variable feature results

Visualize normalized results produced by
[`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md).
The summary view shows feature ranks and significance, while the surface
view reuses
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md)
to draw spatial expression for selected features.

## Usage

``` r
SpatialVariableFeaturePlot(
  srt,
  plot_type = c("summary", "surface", "combined"),
  features = NULL,
  nfeatures = 10,
  score_col = "score",
  assay = NULL,
  layer = NULL,
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
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- plot_type:

  Plot type: `"summary"`, `"surface"`, or `"combined"`.

- features:

  Features to plot. If `NULL`, top features from the stored spatial
  variable feature result are used.

- nfeatures:

  Number of top features used when `features = NULL`.

- score_col:

  Result column used for the summary x-axis.

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

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

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

## Value

A `ggplot` or `patchwork` object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
spatial <- RunSpatialVariableFeatures(
  spatial,
  assay = "Spatial",
  nfeatures = 10,
  verbose = FALSE
)
SpatialVariableFeaturePlot(spatial, plot_type = "summary")
```
