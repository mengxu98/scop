# Plot spatial gradient screening results

Visualize normalized results produced by
[`RunSpatialGradientFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialGradientFeatures.md)

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
  line_fit = c("stored", "lm"),
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object.

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

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image, image.alpha:

  Draw the spatial image beneath spots.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- flip.y:

  Reverse the y axis for metadata coordinates.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- stroke:

  Point border width.

- palette, palcolor:

  Color palette passed to SCOP plotting helpers.

- legend.position:

  Legend position for surface, line, and model plots.

- theme_use:

  Theme name or function.

- theme_args:

  Theme name or function, plus extra theme arguments.

- line_size:

  Size of fitted gradient lines.

- line_alpha:

  Alpha for raw value points.

- line_fit:

  Gradient line source. `"stored"` uses the saved `screening$estimate`
  values produced by the selected backend. `"lm"` draws a fresh linear
  fit from `screening$value`, which is useful for showing a simple
  monotonic trend even when the backend stores a smoothed curve.

- nrow, ncol, byrow:

  Layout controls for multi-feature plots.

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

## Value

A `ggplot` or `patchwork` object.

## Examples

``` r
data(visium_human_pancreas_results_sub)
SpatialGradientPlot(
  visium_human_pancreas_results_sub,
  result_name = "scop_gradient_fixture",
  plot_type = "surface",
  features = rownames(visium_human_pancreas_results_sub)[1:2],
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  pt.size = 1.2
)
```
