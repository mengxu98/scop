# Spatial spot plot

Spatial spot plot

## Usage

``` r
SpatialSpotPlot(
  srt = NULL,
  group.by = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  values = NULL,
  plot_type = c("point", "pie"),
  plot.data = NULL,
  spot.by = NULL,
  color.by = NULL,
  geom = c("point", "jitter"),
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  crop = TRUE,
  coord.cols = c("col", "row"),
  flip.y = FALSE,
  show_axes = FALSE,
  split.by = NULL,
  cells = NULL,
  show_na = FALSE,
  pt.size = NULL,
  pie.radius = NULL,
  pie.radius.scale = 0.45,
  pt.alpha = 0.9,
  stroke = 0.1,
  jitter_width = 0.25,
  jitter_height = 0.25,
  palette = "Spectral",
  palcolor = NULL,
  bg_color = "grey20",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_spatial",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  verbose = TRUE,
  object = NULL,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object. The same object may be supplied as `object =` for
  consistency with spatial plotting APIs.

- group.by:

  Metadata columns to color spots by.

- features:

  Features to color spots by (from `assay`/`layer`).

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- values:

  Spot-level values (named vector, matrix, or data.frame).

- plot_type:

  `"point"` or `"pie"`. Pie uses numeric `group.by` columns, `values`,
  or `"<prefix>_prop_*"`/`"<prefix>_frac_*"` when `group.by` is
  `"<prefix>_dominant_type"`.

- plot.data, spot.by, color.by:

  Long-format data for repeated spatial points (e.g. cell-to-spot
  assignments).

- geom:

  `"point"` or `"jitter"` for `plot.data`.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image, image.alpha:

  Draw the spatial image beneath spots.

- crop:

  Crop the panel to plotted spots.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- flip.y:

  Reverse the y axis for metadata or FOV coordinates without a raster
  image. The default `FALSE` agrees with spatial networks and cell
  boundaries. Raster-backed display coordinates already include their
  flip.

- show_axes:

  Keep axis text, ticks, and grid when `theme_use` is `"theme_spatial"`.

- split.by:

  Metadata column to facet by.

- cells:

  Cell names to include.

- show_na:

  If `TRUE`, color `NA` groups with `bg_color`; if `FALSE`, drop them.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- pie.radius, pie.radius.scale:

  Pie radius. `NULL` estimates from spot spacing, then multiplies by
  `pie.radius.scale`.

- stroke:

  Point border width.

- jitter_width, jitter_height:

  Jitter when `geom = "jitter"`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- bg_color:

  Point border color.

- legend.position, legend.direction, legend.title:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- theme_use:

  Theme name or function.

- theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork::patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- verbose:

  Whether to print messages.

- object:

  Optional alias for `srt`. Supply exactly one of `srt` or `object`.

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.

## Examples

``` r
data(visium_human_pancreas_sub)
SpatialSpotPlot(
  object = visium_human_pancreas_sub,
  group.by = "coda_label"
)


SpatialSpotPlot(
  object = visium_human_pancreas_sub,
  features = rownames(visium_human_pancreas_sub)[1:2],
  layer = "counts"
)
```
