# SpatialEcoTyper spatial plot

Plot SpatialEcoTyper labels on spatial coordinates using the spatial
plotting style.

## Usage

``` r
SpatialEcoTyperSpatialPlot(
  srt,
  group.by = "SpatialEcoTyper_SE",
  x.by = "X",
  y.by = "Y",
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c(x.by, y.by),
  palette = "Paired",
  palcolor = NULL,
  theme_use = "theme_scop",
  ...,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object. The same object may be supplied as `object =` for
  consistency with spatial plotting APIs.

- group.by:

  Metadata column containing SpatialEcoTyper labels.

- x.by, y.by:

  Metadata coordinate columns used when no image coordinates are
  available.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image:

  Whether to draw the selected spatial image.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- palette, palcolor:

  Palette passed to `palette_colors()`.

- theme_use:

  Theme name or function.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.

## Examples

``` r
data(visium_human_pancreas_pair_sub)
# The stored `domain` column is a real PRECAST grouping used only for this
# fast spatial plotting example; it is not a new SpatialEcoTyper inference.
SpatialEcoTyperSpatialPlot(
  visium_human_pancreas_pair_sub,
  group.by = "domain",
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  pt.size = 1.5
)
```
