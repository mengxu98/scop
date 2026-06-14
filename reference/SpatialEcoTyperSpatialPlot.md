# SpatialEcoTyper spatial plot

Plot SpatialEcoTyper labels on spatial coordinates using scop's spatial
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
  palette = "Spectral",
  palcolor = NULL,
  theme_use = "theme_scop",
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Metadata column containing SpatialEcoTyper labels.

- x.by, y.by:

  Metadata coordinate columns used when no image coordinates are
  available.

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- overlay_image:

  Whether to draw the spatial image beneath spots.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- palette, palcolor:

  Palette passed to `palette_colors()`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.
