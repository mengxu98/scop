# Plot cell2location spatial results

Plot cell2location spatial results

## Usage

``` r
Cell2locationPlot(
  object,
  plot_type = c("proportion", "abundance", "dominant", "pie"),
  cell_types = NULL,
  prefix = "Cell2location",
  tool_name = "Cell2location",
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...,
  image.scale = c("lowres", "hires"),
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object returned by
  [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md).

- plot_type:

  Result to draw: normalized proportion, q05 absolute abundance,
  dominant cell type, or a spot-level proportion pie.

- cell_types:

  Optional cell types to display for abundance, proportion, or pie
  plots.

- prefix:

  Metadata prefix used by
  [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md).

- tool_name:

  Name of the `srt@tools` result entry.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image:

  Whether to draw the selected spatial image.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- combine:

  Whether to combine matrix or dominant point maps. If `FALSE`, return a
  named list of plots.

- nrow, ncol, byrow:

  Point-map layout controls.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `ggplot`, `patchwork`, or list of plots.

## Details

Abundance plots use posterior q05 values. Proportion plots normalize
these values across cell types within each spot. Matrix point maps share
a default display range. Explicit cutoff or quantile arguments passed
through `...` retain the scales computed by
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md),
including when returning separate panels. Saved plotting matrices retain
every original spot; filtered observations have `NA` values, while
`cells` records only the modeled spots. Older compact results are
expanded only when their modeled IDs and `dropped_spots` record account
for all missing rows. Unexplained missing or unknown IDs are errors.
