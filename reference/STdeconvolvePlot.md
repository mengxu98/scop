# Plot STdeconvolve topic proportions

Plot stored topic proportions without rerunning the optional backend.
Multi-topic point maps use one shared proportion scale and reserve title
space above the automatic layout.

## Usage

``` r
STdeconvolvePlot(
  srt,
  tool_name = "STdeconvolve",
  topics = NULL,
  prefix = NULL,
  plot_type = c("point", "pie"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object.

- tool_name:

  Result key written to `srt@tools` by
  [`RunSTdeconvolve()`](https://mengxu98.github.io/scop/reference/RunSTdeconvolve.md).

- topics:

  Topic names, topic numbers, or metadata columns to plot. If `NULL`,
  all topics in the stored result are used.

- prefix:

  Metadata prefix used by
  [`RunSTdeconvolve()`](https://mengxu98.github.io/scop/reference/RunSTdeconvolve.md).
  If `NULL`, the prefix is read from the stored result.

- plot_type:

  `"point"` or `"pie"`. Pie uses numeric `group.by` columns, `values`,
  or `"<prefix>_prop_*"`/`"<prefix>_frac_*"` when `group.by` is
  `"<prefix>_dominant_type"`.

- combine:

  Whether to combine point plots. If `FALSE`, a named list of plots is
  returned.

- nrow, ncol, byrow:

  Point-plot layout controls. When both `nrow` and `ncol` are `NULL`, a
  near-square layout with at most three columns is used.

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
data(visium_human_pancreas_results_sub)
STdeconvolvePlot(
  visium_human_pancreas_results_sub,
  topics = 1:2,
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
```
