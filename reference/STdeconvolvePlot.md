# Plot STdeconvolve topic proportions

Plot STdeconvolve topic proportions

## Usage

``` r
STdeconvolvePlot(
  srt,
  topics = NULL,
  prefix = "STdeconvolve",
  plot_type = c("point", "pie"),
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- topics:

  Topic names, topic numbers, or metadata columns to plot. If `NULL`,
  all `"<prefix>_prop_*"` columns are used for point plots.

- prefix:

  Metadata prefix used by
  [`RunSTdeconvolve()`](https://mengxu98.github.io/scop/reference/RunSTdeconvolve.md).

- plot_type:

  Plot type. `"point"` keeps the default spot plot behavior. `"pie"`
  draws spot-level pies from numeric metadata columns supplied to
  `group.by` or from a numeric matrix/data.frame supplied to `values`.
  When `group.by` is a single `"<prefix>_dominant_type"` column,
  matching `"<prefix>_prop_*"` or `"<prefix>_frac_*"` numeric metadata
  columns are used automatically.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.
