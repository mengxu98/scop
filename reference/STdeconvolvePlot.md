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

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
topic_weights <- data.frame(
  STdeconvolve_prop_topic_1 = seq(0.75, 0.20, length.out = ncol(spatial)),
  STdeconvolve_prop_topic_2 = seq(0.20, 0.70, length.out = ncol(spatial)),
  STdeconvolve_prop_topic_3 = 0.10,
  row.names = colnames(spatial)
)
topic_weights <- topic_weights / rowSums(topic_weights)
spatial <- Seurat::AddMetaData(spatial, topic_weights)
spatial$STdeconvolve_dominant_type <- sub(
  "^STdeconvolve_prop_",
  "",
  colnames(topic_weights)[max.col(topic_weights)]
)

STdeconvolvePlot(
  spatial,
  topics = 1:2,
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)

if (requireNamespace("scatterpie", quietly = TRUE)) {
  STdeconvolvePlot(
    spatial,
    plot_type = "pie",
    overlay_image = FALSE,
    coord.cols = c("x", "y")
  )
}
```
