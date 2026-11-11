# Plot positional ATAC region signal

Plot positional ATAC region signal

## Usage

``` r
ATACRegionHeatmap(
  x,
  normalize = TRUE,
  per.base = TRUE,
  groups = NULL,
  log1p = FALSE,
  order.rows = FALSE,
  low = "white",
  high = "#2166AC"
)

ATACRegionPlot(x, normalize = TRUE, per.base = TRUE, groups = NULL)
```

## Arguments

- x:

  Output of
  [`RunATACRegionMatrix()`](https://mengxu98.github.io/scop/reference/RunATACRegionMatrix.md).

- normalize:

  Divide counts by the number of selected cells in each group.

- per.base:

  Divide counts by actual bin widths (including a shorter last bin).
  TRUE makes differently sized bins comparable.

- groups:

  Optional unique group names to plot, in display order.

- log1p:

  Apply log1p to displayed heatmap values only.

- order.rows:

  Order heatmap rows by their total displayed signal across all selected
  groups; the same order is used in every panel.

- low, high:

  Heatmap color endpoints.

## Value

A ggplot with its summarized plotting data. The line plot shows the mean
across supplied regions, not a confidence interval or enrichment test.

## See also

[`RunATACRegionMatrix()`](https://mengxu98.github.io/scop/reference/RunATACRegionMatrix.md)
