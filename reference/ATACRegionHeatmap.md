# Plot positional ATAC region signal

Display raw or explicitly normalized endpoint signal from ATAC regions.

## Usage

``` r
ATACRegionHeatmap(x, normalize = TRUE, per.base = TRUE, groups = NULL,
  log1p = FALSE, order.rows = FALSE, low = "white", high = "#2166AC")
ATACRegionPlot(x, normalize = TRUE, per.base = TRUE, groups = NULL)
```

## Arguments

- x:

  Output of RunATACRegionMatrix.

- normalize:

  Divide by selected cells per group.

- per.base:

  Divide by actual bin widths, including the shorter last bin.

- groups:

  Unique group names in display order; NULL includes all.

- log1p:

  Apply log1p to displayed heatmap values.

- order.rows:

  Order all panels by total signal across selected groups.

- low,high:

  Heatmap color endpoints.

## Value

A ggplot with summarized plotting data. The line plot averages across
regions, not across independent biological replicates; it is not a
confidence interval or enrichment test.

## See also

[`RunATACRegionMatrix`](https://mengxu98.github.io/scop/reference/RunATACRegionMatrix.md)
