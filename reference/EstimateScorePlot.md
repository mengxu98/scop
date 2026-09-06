# ESTIMATE score plots

Visualize ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity
scores from a
[`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
result.

## Usage

``` r
EstimateScorePlot(
  object = NULL,
  score.data = NULL,
  group.by = NULL,
  group.data = NULL,
  plot_type = c("violin", "box", "heatmap", "cor"),
  scores = c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
  add_stat = TRUE,
  ...
)
```

## Arguments

- object:

  Optional
  [`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
  bundle, `SummarizedExperiment`, or `Seurat` object containing ESTIMATE
  results.

- score.data:

  Optional score matrix or data frame with samples in rows.

- group.by:

  Optional grouping column for grouped violin and box plots.

- group.data:

  Optional named vector or data frame containing sample groups.

- plot_type:

  Plot type.

- scores:

  ESTIMATE score columns to plot.

- add_stat:

  Whether to add group comparison labels to violin or box plots.
  Requires `ggpubr`.

- ...:

  Additional plotting arguments.

## Value

A `ggplot` object.
