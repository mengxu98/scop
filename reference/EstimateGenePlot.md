# Gene and ESTIMATE score relationship plots

Compare ESTIMATE scores between target-gene high and low expression
groups, or draw continuous gene-expression correlations with ESTIMATE
scores.

## Usage

``` r
EstimateGenePlot(
  object = NULL,
  score.data = NULL,
  gene.data = NULL,
  features,
  assay = NULL,
  layer = "data",
  plot_type = c("violin", "scatter"),
  split = c("median", "quantile"),
  quantile = 0.5,
  cor_method = c("spearman", "pearson", "kendall"),
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

- gene.data:

  Optional gene expression matrix with genes in rows and samples in
  columns.

- features:

  Target genes to plot.

- assay:

  Assay used for `Seurat` or `SummarizedExperiment` expression.

- layer:

  Assay layer used for `Seurat` expression.

- plot_type:

  Plot type.

- split:

  Split method for high/low groups in violin plots.

- quantile:

  Quantile cutoff used when `split = "quantile"`.

- cor_method:

  Correlation method for scatter plots.

- add_stat:

  Whether to add group comparison labels to violin plots.

- ...:

  Additional plotting arguments.

## Value

A `ggplot` object.
