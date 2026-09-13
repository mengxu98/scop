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

  Assay used for \`Seurat\` or \`SummarizedExperiment\` expression.

- layer:

  Assay layer used for \`Seurat\` expression.

- plot_type:

  Plot type.

- split:

  Split method for high/low groups in violin plots.

- quantile:

  Quantile cutoff used when \`split = "quantile"\`.

- cor_method:

  Correlation method for scatter plots.

- add_stat:

  Whether to add group comparison labels to violin or box plots.
  Requires `ggpubr`.

- ...:

  Additional plotting arguments.

## Value

A \`ggplot\` object.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-13 21:40:34] Start standard processing workflow...
#> ℹ [2026-09-13 21:40:34] Checking a list of <Seurat>...
#> ! [2026-09-13 21:40:34] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-13 21:40:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 21:40:34] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 21:40:34] Use the separate HVF from `srt_list`
#> ℹ [2026-09-13 21:40:34] Number of available HVF: 2000
#> ℹ [2026-09-13 21:40:34] Finished check
#> ℹ [2026-09-13 21:40:34] Perform `ScaleData()`
#> ℹ [2026-09-13 21:40:34] Perform pca linear dimension reduction
#> ℹ [2026-09-13 21:40:34] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-13 21:40:35] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-13 21:40:35] Reorder clusters...
#> ℹ [2026-09-13 21:40:35] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-13 21:40:35] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-13 21:40:41] Standard processing workflow completed
pancreas_sub <- RunESTIMATE(pancreas_sub, group.by = "SubCellType")
#> ℹ [2026-09-13 21:40:41] ESTIMATE input: 15997 genes x 8 samples
#> ℹ [2026-09-13 21:40:41] Use 8011 ESTIMATE common genes for scoring
#> ℹ [2026-09-13 21:40:41] ESTIMATE signature overlap: stromal 79, immune 77
EstimateGenePlot(
  pancreas_sub,
  features = c("GCG", "INS"),
  group.by = "SubCellType"
)
```
