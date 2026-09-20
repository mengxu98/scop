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

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 21:37:36] Start standard processing workflow...
#> ℹ [2026-09-20 21:37:36] Checking a list of <Seurat>...
#> ! [2026-09-20 21:37:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 21:37:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:37:36] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:37:36] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 21:37:36] Number of available HVF: 2000
#> ℹ [2026-09-20 21:37:36] Finished check
#> ℹ [2026-09-20 21:37:36] Perform `ScaleData()`
#> ℹ [2026-09-20 21:37:36] Perform pca linear dimension reduction
#> ℹ [2026-09-20 21:37:36] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 21:37:37] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 21:37:37] Reorder clusters...
#> ℹ [2026-09-20 21:37:37] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 21:37:37] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 21:37:42] Standard processing workflow completed
pancreas_sub <- RunESTIMATE(pancreas_sub, group.by = "SubCellType")
#> ℹ [2026-09-20 21:37:42] ESTIMATE input: 15997 genes x 8 samples
#> ℹ [2026-09-20 21:37:42] Use 8011 ESTIMATE common genes for scoring
#> ℹ [2026-09-20 21:37:42] ESTIMATE signature overlap: stromal 79, immune 77
EstimateScorePlot(pancreas_sub, group.by = "SubCellType")
```
