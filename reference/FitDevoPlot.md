# Plot FitDevo results

Visualize FitDevo developmental potential scores with `scop` dimensional
and grouped summary plots.

## Usage

``` r
FitDevoPlot(
  srt,
  reduction = NULL,
  group.by = NULL,
  score.name = "FitDevo_Score",
  relative.name = "FitDevo_Relative",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
)
```

## Arguments

- srt:

  A `Seurat` object processed by
  [`RunFitDevo()`](https://mengxu98.github.io/scop/reference/RunFitDevo.md).

- reduction:

  Reduction used by
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  and
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md).

- group.by:

  Optional metadata column used for phenotype and score distribution
  plots.

- score.name:

  Metadata column containing the FitDevo score.

- relative.name:

  Metadata column containing the FitDevo relative rank.

- combine:

  Whether to combine plots with `patchwork`.

- nrow, ncol, byrow:

  Layout arguments passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- pt.size, pt.alpha:

  Point size and alpha.

- palette, palcolor:

  Palette arguments for grouped plots.

- theme_use, theme_args:

  Theme arguments passed to `scop` plot helpers.

- ...:

  Additional arguments passed to
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  and
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md).

## Value

A patchwork object or a named list of `ggplot` objects.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:34:18] Start standard processing workflow...
#> ℹ [2026-09-06 21:34:18] Checking a list of <Seurat>...
#> ! [2026-09-06 21:34:18] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:34:18] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:34:18] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:34:19] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:34:19] Number of available HVF: 2000
#> ℹ [2026-09-06 21:34:19] Finished check
#> ℹ [2026-09-06 21:34:19] Perform `ScaleData()`
#> ℹ [2026-09-06 21:34:19] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:34:19] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:34:19] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:34:19] Reorder clusters...
#> ℹ [2026-09-06 21:34:19] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:34:19] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:34:25] Standard processing workflow completed
pancreas_sub <- RunFitDevo(pancreas_sub, verbose = FALSE)
FitDevoPlot(pancreas_sub)
```
