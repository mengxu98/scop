# DEtest Manhattan Plot

Draw a Manhattan-style plot of differential expression results by cell
type.

## Usage

``` r
DEtestManhattanPlot(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  group_palette = "Paired",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  nlabel = 5,
  features_label = NULL,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  manhattan.bg = "white",
  jitter_width = 0.5,
  jitter_height = 0.4,
  aspect.ratio = NULL,
  xlab = NULL,
  ylab = NULL
)
```

## Arguments

- srt:

  An object of class `Seurat` containing the results of differential
  expression analysis.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- test.use:

  A character string specifying the type of statistical test to use.
  Default is `"wilcox"`.

- res:

  A `data.frame` or `data.table` with differential expression results.
  When `res` is provided, `srt` will be ignored. The data.frame must
  contain columns: `gene`, `group1` (factor or character), `avg_log2FC`,
  `p_val_adj`, and optionally `pct.1` and `pct.2` for calculating
  `diff_pct`.

- DE_threshold:

  A character string specifying the threshold for differential
  expression (used to highlight significant genes in all plot types).
  Default is `"avg_log2FC > 0 & p_val_adj < 0.05"`.

- group_palette:

  Palette for cell types (groups) in Manhattan plot. Default is
  `"Paired"`.

- group_palcolor:

  Custom colors for cell types (groups) in Manhattan plot. Default is
  `NULL`.

- pt.size:

  The size of the points. Default is `1`.

- pt.alpha:

  The transparency of the data points. Default is `1`.

- cols.highlight:

  A character string specifying the color for highlighted points.
  Default is `"black"`.

- sizes.highlight:

  The size of the highlighted points. Default is `1`.

- alpha.highlight:

  The transparency of the highlighted points. Default is `1`.

- stroke.highlight:

  The stroke width for the highlighted points. Default is `0.5`.

- nlabel:

  An integer value specifying the number of labeled points per group.
  Default is `5`.

- features_label:

  A character vector specifying the feature labels to plot. Default is
  `NULL`.

- label.fg:

  A character string specifying the color for the labels' foreground.
  Default is `"black"`.

- label.bg:

  A character string specifying the color for the labels' background.
  Default is `"white"`.

- label.bg.r:

  The radius of the rounding of the labels' background. Default is
  `0.1`.

- label.size:

  The size of the labels. Default is `4`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"RdBu"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- theme_use:

  Theme to use for the plot. Default is `"theme_scop"`.

- theme_args:

  A list of additional arguments to pass to the theme function. Default
  is [`list()`](https://rdrr.io/r/base/list.html).

- manhattan.bg:

  Background color for Manhattan plot. Default is `"white"`.

- jitter_width:

  Horizontal jitter range for points in Manhattan plot. Default is
  `0.5`.

- jitter_height:

  Vertical jitter range for points in Manhattan plot. Default is `0.4`.

- aspect.ratio:

  Aspect ratio of the panel. Default is `NULL`.

- xlab:

  A character string specifying the x-axis label.

- ylab:

  A character string specifying the y-axis label.

## See also

[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestRingPlot](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-02-27 15:26:21] Start standard scop workflow...
#> ℹ [2026-02-27 15:26:22] Checking a list of <Seurat>...
#> ! [2026-02-27 15:26:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-02-27 15:26:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-02-27 15:26:24] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-02-27 15:26:24] Use the separate HVF from srt_list
#> ℹ [2026-02-27 15:26:24] Number of available HVF: 2000
#> ℹ [2026-02-27 15:26:25] Finished check
#> ℹ [2026-02-27 15:26:25] Perform `Seurat::ScaleData()`
#> ℹ [2026-02-27 15:26:25] Perform pca linear dimension reduction
#> ℹ [2026-02-27 15:26:26] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-02-27 15:26:26] Reorder clusters...
#> ℹ [2026-02-27 15:26:26] Perform umap nonlinear dimension reduction
#> ℹ [2026-02-27 15:26:26] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-02-27 15:26:29] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-02-27 15:26:32] Run scop standard workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> ℹ [2026-02-27 15:26:33] Data type is log-normalized
#> ℹ [2026-02-27 15:26:33] Start differential expression test
#> ℹ [2026-02-27 15:26:33] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-02-27 15:26:33] Using 1 core
#> ⠙ [2026-02-27 15:26:33] Running for Ductal [1/5] ■■■■■■■                       …
#> ✔ [2026-02-27 15:26:33] Completed 5 tasks in 1.1s
#> 
#> ℹ [2026-02-27 15:26:33] Building results
#> ✔ [2026-02-27 15:26:34] Differential expression test completed
DEtestManhattanPlot(
  pancreas_sub,
  group.by = "CellType"
)
```
