# Volcano Plot

Generate a volcano plot based on differential expression analysis
results.

## Usage

``` r
VolcanoPlot(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  x_metric = "diff_pct",
  palette = "RdBu",
  palcolor = NULL,
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
  aspect.ratio = NULL,
  xlab = x_metric,
  ylab = "-log10(p-adjust)",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
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

- x_metric:

  A character string specifying the metric to use for the x-axis (only
  for volcano plot). Default is `"diff_pct"`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"RdBu"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

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

- aspect.ratio:

  Aspect ratio of the panel. Default is `NULL`.

- xlab:

  A character string specifying the x-axis label.

- ylab:

  A character string specifying the y-axis label.

- theme_use:

  Theme to use for the plot. Default is `"theme_scop"`.

- theme_args:

  A list of additional arguments to pass to the theme function. Default
  is [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Whether to combine multiple plots into one. Default is `TRUE`.

- nrow:

  Number of rows for combined plots. Default is `NULL`.

- ncol:

  Number of columns for combined plots. Default is `NULL`.

- byrow:

  Whether to fill plots by row. Default is `TRUE`.

## See also

[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[DEtestManhattanPlot](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md),
[DEtestRingPlot](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:56:32] Start standard processing workflow...
#> ℹ [2026-04-02 16:56:33] Checking a list of <Seurat>...
#> ! [2026-04-02 16:56:33] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:56:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:56:34] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:56:35] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:56:35] Number of available HVF: 2000
#> ℹ [2026-04-02 16:56:35] Finished check
#> ℹ [2026-04-02 16:56:35] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:56:36] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:56:40] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:56:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:56:40] Reorder clusters...
#> ℹ [2026-04-02 16:56:40] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:56:40] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:56:40] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:56:40] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:56:40] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:56:43] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:56:46] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-04-02 16:56:51] Data type is log-normalized
#> ℹ [2026-04-02 16:56:51] Start differential expression test
#> ℹ [2026-04-02 16:56:51] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-04-02 16:56:51] Using 1 core
#> ⠙ [2026-04-02 16:56:51] Running for Ductal [1/5] ■■■■■■■                       …
#> ✔ [2026-04-02 16:56:51] Completed 5 tasks in 676ms
#> 
#> ℹ [2026-04-02 16:56:51] Building results
#> ! [2026-04-02 16:56:51] Found 5 failed results
#> ℹ [2026-04-02 16:56:52] ✖ Error details:
#> ℹ                       ✖ "Ductal": The total size of the 3 globals exported for future expression (‘FUN()’) is 556.27 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.94 MiB of class ‘function’), ‘data.use’ (4.32 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Ngn3-high-EP": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.18 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.54 MiB of class ‘function’), ‘data.use’ (2.63 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Endocrine": The total size of the 3 globals exported for future expression (‘FUN()’) is 556.21 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.90 MiB of class ‘function’), ‘data.use’ (4.30 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Ngn3-low-EP": The total size of the 3 globals exported for future expression (‘FUN()’) is 553.31 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.99 MiB of class ‘function’), ‘data.use’ (3.33 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Pre-endocrine": The total size of the 3 globals exported for future expression (‘FUN()’) is 554.90 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.03 MiB of class ‘function’), ‘data.use’ (3.87 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> Error in `[.data.frame`(AllMarkers, , "group1"): undefined columns selected
VolcanoPlot(
  pancreas_sub,
  group.by = "CellType",
  ncol = 2
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first

VolcanoPlot(
  pancreas_sub,
  group.by = "CellType",
  DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05",
  ncol = 2
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first

VolcanoPlot(
  pancreas_sub,
  group.by = "CellType",
  x_metric = "avg_log2FC",
  ncol = 2
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first
```
