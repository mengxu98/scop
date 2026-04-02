# Differential Expression Test Plot

Differential Expression Test Plot

## Usage

``` r
DEtestPlot(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  plot_type = c("volcano", "manhattan", "ring"),
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  x_metric = "diff_pct",
  y_metric = c("p_val_adj", "p_val"),
  x_order = c("gene", "index"),
  palette = "RdBu",
  palcolor = NULL,
  group_palette = "Chinese",
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
  aspect.ratio = NULL,
  xlab = NULL,
  ylab = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  manhattan.bg = "white",
  jitter_width = 0.5,
  jitter_height = 0.4,
  tile_height = 0.3,
  tile_gap = 0.1,
  ring_segments = TRUE,
  seed = 11
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

- plot_type:

  Type of plot to create. Options: `"volcano"`, `"manhattan"`, or
  `"ring"`. Default is `"volcano"`.

- DE_threshold:

  A character string specifying the threshold for differential
  expression (used to highlight significant genes in all plot types).
  Default is `"avg_log2FC > 0 & p_val_adj < 0.05"`.

- x_metric:

  A character string specifying the metric to use for the x-axis (only
  for volcano plot). Default is `"diff_pct"`.

- y_metric:

  A character string specifying the metric to use for the y-axis (only
  for Manhattan plot, not used currently). Options: `"p_val"` or
  `"p_val_adj"`. Default is `"p_val_adj"`.

- x_order:

  A character string specifying how to order genes on x-axis (only for
  Manhattan plot, not used currently). Options: `"gene"` (alphabetical
  by gene name) or `"index"` (by data order). Default is `"gene"`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"RdBu"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- group_palette:

  Palette for cell types (groups) in Manhattan plot. Default is
  `"Chinese"`.

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

- manhattan.bg:

  Background color for Manhattan plot. Default is `"white"`.

- jitter_width:

  Horizontal jitter range for points in Manhattan plot. Default is
  `0.5`.

- jitter_height:

  Vertical jitter range for points in Manhattan plot. Default is `0.4`.

- tile_height:

  Height of the cell-type track in ring plot. Default is `0.3`.

- tile_gap:

  Gap between the track and nudged points in ring plot. Default is
  `0.1`.

- ring_segments:

  Whether to draw segment lines between cell types in ring plot. Default
  is `TRUE`.

- seed:

  Random seed for jitter in ring plot. Default is `11`.

## See also

[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestManhattanPlot](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md),
[DEtestRingPlot](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 15:29:57] Start standard processing workflow...
#> ℹ [2026-04-02 15:29:57] Checking a list of <Seurat>...
#> ! [2026-04-02 15:29:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:29:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:29:59] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:29:59] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:29:59] Number of available HVF: 2000
#> ℹ [2026-04-02 15:30:00] Finished check
#> ℹ [2026-04-02 15:30:00] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:30:00] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:30:05] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:30:05] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:30:05. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-12); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-12 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-04-02 15:30:11] Infinite values detected
#> ! [2026-04-02 15:30:11] Data in the 'data' layer is unknown. Please check the data type
#> ℹ [2026-04-02 15:30:11] Start differential expression test
#> ℹ [2026-04-02 15:30:11] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-04-02 15:30:11] Using 1 core
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> ⠙ [2026-04-02 15:30:11] Running for Ductal [1/5] ■■■■■■■                       …
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> ✔ [2026-04-02 15:30:11] Completed 5 tasks in 40ms
#> 
#> ℹ [2026-04-02 15:30:11] Building results
#> ! [2026-04-02 15:30:11] Found 5 failed results
#> ℹ [2026-04-02 15:30:11] ✖ Error details:
#> ℹ                       ✖ "Ductal": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Ngn3-high-EP": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Endocrine": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Ngn3-low-EP": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Pre-endocrine": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Error in `[.data.frame`(AllMarkers, , "group1"): undefined columns selected

DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "volcano",
  ncol = 2
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first

DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "manhattan"
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first

DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "ring"
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first

de_results1 <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
DEtestPlot(
  res = de_results1,
  plot_type = "volcano",
  ncol = 2
)
#> Error in DEtestPlot(res = de_results1, plot_type = "volcano", ncol = 2): argument "srt" is missing, with no default

de_results2 <- Seurat::FindMarkers(
  pancreas_sub,
  group.by = "CellType",
  ident.1 = "Ductal",
  ident.2 = "Endocrine"
)
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
DEtestPlot(
  res = de_results2,
  plot_type = "volcano"
)
#> Error: object 'de_results2' not found

de_results3 <- Seurat::FindAllMarkers(
  pancreas_sub,
  group.by = "CellType"
)
#> Calculating cluster Ductal
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Calculating cluster Ngn3-high-EP
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Calculating cluster Endocrine
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Calculating cluster Ngn3-low-EP
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Calculating cluster Pre-endocrine
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> Warning: No DE genes identified
#> Warning: The following tests were not performed: 
#> Warning: When testing Ductal versus all:
#>  error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Warning: When testing Ngn3-high-EP versus all:
#>  error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Warning: When testing Endocrine versus all:
#>  error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Warning: When testing Ngn3-low-EP versus all:
#>  error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Warning: When testing Pre-endocrine versus all:
#>  error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
DEtestPlot(
  res = de_results3,
  plot_type = "volcano",
  ncol = 2
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Missing required column: "gene". Please provide gene names in a 'gene'
#> column or as row names.
```
