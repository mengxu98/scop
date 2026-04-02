# DEtest Ring Plot

Draw a circular (ring) plot of differential expression results by cell
type.

## Usage

``` r
DEtestRingPlot(
  srt,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
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
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  tile_height = 0.3,
  tile_gap = 0.1,
  jitter_width = 0.5,
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

- DE_threshold:

  A character string specifying the threshold for differential
  expression (used to highlight significant genes in all plot types).
  Default is `"avg_log2FC > 0 & p_val_adj < 0.05"`.

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

- tile_height:

  Height of the cell-type track in ring plot. Default is `0.3`.

- tile_gap:

  Gap between the track and nudged points in ring plot. Default is
  `0.1`.

- jitter_width:

  Horizontal jitter range for points in Manhattan plot. Default is
  `0.5`.

- ring_segments:

  Whether to draw segment lines between cell types in ring plot. Default
  is `TRUE`.

- seed:

  Random seed for jitter in ring plot. Default is `11`.

## See also

[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestManhattanPlot](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 15:30:15] Start standard processing workflow...
#> ℹ [2026-04-02 15:30:16] Checking a list of <Seurat>...
#> ! [2026-04-02 15:30:16] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:30:16] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:30:17] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:30:18] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:30:18] Number of available HVF: 2000
#> ℹ [2026-04-02 15:30:18] Finished check
#> ℹ [2026-04-02 15:30:18] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:30:19] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:30:24] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:30:24] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:30:24. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-13); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-13 (‘future_lapply-1’)
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
#> ! [2026-04-02 15:30:30] Infinite values detected
#> ! [2026-04-02 15:30:30] Data in the 'data' layer is unknown. Please check the data type
#> ℹ [2026-04-02 15:30:30] Start differential expression test
#> ℹ [2026-04-02 15:30:30] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-04-02 15:30:30] Using 1 core
#> Warning: No layers found matching search pattern provided
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘data’ is empty
#> ⠙ [2026-04-02 15:30:30] Running for Ductal [1/5] ■■■■■■■                       …
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
#> ✔ [2026-04-02 15:30:30] Completed 5 tasks in 34ms
#> 
#> ℹ [2026-04-02 15:30:30] Building results
#> ! [2026-04-02 15:30:30] Found 5 failed results
#> ℹ [2026-04-02 15:30:30] ✖ Error details:
#> ℹ                       ✖ "Ductal": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Ngn3-high-EP": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Endocrine": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Ngn3-low-EP": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> ℹ                       ✖ "Pre-endocrine": error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds
#> Error in `[.data.frame`(AllMarkers, , "group1"): undefined columns selected
DEtestRingPlot(
  pancreas_sub,
  group.by = "CellType"
)
#> Error in get_de_data(srt, group.by, test.use, DE_threshold, res): Cannot find the DEtest result for the group "CellType". Perform
#> `RunDEtest()` first
```
