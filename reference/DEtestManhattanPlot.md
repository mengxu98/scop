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
  group_use = NULL,
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
  only.pos = FALSE,
  label.by = c("p_val_adj", "p_val", "diff_pct", "avg_log2FC"),
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  palette = "RdBu",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  manhattan.bg = "white",
  group_track_width = NULL,
  group_track_height = NULL,
  jitter_width = 0.5,
  jitter_height = 0,
  seed = 11,
  aspect.ratio = NULL,
  xlab = NULL,
  ylab = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object or `SummarizedExperiment` object containing the
  results of differential expression analysis.

- group.by:

  Metadata column(s) used to color cells.

- test.use:

  Type of statistical test to use.

- res:

  A `data.frame` or `data.table` with differential expression results.
  When `res` is provided, `srt` will be ignored. The data.frame must
  contain columns: `gene`, `group1` (factor or character), `avg_log2FC`,
  `p_val_adj`, and optionally `pct.1` and `pct.2` for calculating
  `diff_pct`.

- group_use:

  Groups to plot. Default is `NULL` (all groups).

- DE_threshold:

  Threshold for differential expression (used to highlight significant
  genes in all plot types). Default is `"p_val < 0.05"` for sample-level
  methods (`"edgeR"` and `"limma"`). For cell-level volcano plots, it is
  `"abs(avg_log2FC) > 0 & p_val_adj < 0.05"` when `only.pos = FALSE`,
  and `"avg_log2FC > 0 & p_val_adj < 0.05"` otherwise.

- group_palette:

  Palette for cell types (groups) in Manhattan plot.

- group_palcolor:

  Custom colors for cell types (groups) in Manhattan plot.

- pt.size:

  The size of the points.

- pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- cols.highlight:

  Color for highlighted points.

- sizes.highlight:

  The size of the highlighted points.

- alpha.highlight:

  The transparency of the highlighted points.

- stroke.highlight:

  The stroke width for the highlighted points.

- nlabel:

  An integer value specifying the number of labeled points per group.

- features_label:

  Feature labels to plot.

- only.pos:

  Whether to show only positive log2 fold-change results in differential
  expression visualizations.

- label.by:

  Metric used to select automatic labels when `features_label = NULL`.
  Options are `"p_val_adj"`, `"p_val"`, `"diff_pct"`, and
  `"avg_log2FC"`. Smaller p-values are ranked first; `diff_pct` and
  `avg_log2FC` use the strongest positive and negative effects within
  each group.

- label.fg:

  Color for the labels' foreground.

- label.bg:

  Color for the labels' background.

- label.bg.r:

  The radius of the rounding of the labels' background.

- label.size:

  The size of the labels.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).

- palcolor:

  Custom colors used to create a color palette.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- manhattan.bg:

  Background color for Manhattan plot.

- group_track_width:

  Width of the centered cell-type track in Manhattan plot. Default is
  `NULL`, which uses the current automatic width.

- group_track_height:

  Height of the centered cell-type track in Manhattan plot. Default is
  `NULL`, which uses the current automatic height.

- jitter_width:

  Horizontal jitter range for points in Manhattan plot.

- jitter_height:

  Vertical jitter range for points in Manhattan plot.

- seed:

  Random seed for jitter in Manhattan and ring plots.

- aspect.ratio:

  Aspect ratio of the panel.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestRingPlot](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:06:22] Start standard processing workflow...
#> ℹ [2026-09-06 21:06:23] Checking a list of <Seurat>...
#> ! [2026-09-06 21:06:23] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:06:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:06:23] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:06:24] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:06:24] Number of available HVF: 2000
#> ℹ [2026-09-06 21:06:24] Finished check
#> ℹ [2026-09-06 21:06:24] Perform `ScaleData()`
#> ℹ [2026-09-06 21:06:24] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:06:24] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:06:24] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:06:24] Reorder clusters...
#> ℹ [2026-09-06 21:06:24] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:06:24] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:06:30] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> ℹ [2026-09-06 21:06:31] Data type is log-normalized
#> ℹ [2026-09-06 21:06:31] Start differential expression test
#> ℹ [2026-09-06 21:06:31] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-09-06 21:06:31] Using 1 core
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> ⠙ [2026-09-06 21:06:31] Running for Ductal [1/5] ■■          20% | ETA: 49s
#> ⠹ [2026-09-06 21:06:31] Running for Ngn3-high-EP [2/5] ■■■■        40% | ETA: 3…
#> ⠸ [2026-09-06 21:06:31] Running for Endocrine [3/5] ■■■■■■      60% | ETA: 21s
#> ⠼ [2026-09-06 21:06:31] Running for Ngn3-low-EP [4/5] ■■■■■■■■    80% | ETA: 10s
#> ✔ [2026-09-06 21:06:31] Completed 5 tasks in 47.1s
#> 
#> ℹ [2026-09-06 21:06:31] Building results
#> ✔ [2026-09-06 21:07:18] Differential expression test completed
DEtestManhattanPlot(
  pancreas_sub,
  group.by = "CellType"
)
```
