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
  tile_height = 0.3,
  tile_gap = 0.1,
  jitter_width = 0.5,
  ring_segments = TRUE,
  seed = 11
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

- tile_height:

  Height of the cell-type track in ring plot.

- tile_gap:

  Gap between the track and nudged points in ring plot.

- jitter_width:

  Horizontal jitter range for points in Manhattan plot.

- ring_segments:

  Whether to draw segment lines between cell types in ring plot.

- seed:

  Random seed for jitter in Manhattan and ring plots.

## See also

[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md),
[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestManhattanPlot](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:21:24] Start standard processing workflow...
#> ℹ [2026-09-06 21:21:24] Checking a list of <Seurat>...
#> ! [2026-09-06 21:21:24] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:21:24] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:21:24] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:21:25] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:21:25] Number of available HVF: 2000
#> ℹ [2026-09-06 21:21:25] Finished check
#> ℹ [2026-09-06 21:21:25] Perform `ScaleData()`
#> ℹ [2026-09-06 21:21:25] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:21:25] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:21:25] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:21:25] Reorder clusters...
#> ℹ [2026-09-06 21:21:25] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:21:25] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:21:31] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> ℹ [2026-09-06 21:21:32] Data type is log-normalized
#> ℹ [2026-09-06 21:21:32] Start differential expression test
#> ℹ [2026-09-06 21:21:32] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-09-06 21:21:32] Using 1 core
#> ⠙ [2026-09-06 21:21:32] Running for Ductal [1/5] ■■          20% | ETA: 47s
#> ⠹ [2026-09-06 21:21:32] Running for Ngn3-high-EP [2/5] ■■■■        40% | ETA: 3…
#> ⠸ [2026-09-06 21:21:32] Running for Endocrine [3/5] ■■■■■■      60% | ETA: 20s
#> ⠼ [2026-09-06 21:21:32] Running for Ngn3-low-EP [4/5] ■■■■■■■■    80% | ETA:  9s
#> ✔ [2026-09-06 21:21:32] Completed 5 tasks in 45.1s
#> 
#> ℹ [2026-09-06 21:21:32] Building results
#> ✔ [2026-09-06 21:22:17] Differential expression test completed
DEtestRingPlot(
  pancreas_sub,
  group.by = "CellType"
)
```
