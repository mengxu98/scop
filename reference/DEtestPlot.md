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
  group_use = NULL,
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  x_metric = NULL,
  y_metric = c("p_val_adj", "p_val"),
  x_order = c("gene", "index"),
  palette = "RdBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.background = "grey80",
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
  group_track_width = NULL,
  group_track_height = NULL,
  jitter_width = 0.5,
  jitter_height = 0,
  tile_height = 0.3,
  tile_gap = 0.1,
  ring_segments = TRUE,
  seed = 11,
  threshold_method = c("rectangular", "hyperbolic"),
  hyperbola_c = 6,
  annotate_enrichment = FALSE,
  enrich_from = c("Enrichment", "GSEA", "GSVA"),
  enrich_db = NULL,
  enrich_terms = NULL,
  enrich_top_terms = 3,
  enrich_padj_cutoff = 0.05,
  enrich_gsva_score_cutoff = NULL,
  gsva_method = NULL,
  enrich_nlabel = 15
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

- plot_type:

  Type of plot to create. Options: `"volcano"`, `"manhattan"`, or
  `"ring"`.

- group_use:

  Groups to plot. Default is `NULL` (all groups).

- DE_threshold:

  Threshold for differential expression (used to highlight significant
  genes in all plot types). Default is `"p_val < 0.05"` for sample-level
  methods (`"edgeR"` and `"limma"`). For cell-level volcano plots, it is
  `"abs(avg_log2FC) > 0 & p_val_adj < 0.05"` when `only.pos = FALSE`,
  and `"avg_log2FC > 0 & p_val_adj < 0.05"` otherwise.

- x_metric:

  Metric to use for the x-axis (only for volcano plot). Default is
  `NULL`, which uses `"avg_log2FC"` for sample-level methods (`"edgeR"`
  and `"limma"`) and `"diff_pct"` otherwise.

- y_metric:

  Metric to use for the y-axis. Options: `"p_val"` or `"p_val_adj"`.
  Default is `"p_val"` for sample-level methods (`"edgeR"` and
  `"limma"`) and `"p_val_adj"` otherwise.

- x_order:

  How to order genes on x-axis (only for Manhattan plot, not used
  currently). Options: `"gene"` (alphabetical by gene name) or `"index"`
  (by data order).

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).

- palcolor:

  Custom colors used to create a color palette.

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

- cols.background:

  Color for non-DE background points in volcano plots.

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

- aspect.ratio:

  Aspect ratio of the panel.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

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

- tile_height:

  Height of the cell-type track in ring plot.

- tile_gap:

  Gap between the track and nudged points in ring plot.

- ring_segments:

  Whether to draw segment lines between cell types in ring plot.

- seed:

  Random seed for jitter in Manhattan and ring plots.

- threshold_method:

  Volcano significance threshold method. Options are `"rectangular"`
  (legacy DE_threshold) or `"hyperbolic"`
  (`|log2FC * -log10(padj)| > c`).

- hyperbola_c:

  Numeric cutoff `c` for hyperbolic volcano threshold.

- annotate_enrichment:

  Whether to annotate enrichment-hit genes on volcano plots. Enrichment
  results are read from existing results in `srt@tools` only.

- enrich_from:

  Character vector specifying enrichment result source(s) to annotate.
  Options are `"Enrichment"`, `"GSEA"`, `"GSVA"`.

- enrich_db:

  Optional database filter for enrichment annotation, e.g. `"GO_BP"` or
  `"KEGG"`.

- enrich_terms:

  Optional whitelist of enrichment term IDs or names for annotation.

- enrich_top_terms:

  Number of top enriched terms selected per source/group/database.

- enrich_padj_cutoff:

  Adjusted p-value cutoff for `"Enrichment"` and `"GSEA"` annotation.

- enrich_gsva_score_cutoff:

  Optional absolute GSVA score cutoff for `"GSVA"` annotation.

- gsva_method:

  Optional GSVA method filter (e.g. `"gsva"` or `"ssgsea"`) when
  multiple GSVA tool slots exist.

- enrich_nlabel:

  Maximum number of enrichment-derived labels added per group. Labels
  from `features_label` are always retained.

## See also

[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md),
[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[DEtestManhattanPlot](https://mengxu98.github.io/scop/reference/DEtestManhattanPlot.md),
[DEtestRingPlot](https://mengxu98.github.io/scop/reference/DEtestRingPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:07:21] Start standard processing workflow...
#> ℹ [2026-09-06 21:07:22] Checking a list of <Seurat>...
#> ! [2026-09-06 21:07:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:07:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:07:22] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:07:22] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:07:22] Number of available HVF: 2000
#> ℹ [2026-09-06 21:07:22] Finished check
#> ℹ [2026-09-06 21:07:22] Perform `ScaleData()`
#> ℹ [2026-09-06 21:07:22] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:07:22] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:07:23] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:07:23] Reorder clusters...
#> ℹ [2026-09-06 21:07:23] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:07:23] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:07:28] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> ℹ [2026-09-06 21:07:29] Data type is log-normalized
#> ℹ [2026-09-06 21:07:29] Start differential expression test
#> ℹ [2026-09-06 21:07:29] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-09-06 21:07:29] Using 1 core
#> ⠙ [2026-09-06 21:07:29] Running for Ductal [1/5] ■■          20% | ETA: 49s
#> ⠹ [2026-09-06 21:07:29] Running for Ngn3-high-EP [2/5] ■■■■        40% | ETA: 3…
#> ⠸ [2026-09-06 21:07:29] Running for Endocrine [3/5] ■■■■■■      60% | ETA: 20s
#> ⠼ [2026-09-06 21:07:29] Running for Ngn3-low-EP [4/5] ■■■■■■■■    80% | ETA: 10s
#> ✔ [2026-09-06 21:07:29] Completed 5 tasks in 46.2s
#> 
#> ℹ [2026-09-06 21:07:29] Building results
#> ✔ [2026-09-06 21:08:15] Differential expression test completed

DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "volcano",
  ncol = 2
)


DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "volcano",
  group_use = c("Ductal", "Endocrine"),
  ncol = 2
)


DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  test.use = "wilcox",
  plot_type = "volcano",
  x_metric = "avg_log2FC",
  y_metric = "p_val_adj",
  group_use = c("Ductal", "Endocrine"),
  DE_threshold = "abs(avg_log2FC) > 0.25 & p_val_adj < 0.05"
)


DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "volcano",
  threshold_method = "hyperbolic",
  hyperbola_c = 6,
  ncol = 2
)


pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-09-06 21:08:22] Start Enrichment analysis
#> ℹ [2026-09-06 21:08:22] Species: "Mus_musculus"
#> 
#> 
#> ℹ [2026-09-06 21:13:25] Preparing database: GO_BP
#> ℹ [2026-09-06 21:13:56] Convert ID types for the GO_BP database
#> ℹ [2026-09-06 21:13:56] Converted ID types using local annotation package org.Mm.eg.db
#> ℹ [2026-09-06 21:13:57] Permform enrichment...
#> ℹ [2026-09-06 21:13:59] Using 1 core
#> ⠙ [2026-09-06 21:13:59] Running for 1 [1/5] ■■          20% | ETA:  2s
#> ✔ [2026-09-06 21:13:59] Completed 5 tasks in 2.6s
#> 
#> ℹ [2026-09-06 21:13:59] Building results
#> ✔ [2026-09-06 21:14:01] Enrichment analysis done
DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "volcano",
  threshold_method = "hyperbolic",
  hyperbola_c = 6,
  annotate_enrichment = TRUE,
  enrich_from = "Enrichment",
  enrich_db = "GO_BP",
  enrich_top_terms = 3,
  enrich_nlabel = 15,
  ncol = 2
)


DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "manhattan"
)


DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "ring"
)


de_results1 <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
DEtestPlot(
  res = de_results1,
  plot_type = "volcano",
  ncol = 2
)


de_results2 <- Seurat::FindMarkers(
  pancreas_sub,
  group.by = "CellType",
  ident.1 = "Ductal",
  ident.2 = "Endocrine"
)
DEtestPlot(
  res = de_results2,
  plot_type = "volcano"
)


de_results3 <- Seurat::FindAllMarkers(
  pancreas_sub,
  group.by = "CellType"
)
#> Calculating cluster Ductal
#> Calculating cluster Ngn3-high-EP
#> Calculating cluster Endocrine
#> Calculating cluster Ngn3-low-EP
#> Calculating cluster Pre-endocrine
DEtestPlot(
  res = de_results3,
  plot_type = "volcano",
  ncol = 2
)
```
