# Ferris Wheel Plot

Draw a Ferris wheel-style donut plot for up/down counts across pathways
or other categories. The center donut shows the total counts, and outer
donuts show category-level counts with radius scaled by category size.

## Usage

``` r
FerrisWheelPlot(
  data = NULL,
  res = NULL,
  label_col = "pathway",
  up_col = "up",
  down_col = "down",
  total_up = NULL,
  total_down = NULL,
  de_results = NULL,
  term_col = "Description",
  group_col = "Groups",
  count_col = "Count",
  padj_col = "p.adjust",
  up_group = "Up",
  down_group = "Down",
  padj_cutoff = 0.05,
  top_n = 10,
  up_label = "Up-regulated Genes",
  down_label = "Down-regulated Genes",
  up_color = NULL,
  down_color = NULL,
  palette = "Chinese",
  palcolor = NULL,
  center_label = "Significant Genes",
  center_radius = 1.25,
  center_width = 0.28,
  outer_distance = 2.55,
  outer_width = 0.12,
  min_outer_radius = 0.28,
  max_outer_radius = 0.48,
  label_wrap = 26,
  label_case = c("title", "none"),
  label.stroke = 0.2,
  label.stroke.color = "white",
  label.size = 3,
  number.size = 3,
  center.size = 4.5,
  legend.size = 3.5,
  line.color = "grey35",
  line.linewidth = 0.35,
  legend.position = c("bottom", "none"),
  direction = 1,
  start = pi/2,
  theme_use = "theme_scop",
  theme_args = list()
)
```

## Arguments

- data:

  A `data.frame` containing category labels and two non-negative count
  columns. If `NULL`, the table is built from `res`.

- res:

  A result list returned by
  [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  or an enrichment `data.frame`. When provided, `FerrisWheelPlot()`
  summarizes enriched terms into up/down counts internally.

- label_col:

  Column name for category labels. Default is `"pathway"`.

- up_col:

  Column name for up-regulated counts. Default is `"up"`.

- down_col:

  Column name for down-regulated counts. Default is `"down"`.

- total_up, total_down:

  Total counts shown in the center donut. When `NULL`, the corresponding
  column sum is used. For enrichment input, these can be inferred from
  `de_results`.

- de_results:

  Optional differential expression result table used to infer `total_up`
  and `total_down` from `avg_log2FC`.

- term_col, group_col, count_col, padj_col:

  Column names used when `res` is an enrichment result. Defaults match
  [`RunEnrichment()`](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  output.

- up_group, down_group:

  Group labels in enrichment results. Default matches examples that use
  `"Up"` and `"Down"` as `geneID_groups`.

- padj_cutoff:

  Adjusted p-value cutoff used to select enriched terms from `res`.

- top_n:

  Number of enriched terms to show when `res` is provided. Default is
  `10`.

- up_label, down_label:

  Legend labels and center text labels. Defaults are
  `"Up-regulated Genes"` and `"Down-regulated Genes"`.

- up_color, down_color:

  Optional colors for up/down counts. When `NULL`, colors are selected
  from the `"Chinese"` palette used by `scop`.

- palette, palcolor:

  Palette passed to
  [`thisplot::palette_colors()`](https://mengxu98.github.io/thisplot/reference/palette_colors.html)
  when `up_color` or `down_color` is not provided.

- center_label:

  Text shown inside the center donut. Default is `"Significant Genes"`.

- center_radius:

  Radius of the center donut. Default is `1.25`.

- center_width:

  Width of the center donut band. Default is `0.28`.

- outer_distance:

  Distance from plot center to outer donut centers. Default is `2.55`.

- outer_width:

  Width of outer donut bands. Default is `0.12`.

- min_outer_radius, max_outer_radius:

  Minimum and maximum outer donut radii after scaling by category count.

- label_wrap:

  Maximum characters per category label line. Set `NULL` to disable
  wrapping.

- label_case:

  Label case for category names. `"title"` capitalizes words by default;
  `"none"` keeps labels unchanged.

- label.stroke:

  White outline width around text labels. Default is `0.1`. Set `0` to
  disable the outline.

- label.stroke.color:

  Outline color for text labels.

- label.size, number.size, center.size, legend.size:

  Text sizes.

- line.color, line.linewidth:

  Connector line styling.

- legend.position:

  Position of the compact legend. Use `"bottom"` or `"none"`.

- direction:

  Drawing direction for the outer donuts. `1` is clockwise and `-1` is
  counterclockwise.

- start:

  Start angle in radians for the first category. Default is `pi / 2`.

- theme_use:

  Theme function name. Default is `"theme_scop"`.

- theme_args:

  Additional arguments passed to `theme_use`.

## Value

A `ggplot` object.

## See also

[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
[DEtestPlot](https://mengxu98.github.io/scop/reference/DEtestPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-27 19:37:07] Start standard processing workflow...
#> ℹ [2026-06-27 19:37:08] Checking a list of <Seurat>...
#> ! [2026-06-27 19:37:08] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-27 19:37:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 19:37:08] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 19:37:08] Use the separate HVF from `srt_list`
#> ℹ [2026-06-27 19:37:08] Number of available HVF: 2000
#> ℹ [2026-06-27 19:37:08] Finished check
#> ℹ [2026-06-27 19:37:08] Perform `ScaleData()`
#> ℹ [2026-06-27 19:37:08] Perform pca linear dimension reduction
#> ℹ [2026-06-27 19:37:09] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-27 19:37:09] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-27 19:37:09] Reorder clusters...
#> ℹ [2026-06-27 19:37:09] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-27 19:37:09] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-27 19:37:15] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  only.pos = FALSE
)
#> ℹ [2026-06-27 19:37:15] Data type is log-normalized
#> ℹ [2026-06-27 19:37:15] Start differential expression test
#> ℹ [2026-06-27 19:37:15] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-06-27 19:37:15] Using 1 core
#> ⠙ [2026-06-27 19:37:15] Running for Ductal [1/5] ■■          20% | ETA:  0s
#> ✔ [2026-06-27 19:37:15] Completed 5 tasks in 499ms
#> 
#> ℹ [2026-06-27 19:37:15] Building results
#> ✔ [2026-06-27 19:37:16] Differential expression test completed

de_df <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
de_df <- de_df[
  de_df$p_val_adj < 0.05 & abs(de_df$avg_log2FC) > 0.25,
  ,
  drop = FALSE
]
de_df$direction <- ifelse(de_df$avg_log2FC > 0, "Up", "Down")

enrich_out <- RunEnrichment(
  geneID = de_df$gene,
  geneID_groups = de_df$direction,
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ℹ [2026-06-27 19:37:16] Start Enrichment analysis
#> ℹ [2026-06-27 19:37:16] Species: "Mus_musculus"
#> ℹ [2026-06-27 19:37:16] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-27 19:31:20
#> ℹ [2026-06-27 19:37:17] Permform enrichment...
#> ℹ [2026-06-27 19:37:18] Using 1 core
#> ⠙ [2026-06-27 19:37:18] Running for 1 [1/2] ■■■■■       50% | ETA:  1s
#> ✔ [2026-06-27 19:37:18] Completed 2 tasks in 1.2s
#> 
#> ℹ [2026-06-27 19:37:18] Building results
#> ✔ [2026-06-27 19:37:20] Enrichment analysis done
FerrisWheelPlot(
  res = enrich_out,
  de_results = de_df
)
```
