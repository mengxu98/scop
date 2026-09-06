# Plot SCENIC or SCENIC+ regulon activity

Plot RSS and regulon activity from
[`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md)
or
[`RunSCENICPlus()`](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md).
Use
[`SCENICPlusPlot()`](https://mengxu98.github.io/scop/reference/SCENICPlusPlot.md)
for SCENIC+ defaults. Network `palette = "RdYlBu"` uses `"Chinese"`.
`"network"` draws one hub per TF. Network plots have no title. When one
TF has multiple regulons, network legends show the regulon with the
highest RSS; the corresponding cell type and score remain available in
the returned annotation table.

## Usage

``` r
SCENICPlot(
  srt,
  group.by,
  tool_name = "SCENIC",
  assay = "scenic",
  layer = "data",
  plot_type = c("rss_rank", "rss_heatmap", "rss_dotplot", "heatmap_dotplot",
    "activity_heatmap", "activity_violin", "activity_dim", "eregulon_dim",
    "activity_cor_dumbbell", "regulon_size", "network_graph", "network", "egrn",
    "overlap", "target_bar", "coverage"),
  features = NULL,
  reduction = NULL,
  dims = c(1, 2),
  cor.features = NULL,
  cor.feature.labels = NULL,
  cor_method = c("spearman", "pearson", "kendall"),
  p_cutoff = 0.05,
  cor_sort.by = c("difference", "max_abs", "input"),
  cor_cols = c("#56B4E9", "#E83F6F"),
  cor_xlim = NULL,
  cor_label = TRUE,
  top_n = 12,
  regulon_label = c("auto", "regulon", "tf"),
  rss_rank_yscale = c("shared", "free"),
  activity_scale = FALSE,
  rss_scale = FALSE,
  heatmap_show_row_names = FALSE,
  heatmap_show_column_names = FALSE,
  heatmap_cluster_rows = TRUE,
  heatmap_cluster_columns = FALSE,
  heatmap_order = c("cluster", "group", "input"),
  heatmap_row_names_side = "right",
  heatmap_column_names_side = "top",
  heatmap_row_names_rot = 0,
  heatmap_column_names_rot = 45,
  heatmap_border = TRUE,
  heatmap_palette = NULL,
  heatmap_palcolor = NULL,
  heatmap_group_palette = "Chinese",
  heatmap_group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list(),
  palette = "RdYlBu",
  palcolor = NULL,
  compare_expression = FALSE,
  expression_assay = NULL,
  expression_layer = "data",
  expression_scale = TRUE,
  violin_args = list(),
  dim_args = list(),
  atac_assay = NULL,
  extend.upstream = 50000,
  extend.downstream = 50000,
  coverage_args = list(),
  max_targets = 20,
  max_edges = Inf,
  network_layout = c("auto", "star", "hub", "tripartite", "kk", "fr"),
  network_tf = NULL,
  network_include_regions = TRUE,
  label_nodes = c("auto", "tfs", "all", "none"),
  network_label_top_n = 60,
  combine = TRUE,
  ncol = 3,
  return_data = TRUE,
  title = NULL,
  point_color = "#B8C2CC",
  top_color = "#B64342",
  point_size = 1.8,
  point_alpha = 0.85,
  highlight_tf = NULL,
  highlight_color = "#7A0177",
  highlight_point_size = 2,
  highlight_linewidth = 0.5,
  label_size = 3,
  label_max_overlaps = Inf,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object with SCENIC or SCENIC+ results.

- group.by:

  Metadata column for cell groups.

- tool_name:

  `srt@tools` entry. Default `"SCENIC"`.

- assay:

  Fallback assay for regulon activity.

- layer:

  Assay layer for regulon activity.

- plot_type:

  One of `"rss_rank"`, `"rss_heatmap"`, `"rss_dotplot"`,
  `"heatmap_dotplot"`, `"activity_heatmap"`, `"activity_violin"`,
  `"activity_dim"`, `"eregulon_dim"`, `"activity_cor_dumbbell"`,
  `"regulon_size"`, `"network"`, `"network_graph"`, `"egrn"`,
  `"overlap"`, `"target_bar"`, `"coverage"`.

- features:

  TF or regulon names.

- reduction:

  Reduction for `"activity_dim"` / `"eregulon_dim"`.

- dims:

  Two reduction dimensions.

- cor.features:

  Two metadata columns or genes for `"activity_cor_dumbbell"`.

- cor.feature.labels:

  Optional legend labels for `cor.features`.

- cor_method:

  Correlation method for `"activity_cor_dumbbell"`.

- p_cutoff:

  Significance cutoff for `"activity_cor_dumbbell"`.

- cor_sort.by:

  Row order: `"difference"`, `"max_abs"`, or `"input"`.

- cor_cols:

  Two colors for the two correlation targets.

- cor_xlim:

  Optional x-axis limits for the dumbbell plot.

- cor_label:

  Whether to label dumbbell points.

- top_n:

  Top regulons per group.

- regulon_label:

  Label source for RSS rank plots and network legends. `"auto"` keeps
  TF-only labels when a TF has one regulon and shows the regulon name
  when multiple regulons share a TF; `"regulon"` always shows regulon
  names; `"tf"` always shows TF names.

- rss_rank_yscale:

  Use a `"shared"` y-axis across RSS rank panels for direct score
  comparison, or `"free"` for the legacy per-panel scale.

- activity_scale:

  Z-score regulons in `"activity_heatmap"`.

- rss_scale:

  Z-score regulons in `"rss_heatmap"`.

- heatmap_show_row_names, heatmap_show_column_names:

  Heatmap names.

- heatmap_cluster_rows, heatmap_cluster_columns:

  Heatmap clustering.

- heatmap_order:

  `"cluster"`, `"group"`, or `"input"`.

- heatmap_row_names_side, heatmap_column_names_side:

  Name sides.

- heatmap_row_names_rot, heatmap_column_names_rot:

  Name rotation.

- heatmap_border:

  Heatmap borders.

- heatmap_palette, heatmap_palcolor:

  Heatmap palette.

- heatmap_group_palette, heatmap_group_palcolor:

  Group annotation palette.

- heatmap_limits:

  Color-scale limits for RSS/activity heatmaps.

- heatmap_args:

  Extra arguments for
  [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  /
  [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md).

- palette, palcolor:

  Palette for RSS, size, target, and network colors.

- compare_expression:

  Add TF expression next to `"activity_dim"`.

- expression_assay, expression_layer:

  Assay and layer for TF expression.

- expression_scale:

  Z-score TF expression in `"heatmap_dotplot"`.

- violin_args:

  Extra arguments for
  [`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)
  in `"activity_violin"`.

- dim_args:

  Extra arguments for
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md).

- atac_assay:

  Chromatin assay for `"coverage"`.

- extend.upstream, extend.downstream:

  Flanks for `"coverage"`.

- coverage_args:

  Extra arguments for
  [`CoverageTrackPlot()`](https://mengxu98.github.io/scop/reference/CoverageTrackPlot.md).

- max_targets:

  Max targets per TF in network plots.

- max_edges:

  Max edges in `"network_graph"`.

- network_layout:

  `"auto"`, `"star"`, `"kk"`, `"hub"`, `"tripartite"`, `"fr"`.
  `"network"` uses `"star"` for one TF and `"hub"` for several.

- network_tf:

  TFs for `"network"` / `"egrn"`.

- network_include_regions:

  Use SCENIC+ triplets when stored.

- label_nodes:

  `"auto"`, `"tfs"`, `"all"`, or `"none"`.

- network_label_top_n:

  Max TF labels in `"network_graph"`.

- combine:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- ncol:

  Number of columns of the combined plot.

- return_data:

  Return RSS tables with plots.

- title:

  Optional combined-plot title. Ignored by `"network"`,
  `"network_graph"`, and `"egrn"`.

- point_color:

  Rank-plot point color.

- top_color:

  Top-regulon point color.

- point_size:

  Point size.

- point_alpha:

  Rank-plot alpha.

- highlight_tf:

  TFs or regulons to highlight.

- highlight_color:

  Highlight color.

- highlight_point_size:

  Highlight point size.

- highlight_linewidth:

  Highlight line width.

- label_size:

  Top-regulon label size.

- label_max_overlaps:

  Max overlapping labels.

- verbose:

  Print messages.

- ...:

  Extra heatmap arguments such as `width` and `height`.

## Value

Plot, or a list with `rss_matrix`, `rank_table`, `top_table`, `plots`,
and `plot` when `return_data = TRUE`.

## See also

[RunSCENIC](https://mengxu98.github.io/scop/reference/RunSCENIC.md),
[RunSCENICPlus](https://mengxu98.github.io/scop/reference/RunSCENICPlus.md),
[SCENICPlusPlot](https://mengxu98.github.io/scop/reference/SCENICPlusPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunSCENIC(
  pancreas_sub,
  species = "Mus_musculus",
  backend = "cpp",
  work_dir = "test/scenic"
)
scenic_rss <- SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "rss_rank")
example_tfs <- unique(scenic_rss$top_table$TF)[1:3]
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "heatmap_dotplot")
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "network", features = example_tfs)
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "network_graph", features = example_tfs)
SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "egrn", features = example_tfs)
SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "coverage", features = example_tfs)
} # }
```
