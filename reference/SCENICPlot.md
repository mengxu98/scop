# Plot top regulon specificity scores from SCENIC results

Calculate regulon specificity score (RSS) from SCENIC regulon activity
and plot the top regulons for each group.

## Usage

``` r
SCENICPlot(
  srt,
  group.by,
  tool_name = "SCENIC",
  assay = "scenic",
  layer = "data",
  plot_type = c("rss_rank", "rss_heatmap", "rss_dotplot", "activity_heatmap",
    "activity_violin", "activity_dim", "regulon_size", "network_graph", "network",
    "target_bar"),
  features = NULL,
  reduction = NULL,
  dims = c(1, 2),
  top_n = 12,
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
  max_targets = 20,
  max_edges = Inf,
  network_layout = c("fr", "nicely", "kk", "lgl", "drl"),
  network_tf = NULL,
  label_nodes = c("tfs", "all", "none"),
  network_label_top_n = 60,
  combine = TRUE,
  ncol = 3,
  return_data = TRUE,
  title = NULL,
  point_color = "#1F77B4",
  top_color = "#DC050C",
  point_size = 2,
  point_alpha = 0.5,
  highlight_tf = NULL,
  highlight_color = "#7A0177",
  highlight_point_size = 2,
  highlight_linewidth = 0.5,
  label_size = 3,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object containing SCENIC results from
  [`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md).

- group.by:

  Metadata column used as the cell group annotation.

- tool_name:

  Name of the `srt@tools` entry storing SCENIC results.

- assay:

  Assay used as a fallback source of regulon activity.

- layer:

  Assay layer used as a fallback source of regulon activity.

- plot_type:

  Plot type. `"rss_rank"` keeps the original regulon RSS rank plot.
  Other options summarize RSS, regulon activity, regulon sizes, or
  TF-target subnetworks.

- features:

  Optional TF/regulon names used by activity, network, and target plots.
  Values can match either `"Sox9"` or `"Sox9(+)"`.

- reduction:

  Dimensional reduction used when `plot_type = "activity_dim"`. If
  `NULL`, a UMAP/tSNE/PCA-like reduction is selected when available.

- dims:

  Two reduction dimensions used when `plot_type = "activity_dim"`.

- top_n:

  Number of top regulons labeled for each group.

- activity_scale:

  Whether to z-score each regulon across groups in
  `plot_type = "activity_heatmap"`. The default is `FALSE` so that the
  heatmap shows mean regulon activity and does not collapse constant
  regulons to zero.

- rss_scale:

  Whether to z-score each regulon across groups in
  `plot_type = "rss_heatmap"`. Use `rss_scale = TRUE`,
  `activity_scale = TRUE`, and the same `heatmap_limits` value when RSS
  and activity heatmaps should use a comparable row-wise relative scale.

- heatmap_show_row_names, heatmap_show_column_names:

  Whether to show row and column names in `plot_type = "rss_heatmap"`
  and `plot_type = "activity_heatmap"`.

- heatmap_cluster_rows, heatmap_cluster_columns:

  Whether to cluster rows and columns in SCENIC heatmaps.

- heatmap_order:

  Row ordering strategy for `plot_type = "rss_heatmap"` and
  `plot_type = "activity_heatmap"`. `"cluster"` keeps the existing
  dendrogram-based order, `"group"` groups regulons by the group where
  each regulon reaches its maximum heatmap value, and `"input"` keeps
  the resolved feature order. `"group"` and `"input"` disable row
  clustering so the chosen order is preserved.

- heatmap_row_names_side, heatmap_column_names_side:

  Sides used for row and column names in SCENIC heatmaps.

- heatmap_row_names_rot, heatmap_column_names_rot:

  Rotation angles for row and column names in SCENIC heatmaps.

- heatmap_border:

  Whether to draw heatmap borders in SCENIC heatmaps.

- heatmap_palette, heatmap_palcolor:

  Palette passed to
  [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  or
  [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  for SCENIC heatmaps. If `heatmap_palette = NULL`, a sensible default
  is selected for each heatmap type.

- heatmap_group_palette, heatmap_group_palcolor:

  Group annotation palette passed to
  [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  or
  [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  for SCENIC heatmaps.

- heatmap_limits:

  Optional two-length numeric vector used as the color scale limits for
  `plot_type = "rss_heatmap"` and `plot_type = "activity_heatmap"`. For
  example, `c(-2, 2)` fixes both z-score heatmaps to the same legend
  range.

- heatmap_args:

  Additional arguments passed to
  [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  for `plot_type = "activity_heatmap"` or
  [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  for `plot_type = "rss_heatmap"`.

- max_targets:

  Maximum number of target genes shown per TF/regulon in network-style
  plots.

- max_edges:

  Maximum number of TF-target edges shown in global network plots. Edges
  are ranked by absolute weight when a weight column is present.

- network_layout:

  Graph layout used by network plots. `"fr"` matches the force-directed
  Fruchterman-Reingold layout used in Pando examples.

- network_tf:

  Optional TF names used when `plot_type = "network"`. If `NULL`,
  `features`, `highlight_tf`, or the top RSS regulons are used.

- label_nodes:

  Which nodes to label in network plots.

- network_label_top_n:

  Maximum number of high-degree TF nodes labeled in
  `plot_type = "network_graph"` when `label_nodes = "tfs"`.

- combine:

  Whether to combine group plots with
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- ncol:

  Number of columns used when `combine = TRUE`.

- return_data:

  Whether to return RSS matrices and ranking tables together with plots.
  If `FALSE`, only the plot object or plot list is returned.

- title:

  Optional title added to the combined plot.

- point_color:

  Color for all regulon rank points.

- top_color:

  Color for top regulon rank points.

- point_size:

  Point size.

- point_alpha:

  Alpha for all regulon rank points.

- highlight_tf:

  Optional TF or regulon names to highlight in every group plot. Values
  can match either `TF` or `regulon`, for example `"Sox9"` or
  `"Sox9(+)"`.

- highlight_color:

  Color for highlighted TF or regulon points and rank lines.

- highlight_point_size:

  Point size for highlighted TFs or regulons.

- highlight_linewidth:

  Line width for highlighted TF or regulon rank lines.

- label_size:

  Text size for top regulon labels.

- verbose:

  Whether to print messages.

- ...:

  Additional arguments passed directly to the underlying
  [`GroupHeatmap()`](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)
  or
  [`FeatureHeatmap()`](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
  call when `plot_type` is `"activity_heatmap"` or `"rss_heatmap"`. For
  example, `width` and `height` can be supplied directly.

## Value

A list containing `rss_matrix`, `rank_table`, `top_table`, `plots`, and
`plot` when `return_data = TRUE`; otherwise a plot object or list of
plots.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSCENIC(
  pancreas_sub,
  species = "Mus_musculus"
)

scenic_rss <- SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "rss_rank"
)
scenic_rss$plot
example_regulons <- unique(scenic_rss$top_table$regulon)[1:2]
example_tfs <- unique(scenic_rss$top_table$TF)[1:2]

SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "rss_heatmap")
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "rss_heatmap",
  width = 2,
  height = 3
)
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "rss_dotplot")
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "activity_heatmap")
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "rss_heatmap",
  heatmap_order = "group",
  heatmap_cluster_columns = FALSE
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "rss_heatmap",
  rss_scale = TRUE,
  heatmap_order = "group",
  heatmap_limits = c(-2, 2)
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "activity_heatmap",
  activity_scale = TRUE,
  heatmap_limits = c(-2, 2)
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "activity_violin",
  features = example_regulons
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "activity_dim",
  features = example_regulons
)
SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "regulon_size")
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "network_graph",
  max_targets = 10,
  max_edges = 500,
  label_nodes = "tfs"
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "network",
  network_tf = example_tfs,
  max_targets = 30
)
SCENICPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "target_bar",
  features = example_regulons,
  max_targets = 20
)
} # }
```
