# Plot scFEA module flux heatmap

Plot scFEA module flux heatmap

## Usage

``` r
scFEAHeatmap(
  srt,
  assay = "scFEAflux",
  layer = "data",
  group.by = NULL,
  features = NULL,
  modules = NULL,
  label_by = c("module", "reaction", "module_reaction"),
  add_sm_anno = TRUE,
  scale_rows = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  sm_anno_label_rot = 0,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_limit = 2,
  heatmap_column_width = NULL,
  column_names_rot = 90,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  feature_split_palette = "simspec",
  feature_split_palcolor = NULL,
  border = TRUE,
  use_raster = TRUE,
  raster_by_magick = FALSE,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 11),
  row_title_gp = grid::gpar(fontsize = 7.5),
  legend_title_gp = grid::gpar(fontsize = 8),
  legend_labels_gp = grid::gpar(fontsize = 7),
  row_gap = grid::unit(1.2, "mm"),
  sm_anno_size = grid::unit(2.5, "mm"),
  group_anno_size = grid::unit(3, "mm"),
  mark_features = NULL,
  mark_label_by = NULL,
  mark_labels_gp = NULL,
  mark_annotation_width = NULL,
  mark_link_width = grid::unit(5, "mm"),
  name = NULL,
  ht_params = list(),
  width = NULL,
  height = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by \[RunscFEA()\].

- assay:

  Flux assay name.

- layer:

  Flux assay layer.

- group.by:

  Metadata column used to aggregate cells.

- features:

  Optional scFEA modules to plot. Supports Seurat-safe feature ids
  (\`"M-1"\`), scFEA ids (\`"M_1"\`), labels (\`"M1"\`), or reaction
  labels.

- modules:

  Optional scFEA modules to plot. This is a user-facing alias for
  \`features\` that also accepts combined labels such as \`"M150: PRPP
  -\> UMP"\`.

- label_by:

  Row label mode: module id, reaction, or both.

- add_sm_anno:

  Whether to add \`SM_anno\` row annotation and split rows by pathway
  class.

- scale_rows:

  Whether to z-score each row after group aggregation.

- cluster_rows, cluster_columns:

  Passed to \[ComplexHeatmap::Heatmap()\].

- sm_anno_label_rot:

  Rotation angle for \`SM_anno\` row-split labels.

- show_row_names:

  Whether to show row labels.

- show_column_names:

  Whether to show aggregated group labels.

- heatmap_limit:

  Numeric clipping limit for row z-scores.

- heatmap_column_width:

  Width for each aggregated group column. Numeric values are interpreted
  as millimeters.

- column_names_rot:

  Rotation angle for heatmap column labels.

- heatmap_palette, heatmap_palcolor:

  Continuous heatmap palette passed to \`palette_colors()\`.

- group_palette, group_palcolor:

  Palette for the top group annotation.

- feature_split_palette, feature_split_palcolor:

  Palette for \`SM_anno\` row annotation.

- border:

  Whether to draw borders around heatmap cells and annotations.

- use_raster, raster_by_magick:

  Raster settings passed to \[ComplexHeatmap::Heatmap()\].

- column_names_gp, column_title_gp, row_title_gp:

  Font settings passed to \[grid::gpar()\]-aware ComplexHeatmap
  arguments.

- legend_title_gp, legend_labels_gp:

  Font settings for heatmap and annotation legends.

- row_gap:

  Gap between \`SM_anno\` row splits.

- sm_anno_size:

  Size of the left \`SM_anno\` annotation strip.

- group_anno_size:

  Size of the top group annotation strip.

- mark_features:

  Optional scFEA modules to mark with \[ComplexHeatmap::anno_mark()\].
  Accepts the same module ids and labels as \`features\`.

- mark_label_by:

  Label mode for marked rows. Defaults to \`label_by\`.

- mark_labels_gp:

  Font settings for marked row labels. If \`NULL\`, a compact default is
  chosen from \`mark_label_by\`.

- mark_annotation_width:

  Width of the right mark annotation. Numeric values are interpreted as
  millimeters.

- mark_link_width:

  Link width used by \[ComplexHeatmap::anno_mark()\].

- name:

  Heatmap legend title.

- ht_params:

  Additional parameters passed to \[ComplexHeatmap::Heatmap()\],
  overriding defaults when names overlap.

- width, height:

  Suggested export size. Stored as attributes on the returned heatmap
  object.

- verbose:

  Whether to print messages.

## Value

A \`ComplexHeatmap\` heatmap object.
