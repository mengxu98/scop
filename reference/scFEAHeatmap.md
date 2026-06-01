# Plot scFEA module flux heatmap

Aggregates scFEA module flux by Seurat metadata groups and draws a
heatmap with optional M168 supermodule annotation.

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
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  sm_anno_label_rot = 0,
  show_row_names = FALSE,
  heatmap_limit = 2,
  heatmap_column_width = grid::unit(10, "mm"),
  column_names_rot = 45,
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
  width = NULL,
  height = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by
  [`RunscFEA()`](https://mengxu98.github.io/scop/reference/RunscFEA.md).

- assay:

  Flux assay name.

- layer:

  Flux assay layer.

- group.by:

  Metadata column used to aggregate cells. If `NULL`, all cells are
  averaged together.

- features:

  Optional scFEA modules to plot. Supports Seurat-safe feature ids such
  as `"M-1"`, scFEA ids such as `"M_1"`, labels such as `"M1"`, or
  reaction labels.

- modules:

  Optional scFEA modules to plot. This is a user-facing alias for
  `features` that also accepts combined labels such as
  `"M150: PRPP -> UMP"`.

- label_by:

  Row label mode: module id, reaction, or both.

- add_sm_anno:

  Whether to add `SM_anno` row annotation and split rows by pathway
  class.

- scale_rows:

  Whether to z-score each row after group aggregation.

- cluster_rows, cluster_columns:

  Passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- sm_anno_label_rot:

  Rotation angle for `SM_anno` row-split labels.

- show_row_names:

  Whether to show row labels.

- heatmap_limit:

  Numeric clipping limit for row z-scores.

- heatmap_column_width:

  Width for each aggregated group column. Numeric values are interpreted
  as millimeters.

- column_names_rot:

  Rotation angle for heatmap column labels.

- column_names_gp, column_title_gp, row_title_gp:

  Font settings passed to
  [`grid::gpar()`](https://rdrr.io/r/grid/gpar.html)-aware
  ComplexHeatmap arguments.

- legend_title_gp, legend_labels_gp:

  Font settings for heatmap and annotation legends.

- row_gap:

  Gap between `SM_anno` row splits.

- sm_anno_size:

  Size of the left `SM_anno` annotation strip.

- group_anno_size:

  Size of the top group annotation strip.

- mark_features:

  Optional scFEA modules to mark with
  [`ComplexHeatmap::anno_mark()`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_mark.html).
  Accepts the same module ids and labels as `features`.

- mark_label_by:

  Label mode for marked rows. Defaults to `label_by`.

- mark_labels_gp:

  Font settings for marked row labels. If `NULL`, a compact default is
  chosen from `mark_label_by`.

- mark_annotation_width:

  Width of the right mark annotation. Numeric values are interpreted as
  millimeters.

- mark_link_width:

  Link width used by
  [`ComplexHeatmap::anno_mark()`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_mark.html).

- name:

  Heatmap legend title.

- width, height:

  Suggested export size. Stored as attributes on the returned heatmap
  object.

- verbose:

  Whether to print messages.

## Value

A ComplexHeatmap heatmap object.
