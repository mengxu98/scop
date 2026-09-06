# Cell Dimensional Plot

Color cells on a dimensionality reduction by metadata groups, with
optional labels, marks, density, lineages, PAGA, and RNA velocity
layers.

## Usage

``` r
CellDimPlot(
  srt,
  group.by,
  label.by = NULL,
  mark.by = NULL,
  legend.by = NULL,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  show_na = FALSE,
  show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey80",
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  label_insitu = FALSE,
  label_repel = FALSE,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black",
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  add_mark = FALSE,
  mark_type = c("hull", "ellipse", "rect", "circle"),
  mark_expand = grid::unit(3, "mm"),
  mark_alpha = 0.1,
  mark_linetype = 1,
  mark_linewidth = 0.5,
  mark_border = NULL,
  mark_palette = palette,
  mark_palcolor = NULL,
  add_grid = FALSE,
  grid_n = 12,
  grid_color = "black",
  grid_size = 0.25,
  grid_alpha = 0.35,
  lineages = NULL,
  lineages_trim = c(0.01, 0.99),
  lineages_span = 0.75,
  lineages_palette = "Dark2",
  lineages_palcolor = NULL,
  lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
  lineages_linewidth = 1,
  lineages_line_bg = "white",
  lineages_line_bg_stroke = 0.5,
  lineages_whiskers = FALSE,
  lineages_whiskers_linewidth = 0.5,
  lineages_whiskers_alpha = 0.5,
  stat.by = NULL,
  stat_type = "percent",
  stat_plot_type = "pie",
  stat_plot_position = c("stack", "dodge"),
  stat_plot_size = 0.2,
  stat_plot_palette = "Set1",
  stat_palcolor = NULL,
  stat_plot_alpha = 1,
  stat_plot_label = FALSE,
  stat_plot_label_size = 3,
  graph = NULL,
  edge_size = c(0.05, 0.5),
  edge_alpha = 0.1,
  edge_color = "grey40",
  paga = NULL,
  paga_type = "connectivities",
  paga_node_size = 4,
  paga_edge_threshold = 0.01,
  paga_edge_size = c(0.2, 1),
  paga_edge_color = "grey40",
  paga_edge_alpha = 0.5,
  paga_transition_threshold = 0.01,
  paga_transition_size = c(0.2, 1),
  paga_transition_color = "black",
  paga_transition_alpha = 1,
  paga_show_transition = FALSE,
  velocity = NULL,
  velocity_plot_type = "raw",
  velocity_n_neighbors = ceiling(ncol(srt@assays[[1]])/50),
  velocity_density = 1,
  velocity_smooth = 0.5,
  velocity_scale = 1,
  velocity_min_mass = 1,
  velocity_cutoff_perc = 5,
  velocity_arrow_color = "black",
  velocity_arrow_angle = 20,
  streamline_L = 5,
  streamline_minL = 1,
  streamline_res = 1,
  streamline_n = 15,
  streamline_width = c(0, 0.8),
  streamline_alpha = 1,
  streamline_color = NULL,
  streamline_palette = "RdYlBu",
  streamline_palcolor = NULL,
  streamline_bg_color = "white",
  streamline_bg_stroke = 0.5,
  hex = FALSE,
  hex.linewidth = 0.5,
  hex.count = TRUE,
  hex.bins = 50,
  hex.binwidth = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- label.by:

  Metadata column for labels. `NULL` uses `group.by`.

- mark.by:

  Metadata column for marks. `NULL` uses `legend.by` (nested legend) or
  `group.by`.

- legend.by:

  Parent group for a nested legend. `NULL` uses a standard legend.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Length-2 vector of dimensions to plot.

- split.by:

  Metadata column to facet by.

- cells:

  Cell names to include.

- show_na:

  If `TRUE`, color `NA` groups with `bg_color`; if `FALSE`, drop them.

- show_stat:

  Show cell-count statistics on the plot.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- bg_color:

  Color for background / `NA` points.

- label, label_insitu, label.size, label.fg, label.bg, label.bg.r:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label_repel, label_repulsion, label_point_size, label_point_color,
  label_segment_color:

  Repel labels away from cluster centers.

- cells.highlight, cols.highlight, sizes.highlight, alpha.highlight,
  stroke.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- add_density, density_color, density_filled, density_filled_palette,
  density_filled_palcolor:

  Density contours; `density_filled` draws filled bands.

- add_mark, mark_type, mark_expand, mark_alpha, mark_linetype,
  mark_linewidth, mark_border:

  Marks around groups. `mark_type` is `"hull"`, `"ellipse"`, `"rect"`,
  or `"circle"`. `mark_border = NULL` uses group colors.

- mark_palette, mark_palcolor:

  Palette for `mark.by` groups and nested-legend headers.

- add_grid, grid_n, grid_color, grid_size, grid_alpha:

  Background point grid (useful for atlas-style blank axes).

- lineages, lineages_trim, lineages_span, lineages_palette,
  lineages_palcolor, lineages_arrow, lineages_linewidth,
  lineages_line_bg, lineages_line_bg_stroke, lineages_whiskers,
  lineages_whiskers_linewidth, lineages_whiskers_alpha:

  Pseudotime lineages as
  [stats::loess](https://rdrr.io/r/stats/loess.html) curves. See
  [grid::arrow](https://rdrr.io/r/grid/arrow.html) for `lineages_arrow`.

- stat.by, stat_type, stat_plot_type, stat_plot_size, stat_plot_palette,
  stat_palcolor, stat_plot_position, stat_plot_alpha, stat_plot_label,
  stat_plot_label_size:

  Inset composition plot. `stat_type` is `"percent"` or `"count"`.

- graph, edge_size, edge_alpha, edge_color:

  Neighbor-graph edges.

- paga, paga_type, paga_node_size, paga_edge_threshold, paga_edge_size,
  paga_edge_color, paga_edge_alpha, paga_show_transition,
  paga_transition_threshold, paga_transition_size,
  paga_transition_color, paga_transition_alpha:

  PAGA graph layer. `paga_type` is `"connectivities"` or
  `"connectivities_tree"`.

- velocity, velocity_plot_type, velocity_n_neighbors, velocity_density,
  velocity_smooth, velocity_scale, velocity_min_mass,
  velocity_cutoff_perc, velocity_arrow_color, velocity_arrow_angle:

  RNA velocity layer.

- streamline_L, streamline_minL, streamline_res, streamline_n,
  streamline_width, streamline_alpha, streamline_color,
  streamline_palette, streamline_palcolor, streamline_bg_color,
  streamline_bg_stroke:

  Streamline appearance for velocity plots.

- hex, hex.count, hex.bins, hex.binwidth, hex.linewidth:

  Hexagonal bins instead of points.

- raster, raster.dpi:

  Rasterize points. `raster = NULL` rasterizes when there are more than
  100,000 cells.

- aspect.ratio:

  Panel aspect ratio.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend.position, legend.direction, legend.title:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- force:

  Draw even when a grouping has more than 100 levels.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## Value

A `ggplot`, `patchwork`, or list of `ggplot` objects.

## See also

[CellDimPlot3D](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[FeatureDimPlot3D](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:01:17] Start standard processing workflow...
#> ℹ [2026-09-06 21:01:18] Checking a list of <Seurat>...
#> ! [2026-09-06 21:01:18] Data 1/1 of the `srt_list` is "unknown"
#> Warning: Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:01:18] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:01:18] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:01:18] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:01:18] Number of available HVF: 2000
#> ℹ [2026-09-06 21:01:18] Finished check
#> ℹ [2026-09-06 21:01:18] Perform `ScaleData()`
#> ℹ [2026-09-06 21:01:18] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:01:18] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:01:19] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:01:19] Reorder clusters...
#> ℹ [2026-09-06 21:01:19] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:01:19] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:01:24] Standard processing workflow completed
p1 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
p1


thisplot::panel_fix(
  p1,
  height = 2,
  raster = TRUE,
  dpi = 300
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  theme_use = "theme_blank"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  theme_use = ggplot2::theme_classic,
  theme_args = list(base_size = 16)
)


# Highlight cells
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  cells.highlight = colnames(
    pancreas_sub
  )[pancreas_sub$SubCellType == "Epsilon"]
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  split.by = "Phase",
  reduction = "UMAP",
  cells.highlight = TRUE,
  theme_use = "theme_blank",
  legend.position = "none"
)


# Add group labels
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label.fg = "orange",
  label.bg = "red",
  label.size = 5
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  label_segment_color = "red"
)


# Add various shape of marks
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_expand = grid::unit(1, "mm")
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_alpha = 0.3
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_linetype = 2
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "ellipse"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "rect"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "circle"
)


# Add a density layer
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE,
  density_filled = TRUE
)
#> Warning: Removed 396 rows containing missing values or values outside the scale range
#> (`geom_raster()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE,
  density_filled = TRUE,
  density_filled_palette = "Blues",
  cells.highlight = TRUE
)
#> Warning: Removed 396 rows containing missing values or values outside the scale range
#> (`geom_raster()`).


# Add statistical charts
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase"
)


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase",
  stat_plot_type = "ring",
  stat_plot_label = TRUE,
  stat_plot_size = 0.15
)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase",
  stat_plot_type = "bar",
  stat_type = "count",
  stat_plot_position = "dodge"
)


# Chane the plot type from point to the hexagonal bin
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE,
  hex.bins = 20
)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE,
  hex.count = FALSE
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


# Show neighbors graphs on the plot
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  graph = "Standardpca_SNN"
)


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  graph = "Standardpca_SNN",
  edge_color = "grey80"
)


# Show lineages based on the pseudotime
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  show_plot = FALSE
)

FeatureDimPlot(
  pancreas_sub,
  features = paste0("Lineage", 1:2),
  reduction = "UMAP"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2)
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_whiskers = TRUE
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1
)


# Show PAGA results on the plot
pancreas_sub <- RunPAGA(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  backend = "cpp",
  return_seurat = TRUE
)
#> ℹ [2026-09-06 21:01:42] Running PAGA with BiocNeighbors using 29 neighbors
#> ✔ [2026-09-06 21:01:42] PAGA cpp backend completed

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@tools[["PAGA"]]
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@tools[["PAGA"]],
  paga_type = "connectivities_tree"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  label = TRUE,
  label_repel = TRUE,
  label_insitu = TRUE,
  label_segment_color = "transparent",
  paga = pancreas_sub@tools[["PAGA"]],
  paga_edge_threshold = 0.1,
  paga_edge_color = "black",
  paga_edge_alpha = 1,
  legend.position = "none",
  theme_use = "theme_blank"
)


# Show RNA velocity results on the plot
pancreas_sub <- RunSCVELO(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  mode = "stochastic",
  backend = "cpp",
  show_plot = FALSE,
  return_seurat = TRUE
)
#> ℹ [2026-09-06 21:01:43] Running scanpy-compatible preprocessing (15998 features -> filter + normalize)...
#> ℹ [2026-09-06 21:01:46] Running scVelo "stochastic" mode with `backend = 'cpp'` (9699 features)
#> ✔ [2026-09-06 21:01:49] scVelo "stochastic" mode completed
#> ✔ [2026-09-06 21:01:49] scVelo cpp backend completed

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  velocity = "stochastic"
)
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_segment()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "grid"
)
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_segment()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "grid",
  velocity_scale = 1.5
)
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_segment()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "stream"
)


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  label = TRUE,
  label_insitu = TRUE,
  velocity = "stochastic",
  velocity_plot_type = "stream",
  velocity_arrow_color = "yellow",
  velocity_density = 2,
  velocity_smooth = 1,
  streamline_n = 20,
  streamline_color = "black",
  legend.position = "none",
  theme_use = "theme_blank"
)
```
