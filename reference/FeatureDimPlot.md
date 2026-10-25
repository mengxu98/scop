# Feature values on a 2D reduction

Color cells on a dimensionality reduction by feature values (assay genes
or numeric metadata).

## Usage

``` r
FeatureDimPlot(
  object,
  features,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
  palette = ifelse(isTRUE(compare_features), "Set1", "Spectral"),
  palcolor = NULL,
  pt.size = NULL,
  pt.alpha = 1,
  bg_cutoff = 0,
  bg_color = "grey80",
  keep_scale = "feature",
  lower_quantile = 0,
  upper_quantile = 0.99,
  lower_cutoff = NULL,
  upper_cutoff = NULL,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  calculate_coexp = FALSE,
  compare_features = FALSE,
  color_blend_mode = c("blend", "average", "screen", "multiply"),
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  label_insitu = TRUE,
  label_repel = FALSE,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black",
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
  graph = NULL,
  edge_size = c(0.05, 0.5),
  edge_alpha = 0.1,
  edge_color = "grey40",
  hex = FALSE,
  hex.linewidth = 0.5,
  hex.color = "grey90",
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
  verbose = TRUE,
  srt = NULL,
  title.face = NULL,
  title.color = NULL,
  show_points = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- features:

  Features to plot: a character vector or a named list of assay gene
  names or numeric metadata columns.

- reduction:

  Dimensionality reduction to use. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Length-2 vector of dimensions to plot.

- split.by:

  Metadata column to split the analysis or plot by.

- cells:

  Cell names to include.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- show_stat:

  Show cell-count statistics on the plot.

- palette:

  Color palette name.

- palcolor:

  Custom colors used to create a color palette.

- pt.size:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- bg_cutoff:

  Values below this cutoff are colored with `bg_color`.

- bg_color:

  Color for background / `NA` points.

- keep_scale:

  Color scale across panels:

  - `NULL`: each panel is scaled independently (not comparable).

  - `"feature"`: scale each feature across `split.by` panels.

  - `"all"`: one scale for every feature and panel.

- lower_quantile, upper_quantile, lower_cutoff, upper_cutoff:

  Per-feature quantile or absolute cutoffs for the color scale.

- add_density, density_color, density_filled, density_filled_palette,
  density_filled_palcolor:

  Add a smooth feature-expression layer. For unblended features, density
  contours are colored by local Gaussian-weighted mean feature values,
  including zero values. Contour shape represents cell density; color
  represents local mean expression and may vary along one contour.
  `density_filled = TRUE` draws the local mean expression surface with
  fine density contours instead. With points visible, the original
  point-expression scale is preserved and local means use that scale.
  Without points, the scale is based on the local means. Cutoffs and
  `keep_scale` apply to the selected scale. `bg_cutoff` only affects
  cell marks. When neither fill-color option is supplied, expression
  surfaces inherit `palette` and `palcolor`; explicit fill options
  override this default. For compatibility, explicitly supplying
  `density_color` retains the old fixed-color cell-density layer.
  Blended multi-feature plots using `compare_features = TRUE` also
  retain the original cell-density layer, including filled density bands
  and their original default colors. Use
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  for cell-distribution density layers.

- cells.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- cols.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- sizes.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- alpha.highlight, stroke.highlight:

  Transparency and stroke width of the highlighted points.

- calculate_coexp:

  Plot the geometric mean of the features.

- compare_features:

  Show multiple features on one plot.

- color_blend_mode:

  Blend mode when `compare_features = TRUE`.

- label:

  Label high-expressing cells.

- label.size:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.fg:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.bg:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.bg.r:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label_insitu:

  Use feature names instead of numbers when `compare_features = TRUE`.

- label_repel, label_repulsion, label_point_size, label_point_color,
  label_segment_color:

  Repel labels away from cluster centers.

- lineages:

  Pseudotime columns to use. If `NULL`, lineage-like pseudotime columns
  such as `Lineage1`, `Lineage2`, `prefix_Lineage1`, or `pseudotime` are
  detected and merged into one global pseudotime for a single panel. Use
  `"all"` to plot each detected lineage in separate panels.

- lineages_trim, lineages_span, lineages_palcolor, lineages_arrow,
  lineages_linewidth, lineages_line_bg, lineages_line_bg_stroke,
  lineages_whiskers, lineages_whiskers_linewidth,
  lineages_whiskers_alpha:

  Pseudotime lineages as
  [stats::loess](https://rdrr.io/r/stats/loess.html) curves. See
  [grid::arrow](https://rdrr.io/r/grid/arrow.html) for `lineages_arrow`.

- lineages_palette:

  Color palette used for lineage groups.

- graph:

  Neighbor-graph edges.

- edge_size, edge_alpha, edge_color:

  Neighbor-graph edges.

- hex:

  Whether to use hexagonal binning for 2D assignments.

- hex.linewidth:

  Line width of six-point hexagons.

- hex.color:

  Border color of hexagonal bins.

- hex.bins:

  Number of hexagonal bins.

- hex.binwidth:

  Width of hexagonal bins.

- raster, raster.dpi:

  Rasterize points. `raster = NULL` rasterizes when there are more than
  100,000 cells.

- aspect.ratio:

  Panel aspect ratio.

- title:

  Plot title. `NULL` hides the title for merged/single panels. When
  multiple lineages are plotted and `title` is `NULL`, each panel is
  titled with its lineage column.

- subtitle:

  Plot subtitle.

- xlab:

  Plot labels.

- ylab:

  Plot labels.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- legend.title:

  Title of the legend. `NULL` uses the group name.

- theme_use:

  Theme name or function.

- theme_args:

  Theme name or function, plus extra theme arguments.

- combine:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- nrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- ncol:

  Number of columns of the combined plot.

- byrow:

  Whether to fill matrices/plots by row.

- force:

  Draw even when more than 100 features are requested.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- title.face:

  Font face of individual feature headings. `NULL` (default) uses italic
  for assay genes and plain for numeric metadata or embeddings. Set to
  `"plain"`, `"bold"`, `"italic"`, or `"bold.italic"` to override.

- title.color:

  Color of individual feature headings: `NULL` inherits the theme, an
  unnamed color applies to every feature, or a named color vector maps
  feature names to colors. Unspecified features inherit the theme. Names
  refer to features, not group names or panel order. These heading
  options do not affect the custom `title`, split labels, expression
  palette, or blended panels (`compare_features = TRUE`).

- show_points:

  Draw cell points or hexagonal bins and highlighted cell points. Set to
  `FALSE` for contour-only or surface-only plots. Labels and graph or
  lineage overlays are controlled separately. Defaults to `TRUE`.

## Details

Expression-colored contours use a 240 by 240 grid, Gaussian bandwidths
from [`stats::bw.nrd()`](https://rdrr.io/r/stats/bandwidth.html), and
nine density levels from 8% to 92% of peak density. Major contours use a
wider stroke every third level. Grid limits extend three bandwidths
beyond the finite coordinates to avoid clipping outer contours. Local
means are interpolated onto contour vertices; filled surfaces omit grid
points below 2% of peak kernel support, fade in between 2% and 3.5%
support, and interpolate raster pixels for smooth edges. Cells with
non-finite coordinates or feature values are omitted from this
calculation. Panels with fewer than three finite cells or no variation
in either coordinate return no contours with a warning. No additional
public plotting function or dependency is required.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-10-25 20:54:23] Start standard processing workflow...
#> ℹ [2026-10-25 20:54:23] Checking a list of <Seurat>...
#> ! [2026-10-25 20:54:23] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 20:54:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-25 20:54:23] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-25 20:54:23] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 20:54:23] Number of available HVF: 2000
#> ℹ [2026-10-25 20:54:23] Finished check
#> ℹ [2026-10-25 20:54:23] Perform `ScaleData()`
#> ℹ [2026-10-25 20:54:23] Perform pca linear dimension reduction
#> ℹ [2026-10-25 20:54:24] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-25 20:54:24] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-25 20:54:24] Reorder clusters...
#> ℹ [2026-10-25 20:54:24] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 20:54:24] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-25 20:54:31] Standard processing workflow completed
FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score", reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  bg_cutoff = -Inf
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  theme_use = "theme_blank"
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  theme_use = ggplot2::theme_classic,
  theme_args = list(base_size = 16)
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP"
) |> thisplot::panel_fix(
  height = 2,
  raster = TRUE,
  dpi = 30
)


# Label and highlight cell points
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  label = TRUE,
  cells.highlight = colnames(
    pancreas_sub
  )[pancreas_sub$SubCellType == "Delta"]
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  split.by = "Phase",
  reduction = "UMAP",
  cells.highlight = TRUE,
  theme_use = "theme_blank"
)


FeatureDimPlot(
  pancreas_sub, features = c("Rbp4", "Sst"), reduction = "UMAP",
  title.color = c(Rbp4 = "#2A6F97", Sst = "#B14B28")
)


# Add a density layer
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  add_density = TRUE,
  show_points = FALSE,
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  palcolor = c("#345C96", "#349AA4", "#91C6B0", "#DDBF63",
               "#D77A52", "#AC4968", "#753858")
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  label = TRUE,
  add_density = TRUE,
  density_filled = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  add_density = TRUE,
  show_points = FALSE,
  density_filled = TRUE,
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  density_filled_palcolor = c("#F2F7F5", "#C5DFD8", "#89BDB5",
                             "#4A969E", "#286B87", "#253D63")
)


# Chane the plot type from point to the hexagonal bin
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  hex = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  hex = TRUE,
  hex.bins = 20
)


# Show lineages on the plot based on the pseudotime
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2"
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2",
  lineages_whiskers = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2",
  lineages_span = 0.1
)


# Input a named feature list
markers <- list(
  "Ductal" = c("Sox9", "Anxa2", "Bicc1"),
  "EPs" = c("Neurog3", "Hes6"),
  "Pre-endocrine" = c("Fev", "Neurod1"),
  "Endocrine" = c("Rbp4", "Pyy"),
  "Beta" = "Ins1",
  "Alpha" = "Gcg",
  "Delta" = "Sst",
  "Epsilon" = "Ghrl"
)
FeatureDimPlot(
  pancreas_sub,
  features = markers,
  reduction = "UMAP",
  theme_use = "theme_blank",
  theme_args = list(
    plot.subtitle = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 8)
  )
)


# Plot multiple features with different scales
endocrine_markers <- c(
  "Beta" = "Ins1",
  "Alpha" = "Gcg",
  "Delta" = "Sst",
  "Epsilon" = "Ghrl"
)
FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  lower_quantile = 0,
  upper_quantile = 0.8
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  lower_cutoff = 1,
  upper_cutoff = 4
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  keep_scale = "all"
)


FeatureDimPlot(
  pancreas_sub,
  c("Delta" = "Sst", "Epsilon" = "Ghrl"),
  split.by = "Phase",
  reduction = "UMAP",
  keep_scale = "feature"
)


# Plot multiple features on one picture
FeatureDimPlot(
  pancreas_sub,
  features = endocrine_markers,
  pt.size = 1,
  compare_features = TRUE,
  color_blend_mode = "blend",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "blend",
  title = "blend",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "average",
  title = "average",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "screen",
  title = "screen",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "multiply",
  title = "multiply",
  label = TRUE
)
```
