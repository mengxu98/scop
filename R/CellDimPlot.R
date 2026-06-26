#' @title Cell Dimensional Plot
#'
#' @description
#' Visualize cell groups on a 2-dimensional reduction plot.
#' Plotting cell points on a reduced 2D plane and coloring according to the groups.
#'
#' @md
#' @inheritParams standard_scop
#' @param group.by Name of one or more meta.data columns to group (color) cells by.
#' @param label.by Name of a meta.data column used to place group labels. If
#' `NULL`, labels use `group.by`.
#' @param mark.by Name of a meta.data column used to draw group marks. If
#' `NULL`, marks use `legend.by` when a nested legend is requested, otherwise
#' marks use `group.by`.
#' @param legend.by Name of a meta.data column used as the parent group in a
#' nested legend. The `group.by` levels are shown under each `legend.by` level.
#' If `NULL`, the standard legend is used.
#' @param reduction Which dimensionality reduction to use.
#' If not specified, will use the reduction returned by [DefaultReduction].
#' @param split.by Name of a column in meta.data column to split plot by.
#' Default is `NULL`.
#' @param palette Color palette name.
#' Available palettes can be found in [thisplot::show_palettes].
#' Default is `"Chinese"`.
#' @param palcolor Custom colors used to create a color palette.
#' Default is `NULL`.
#' @param bg_color Color value for background(NA) points.
#' @param pt.size The size of the points in the plot.
#' @param pt.alpha The transparency of the data points.
#' Default is `1`.
#' @param cells.highlight A logical or character vector specifying the cells to highlight in the plot.
#' If `TRUE`, all cells are highlighted. If `FALSE`, no cells are highlighted.
#' Default is `NULL`.
#' @param cols.highlight Color used to highlight the cells.
#' @param sizes.highlight Size of highlighted cell points.
#' @param alpha.highlight Transparency of highlighted cell points.
#' @param stroke.highlight Border width of highlighted cell points.
#' @param legend.position The position of legends,
#' one of `"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`.
#' Default is `"right"`.
#' @param legend.direction The direction of the legend in the plot.
#' Can be one of `"vertical"` or `"horizontal"`.
#' @param legend.title Title for the legend. Default is `NULL`, which uses the group name.
#' @param combine Combine plots into a single `patchwork` object.
#' If `FALSE`, return a list of ggplot objects.
#' @param nrow Number of rows in the combined plot.
#' Default is `NULL`, which means determined automatically based on the number of plots.
#' @param ncol Number of columns in the combined plot.
#' Default is `NULL`, which means determined automatically based on the number of plots.
#' @param byrow Whether to arrange the plots by row in the combined plot.
#' Default is `TRUE`.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param show_na Whether to assign a color from the color palette to NA group.
#' If `TRUE`, cell points with NA level will be colored by `bg_color`.
#' If `FALSE`, cell points with NA level will be removed from the plot.
#' @param show_stat Whether to show statistical information on the plot.
#' @param label Whether to label the cell groups.
#' @param label_insitu Whether to place the raw labels (group names) in the center of the cells with the corresponding group.
#' Default is `FALSE`, which using numbers instead of raw labels.
#' @param label.size Size of labels.
#' @param label.fg Foreground color of label.
#' @param label.bg Background color of label.
#' @param label.bg.r Background ratio of label.
#' @param label_repel Logical value indicating whether the label is repel away from the center points.
#' @param label_repulsion Force of repulsion between overlapping text labels. Default is `20`.
#' @param label_point_size Size of the center points.
#' @param label_point_color Color of the center points.
#' @param label_segment_color Color of the line segment for labels.
#' @param add_density Whether to add a density layer on the plot.
#' @param density_color Color of the density contours lines.
#' @param density_filled Whether to add filled contour bands instead of contour lines.
#' @param density_filled_palette Color palette used to fill contour bands.
#' @param density_filled_palcolor Custom colors used to fill contour bands.
#' @param add_mark Whether to add marks around cell groups. Default is `FALSE`.
#' @param mark_type Type of mark to add around cell groups.
#' One of "hull", "ellipse", "rect", or "circle". Default is `"hull"`.
#' @param mark_expand Expansion of the mark around the cell group.
#' Default is `grid::unit(3, "mm")`.
#' @param mark_alpha Transparency of the mark.
#' Default is `0.1`.
#' @param mark_linetype Line type of the mark border.
#' Default is `1` (solid line).
#' @param mark_linewidth Line width of the mark border.
#' @param mark_border Fixed border color for marks. If `NULL`, mark borders use
#' the group colors.
#' @param mark_palette Color palette name for `mark.by` groups and nested legend
#' parent headers. Defaults to `palette`.
#' @param mark_palcolor Custom colors for `mark.by` groups and nested legend
#' parent headers.
#' @param add_grid Whether to add a background point grid on the reduction panel.
#' This is useful for atlas-style panels with blank axes.
#' @param grid_n Number of grid points along each axis when `add_grid = TRUE`.
#' @param grid_color,grid_size,grid_alpha Color, size, and alpha of the
#' background grid points.
#' @param lineages Lineages/pseudotime to add to the plot.
#' If specified, curves will be fitted using [stats::loess] method.
#' @param lineages_trim Trim the leading and the trailing data in the lineages.
#' @param lineages_span The parameter α which controls the degree of smoothing in [stats::loess] method.
#' @param lineages_palette Color palette used for lineages.
#' @param lineages_palcolor Custom colors used for lineages.
#' @param lineages_arrow Set arrows of the lineages. See [grid::arrow].
#' @param lineages_linewidth Width of fitted curve lines for lineages.
#' @param lineages_line_bg Background color of curve lines for lineages.
#' @param lineages_line_bg_stroke Border width of curve lines background.
#' @param lineages_whiskers Whether to add whiskers for lineages.
#' @param lineages_whiskers_linewidth Width of whiskers for lineages.
#' @param lineages_whiskers_alpha Transparency of whiskers for lineages.
#' @param stat.by The name of a metadata column to stat.
#' @param stat_type Set stat types ("percent" or "count").
#' @param stat_plot_type Set the statistical plot type.
#' @param stat_plot_size Set the statistical plot size. Default is `0.2`.
#' @param stat_plot_palette Color palette used in statistical plot.
#' @param stat_palcolor Custom colors used in statistical plot
#' @param stat_plot_position Position adjustment in statistical plot.
#' @param stat_plot_alpha Transparency of the statistical plot.
#' @param stat_plot_label Whether to add labels in the statistical plot.
#' @param stat_plot_label_size Label size in the statistical plot.
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param edge_size Size of edges.
#' @param edge_alpha Transparency of edges.
#' @param edge_color Color of edges.
#' @param paga Specify the calculated paga results to add a PAGA graph layer to the plot.
#' @param paga_type PAGA plot type. "connectivities" or "connectivities_tree".
#' @param paga_node_size Size of the nodes in PAGA plot.
#' @param paga_edge_threshold Threshold of edge connectivities in PAGA plot.
#' @param paga_edge_size Size of edges in PAGA plot.
#' @param paga_edge_color Color of edges in PAGA plot.
#' @param paga_edge_alpha Transparency of edges in PAGA plot.
#' @param paga_show_transition Whether to show transitions between edges.
#' @param paga_transition_threshold Threshold of transition edges in PAGA plot.
#' @param paga_transition_size Size of transition edges in PAGA plot.
#' @param paga_transition_color Color of transition edges in PAGA plot.
#' @param paga_transition_alpha Transparency of transition edges in PAGA plot.
#' @param velocity Specify the calculated RNA velocity mode to add a velocity layer to the plot.
#' @param velocity_plot_type Set the velocity plot type.
#' @param velocity_n_neighbors Set the number of neighbors used in velocity plot.
#' @param velocity_density Set the density value used in velocity plot.
#' @param velocity_smooth Set the smooth value used in velocity plot.
#' @param velocity_scale Set the scale value used in velocity plot.
#' @param velocity_min_mass Set the min_mass value used in velocity plot.
#' @param velocity_cutoff_perc Set the cutoff_perc value used in velocity plot.
#' @param velocity_arrow_color Color of arrows in velocity plot.
#' @param velocity_arrow_angle Angle of arrows in velocity plot.
#' @param streamline_L Typical length of a streamline in x and y units
#' @param streamline_minL Minimum length of segments to show.
#' @param streamline_res Resolution parameter (higher numbers increases the resolution).
#' @param streamline_n Number of points to draw.
#' @param streamline_width Size of streamline.
#' @param streamline_alpha Transparency of streamline.
#' @param streamline_color Color of streamline.
#' @param streamline_palette Color palette used for streamline.
#' @param streamline_palcolor Custom colors used for streamline.
#' @param streamline_bg_color Background color of streamline.
#' @param streamline_bg_stroke Border width of streamline background.
#' @param hex Whether to chane the plot type from point to the hexagonal bin.
#' @param hex.count Whether show cell counts in each hexagonal bin.
#' @param hex.bins Number of hexagonal bins.
#' @param hex.binwidth Hexagonal bin width.
#' @param hex.linewidth Border width of hexagonal bins.
#' @param raster Convert points to raster format.
#' Default is `NULL`, which automatically rasterizes if plotting more than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots.
#' Default is `c(512, 512)`.
#' @param theme_use Theme used. Can be a character string or a theme function.
#' Default is `"theme_scop"`.
#' @param aspect.ratio Aspect ratio of the panel. Default is `1`.
#' @param title The text for the title.
#' Default is `NULL`.
#' @param subtitle The text for the subtitle for the plot which will be displayed below the title.
#' Default is `NULL`.
#' @param xlab The x-axis label of the plot.
#' Default is `NULL`.
#' @param ylab The y-axis label of the plot.
#' Default is `NULL`.
#' @param force Whether to force drawing regardless of maximum levels in any cell group is greater than 100.
#' Default is `FALSE`.
#' @param cells A character vector of cell names to use.
#' @param theme_args Other arguments passed to the `theme_use`.
#' Default is `list()`.
#'
#' @seealso
#' [CellDimPlot3D], [FeatureDimPlot], [FeatureDimPlot3D]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' p1 <- CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' p1
#'
#' thisplot::panel_fix(
#'   p1,
#'   height = 2,
#'   raster = TRUE,
#'   dpi = 300
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   theme_use = "theme_blank"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   theme_use = ggplot2::theme_classic,
#'   theme_args = list(base_size = 16)
#' )
#'
#' # Highlight cells
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   cells.highlight = colnames(
#'     pancreas_sub
#'   )[pancreas_sub$SubCellType == "Epsilon"]
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   reduction = "UMAP",
#'   cells.highlight = TRUE,
#'   theme_use = "theme_blank",
#'   legend.position = "none"
#' )
#'
#' # Add group labels
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE,
#'   label.fg = "orange",
#'   label.bg = "red",
#'   label.size = 5
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE,
#'   label_insitu = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE,
#'   label_insitu = TRUE,
#'   label_repel = TRUE,
#'   label_segment_color = "red"
#' )
#'
#' # Add various shape of marks
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_expand = grid::unit(1, "mm")
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_alpha = 0.3
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_linetype = 2
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_type = "ellipse"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_type = "rect"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_type = "circle"
#' )
#'
#' # Add a density layer
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_density = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_density = TRUE,
#'   density_filled = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_density = TRUE,
#'   density_filled = TRUE,
#'   density_filled_palette = "Blues",
#'   cells.highlight = TRUE
#' )
#'
#' # Add statistical charts
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   stat.by = "Phase"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   stat.by = "Phase",
#'   stat_plot_type = "ring",
#'   stat_plot_label = TRUE,
#'   stat_plot_size = 0.15
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   stat.by = "Phase",
#'   stat_plot_type = "bar",
#'   stat_type = "count",
#'   stat_plot_position = "dodge"
#' )
#'
#' # Chane the plot type from point to the hexagonal bin
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   hex = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   hex = TRUE,
#'   hex.bins = 20
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   hex = TRUE,
#'   hex.count = FALSE
#' )
#'
#' # Show neighbors graphs on the plot
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   graph = "Standardpca_SNN"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   graph = "Standardpca_SNN",
#'   edge_color = "grey80"
#' )
#'
#' # Show lineages based on the pseudotime
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   show_plot = FALSE
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = paste0("Lineage", 1:2),
#'   reduction = "UMAP"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2)
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_whiskers = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_span = 0.1
#' )
#'
#' # Show PAGA results on the plot
#' pancreas_sub <- RunPAGA(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   backend = "cpp",
#'   return_seurat = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga,
#'   paga_type = "connectivities_tree"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   label = TRUE,
#'   label_repel = TRUE,
#'   label_insitu = TRUE,
#'   label_segment_color = "transparent",
#'   paga = pancreas_sub@misc$paga,
#'   paga_edge_threshold = 0.1,
#'   paga_edge_color = "black",
#'   paga_edge_alpha = 1,
#'   legend.position = "none",
#'   theme_use = "theme_blank"
#' )
#'
#' # Show RNA velocity results on the plot
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   mode = "stochastic",
#'   backend = "cpp",
#'   return_seurat = TRUE
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   velocity = "stochastic"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   velocity = "stochastic",
#'   velocity_plot_type = "grid"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   velocity = "stochastic",
#'   velocity_plot_type = "grid",
#'   velocity_scale = 1.5
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   velocity = "stochastic",
#'   velocity_plot_type = "stream"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   label = TRUE,
#'   label_insitu = TRUE,
#'   velocity = "stochastic",
#'   velocity_plot_type = "stream",
#'   velocity_arrow_color = "yellow",
#'   velocity_density = 2,
#'   velocity_smooth = 1,
#'   streamline_n = 20,
#'   streamline_color = "black",
#'   legend.position = "none",
#'   theme_use = "theme_blank"
#' )
CellDimPlot <- function(
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
  velocity_n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50),
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
) {
  set.seed(seed)
  mark_type <- match.arg(mark_type)
  if (!is.null(label.by) && length(label.by) != 1L) {
    log_message(
      "{.arg label.by} must be a single meta.data column name",
      message_type = "error"
    )
  }
  if (!is.null(mark.by) && length(mark.by) != 1L) {
    log_message(
      "{.arg mark.by} must be a single meta.data column name",
      message_type = "error"
    )
  }
  if (!is.null(legend.by) && length(legend.by) != 1L) {
    log_message(
      "{.arg legend.by} must be a single meta.data column name",
      message_type = "error"
    )
  }

  check_r("ggnewscale", verbose = FALSE)
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in unique(c(group.by, label.by, mark.by, legend.by, split.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      log_message(
        "{.val {i}} is not in the meta.data of srt object",
        message_type = "error"
      )
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(
        srt@meta.data[[i]],
        levels = unique(srt@meta.data[[i]])
      )
    }
    if (isTRUE(show_na) && any(is.na(srt@meta.data[[i]]))) {
      raw_levels <- unique(c(levels(srt@meta.data[[i]]), "NA"))
      srt@meta.data[[i]] <- as.character(srt@meta.data[[i]])
      srt@meta.data[[i]][is.na(srt@meta.data[[i]])] <- "NA"
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = raw_levels)
    }
  }
  for (l in lineages) {
    if (!l %in% colnames(srt@meta.data)) {
      log_message(
        "Lineage {.val {l}} is not in the meta.data of srt object",
        message_type = "error"
      )
    }
  }
  if (!is.null(graph) && !graph %in% names(srt@graphs)) {
    log_message(
      "Graph {.val {graph}} is not exist in the srt object",
      message_type = "error"
    )
  }
  if (!is.null(graph)) {
    graph <- srt@graphs[[graph]]
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }
  if (!is.null(cells.highlight) && isFALSE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "No cells in '{.arg cells.highlight}' found in srt",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "Some cells in '{.arg cells.highlight}' not found in srt",
        message_type = "warning",
        verbose = verbose
      )
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }

  dat_meta <- srt@meta.data[
    ,
    unique(c(group.by, label.by, mark.by, legend.by, split.by)),
    drop = FALSE
  ]
  nlev <- sapply(dat_meta, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && isFALSE(force)) {
    log_message(
      "The following variables have more than 100 levels: {.val {names(nlev)}}",
      message_type = "warning",
      verbose = verbose
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(invisible(NULL))
    }
  }

  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_use <- cbind(dat_dim, dat_meta[row.names(dat_dim), , drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) > 1e5)
  if (isTRUE(raster)) {
    check_r("scattermore", verbose = FALSE)
  }
  if (!is.null(raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(raster.dpi) != 2) {
      log_message(
        "'{.arg raster.dpi}' must be a two-length numeric vector",
        message_type = "error"
      )
    }
  }
  if (!is.null(stat.by)) {
    subplots <- CellStatPlot(
      srt,
      cells = cells,
      stat.by = stat.by,
      group.by = group.by,
      split.by = split.by,
      stat_type = stat_type,
      plot_type = stat_plot_type,
      position = stat_plot_position,
      palette = stat_plot_palette,
      palcolor = stat_palcolor,
      alpha = stat_plot_alpha,
      label = stat_plot_label,
      label.size = stat_plot_label_size,
      legend.position = "bottom",
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      individual = TRUE,
      combine = FALSE
    )
  }

  if (!is.null(lineages)) {
    lineages_layers <- LineagePlot(
      srt,
      cells = cells,
      lineages = lineages,
      reduction = reduction,
      dims = dims,
      trim = lineages_trim,
      span = lineages_span,
      palette = lineages_palette,
      palcolor = lineages_palcolor,
      lineages_arrow = lineages_arrow,
      linewidth = lineages_linewidth,
      line_bg = lineages_line_bg,
      line_bg_stroke = lineages_line_bg_stroke,
      whiskers = lineages_whiskers,
      whiskers_linewidth = lineages_whiskers_linewidth,
      whiskers_alpha = lineages_whiskers_alpha,
      aspect.ratio = aspect.ratio,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = "bottom",
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[
      !names(lineages_layers) %in% c("lab_layer", "theme_layer")
    ]
  }

  if (!is.null(paga)) {
    if (split.by != "All.groups") {
      log_message(
        "paga can only plot on the non-split data",
        message_type = "error"
      )
    }
    paga_layers <- PAGAPlot(
      srt,
      cells = cells,
      paga = paga,
      type = paga_type,
      reduction = reduction,
      dims = dims,
      node_palette = palette,
      node_palcolor = palcolor,
      node_size = paga_node_size,
      edge_threshold = paga_edge_threshold,
      edge_size = paga_edge_size,
      edge_color = paga_edge_color,
      edge_alpha = paga_edge_alpha,
      transition_threshold = paga_transition_threshold,
      transition_size = paga_transition_size,
      transition_color = paga_transition_color,
      transition_alpha = paga_transition_alpha,
      show_transition = paga_show_transition,
      aspect.ratio = aspect.ratio,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = "bottom",
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      return_layer = TRUE
    )
    paga_layers <- paga_layers[
      !names(paga_layers) %in% c("lab_layer", "theme_layer")
    ]
  }

  if (!is.null(velocity)) {
    if (split.by != "All.groups") {
      log_message(
        "velocity can only plot on the non-split data",
        message_type = "error"
      )
    }
    velocity_layers <- VelocityPlot(
      srt,
      cells = cells,
      reduction = reduction,
      dims = dims,
      velocity = velocity,
      plot_type = velocity_plot_type,
      group.by = group.by,
      group_palette = palette,
      group_palcolor = palcolor,
      n_neighbors = velocity_n_neighbors,
      density = velocity_density,
      smooth = velocity_smooth,
      scale = velocity_scale,
      min_mass = velocity_min_mass,
      cutoff_perc = velocity_cutoff_perc,
      arrow_color = velocity_arrow_color,
      arrow_angle = velocity_arrow_angle,
      streamline_L = streamline_L,
      streamline_minL = streamline_minL,
      streamline_res = streamline_res,
      streamline_n = streamline_n,
      streamline_width = streamline_width,
      streamline_alpha = streamline_alpha,
      streamline_color = streamline_color,
      streamline_palette = streamline_palette,
      streamline_palcolor = streamline_palcolor,
      streamline_bg_color = streamline_bg_color,
      streamline_bg_stroke = streamline_bg_stroke,
      aspect.ratio = aspect.ratio,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = "bottom",
      legend.direction = legend.direction,
      theme_use = theme_void,
      theme_args = theme_args,
      return_layer = TRUE
    )
    velocity_layers <- velocity_layers[
      !names(velocity_layers) %in% c("lab_layer", "theme_layer")
    ]
  }

  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  comb <- expand.grid(
    split = levels(dat_use[[split.by]]),
    group = group.by,
    stringsAsFactors = FALSE
  )
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["group"]])
  plist <- lapply(
    stats::setNames(rownames(comb), rownames(comb)), function(i) {
      g <- comb[i, "group"]
      s <- comb[i, "split"]
      legend_by_use <- legend.by
      label_by_use <- label.by %||% g
      mark_by_use <- mark.by %||% legend_by_use %||% g
      colors <- cell_dim_palette_colors(
        levels(dat_use[[g]]),
        palette = palette,
        palcolor = palcolor,
        NA_keep = TRUE
      )

      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      if (isFALSE(show_na)) {
        dat <- dat[!is.na(dat[[g]]), , drop = FALSE]
      }
      legend_list <- list()
      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
      cells_highlight_use <- cells.highlight
      if (isTRUE(cells_highlight_use)) {
        cells_highlight_use <- rownames(dat)[!is.na(dat[[g]])]
      }

      if (isTRUE(label_insitu)) {
        if (isTRUE(show_stat)) {
          label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
        } else {
          label_use <- paste0(names(labels_tb))
        }
      } else {
        if (isTRUE(label)) {
          if (isTRUE(show_stat)) {
            label_use <- paste0(
              seq_along(labels_tb),
              ": ",
              names(labels_tb),
              "(",
              labels_tb,
              ")"
            )
          } else {
            label_use <- paste0(seq_along(labels_tb), ": ", names(labels_tb))
          }
        } else {
          if (isTRUE(show_stat)) {
            label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
          } else {
            label_use <- paste0(names(labels_tb))
          }
        }
      }

      dat[["x"]] <- dat[[paste0(reduction_key, dims[1])]]
      dat[["y"]] <- dat[[paste0(reduction_key, dims[2])]]
      dat[["group.by"]] <- dat[[g]]
      dat[["label.by"]] <- dat[[label_by_use]]
      dat[["mark.by"]] <- dat[[mark_by_use]]
      if (!is.null(legend_by_use)) {
        dat[["legend.by"]] <- dat[[legend_by_use]]
      }
      dat[, "split.by"] <- s
      dat <- dat[
        order(dat[, "group.by"], decreasing = FALSE, na.last = FALSE), ,
        drop = FALSE
      ]
      naindex <- which(is.na(dat[, "group.by"]))
      naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
      dat <- dat[
        c(1:naindex, sample((min(naindex + 1, nrow(dat))):nrow(dat))), ,
        drop = FALSE
      ]
      if (isTRUE(show_stat)) {
        subtitle_use <- subtitle %||%
          paste0(s, " nCells:", sum(!is.na(dat[["group.by"]])))
      } else {
        subtitle_use <- subtitle
      }

      if (isTRUE(add_mark)) {
        mark_fun <- switch(
          EXPR = mark_type,
          "ellipse" = "geom_mark_ellipse",
          "hull" = "geom_mark_hull",
          "rect" = "geom_mark_rect",
          "circle" = "geom_mark_circle"
        )
        mark_mapping <- aes(
          x = .data[["x"]],
          y = .data[["y"]],
          fill = .data[["mark.by"]]
        )
        mark_colors <- cell_dim_palette_colors(
          levels(droplevels(dat[["mark.by"]])),
          palette = mark_palette,
          palcolor = mark_palcolor,
          NA_keep = TRUE
        )
        mark_scales <- list(
          scale_fill_manual(values = mark_colors),
          ggnewscale::new_scale_fill()
        )
        if (is.null(mark_border)) {
          mark_mapping <- aes(
            x = .data[["x"]],
            y = .data[["y"]],
            color = .data[["mark.by"]],
            fill = .data[["mark.by"]]
          )
          mark_scales <- c(
            list(scale_color_manual(values = mark_colors)),
            mark_scales,
            list(ggnewscale::new_scale_color())
          )
        }
        mark_args <- list(
          data = dat[!is.na(dat[["mark.by"]]), , drop = FALSE],
          mapping = mark_mapping,
          expand = mark_expand,
          alpha = mark_alpha,
          linetype = mark_linetype,
          linewidth = mark_linewidth,
          con.size = mark_linewidth,
          show.legend = FALSE,
          inherit.aes = FALSE
        )
        if (!is.null(mark_border)) {
          mark_args[["colour"]] <- mark_border
        }
        mark <- c(list(do.call(mark_fun, mark_args)), mark_scales)
      } else {
        mark <- NULL
      }

      if (!is.null(graph)) {
        net_mat <- as_matrix(graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- reshape2::melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
        net_df[, "value"] <- as.numeric(net_df[, "value"])
        net_df[, "Var1"] <- as.character(net_df[, "Var1"])
        net_df[, "Var2"] <- as.character(net_df[, "Var2"])
        net_df[, "x"] <- dat[net_df[, "Var1"], "x"]
        net_df[, "y"] <- dat[net_df[, "Var1"], "y"]
        net_df[, "xend"] <- dat[net_df[, "Var2"], "x"]
        net_df[, "yend"] <- dat[net_df[, "Var2"], "y"]
        net <- list(
          geom_segment(
            data = net_df,
            mapping = aes(
              x = x,
              y = y,
              xend = xend,
              yend = yend,
              linewidth = value
            ),
            color = edge_color,
            alpha = edge_alpha,
            show.legend = FALSE
          ),
          scale_linewidth_continuous(range = edge_size)
        )
      } else {
        net <- NULL
      }

      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_color <- palette_colors(
            palette = density_filled_palette,
            palcolor = density_filled_palcolor
          )
          density <- list(
            stat_density_2d(
              geom = "raster",
              aes(x = .data[["x"]], y = .data[["y"]], fill = after_stat(density)),
              contour = FALSE,
              inherit.aes = FALSE,
              show.legend = FALSE
            ),
            scale_fill_gradientn(name = "Density", colours = filled_color),
            ggnewscale::new_scale_fill()
          )
        } else {
          density <- geom_density_2d(
            aes(x = .data[["x"]], y = .data[["y"]]),
            color = density_color,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
        }
      } else {
        density <- NULL
      }

      grid_layer <- cell_dim_grid_layer(
        dat = dat,
        add_grid = add_grid,
        grid_n = grid_n,
        grid_color = grid_color,
        grid_size = grid_size,
        grid_alpha = grid_alpha
      )

      p <- ggplot(dat) +
        grid_layer +
        mark +
        net +
        density +
        labs(title = title, subtitle = subtitle_use, x = xlab, y = ylab) +
        scale_x_continuous(
          limits = c(
            min(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE),
            max(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE)
          )
        ) +
        scale_y_continuous(
          limits = c(
            min(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE),
            max(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE)
          )
        ) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )

      if (split.by != "All.groups") {
        p <- p + facet_grid(. ~ split.by)
      }

      if (isTRUE(raster)) {
        p <- p +
          scattermore::geom_scattermore(
            data = dat[is.na(dat[, "group.by"]), , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]]),
            color = bg_color,
            pointsize = ceiling(pt.size),
            alpha = pt.alpha,
            pixels = raster.dpi
          ) +
          scattermore::geom_scattermore(
            data = dat[!is.na(dat[, "group.by"]), , drop = FALSE],
            mapping = aes(
              x = .data[["x"]],
              y = .data[["y"]],
              color = .data[["group.by"]]
            ),
            pointsize = ceiling(pt.size),
            alpha = pt.alpha,
            pixels = raster.dpi
          )
      } else if (isTRUE(hex)) {
        check_r("hexbin", verbose = FALSE)
        if (isTRUE(hex.count)) {
          p <- p +
            geom_hex(
              mapping = aes(
                x = .data[["x"]],
                y = .data[["y"]],
                fill = .data[["group.by"]],
                color = .data[["group.by"]],
                alpha = after_stat(count)
              ),
              linewidth = hex.linewidth,
              bins = hex.bins,
              binwidth = hex.binwidth
            )
        } else {
          p <- p +
            geom_hex(
              mapping = aes(
                x = .data[["x"]],
                y = .data[["y"]],
                fill = .data[["group.by"]],
                color = .data[["group.by"]]
              ),
              linewidth = hex.linewidth,
              bins = hex.bins,
              binwidth = hex.binwidth
            )
        }
      } else {
        p <- p +
          geom_point(
            mapping = aes(
              x = .data[["x"]],
              y = .data[["y"]],
              color = .data[["group.by"]]
            ),
            shape = "circle",
            size = pt.size,
            alpha = pt.alpha
          )
      }

      if (!is.null(cells_highlight_use) && isFALSE(hex)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells_highlight_use)
        if (nrow(cell_df) > 0) {
          if (isTRUE(raster)) {
            p <- p +
              scattermore::geom_scattermore(
                data = cell_df,
                aes(x = .data[["x"]], y = .data[["y"]]),
                color = cols.highlight,
                pointsize = floor(sizes.highlight) + stroke.highlight,
                alpha = alpha.highlight,
                pixels = raster.dpi
              ) +
              scattermore::geom_scattermore(
                data = cell_df,
                aes(
                  x = .data[["x"]],
                  y = .data[["y"]],
                  color = .data[["group.by"]]
                ),
                pointsize = floor(sizes.highlight),
                alpha = alpha.highlight,
                pixels = raster.dpi
              )
          } else {
            p <- p +
              geom_point(
                data = cell_df,
                aes(x = .data[["x"]], y = .data[["y"]]),
                shape = "circle",
                color = cols.highlight,
                size = sizes.highlight + stroke.highlight,
                alpha = alpha.highlight
              ) +
              geom_point(
                data = cell_df,
                aes(
                  x = .data[["x"]],
                  y = .data[["y"]],
                  color = .data[["group.by"]]
                ),
                shape = "circle",
                size = sizes.highlight,
                alpha = alpha.highlight
              )
          }
        }
      }
      legend_title_use <- if (is.null(legend.title)) paste0(g, ":") else legend.title
      color_guide <- if (is.null(legend.by)) {
        guide_legend(
          title.hjust = 0,
          order = 1,
          override.aes = list(shape = "circle", size = 4, alpha = 1)
        )
      } else {
        "none"
      }
      p <- p +
        scale_color_manual(
          name = legend_title_use,
          values = colors[names(labels_tb)],
          labels = label_use,
          na.value = bg_color,
          guide = color_guide
        )
      if (isTRUE(hex)) {
        fill_guide <- if (is.null(legend.by)) {
          guide_legend(title.hjust = 0, order = 1)
        } else {
          "none"
        }
        p <- p +
          scale_fill_manual(
            name = legend_title_use,
            values = colors[names(labels_tb)],
            labels = label_use,
            na.value = bg_color,
            guide = fill_guide
          )
      }
      p_base <- p
      nested_legend <- NULL
      if (!is.null(legend.by)) {
        nested_legend <- cell_dim_nested_legend_grob(
          dat = dat[!is.na(dat[["group.by"]]) & !is.na(dat[["legend.by"]]), , drop = FALSE],
          colors = colors[names(labels_tb)],
          parent_colors = cell_dim_palette_colors(
            levels(droplevels(dat[["legend.by"]])),
            palette = mark_palette,
            palcolor = mark_palcolor,
            NA_keep = TRUE
          ),
          title = if (is.null(legend.title)) paste0(legend.by, ":") else legend.title,
          legend_direction = legend.direction,
          legend_position = legend.position
        )
      }

      if (!is.null(stat.by)) {
        coor_df <- stats::aggregate(
          p$data[, c("x", "y")],
          by = list(p$data[["group.by"]]),
          FUN = stats::median
        )
        colnames(coor_df)[1] <- "group"
        x_range <- diff(layer_scales(p)$x$range$range)
        y_range <- diff(layer_scales(p)$y$range$range)
        stat_plot <- subplots[paste0(g, ":", levels(dat[, "group.by"]), ":", s)]
        names(stat_plot) <- levels(dat[, "group.by"])

        stat_plot_list <- list()
        for (i in seq_len(nrow(coor_df))) {
          stat_plot_list[[i]] <- ggplot2::annotation_custom(
            as_grob(
              stat_plot[[coor_df[i, "group"]]] +
                theme_void() +
                theme(legend.position = "none")
            ),
            xmin = coor_df[i, "x"] - x_range * stat_plot_size / 2,
            ymin = coor_df[i, "y"] - y_range * stat_plot_size / 2,
            xmax = coor_df[i, "x"] + x_range * stat_plot_size / 2,
            ymax = coor_df[i, "y"] + y_range * stat_plot_size / 2
          )
        }
        p <- p + stat_plot_list
        legend_list[["stat.by"]] <- get_legend(
          stat_plot[[coor_df[i, "group"]]] +
            theme(legend.position = "bottom")
        )
      }

      if (!is.null(lineages)) {
        lineages_layers <- c(list(ggnewscale::new_scale_color()), lineages_layers)
        suppressMessages({
          legend_list[["lineages"]] <- get_legend(
            ggplot() +
              lineages_layers +
              theme_scop(
                legend.position = "bottom",
                legend.direction = legend.direction
              )
          )
        })
        p <- suppressWarnings({
          p + lineages_layers + theme(legend.position = "none")
        })
        if (is.null(legend_list[["lineages"]])) {
          legend_list["lineages"] <- list(NULL)
        }
      }

      if (!is.null(paga)) {
        paga_layers <- c(list(ggnewscale::new_scale_color()), paga_layers)
        if (g != paga$groups) {
          suppressMessages({
            legend_list[["paga"]] <- get_legend(
              ggplot() +
                paga_layers +
                theme_scop(
                  legend.position = "bottom",
                  legend.direction = legend.direction
                )
            )
          })
        }
        p <- suppressWarnings({
          p + paga_layers + theme(legend.position = "none")
        })
        if (is.null(legend_list[["paga"]])) {
          legend_list["paga"] <- list(NULL)
        }
      }

      if (!is.null(velocity)) {
        velocity_layers <- c(
          list(ggnewscale::new_scale_color()),
          list(ggnewscale::new_scale("size")),
          velocity_layers
        )
        if (velocity_plot_type == "stream" && is.null(streamline_color)) {
          suppressMessages({
            legend_list[["velocity"]] <- get_legend(
              ggplot() +
                velocity_layers +
                theme_scop(
                  legend.position = "bottom",
                  legend.direction = legend.direction
                )
            )
          })
        }
        p <- suppressWarnings({
          p + velocity_layers + theme(legend.position = "none")
        })
        if (is.null(legend_list[["velocity"]])) {
          legend_list["velocity"] <- list(NULL)
        }
      }

      if (isTRUE(label)) {
        label_df <- stats::aggregate(
          p$data[, c("x", "y")],
          by = list(p$data[["label.by"]]),
          FUN = stats::median
        )
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
        if (isFALSE(label_insitu)) {
          label_df[, "label"] <- seq_len(nrow(label_df))
        }
        if (isTRUE(label_repel)) {
          p <- p +
            geom_point(
              data = label_df,
              mapping = aes(x = .data[["x"]], y = .data[["y"]]),
              shape = "circle",
              color = label_point_color,
              size = label_point_size
            ) +
            ggrepel::geom_text_repel(
              data = label_df,
              aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
              fontface = "bold",
              min.segment.length = 0,
              segment.color = label_segment_color,
              point.size = label_point_size,
              max.overlaps = 100,
              force = label_repulsion,
              color = label.fg,
              bg.color = label.bg,
              bg.r = label.bg.r,
              size = label.size,
              inherit.aes = FALSE
            )
        } else {
          p <- p +
            ggrepel::geom_text_repel(
              data = label_df,
              aes(
                x = .data[["x"]],
                y = .data[["y"]],
                label = .data[["label"]]
              ),
              fontface = "bold",
              min.segment.length = 0,
              segment.color = label_segment_color,
              point.size = NA,
              max.overlaps = 100,
              force = 0,
              color = label.fg,
              bg.color = label.bg,
              bg.r = label.bg.r,
              size = label.size,
              inherit.aes = FALSE
            )
        }
      }

      if (length(legend_list) > 0) {
        legend_list <- legend_list[!sapply(legend_list, is.null)]
        if (is.null(nested_legend)) {
          legend_base <- get_legend(
            p_base +
              theme_scop(
                legend.position = "bottom",
                legend.direction = legend.direction
              )
          )
        } else {
          legend_base <- nested_legend
        }
        if (legend.direction == "vertical") {
          legend <- do.call(cbind, c(list(base = legend_base), legend_list))
        } else {
          legend <- do.call(rbind, c(list(base = legend_base), legend_list))
        }
        gtable <- as_grob(p + theme(legend.position = "none"))
        gtable <- add_grob(gtable, legend, legend.position)
        p <- patchwork::wrap_plots(gtable)
      } else if (!is.null(nested_legend)) {
        gtable <- as_grob(p + theme(legend.position = "none"))
        gtable <- add_grob(
          gtable,
          nested_legend,
          legend.position,
          space = attr(nested_legend, "legend_space") %||% NULL,
          clip = "off"
        )
        p <- patchwork::wrap_plots(gtable)
      }

      return(p)
    }
  )

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(
        plotlist = plist,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      )
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

cell_dim_grid_layer <- function(
  dat,
  add_grid = FALSE,
  grid_n = 12,
  grid_color = "black",
  grid_size = 0.25,
  grid_alpha = 0.35
) {
  if (!isTRUE(add_grid)) {
    return(NULL)
  }
  if (!is.numeric(grid_n) || length(grid_n) != 1L || is.na(grid_n) || grid_n < 2) {
    log_message(
      "{.arg grid_n} must be a single numeric value greater than or equal to 2",
      message_type = "error"
    )
  }
  x_range <- range(dat[["x"]], na.rm = TRUE)
  y_range <- range(dat[["y"]], na.rm = TRUE)
  if (!all(is.finite(c(x_range, y_range)))) {
    return(NULL)
  }
  grid_df <- expand.grid(
    x = seq(x_range[1], x_range[2], length.out = grid_n),
    y = seq(y_range[1], y_range[2], length.out = grid_n)
  )
  ggplot2::geom_point(
    data = grid_df,
    mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
    inherit.aes = FALSE,
    shape = 16,
    color = grid_color,
    size = grid_size,
    alpha = grid_alpha,
    show.legend = FALSE
  )
}

cell_dim_palette_colors <- function(
  x,
  palette = "Chinese",
  palcolor = NULL,
  ...
) {
  if (!is.null(palcolor) && !is.null(names(palcolor))) {
    x_levels <- if (is.factor(x)) {
      levels(x)
    } else {
      unique(as.character(x))
    }
    x_levels <- x_levels[!is.na(x_levels) & nzchar(x_levels)]
    if (length(x_levels) > 0L && all(x_levels %in% names(palcolor))) {
      palcolor <- palcolor[x_levels]
    }
  }
  palette_colors(
    x = x,
    palette = palette,
    palcolor = palcolor,
    ...
  )
}

cell_dim_nested_legend_grob <- function(
  dat,
  colors,
  parent_colors = NULL,
  title = NULL,
  legend_direction = "horizontal",
  legend_position = "bottom"
) {
  if (nrow(dat) == 0L) {
    return(NULL)
  }
  legend_data <- cell_dim_nested_legend_data(
    dat = dat,
    colors = colors,
    parent_colors = parent_colors
  )
  if (is.null(legend_data)) {
    return(NULL)
  }
  parent_levels <- legend_data[["parent_levels"]]
  parent_colors <- legend_data[["parent_colors"]]
  nested_items <- legend_data[["nested_items"]]
  side_legend <- legend_position %in% c("left", "right")
  dense_side_legend <- isTRUE(side_legend) && length(parent_levels) > 8L
  parent_text_size <- if (isTRUE(dense_side_legend)) 2.45 else 3
  child_text_size <- if (isTRUE(dense_side_legend)) 2.45 else 3
  child_point_size <- if (isTRUE(dense_side_legend)) 1.8 else 2.2

  if (isTRUE(side_legend)) {
    item_step <- 0.9
    group_gap <- 1.1
    max_column_depth <- if (isTRUE(dense_side_legend)) 35 else 18
    label_values <- c(parent_levels, unlist(nested_items, use.names = FALSE))
    max_label_chars <- max(nchar(label_values), na.rm = TRUE)
    column_width <- max(2.8, min(4.8, 0.18 * max_label_chars + 1.35))

    side_layout <- data.frame(
      parent = character(),
      column = integer(),
      header_y = numeric(),
      stringsAsFactors = FALSE
    )
    column_i <- 1L
    y_cursor <- 0
    for (parent_i in parent_levels) {
      n_items <- length(nested_items[[parent_i]])
      group_depth <- item_step * n_items + group_gap
      if (nrow(side_layout) > 0L && abs(y_cursor - group_depth) > max_column_depth) {
        column_i <- column_i + 1L
        y_cursor <- 0
      }
      side_layout <- rbind(
        side_layout,
        data.frame(
          parent = parent_i,
          column = column_i,
          header_y = y_cursor,
          stringsAsFactors = FALSE
        )
      )
      y_cursor <- y_cursor - group_depth
    }
    side_layout[["base_x"]] <- 0.34 + (side_layout[["column"]] - 1L) * column_width
    header_df <- do.call(rbind, lapply(seq_len(nrow(side_layout)), function(i) {
      parent_i <- side_layout[["parent"]][i]
      data.frame(
        parent = parent_i,
        x = side_layout[["base_x"]][i],
        y = side_layout[["header_y"]][i],
        fill = unname(parent_colors[parent_i]),
        stringsAsFactors = FALSE
      )
    }))
    legend_items <- do.call(rbind, lapply(seq_len(nrow(side_layout)), function(i) {
      parent_i <- side_layout[["parent"]][i]
      subgroup_i <- nested_items[[parent_i]]
      data.frame(
        parent = parent_i,
        subgroup = subgroup_i,
        x = side_layout[["base_x"]][i] + 0.58,
        y = side_layout[["header_y"]][i] - seq_along(subgroup_i) * item_step,
        color = unname(colors[subgroup_i]),
        stringsAsFactors = FALSE
      )
    }))
    n_columns <- max(side_layout[["column"]])
    x_limits <- c(0, 0.34 + n_columns * column_width)
    y_limits <- c(-max_column_depth - 0.2, 0.22)
  } else {
    label_values <- c(parent_levels, unlist(nested_items, use.names = FALSE))
    max_label_chars <- max(nchar(label_values), na.rm = TRUE)
    column_width <- max(2, min(3.5, 0.16 * max_label_chars + 1.15))
    parent_x <- 0.35 + (seq_along(parent_levels) - 1) * column_width
    names(parent_x) <- parent_levels
    header_df <- data.frame(
      parent = parent_levels,
      x = unname(parent_x[parent_levels]),
      y = 0,
      fill = unname(parent_colors[parent_levels]),
      stringsAsFactors = FALSE
    )
    legend_items <- do.call(rbind, lapply(parent_levels, function(parent_i) {
      subgroup_i <- nested_items[[parent_i]]
      data.frame(
        parent = parent_i,
        subgroup = subgroup_i,
        x = unname(parent_x[parent_i]) + 0.58,
        y = -seq_along(subgroup_i),
        color = unname(colors[subgroup_i]),
        stringsAsFactors = FALSE
      )
    }))
    x_limits <- c(0, max(parent_x) + column_width * 0.55)
    y_limits <- c(min(legend_items[["y"]], -1) - 0.35, 0.45)
  }
  legend_plot <- ggplot() +
    geom_label(
      data = header_df,
      mapping = aes(
        x = .data[["x"]],
        y = .data[["y"]],
        label = .data[["parent"]],
        fill = .data[["parent"]]
      ),
      color = "white",
      fontface = "bold",
      label.size = 0,
      label.r = grid::unit(0, "pt"),
      label.padding = grid::unit(1.8, "pt"),
      hjust = 0,
      size = parent_text_size,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = legend_items,
      mapping = aes(x = .data[["x"]] - 0.26, y = .data[["y"]]),
      color = legend_items[["color"]],
      size = child_point_size,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = legend_items,
      mapping = aes(x = .data[["x"]] - 0.08, y = .data[["y"]], label = .data[["subgroup"]]),
      hjust = 0,
      size = child_text_size,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = stats::setNames(header_df[["fill"]], header_df[["parent"]])) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits,
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = grid::unit(c(1, 2, 1, 2), "pt"),
      legend.position = "none"
    )

  legend_grob <- as_grob(legend_plot)
  if (isTRUE(side_legend)) {
    attr(legend_grob, "legend_space") <- grid::unit(
      max(1.85, min(7.2, 0.82 * diff(x_limits))),
      "in"
    )
  } else {
    attr(legend_grob, "legend_space") <- grid::unit(
      max(0.42, 0.16 * abs(min(y_limits)) + 0.2),
      "in"
    )
  }
  legend_grob
}

cell_dim_nested_legend_data <- function(
  dat,
  colors,
  parent_colors = NULL
) {
  group_levels <- names(colors)
  group_levels <- group_levels[group_levels %in% as.character(dat[["group.by"]])]
  if (length(group_levels) == 0L) {
    return(NULL)
  }
  parent_levels <- if (is.factor(dat[["legend.by"]])) {
    levels(droplevels(dat[["legend.by"]]))
  } else {
    unique(as.character(dat[["legend.by"]]))
  }
  parent_levels <- parent_levels[nzchar(parent_levels)]
  if (length(parent_levels) == 0L) {
    return(NULL)
  }
  nested_items <- lapply(parent_levels, function(parent_i) {
    subgroup_i <- group_levels[
      group_levels %in% as.character(dat[dat[["legend.by"]] == parent_i, "group.by"])
    ]
    subgroup_i
  })
  names(nested_items) <- parent_levels
  parent_levels <- parent_levels[lengths(nested_items) > 0L]
  nested_items <- nested_items[parent_levels]
  if (length(parent_levels) == 0L) {
    return(NULL)
  }

  fallback_parent_colors <- vapply(parent_levels, function(parent_i) {
    subgroup_i <- nested_items[[parent_i]]
    unname(colors[subgroup_i[1]])
  }, character(1))
  if (is.null(parent_colors)) {
    parent_colors <- fallback_parent_colors
  } else {
    parent_colors <- parent_colors[parent_levels]
    missing_parent_colors <- is.na(parent_colors) | !nzchar(parent_colors)
    parent_colors[missing_parent_colors] <- fallback_parent_colors[missing_parent_colors]
  }

  list(
    group_levels = group_levels,
    parent_levels = parent_levels,
    nested_items = nested_items,
    parent_colors = parent_colors
  )
}

#' @title 3D-Dimensional reduction plot for cell classification visualization.
#'
#' @description
#' Plotting cell points on a reduced 3D space and coloring according to the groups of the cells.
#'
#' @md
#' @inheritParams CellDimPlot
#' @param dims Dimensions to plot, must be a three-length numeric vector specifying x-, y- and z-dimensions
#' @param axis_labs A character vector of length 3 indicating the labels for the axes.
#' @param span The span of the loess smoother for lineages line.
#' @param shape.highlight Shape of the cell to highlight.
#' See \href{https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol}{scattergl-marker-symbol}
#' @param width Width in pixels, defaults to automatic sizing.
#' @param height Height in pixels, defaults to automatic sizing.
#' @param save The name of the file to save the plot to. Must end in ".html".
#'
#' @seealso
#' [CellDimPlot], [FeatureDimPlot3D]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   nonlinear_reduction_dims = 3
#' )
#' CellDimPlot3D(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D"
#' )
#'
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D",
#'   show_plot = FALSE
#' )
#' CellDimPlot3D(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D",
#'   lineages = "Lineage1"
#' )
CellDimPlot3D <- function(
  srt,
  group.by,
  reduction = NULL,
  dims = c(1, 2, 3),
  axis_labs = NULL,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey80",
  pt.size = 1.5,
  cells.highlight = NULL,
  cols.highlight = "black",
  shape.highlight = "circle-open",
  sizes.highlight = 2,
  lineages = NULL,
  lineages_palette = "Dark2",
  span = 0.75,
  width = NULL,
  height = NULL,
  save = NULL,
  force = FALSE,
  verbose = TRUE
) {
  bg_color <- col2hex(bg_color)
  cols.highlight <- col2hex(cols.highlight)

  for (i in c(group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      log_message(
        paste0(i, " is not in the meta.data of srt object."),
        message_type = "error"
      )
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(
        srt@meta.data[[i]],
        levels = unique(srt@meta.data[[i]])
      )
    }
  }
  for (l in lineages) {
    if (!l %in% colnames(srt@meta.data)) {
      log_message(
        paste0(l, " is not in the meta.data of srt object."),
        message_type = "error"
      )
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt, min_dim = 3)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction, min_dim = 3)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      paste0(reduction, " is not in the srt reduction names."),
      message_type = "error"
    )
  }
  if (ncol(srt@reductions[[reduction]]@cell.embeddings) < 3) {
    log_message(
      "Reduction must be in three dimensions or higher.",
      message_type = "error"
    )
  }
  if (!is.null(cells.highlight) && isFALSE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "No cells in 'cells.highlight' found in srt.",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "Some cells in 'cells.highlight' not found in srt.",
        message_type = "warning",
        verbose = verbose
      )
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  reduction_key <- srt@reductions[[reduction]]@key
  if (is.null(axis_labs) || length(axis_labs) != 3) {
    xlab <- paste0(reduction_key, dims[1])
    ylab <- paste0(reduction_key, dims[2])
    zlab <- paste0(reduction_key, dims[3])
  } else {
    xlab <- axis_labs[1]
    ylab <- axis_labs[2]
    zlab <- axis_labs[3]
  }
  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    check_r("htmlwidgets", verbose = FALSE)
    if (!grepl(".html$", save)) {
      log_message(
        "'save' must be a string with .html as a suffix.",
        message_type = "error"
      )
    }
  }

  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_use <- cbind(
    dat_dim[colnames(srt@assays[[1]]), , drop = FALSE],
    srt@meta.data[colnames(srt@assays[[1]]), , drop = FALSE]
  )
  nlev <- sapply(dat_use[, group.by, drop = FALSE], nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && isFALSE(force)) {
    log_message(
      paste0(
        "The following variables have more than 100 levels: ",
        paste(names(nlev), collapse = ",")
      ),
      message_type = "warning",
      verbose = verbose
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(invisible(NULL))
    }
  }
  if (!is.null(lineages)) {
    dat_lineages <- srt@meta.data[, unique(lineages), drop = FALSE]
    dat_use <- cbind(dat_use, dat_lineages[row.names(dat_use), , drop = FALSE])
  }
  dat_use[["group.by"]] <- dat_use[[group.by]]
  if (any(is.na(dat_use[[group.by]]))) {
    n <- as.character(dat_use[[group.by]])
    n[is.na(n)] <- "NA"
    dat_use[[group.by]] <- factor(
      n,
      levels = c(levels(dat_use[[group.by]]), "NA")
    )
  }

  dat_use[["color"]] <- dat_use[[group.by]]
  colors <- palette_colors(
    dat_use[["group.by"]],
    palette = palette,
    palcolor = palcolor,
    NA_color = bg_color,
    NA_keep = TRUE
  )

  dat_use[[paste0(reduction_key, dims[1], "All_cells")]] <- dat_use[[paste0(
    reduction_key,
    dims[1]
  )]]
  dat_use[[paste0(reduction_key, dims[2], "All_cells")]] <- dat_use[[paste0(
    reduction_key,
    dims[2]
  )]]
  dat_use[[paste0(reduction_key, dims[3], "All_cells")]] <- dat_use[[paste0(
    reduction_key,
    dims[3]
  )]]
  cells_highlight_use <- cells.highlight
  if (isTRUE(cells_highlight_use)) {
    cells_highlight_use <- rownames(dat_use)[dat_use[[group.by]] != "NA"]
  }
  if (!is.null(cells_highlight_use)) {
    cells_highlight_use <- cells_highlight_use[
      cells_highlight_use %in% rownames(dat_use)
    ]
    dat_use_highlight <- dat_use[cells_highlight_use, , drop = FALSE]
  }

  p <- plotly::plot_ly(data = dat_use, width = width, height = height)
  p <- plotly::add_trace(
    p = p,
    data = dat_use,
    x = dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
    y = dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
    z = dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
    text = paste0(
      "Cell:",
      rownames(dat_use),
      "\ngroup.by:",
      dat_use[["group.by"]],
      "\ncolor:",
      dat_use[["color"]]
    ),
    type = "scatter3d",
    mode = "markers",
    color = dat_use[[group.by]],
    colors = colors,
    marker = list(size = pt.size),
    showlegend = TRUE,
    visible = TRUE
  )
  if (!is.null(cells_highlight_use)) {
    p <- plotly::add_trace(
      p = p,
      x = dat_use_highlight[[paste0(reduction_key, dims[1], "All_cells")]],
      y = dat_use_highlight[[paste0(reduction_key, dims[2], "All_cells")]],
      z = dat_use_highlight[[paste0(reduction_key, dims[3], "All_cells")]],
      text = paste0(
        "Cell:",
        rownames(dat_use_highlight),
        "\ngroup.by:",
        dat_use_highlight[["group.by"]],
        "\ncolor:",
        dat_use_highlight[["color"]]
      ),
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = sizes.highlight,
        color = cols.highlight,
        symbol = shape.highlight
      ),
      name = "highlight",
      showlegend = TRUE,
      visible = TRUE
    )
  }
  if (!is.null(lineages)) {
    for (l in lineages) {
      dat_sub <- dat_use[!is.na(dat_use[[l]]), , drop = FALSE]
      dat_sub <- dat_sub[order(dat_sub[[l]]), , drop = FALSE]

      xlo <- stats::loess(
        stats::formula(
          paste(
            paste0(reduction_key, dims[1], "All_cells"),
            l,
            sep = "~"
          )
        ),
        data = dat_sub,
        span = span,
        degree = 2
      )
      ylo <- stats::loess(
        stats::formula(
          paste(
            paste0(reduction_key, dims[2], "All_cells"),
            l,
            sep = "~"
          )
        ),
        data = dat_sub,
        span = span,
        degree = 2
      )
      zlo <- stats::loess(
        stats::formula(
          paste(
            paste0(reduction_key, dims[3], "All_cells"),
            l,
            sep = "~"
          )
        ),
        data = dat_sub,
        span = span,
        degree = 2
      )
      dat_smooth <- data.frame(x = xlo$fitted, y = ylo$fitted, z = zlo$fitted)
      dat_smooth <- dat_smooth[
        dat_smooth[["x"]] <=
          max(
            dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
            na.rm = TRUE
          ) &
          dat_smooth[["x"]] >=
            min(
              dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
              na.rm = TRUE
            ), ,
        drop = FALSE
      ]
      dat_smooth <- dat_smooth[
        dat_smooth[["y"]] <=
          max(
            dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
            na.rm = TRUE
          ) &
          dat_smooth[["y"]] >=
            min(
              dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
              na.rm = TRUE
            ), ,
        drop = FALSE
      ]
      dat_smooth <- dat_smooth[
        dat_smooth[["z"]] <=
          max(
            dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
            na.rm = TRUE
          ) &
          dat_smooth[["z"]] >=
            min(
              dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
              na.rm = TRUE
            ), ,
        drop = FALSE
      ]
      dat_smooth <- unique(stats::na.omit(dat_smooth))
      p <- plotly::add_trace(
        p = p,
        x = dat_smooth[, "x"],
        y = dat_smooth[, "y"],
        z = dat_smooth[, "z"],
        text = paste0(
          "Lineage:",
          l
        ),
        type = "scatter3d",
        mode = "lines",
        line = list(
          width = 6,
          color = palette_colors(x = lineages, palette = lineages_palette)[l],
          reverscale = FALSE
        ),
        name = l,
        showlegend = TRUE,
        visible = TRUE
      )
    }
  }

  p <- plotly::layout(
    p = p,
    title = list(
      text = paste0("Total", " (nCells:", nrow(dat_use), ")"),
      font = list(size = 16, color = "black"),
      y = 0.95
    ),
    font = list(size = 12, color = "black"),
    showlegend = TRUE,
    legend = list(
      itemsizing = "constant",
      y = 0.5,
      x = 1,
      xanchor = "left",
      alpha = 1
    ),
    scene = list(
      xaxis = list(
        title = xlab,
        range = c(
          min(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE),
          max(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE)
        )
      ),
      yaxis = list(
        title = ylab,
        range = c(
          min(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE),
          max(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE)
        )
      ),
      zaxis = list(
        title = zlab,
        range = c(
          min(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE),
          max(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE)
        )
      ),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    autosize = FALSE
  )

  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    htmlwidgets::saveWidget(
      widget = plotly::as_widget(p),
      file = save
    )
    unlink(gsub("\\.html", "_files", save), recursive = TRUE)
  }

  return(p)
}
