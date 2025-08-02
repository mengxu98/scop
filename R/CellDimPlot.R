#' Visualize cell groups on a 2-dimensional reduction plot
#'
#' Plotting cell points on a reduced 2D plane and coloring according to the groups.
#'
#' @md
#' @param srt A Seurat object.
#' @param group.by Name of one or more meta.data columns to group (color) cells by (for example, orig.ident).
#' @param reduction Which dimensionality reduction to use.
#' If not specified, will use the reduction returned by \link{DefaultReduction}.
#' @param split.by Name of a column in meta.data column to split plot by.
#' @param palette Name of a color palette name collected in scop. Default is "Paired".
#' @param palcolor Custom colors used to create a color palette.
#' @param bg_color Color value for background(NA) points.
#' @param pt.size Point size.
#' @param pt.alpha Point transparency.
#' @param cells.highlight A vector of cell names to highlight.
#' @param cols.highlight Color used to highlight the cells.
#' @param sizes.highlight Size of highlighted cell points.
#' @param alpha.highlight Transparency of highlighted cell points.
#' @param stroke.highlight Border width of highlighted cell points.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top").
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical")
#' @param combine Combine plots into a single \code{patchwork} object.
#' If \code{FALSE}, return a list of ggplot objects.
#' @param nrow Number of rows in the combined plot.
#' @param ncol Number of columns in the combined plot.
#' @param byrow Logical value indicating if the plots should be arrange by row (default) or by column.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param show_na Whether to assign a color from the color palette to NA group.
#' If \code{FALSE}, cell points with NA level will colored by \code{bg_color}.
#' @param show_stat Whether to show statistical information on the plot.
#' @param label Whether to label the cell groups.
#' @param label_insitu Whether to place the raw labels (group names) in the center of the cells with the corresponding group.
#' Default is \code{FALSE}, which using numbers instead of raw labels.
#' @param label.size Size of labels.
#' @param label.fg Foreground color of label.
#' @param label.bg Background color of label.
#' @param label.bg.r Background ratio of label.
#' @param label_repel Logical value indicating whether the label is repel away from the center points.
#' @param label_repulsion Force of repulsion between overlapping text labels. Defaults to 20.
#' @param label_point_size Size of the center points.
#' @param label_point_color Color of the center points.
#' @param label_segment_color Color of the line segment for labels.
#' @param add_density Whether to add a density layer on the plot.
#' @param density_color Color of the density contours lines.
#' @param density_filled Whether to add filled contour bands instead of contour lines.
#' @param density_filled_palette Color palette used to fill contour bands.
#' @param density_filled_palcolor Custom colors used to fill contour bands.
#' @param add_mark Whether to add marks around cell groups. Default is \code{FALSE}.
#' @param mark_type Type of mark to add around cell groups.
#' One of "hull", "ellipse", "rect", or "circle". Default is "hull".
#' @param mark_expand Expansion of the mark around the cell group.
#' Default is \code{grid::unit(3, "mm")}.
#' @param mark_alpha Transparency of the mark.
#' Default is 0.1.
#' @param mark_linetype Line type of the mark border.
#' Default is 1 (solid line).
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
#' @param stat_plot_size Set the statistical plot size. Defaults to 0.1
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
#' @param raster Convert points to raster format, default is NULL which automatically rasterizes if plotting more than 100,000 cells
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param theme_use Theme used. Can be a character string or a theme function.
#' For example, \code{"theme_blank"} or \code{ggplot2::theme_classic}.
#' @param aspect.ratio Aspect ratio of the panel.
#' @param title The text for the title.
#' @param subtitle The text for the subtitle for the plot which will be displayed below the title.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param force Whether to force drawing regardless of maximum levels in any cell group is greater than 100.
#' @param cells Subset cells to plot.
#' @param theme_args Other arguments passed to the \code{theme_use}.
#' @param seed Random seed set for reproducibility
#'
#' @seealso \link{FeatureDimPlot}
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' p1 <- CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' p1
#'
#' panel_fix(
#'   p1,
#'   height = 2,
#'   raster = TRUE,
#'   dpi = 30
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   theme_use = "theme_blank"
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE,
#'   label.fg = "orange",
#'   label.bg = "red",
#'   label.size = 5
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   label = TRUE,
#'   label_insitu = TRUE
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_expand = grid::unit(1, "mm")
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_alpha = 0.3
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_linetype = 2
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_type = "ellipse"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_mark = TRUE,
#'   mark_type = "rect"
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   add_density = TRUE,
#'   density_filled = TRUE
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   stat.by = "Phase",
#'   stat_plot_type = "ring",
#'   stat_plot_label = TRUE,
#'   stat_plot_size = 0.15
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   hex = TRUE,
#'   hex.bins = 20
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   hex = TRUE,
#'   hex.count = FALSE
#' )
#'
#' # Show neighbors graphs on the plot
#' pancreas_sub <- standard_scop(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   graph = "Standardpca_SNN"
#' )
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
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = paste0("Lineage", 1:2),
#'   reduction = "UMAP"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_whiskers = TRUE
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_span = 0.1
#' )
#'
#' \dontrun{
#' # Show PAGA results on the plot
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   return_seurat = TRUE
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga,
#'   paga_type = "connectivities_tree"
#' )
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
#'   srt = pancreas_sub,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP",
#'   mode = "stochastic",
#'   return_seurat = TRUE
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   paga = pancreas_sub@misc$paga,
#'   paga_show_transition = TRUE
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = NA,
#'   velocity = "stochastic"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   velocity = "stochastic",
#'   velocity_plot_type = "grid"
#' )
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
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   pt.size = 5,
#'   pt.alpha = 0.2,
#'   velocity = "stochastic",
#'   velocity_plot_type = "stream"
#' )
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
#' }
CellDimPlot <- function(
    srt,
    group.by,
    reduction = NULL,
    dims = c(1, 2),
    split.by = NULL,
    cells = NULL,
    show_na = FALSE,
    show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
    pt.size = NULL,
    pt.alpha = 1,
    palette = "Paired",
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
    stat_plot_size = 0.15,
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
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    force = FALSE,
    seed = 11) {
  set.seed(seed)
  mark_type <- match.arg(mark_type)

  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in unique(c(group.by, split.by))) {
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
        paste0(l, " is not in the meta.data of srt object."),
        message_type = "error"
      )
    }
  }
  if (!is.null(graph) && !graph %in% names(srt@graphs)) {
    log_message(
      paste0("Graph ", graph, " is not exist in the srt object."),
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
      paste0(reduction, " is not in the srt reduction names."),
      message_type = "error"
    )
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "No cells in 'cells.highlight' found in srt.",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "Some cells in 'cells.highlight' not found in srt.",
        message_type = "warning"
      )
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }

  dat_meta <- srt@meta.data[, unique(c(group.by, split.by)), drop = FALSE]
  nlev <- sapply(dat_meta, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    log_message(
      paste0(
        "The following variables have more than 100 levels: ",
        paste(names(nlev), collapse = ",")
      ),
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
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
    check_r("scattermore")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      log_message(
        "'raster.dpi' must be a two-length numeric vector",
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
      group_by = group.by,
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
      colors <- palette_scop(
        levels(dat_use[[g]]),
        palette = palette,
        palcolor = palcolor,
        NA_keep = TRUE
      )

      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      legend_list <- list()
      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
      cells.highlight_use <- cells.highlight
      if (isTRUE(cells.highlight_use)) {
        cells.highlight_use <- rownames(dat)[!is.na(dat[[g]])]
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
        mark_fun <- switch(mark_type,
          "ellipse" = "geom_mark_ellipse",
          "hull" = "geom_mark_hull",
          "rect" = "geom_mark_rect",
          "circle" = "geom_mark_circle"
        )
        mark <- list(
          do.call(
            mark_fun,
            list(
              data = dat[!is.na(dat[["group.by"]]), , drop = FALSE],
              mapping = aes(
                x = .data[["x"]],
                y = .data[["y"]],
                color = .data[["group.by"]],
                fill = .data[["group.by"]]
              ),
              expand = mark_expand,
              alpha = mark_alpha,
              linetype = mark_linetype,
              show.legend = FALSE,
              inherit.aes = FALSE
            ),
          ),
          scale_fill_manual(values = colors[names(labels_tb)]),
          scale_color_manual(values = colors[names(labels_tb)]),
          ggnewscale::new_scale_fill(),
          ggnewscale::new_scale_color()
        )
      } else {
        mark <- NULL
      }

      if (!is.null(graph)) {
        net_mat <- Matrix::as.matrix(graph)[rownames(dat), rownames(dat)]
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
          filled_color <- palette_scop(
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

      p <- ggplot(dat) +
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
        check_r("hexbin")
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
            size = pt.size,
            alpha = pt.alpha
          )
      }

      if (!is.null(cells.highlight_use) && !isTRUE(hex)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
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
                size = sizes.highlight,
                alpha = alpha.highlight
              )
          }
        }
      }
      p <- p +
        scale_color_manual(
          name = paste0(g, ":"),
          values = colors[names(labels_tb)],
          labels = label_use,
          na.value = bg_color,
          guide = guide_legend(
            title.hjust = 0,
            order = 1,
            override.aes = list(size = 4, alpha = 1)
          )
        ) +
        scale_fill_manual(
          name = paste0(g, ":"),
          values = colors[names(labels_tb)],
          labels = label_use,
          na.value = bg_color,
          guide = guide_legend(
            title.hjust = 0,
            order = 1
          )
        )
      p_base <- p

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
        if (velocity_plot_type != "raw") {
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
          by = list(p$data[["group.by"]]),
          FUN = stats::median
        )
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
        if (!isTRUE(label_insitu)) {
          label_df[, "label"] <- seq_len(nrow(label_df))
        }
        if (isTRUE(label_repel)) {
          p <- p +
            geom_point(
              data = label_df,
              mapping = aes(x = .data[["x"]], y = .data[["y"]]),
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
        legend_base <- get_legend(
          p_base +
            theme_scop(
              legend.position = "bottom",
              legend.direction = legend.direction
            )
        )
        if (legend.direction == "vertical") {
          legend <- do.call(cbind, c(list(base = legend_base), legend_list))
        } else {
          legend <- do.call(rbind, c(list(base = legend_base), legend_list))
        }
        gtable <- as_grob(p + theme(legend.position = "none"))
        gtable <- add_grob(gtable, legend, legend.position)
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

#' 3D-Dimensional reduction plot for cell classification visualization.
#'
#' Plotting cell points on a reduced 3D space and coloring according to the groups of the cells.
#'
#' @inheritParams CellDimPlot
#' @param dims Dimensions to plot, must be a three-length numeric vector specifying x-, y- and z-dimensions
#' @param axis_labs A character vector of length 3 indicating the labels for the axes.
#' @param span A numeric value specifying the span of the loess smoother for lineages line.
#' @param shape.highlight Shape of the cell to highlight. See \href{https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol}{scattergl-marker-symbol}
#' @param width Width in pixels, defaults to automatic sizing.
#' @param height Height in pixels, defaults to automatic sizing.
#' @param save The name of the file to save the plot to. Must end in ".html".
#' @seealso \link{CellDimPlot}, \link{FeatureDimPlot3D}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' CellDimPlot3D(
#'   srt = pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D"
#' )
#'
#' pancreas_sub <- RunSlingshot(
#'   srt = pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D",
#'   show_plot = FALSE
#' )
#' CellDimPlot3D(
#'   srt = pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D",
#'   lineages = "Lineage1"
#' )
#' }
CellDimPlot3D <- function(
    srt,
    group.by,
    reduction = NULL,
    dims = c(1, 2, 3),
    axis_labs = NULL,
    palette = "Paired",
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
    force = FALSE) {
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
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "No cells in 'cells.highlight' found in srt.",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      log_message(
        "Some cells in 'cells.highlight' not found in srt.",
        message_type = "warning"
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
    check_r("htmlwidgets")
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
  if (length(nlev) > 0 && !isTRUE(force)) {
    log_message(
      paste0(
        "The following variables have more than 100 levels: ",
        paste(names(nlev), collapse = ",")
      ),
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
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
  colors <- palette_scop(
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
  cells.highlight_use <- cells.highlight
  if (isTRUE(cells.highlight_use)) {
    cells.highlight_use <- rownames(dat_use)[dat_use[[group.by]] != "NA"]
  }
  if (!is.null(cells.highlight_use)) {
    cells.highlight_use <- cells.highlight_use[
      cells.highlight_use %in% rownames(dat_use)
    ]
    dat_use_highlight <- dat_use[cells.highlight_use, , drop = FALSE]
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
  if (!is.null(cells.highlight_use)) {
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
          color = palette_scop(x = lineages, palette = lineages_palette)[l],
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
