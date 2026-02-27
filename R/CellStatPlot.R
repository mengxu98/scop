#' @title Statistical plot of cells
#'
#' @md
#' @inheritParams FeatureStatPlot
#' @param NA_color The color to use for missing values.
#' @param NA_stat Whether to include missing values in the plot.
#' Default is `TRUE`.
#' @param stat_level The level(s) of the variable(s) specified in `stat.by` to include in the plot.
#' Default is `NULL`.
#' @param stat_type The type of statistic to compute for the plot.
#' Can be one of `"percent"` or `"count"`.
#' @param position The position adjustment for the plot.
#' Can be one of `"stack"` or `"dodge"`.
#' @param label Whether to add labels on the plot.
#' Default is `FALSE`.
#' @param label.size The size of the labels.
#' @param label.fg The foreground color of the labels.
#' @param label.bg The background color of the labels.
#' @param label.bg.r The radius of the rounded corners of the label background.
#'
#' @seealso [FeatureStatPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' p1 <- CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   label = TRUE
#' )
#' p1
#'
#' thisplot::panel_fix(p1, height = 2, width = 3)
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   stat_type = "count",
#'   position = "dodge",
#'   label = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   palette = "Set1",
#'   stat_type = "count",
#'   position = "dodge"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "bar"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "rose"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "ring"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "pie"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   plot_type = "dot"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "bar"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "rose"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "ring"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "area"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "dot"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "trend"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   plot_type = "bar",
#'   individual = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "bar"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "rose"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "ring"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "area"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "dot"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "trend"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "bar",
#'   position = "dodge",
#'   label = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "rose",
#'   position = "dodge",
#'   label = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = "Phase",
#'   group.by = "CellType",
#'   stat_type = "count",
#'   plot_type = "ring",
#'   position = "dodge",
#'   label = TRUE
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "sankey"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "chord"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c("CellType", "Phase"),
#'   plot_type = "venn",
#'   stat_level = list(
#'     CellType = c("Ductal", "Ngn3-low-EP"),
#'     Phase = "S"
#'   )
#' )
#'
#' pancreas_sub$Progenitor <- pancreas_sub$CellType %in% c("Ngn3-low-EP", "Ngn3-high-EP")
#' pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
#' pancreas_sub$Fancb_Expressed <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )["Fancb", ] > 0
#' pancreas_sub$Dlg3_Expressed <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )["Dlg3", ] > 0
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
#'   ),
#'   plot_type = "venn",
#'   stat_level = "TRUE"
#' )
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
#'   ),
#'   plot_type = "upset",
#'   stat_level = "TRUE"
#' )
#'
#' sum(
#'   pancreas_sub$Progenitor == "FALSE" &
#'     pancreas_sub$G2M == "FALSE" &
#'     pancreas_sub$Fancb_Expressed == "TRUE" &
#'     pancreas_sub$Dlg3_Expressed == "FALSE"
#' )
CellStatPlot <- function(
    srt,
    stat.by,
    group.by = NULL,
    split.by = NULL,
    bg.by = NULL,
    cells = NULL,
    flip = FALSE,
    NA_color = "grey",
    NA_stat = TRUE,
    keep_empty = FALSE,
    individual = FALSE,
    stat_level = NULL,
    plot_type = c(
      "bar",
      "rose",
      "ring",
      "pie",
      "trend",
      "area",
      "dot",
      "sankey",
      "chord",
      "venn",
      "upset"
    ),
    stat_type = c("percent", "count"),
    position = c("stack", "dodge"),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    bg_palette = "Paired",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    label = FALSE,
    label.size = 3.5,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    aspect.ratio = NULL,
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
  cells <- cells %||% colnames(srt@assays[[1]])
  meta_data <- srt@meta.data[cells, , drop = FALSE]

  plot <- StatPlot(
    meta_data,
    stat.by = stat.by,
    group.by = group.by,
    split.by = split.by,
    bg.by = bg.by,
    flip = flip,
    NA_color = NA_color,
    NA_stat = NA_stat,
    keep_empty = keep_empty,
    individual = individual,
    stat_level = stat_level,
    plot_type = plot_type,
    stat_type = stat_type,
    position = position,
    palette = palette,
    palcolor = palcolor,
    alpha = alpha,
    bg_palette = bg_palette,
    bg_palcolor = bg_palcolor,
    bg_alpha = bg_alpha,
    label = label,
    label.size = label.size,
    label.fg = label.fg,
    label.bg = label.bg,
    label.bg.r = label.bg.r,
    aspect.ratio = aspect.ratio,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    force = force,
    seed = seed
  )

  return(plot)
}
