#' Statistical plot of features
#'
#' This function generates a statistical plot for features.
#'
#' @param srt A Seurat object.
#' @param stat.by A character vector specifying the features to plot.
#' @param group.by A character vector specifying the groups to group by. Default is NULL.
#' @param split.by A character vector specifying the variable to split the plot by. Default is NULL.
#' @param plot.by A character vector specifying how to plot the data, by group or feature. Possible values are "group", "feature". Default is "group".
#' @param bg.by A character vector specifying the variable to use as the background color. Default is NULL.
#' @param fill.by A string specifying what to fill the plot by. Possible values are "group", "feature", or "expression". Default is "group".
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#' @param layer A string specifying which layer of the Seurat object to use. Default is "data".
#' @param assay A string specifying which assay to use. Default is NULL.
#' @param keep_empty A logical indicating whether to keep empty levels in the plot. Default is FALSE.
#' @param individual A logical indicating whether to create individual plots for each group. Default is FALSE.
#' @param plot_type A string specifying the type of plot to create. Possible values are "violin", "box", "bar", "dot", or "col". Default is "violin".
#' @param palette A string specifying the color palette to use for filling. Default is "Paired".
#' @param palcolor A character vector specifying specific colors to use for filling. Default is NULL.
#' @param alpha A numeric value specifying the transparency of the plot. Default is 1.
#' @param bg_palette A string specifying the color palette to use for the background. Default is "Paired".
#' @param bg_palcolor A character vector specifying specific colors to use for the background. Default is NULL.
#' @param bg_alpha A numeric value specifying the transparency of the background. Default is 0.2.
#' @param add_box A logical indicating whether to add a box plot to the plot. Default is FALSE.
#' @param box_color A string specifying the color of the box plot. Default is "black".
#' @param box_width A numeric value specifying the width of the box plot. Default is 0.1.
#' @param box_ptsize A numeric value specifying the size of the points of the box plot. Default is 2.
#' @param add_point A logical indicating whether to add individual data points to the plot. Default is FALSE.
#' @param pt.color A string specifying the color of the data points. Default is "grey30".
#' @param pt.size A numeric value specifying the size of the data points. If NULL, the size is automatically determined. Default is NULL.
#' @param pt.alpha A numeric value specifying the transparency of the data points. Default is 1.
#' @param jitter.width A numeric value specifying the width of the jitter. Default is 0.5.
#' @param jitter.height A numeric value specifying the height of the jitter. Default is 0.1.
#' @param add_trend A logical indicating whether to add a trend line to the plot. Default is FALSE.
#' @param trend_color A string specifying the color of the trend line. Default is "black".
#' @param trend_linewidth A numeric value specifying the width of the trend line. Default is 1.
#' @param trend_ptsize A numeric value specifying the size of the points of the trend line. Default is 2.
#' @param add_stat A string specifying which statistical summary to add to the plot. Possible values are "none", "mean", or "median". Default is "none".
#' @param stat_color A string specifying the color of the statistical summary. Default is "black".
#' @param stat_size A numeric value specifying the size of the statistical summary. Default is 1.
#' @param stat_stroke A numeric value specifying the stroke width of the statistical summary. Default is 1.
#' @param stat_shape A numeric value specifying the shape of the statistical summary. Default is 25.
#' @param add_line A numeric value specifying the y-intercept for adding a horizontal line. Default is NULL.
#' @param line_color A string specifying the color of the horizontal line. Default is "red".
#' @param line_size A numeric value specifying the width of the horizontal line. Default is 1.
#' @param line_type A numeric value specifying the type of the horizontal line. Default is 1.
#' @param cells.highlight A logical or character vector specifying the cells to highlight in the plot. If TRUE, all cells are highlighted. If FALSE, no cells are highlighted. Default is NULL.
#' @param cols.highlight A string specifying the color of the highlighted cells. Default is "red".
#' @param sizes.highlight A numeric value specifying the size of the highlighted cells. Default is 1.
#' @param alpha.highlight A numeric value specifying the transparency of the highlighted cells. Default is 1.
#' @param calculate_coexp A logical indicating whether to calculate co-expression values. Default is FALSE.
#' @param same.y.lims A logical indicating whether to use the same y-axis limits for all plots. Default is FALSE.
#' @param y.min A numeric or character value specifying the minimum y-axis limit. If a character value is provided, it must be of the form "qN" where N is a number between 0 and 100 (inclusive) representing the quantile to use for the limit. Default is NULL.
#' @param y.max A numeric or character value specifying the maximum y-axis limit. If a character value is provided, it must be of the form "qN" where N is a number between 0 and 100 (inclusive) representing the quantile to use for the limit. Default is NULL.
#' @param y.trans A string specifying the transformation to apply to the y-axis. Possible values are "identity" or "log2". Default is "identity".
#' @param y.nbreaks An integer specifying the number of breaks to use for the y-axis. Default is 5.
#' @param sort A logical or character value specifying whether to sort the groups on the x-axis. If TRUE, groups are sorted in increasing order. If FALSE, groups are not sorted. If "increasing", groups are sorted in increasing order. If "decreasing", groups are sorted in decreasing order. Default is FALSE.
#' @param stack A logical specifying whether to stack the plots on top of each other. Default is FALSE.
#' @param flip A logical specifying whether to flip the plot vertically. Default is FALSE.
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the groups of interest, to be compared.
#' @param ref_group A string specifying the reference group for pairwise comparisons. Default is NULL.
#' @param pairwise_method A string specifying the method to use for pairwise comparisons. Default is "wilcox.test".
#' @param multiplegroup_comparisons A logical indicating whether to add multiple group comparisons to the plot. Default is FALSE.
#' @param multiple_method A string specifying the method to use for multiple group comparisons. Default is "kruskal.test".
#' @param sig_label A string specifying the label to use for significant comparisons. Possible values are "p.signif" or "p.format". Default is "p.format".
#' @param sig_labelsize A numeric value specifying the size of the significant comparison labels. Default is 3.5.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot. Default is NULL.
#' @param title A string specifying the title of the plot. Default is NULL.
#' @param subtitle A string specifying the subtitle of the plot. Default is NULL.
#' @param xlab A string specifying the label of the x-axis. Default is NULL.
#' @param ylab A string specifying the label of the y-axis. Default is "Expression level".
#' @param legend.position A string specifying the position of the legend. Possible values are "right", "left", "top", "bottom", or "none". Default is "right".
#' @param legend.direction A string specifying the direction of the legend. Possible values are "vertical" or "horizontal". Default is "vertical".
#' @param theme_use A string specifying the theme to use for the plot. Default is "theme_scop".
#' @param theme_args A list of arguments to pass to the theme function. Default is an empty list.
#' @param combine A logical indicating whether to combine the individual plots into a single plot. Default is TRUE.
#' @param nrow An integer specifying the number of rows for the combined plot. Default is NULL.
#' @param ncol An integer specifying the number of columns for the combined plot. Default is NULL.
#' @param byrow A logical specifying whether to fill the combined plot by row or by column. Default is TRUE.
#' @param force A logical indicating whether to force the plot creation even if there are more than 100 levels in a variable. Default is FALSE.
#' @param seed An integer specifying the random seed to use for generating jitter. Default is 11.
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType"
#' ) %>% panel_fix(height = 1, width = 2)
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "box"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "bar"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "dot"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "col"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_box = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_point = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_trend = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_stat = "mean"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_line = 0.2,
#'   line_type = 2
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   add_box = TRUE,
#'   add_trend = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   comparisons = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   fill.by = "expression",
#'   palette = "Blues",
#'   same.y.lims = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   multiplegroup_comparisons = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta"))
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")),
#'   sig_label = "p.format"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   add_box = TRUE, stack = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   legend.position = "top",
#'   legend.direction = "horizontal",
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   stack = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   fill.by = "feature",
#'   plot_type = "box",
#'   group.by = "SubCellType",
#'   bg.by = "CellType", stack = TRUE, flip = TRUE
#' ) %>% panel_fix_overall(
#'   width = 8, height = 5
#' )
#' # As the plot is created by combining,
#' # we can adjust the overall height and width directly.
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "group"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   multiplegroup_comparisons = TRUE,
#'   sig_label = "p.format",
#'   sig_labelsize = 4
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   comparisons = list(c("Neurog3", "Rbp4"), c("Rbp4", "Ins1")),
#'   stack = TRUE
#' )
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ), group.by = "SubCellType",
#'   plot.by = "feature",
#'   stack = TRUE
#' )
#'
#' library(Matrix)
#' data <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "data"
#' )
#' pancreas_sub <- SeuratObject::SetAssayData(
#'   object = pancreas_sub,
#'   layer = "scale.data",
#'   assay = "RNA",
#'   new.data = Matrix::as.matrix(data / rowMeans(data))
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4"),
#'   group.by = "CellType",
#'   layer = "scale.data",
#'   ylab = "FoldChange",
#'   same.y.lims = TRUE,
#'   y.max = 4
#' )
FeatureStatPlot <- function(
    srt,
    stat.by,
    group.by = NULL,
    split.by = NULL,
    bg.by = NULL,
    plot.by = c("group", "feature"),
    fill.by = c("group", "feature", "expression"),
    cells = NULL,
    layer = "data",
    assay = NULL,
    keep_empty = FALSE,
    individual = FALSE,
    plot_type = c("violin", "box", "bar", "dot", "col"),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    bg_palette = "Paired",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    add_box = FALSE,
    box_color = "black",
    box_width = 0.1,
    box_ptsize = 2,
    add_point = FALSE,
    pt.color = "grey30",
    pt.size = NULL,
    pt.alpha = 1,
    jitter.width = 0.4,
    jitter.height = 0.1,
    add_trend = FALSE,
    trend_color = "black",
    trend_linewidth = 1,
    trend_ptsize = 2,
    add_stat = c("none", "mean", "median"),
    stat_color = "black",
    stat_size = 1,
    stat_stroke = 1,
    stat_shape = 25,
    add_line = NULL,
    line_color = "red",
    line_size = 1,
    line_type = 1,
    cells.highlight = NULL,
    cols.highlight = "red",
    sizes.highlight = 1,
    alpha.highlight = 1,
    calculate_coexp = FALSE,
    same.y.lims = FALSE,
    y.min = NULL,
    y.max = NULL,
    y.trans = "identity",
    y.nbreaks = 5,
    sort = FALSE,
    stack = FALSE,
    flip = FALSE,
    comparisons = NULL,
    ref_group = NULL,
    pairwise_method = "wilcox.test",
    multiplegroup_comparisons = FALSE,
    multiple_method = "kruskal.test",
    sig_label = c("p.signif", "p.format"),
    sig_labelsize = 3.5,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = "Expression level",
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
  if (is.null(group.by)) {
    group.by <- "All.groups" # avoid having the same name with split.by. split.by will be All.groups by default
    xlab <- "All groups"
    srt[[group.by]] <- factor("All groups")
  }

  meta.data <- srt@meta.data
  meta.data[["cells"]] <- rownames(meta.data)
  assay <- assay %||% DefaultAssay(srt)
  exp.data <- GetAssayData5(srt, assay = assay, layer = layer)
  plot.by <- match.arg(plot.by)

  if (plot.by == "feature") {
    if (length(group.by) > 1) {
      stop(
        "The 'group.by' must have a length of 1 when 'plot.by' is set to 'feature'"
      )
    }
    if (!is.null(bg.by)) {
      message("'bg.by' is invalid when plot.by is set to 'feature'")
    }
    message("Setting 'group.by' to 'Features' as 'plot.by' is set to 'feature'")
    srt@assays[setdiff(names(srt@assays), assay)] <- NULL
    meta.reshape <- SeuratObject::FetchData(
      srt,
      vars = c(stat.by, group.by, split.by),
      cells = cells %||% rownames(meta.data),
      layer = layer
    )
    meta.reshape[["cells"]] <- rownames(meta.reshape)
    meta.reshape <- reshape2::melt(
      meta.reshape,
      measure.vars = stat.by,
      variable.name = "Features",
      value.name = "Stat.by"
    )
    rownames(meta.reshape) <- paste0(
      meta.reshape[["cells"]],
      "-",
      meta.reshape[["Features"]]
    )
    exp.data <- matrix(
      0,
      nrow = 1,
      ncol = nrow(meta.reshape),
      dimnames = list("Stat.by", rownames(meta.reshape))
    )
    plist <- list()
    for (g in unique(meta.reshape[[group.by]])) {
      if (length(rownames(meta.reshape)[meta.reshape[[group.by]] == g]) > 0) {
        meta.use <- meta.reshape
        meta.use[[group.by]] <- NULL
        colnames(meta.use)[colnames(meta.use) == "Stat.by"] <- g
        p <- ExpressionStatPlot(
          exp.data = exp.data,
          meta.data = meta.use,
          stat.by = g,
          group.by = "Features",
          split.by = split.by,
          bg.by = NULL,
          plot.by = "group",
          fill.by = fill.by,
          cells = rownames(meta.reshape)[meta.reshape[[group.by]] == g],
          keep_empty = keep_empty,
          individual = individual,
          plot_type = plot_type,
          palette = palette,
          palcolor = palcolor,
          alpha = alpha,
          bg_palette = bg_palette,
          bg_palcolor = bg_palcolor,
          bg_alpha = bg_alpha,
          add_box = add_box,
          box_color = box_color,
          box_width = box_width,
          box_ptsize = box_ptsize,
          add_point = add_point,
          pt.color = pt.color,
          pt.size = pt.size,
          pt.alpha = pt.alpha,
          jitter.width = jitter.width,
          jitter.height = jitter.height,
          add_trend = add_trend,
          trend_color = trend_color,
          trend_linewidth = trend_linewidth,
          trend_ptsize = trend_ptsize,
          add_stat = add_stat,
          stat_color = stat_color,
          stat_size = stat_size,
          stat_stroke = stat_stroke,
          stat_shape = stat_shape,
          add_line = add_line,
          line_color = line_color,
          line_size = line_size,
          line_type = line_type,
          cells.highlight = cells.highlight,
          cols.highlight = cols.highlight,
          sizes.highlight = sizes.highlight,
          alpha.highlight = alpha.highlight,
          calculate_coexp = calculate_coexp,
          same.y.lims = same.y.lims,
          y.min = y.min,
          y.max = y.max,
          y.trans = y.trans,
          y.nbreaks = y.nbreaks,
          sort = sort,
          stack = stack,
          flip = flip,
          comparisons = comparisons,
          ref_group = ref_group,
          pairwise_method = pairwise_method,
          multiplegroup_comparisons = multiplegroup_comparisons,
          multiple_method = multiple_method,
          sig_label = sig_label,
          sig_labelsize = sig_labelsize,
          aspect.ratio = aspect.ratio,
          title = title,
          subtitle = subtitle,
          xlab = xlab,
          ylab = ylab,
          legend.position = legend.position,
          legend.direction = legend.direction,
          theme_use = theme_use,
          theme_args = theme_args,
          force = force,
          seed = seed
        )
        plist <- append(plist, p)
      }
    }
    group.by <- "Features"
  } else {
    plist <- ExpressionStatPlot(
      exp.data = exp.data,
      meta.data = meta.data,
      stat.by = stat.by,
      group.by = group.by,
      split.by = split.by,
      bg.by = bg.by,
      plot.by = "group",
      fill.by = fill.by,
      cells = cells,
      keep_empty = keep_empty,
      individual = individual,
      plot_type = plot_type,
      palette = palette,
      palcolor = palcolor,
      alpha = alpha,
      bg_palette = bg_palette,
      bg_palcolor = bg_palcolor,
      bg_alpha = bg_alpha,
      add_box = add_box,
      box_color = box_color,
      box_width = box_width,
      box_ptsize = box_ptsize,
      add_point = add_point,
      pt.color = pt.color,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      jitter.width = jitter.width,
      jitter.height = jitter.height,
      add_trend = add_trend,
      trend_color = trend_color,
      trend_linewidth = trend_linewidth,
      trend_ptsize = trend_ptsize,
      add_stat = add_stat,
      stat_color = stat_color,
      stat_size = stat_size,
      stat_stroke = stat_stroke,
      stat_shape = stat_shape,
      add_line = add_line,
      line_color = line_color,
      line_size = line_size,
      line_type = line_type,
      cells.highlight = cells.highlight,
      cols.highlight = cols.highlight,
      sizes.highlight = sizes.highlight,
      alpha.highlight = alpha.highlight,
      calculate_coexp = calculate_coexp,
      same.y.lims = same.y.lims,
      y.min = y.min,
      y.max = y.max,
      y.trans = y.trans,
      y.nbreaks = y.nbreaks,
      sort = sort,
      stack = stack,
      flip = flip,
      comparisons = comparisons,
      ref_group = ref_group,
      pairwise_method = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons,
      multiple_method = multiple_method,
      sig_label = sig_label,
      sig_labelsize = sig_labelsize,
      aspect.ratio = aspect.ratio,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      force = force,
      seed = seed
    )
  }

  plist_stack <- list()
  if (isTRUE(stack) && length(stat.by) > 1 && isFALSE(individual)) {
    for (g in group.by) {
      plist_g <- plist[
        sapply(strsplit(names(plist), ":"), function(x) x[2]) == g
      ]
      legend <- get_legend(plist_g[[1]])
      if (isTRUE(flip)) {
        lab <- grid::textGrob(
          label = ifelse(is.null(ylab), "Expression level", ylab),
          hjust = 0.5
        )
        plist_g <- lapply(seq_along(plist_g), FUN = function(i) {
          p <- plist_g[[i]]
          if (i != 1) {
            suppressWarnings(
              p <- p +
                theme(
                  legend.position = "none",
                  panel.grid = element_blank(),
                  plot.title = element_blank(),
                  plot.subtitle = element_blank(),
                  axis.title = element_blank(),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(vjust = c(1, 0)),
                  axis.ticks.length.y = grid::unit(0, "pt"),
                  plot.margin = grid::unit(c(0, -0.5, 0, 0), "mm")
                )
            )
          } else {
            suppressWarnings(
              p <- p +
                theme(
                  legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(vjust = c(1, 0)),
                  axis.ticks.length.y = grid::unit(0, "pt"),
                  plot.margin = grid::unit(c(0, -0.5, 0, 0), "mm")
                )
            )
          }
          return(as_grob(p))
        })
        gtable <- do.call(cbind, plist_g)
        gtable <- add_grob(gtable, lab, "bottom", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      } else {
        lab <- grid::textGrob(
          label = ifelse(is.null(ylab), "Expression level", ylab),
          rot = 90,
          hjust = 0.5
        )
        plist_g <- lapply(seq_along(plist_g), FUN = function(i) {
          p <- plist_g[[i]]
          if (i != length(plist_g)) {
            suppressWarnings(
              p <- p +
                theme(
                  legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(vjust = c(0, 1)),
                  axis.ticks.length.x = grid::unit(0, "pt"),
                  plot.margin = grid::unit(c(-0.5, 0, 0, 0), "mm")
                )
            )
            if (i == 1) {
              p <- p +
                theme(
                  plot.title = element_blank(),
                  plot.subtitle = element_blank()
                )
            }
          } else {
            suppressWarnings(
              p <- p +
                theme(
                  legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(vjust = c(0, 1)),
                  axis.ticks.length.x = grid::unit(0, "pt"),
                  plot.margin = grid::unit(c(-0.5, 0, 0, 0), "mm")
                )
            )
          }
          return(as_grob(p))
        })
        gtable <- do.call(rbind, plist_g)
        gtable <- add_grob(gtable, lab, "left", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      }
      gtable <- gtable::gtable_add_padding(
        gtable,
        grid::unit(c(1, 1, 1, 1), units = "cm")
      )
      plot <- patchwork::wrap_plots(gtable)
      plist_stack[[g]] <- plot
    }
  }

  if (length(plist_stack) > 0) {
    plist <- plist_stack
  }
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
