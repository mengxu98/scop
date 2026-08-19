#' @title Feature statistical plots
#'
#' @md
#' @inheritParams scop-params
#' @inheritParams CellDimPlot
#' @inheritParams FeatureDimPlot
#' @param stat.by Features to plot.
#' @param plot.by `"group"` or `"feature"`.
#' @param bg.by Metadata column used as background color.
#' @param fill.by `"group"`, `"feature"`, or `"expression"`.
#' @param cells Cell names to include.
#' @param keep_empty Keep empty factor levels.
#' @param individual One plot per group.
#' @param plot_type `"violin"`, `"box"`, `"bar"`, `"dot"`, or `"col"`.
#' @param alpha Plot transparency.
#' @param bg_palette,bg_palcolor,bg_alpha Background palette and transparency.
#' @param add_box,box_color,box_width,box_ptsize Overlay boxplot.
#' @param add_point,pt.color,jitter.width,jitter.height Overlay jittered points.
#' @param add_trend,trend_color,trend_linewidth,trend_ptsize Overlay trend line.
#' @param add_stat,stat_color,stat_size,stat_stroke,stat_shape
#' Summary statistic (`"none"`, `"mean"`, or `"median"`).
#' @param add_line,line_color,line_size,line_type Horizontal line at this y-intercept.
#' @param cols.highlight,sizes.highlight,alpha.highlight Highlighted cells.
#' @param calculate_coexp Plot the geometric mean of `stat.by`.
#' @param same.y.lims Share y-axis limits across panels.
#' @param y.min,y.max Y-axis limits. Character values like `"q5"` use that quantile.
#' @param y.trans,y.nbreaks Y-axis transform (`"identity"` or `"log2"`) and breaks.
#' @param sort Sort groups on the x-axis: `FALSE`, `TRUE`/`"increasing"`, or `"decreasing"`.
#' @param stack Stack plots vertically.
#' @param flip Flip x and y.
#' @param comparisons,ref_group Pairwise comparisons (length-2 name or index vectors).
#' @param auto_comparison Compare the highest-median group (or `ref_group`) against
#' the others. Requires `split.by = NULL`.
#' @param pairwise_method,multiplegroup_comparisons,multiple_method Comparison tests.
#' @param sig_label,sig_labelsize Significance labels (`"p.signif"` or `"p.format"`).
#' @param aspect.ratio Panel aspect ratio.
#' @param x_text_angle Rotation of x-axis labels.
#' @param ylab Y-axis label.
#' @param grid_major,grid_major_colour,grid_major_linetype,grid_major_linewidth
#' Major panel grid lines.
#' @param ... Additional arguments passed to the plotting helpers.
#'
#' @seealso
#' [CellStatPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType"
#' ) |> thisplot::panel_fix(height = 1, width = 2)
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "box"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   plot_type = "bar"
#' )
#'
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
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_box = TRUE
#' )
#'
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
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_stat = "mean"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   add_line = 0.2,
#'   line_type = 2
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   add_box = TRUE,
#'   add_trend = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("G2M_score", "Fev"),
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   comparisons = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   fill.by = "expression",
#'   palette = "Blues",
#'   same.y.lims = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   multiplegroup_comparisons = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   auto_comparison = TRUE,
#'   sig_label = "p.signif"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta"))
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")),
#'   sig_label = "p.format"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = "Fev",
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   comparisons = TRUE
#' ) + FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = "Fev",
#'   group.by = "SubCellType",
#'   split.by = "Phase",
#'   comparisons = TRUE,
#'   y.max = 5
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Rbp4", "Pyy"),
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   add_box = TRUE, stack = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     # Ductal
#'     "Sox9", "Anxa2", "Bicc1",
#'     # EPs
#'     "Neurog3", "Hes6",
#'     # Pre-endocrine
#'     "Fev", "Neurod1",
#'     # Endocrine
#'     "Rbp4", "Pyy",
#'     # Beta, Alpha, Delta, Epsilon
#'     "Ins1", "Gcg", "Sst", "Ghrl"
#'   ),
#'   legend.position = "top",
#'   legend.direction = "horizontal",
#'   group.by = "SubCellType",
#'   bg.by = "CellType",
#'   stack = TRUE
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     # Ductal
#'     "Sox9", "Anxa2", "Bicc1",
#'     # EPs
#'     "Neurog3", "Hes6",
#'     # Pre-endocrine
#'     "Fev", "Neurod1",
#'     # Endocrine
#'     "Rbp4", "Pyy",
#'     # Beta, Alpha, Delta, Epsilon
#'     "Ins1", "Gcg", "Sst", "Ghrl"
#'   ),
#'   fill.by = "feature",
#'   plot_type = "box",
#'   group.by = "SubCellType",
#'   bg.by = "CellType", stack = TRUE, flip = TRUE
#' ) |> thisplot::panel_fix_overall(
#'   width = 8, height = 5
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "group"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   multiplegroup_comparisons = TRUE,
#'   sig_label = "p.format",
#'   sig_labelsize = 4
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   comparisons = list(
#'     c("Neurog3", "Rbp4"),
#'     c("Rbp4", "Ins1")
#'   ),
#'   stack = TRUE
#' )
#'
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c(
#'     # Ductal
#'     "Sox9", "Anxa2", "Bicc1",
#'     # EPs
#'     "Neurog3", "Hes6",
#'     # Pre-endocrine
#'     "Fev", "Neurod1",
#'     # Endocrine
#'     "Rbp4", "Pyy",
#'     # Beta, Alpha, Delta, Epsilon
#'     "Ins1", "Gcg", "Sst", "Ghrl"
#'   ),
#'   group.by = "SubCellType",
#'   plot.by = "feature",
#'   stack = TRUE
#' )
#'
#' data <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "data"
#' )
#' pancreas_sub <- SeuratObject::SetAssayData(
#'   object = pancreas_sub,
#'   layer = "scale.data",
#'   assay = "RNA",
#'   new.data = data / Matrix::rowMeans(data)
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
  palette = "Chinese",
  palcolor = NULL,
  alpha = 1,
  bg_palette = "Chinese",
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
  auto_comparison = FALSE,
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
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11,
  ...,
  x_text_angle = 45,
  verbose = TRUE
) {
  if (is.null(group.by)) {
    group.by <- "All.groups"
    xlab <- "All groups"
    srt[[group.by]] <- factor("All groups")
  }

  meta.data <- srt@meta.data
  meta.data[["cells"]] <- rownames(meta.data)
  assay <- assay %||% DefaultAssay(srt)
  exp.data <- GetAssayData5(
    srt,
    assay = assay,
    layer = layer
  )
  plot.by <- match.arg(plot.by)
  expression_stat_args <- list(
    split.by = split.by,
    plot.by = "group",
    fill.by = fill.by,
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
    auto_comparison = auto_comparison,
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
    legend.title = legend.title,
    theme_use = theme_use,
    theme_args = theme_args,
    grid_major = grid_major,
    grid_major_colour = grid_major_colour,
    grid_major_linetype = grid_major_linetype,
    grid_major_linewidth = grid_major_linewidth,
    force = force,
    seed = seed,
    x_text_angle = x_text_angle
  )

  if (plot.by == "feature") {
    if (length(group.by) > 1) {
      log_message(
        "{.arg group.by} must have a length of 1 when {.arg plot.by} is set to {.val feature}",
        message_type = "error"
      )
    }
    if (!is.null(bg.by)) {
      log_message(
        "{.arg bg.by} is invalid when {.arg plot.by} is set to {.val feature}",
        message_type = "warning",
        verbose = verbose
      )
    }
    log_message(
      "Setting {.arg group.by} to {.val Features} as {.arg plot.by} is set to {.val feature}",
      verbose = verbose
    )
    srt@assays[setdiff(names(srt@assays), assay)] <- NULL
    meta_reshape <- SeuratObject::FetchData(
      srt,
      vars = c(stat.by, group.by, split.by),
      cells = cells %||% rownames(meta.data),
      layer = layer
    )
    meta_reshape[["cells"]] <- rownames(meta_reshape)
    meta_reshape <- reshape2::melt(
      meta_reshape,
      measure.vars = stat.by,
      variable.name = "Features",
      value.name = "Stat.by"
    )
    rownames(meta_reshape) <- paste0(
      meta_reshape[["cells"]],
      "-",
      meta_reshape[["Features"]]
    )
    exp.data <- matrix(
      0,
      nrow = 1,
      ncol = nrow(meta_reshape),
      dimnames = list("Stat.by", rownames(meta_reshape))
    )
    plist <- list()
    for (g in unique(meta_reshape[[group.by]])) {
      cells_g <- rownames(meta_reshape)[meta_reshape[[group.by]] == g]
      if (length(cells_g) > 0) {
        meta_use <- meta_reshape
        meta_use[[group.by]] <- NULL
        colnames(meta_use)[colnames(meta_use) == "Stat.by"] <- g
        p <- do.call(
          ExpressionStatPlot,
          c(
            list(
              exp.data = exp.data,
              meta.data = meta_use,
              stat.by = g,
              group.by = "Features",
              bg.by = NULL,
              cells = cells_g
            ),
            expression_stat_args
          )
        )
        plist <- append(plist, p)
      }
    }
    group.by <- "Features"
  } else {
    plist <- do.call(
      ExpressionStatPlot,
      c(
        list(
          exp.data = exp.data,
          meta.data = meta.data,
          stat.by = stat.by,
          group.by = group.by,
          bg.by = bg.by,
          cells = cells
        ),
        expression_stat_args
      )
    )
  }

  plist_stack <- list()
  if (isTRUE(stack) && length(stat.by) > 1 && isFALSE(individual)) {
    theme_stack <- tryCatch(
      do.call(theme_use, theme_args),
      error = function(e) NULL
    )
    `%||%` <- function(x, y) if (is.null(x)) y else x
    element_text_to_gpar <- function(el) {
      if (is.null(el) || !inherits(el, "element_text")) {
        return(NULL)
      }
      col_use <- el$colour %||% el$color
      grid::gpar(
        fontsize = el$size,
        col = col_use,
        fontfamily = el$family,
        fontface = el$face
      )
    }
    axis_title_y_gp <- element_text_to_gpar(
      if (!is.null(theme_stack)) {
        ggplot2::calc_element("axis.title.y", theme_stack)
      } else {
        NULL
      }
    )
    axis_title_x_gp <- element_text_to_gpar(
      if (!is.null(theme_stack)) {
        ggplot2::calc_element("axis.title.x", theme_stack)
      } else {
        NULL
      }
    )
    plot_title_el <- NULL
    if (
      !is.null(theme_args[["plot.title"]]) &&
        inherits(theme_args[["plot.title"]], "element_text")
    ) {
      plot_title_el <- theme_args[["plot.title"]]
    } else if (
      !is.null(theme_args[["title"]]) &&
        inherits(theme_args[["title"]], "element_text")
    ) {
      plot_title_el <- theme_args[["title"]]
    } else if (!is.null(theme_stack)) {
      plot_title_el <- ggplot2::calc_element("plot.title", theme_stack)
    }
    plot_title_gp <- element_text_to_gpar(plot_title_el)

    for (g in group.by) {
      if (is.null(names(plist))) {
        plist_g <- plist
      } else {
        plist_g <- plist[
          sapply(strsplit(names(plist), ":"), function(x) x[2]) == g
        ]
        if (length(plist_g) == 0) {
          plist_g <- plist
        }
      }
      legend <- get_legend(plist_g[[1]])
      if (isTRUE(flip)) {
        lab <- grid::textGrob(
          label = ifelse(is.null(ylab), "Expression level", ylab),
          hjust = 0.5,
          gp = axis_title_x_gp
        )
        plist_g <- lapply(
          seq_along(plist_g),
          FUN = function(i) {
            p <- plist_g[[i]]
            if (i != 1) {
              suppressWarnings(
                p <- p +
                  theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    plot.title = element_blank(),
                    plot.subtitle = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
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
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(vjust = c(1, 0)),
                    axis.ticks.length.y = grid::unit(0, "pt"),
                    plot.margin = grid::unit(c(0, -0.5, 0, 0), "mm")
                  )
              )
            }
            p <- p +
              theme(
                plot.title = element_blank(),
                plot.subtitle = element_blank()
              )
            return(as_grob(p))
          }
        )
        gtable <- do.call(cbind, plist_g)
        gtable <- add_grob(gtable, lab, "bottom", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      } else {
        lab <- grid::textGrob(
          label = ifelse(is.null(ylab), "Expression level", ylab),
          rot = 90,
          hjust = 0.5,
          gp = axis_title_y_gp
        )
        plist_g <- lapply(
          seq_along(plist_g),
          FUN = function(i) {
            p <- plist_g[[i]]
            if (i != length(plist_g)) {
              suppressWarnings(
                p <- p +
                  theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(vjust = c(0, 1)),
                    axis.ticks.length.x = grid::unit(0, "pt"),
                    plot.margin = grid::unit(c(-0.5, 0, 0, 0), "mm")
                  )
              )
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
            p <- p +
              theme(
                plot.title = element_blank(),
                plot.subtitle = element_blank()
              )
            return(as_grob(p))
          }
        )
        gtable <- do.call(rbind, plist_g)
        gtable <- add_grob(gtable, lab, "left", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      }
      if (!is.null(title)) {
        title_grob <- grid::textGrob(
          title,
          x = 0,
          hjust = 0,
          gp = plot_title_gp
        )
        title_height <- grid::grobHeight(title_grob) + grid::unit(0.5, "lines")
        gtable <- add_grob(
          gtable,
          title_grob,
          "top",
          title_height,
          clip = "off"
        )
      }
      check_r("gtable", verbose = FALSE)
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
  combine_plot_list(
    plist,
    combine = combine, nrow = nrow, ncol = ncol, byrow = byrow
  )
}
