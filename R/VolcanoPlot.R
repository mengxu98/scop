#' VolcanoPlot
#'
#' Generate a volcano plot based on differential expression analysis results.
#'
#' @md
#' @param srt An object of class `Seurat` containing the results of differential expression analysis.
#' @param group_by A character vector specifying the column in `srt` to group the samples by.
#' Default is `NULL`.
#' @param test.use A character string specifying the type of statistical test to use.
#' Default is "wilcox".
#' @param DE_threshold A character string specifying the threshold for differential expression.
#' Default is "avg_log2FC > 0 & p_val_adj < 0.05".
#' @param x_metric A character string specifying the metric to use for the x-axis.
#' Default is "diff_pct".
#' @param palette A character string specifying the color palette to use for the plot.
#' Default is "RdBu".
#' @param palcolor A character string specifying the color for the palette.
#' Default is `NULL`.
#' @param pt.size A numeric value specifying the size of the points.
#' Default is 1.
#' @param pt.alpha A numeric value specifying the transparency of the points.
#' Default is 1.
#' @param cols.highlight A character string specifying the color for highlighted points.
#' Default is "black".
#' @param sizes.highlight A numeric value specifying the size of the highlighted points.
#' Default is 1.
#' @param alpha.highlight A numeric value specifying the transparency of the highlighted points.
#' Default is 1.
#' @param stroke.highlight A numeric value specifying the stroke width for the highlighted points.
#' Default is 0.5.
#' @param nlabel An integer value specifying the number of labeled points per group.
#' Default is 5.
#' @param features_label A character vector specifying the feature labels to plot.
#' Default is `NULL`.
#' @param label.fg A character string specifying the color for the labels' foreground.
#' Default is "black".
#' @param label.bg A character string specifying the color for the labels' background.
#' Default is "white".
#' @param label.bg.r A numeric value specifying the radius of the rounding of the labels' background.
#' Default is 0.1.
#' @param label.size A numeric value specifying the size of the labels.
#' Default is 4.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' Default is `NULL`.
#' @param xlab A character string specifying the x-axis label.
#' Default is the value of `x_metric`.
#' @param ylab A character string specifying the y-axis label.
#' Default is "-log10(p-adjust)".
#' @param theme_use A character string specifying the theme to use for the plot.
#' Default is "theme_scop".
#' @param theme_args A list of theme arguments to pass to the `theme_use` function.
#' Default is an empty list.
#' @param combine A logical value indicating whether to combine the plots for each group into a single plot.
#' Default is `TRUE`.
#' @param nrow An integer value specifying the number of rows in the combined plot.
#' Default is `NULL`.
#' @param ncol An integer value specifying the number of columns in the combined plot.
#' Default is `NULL`.
#' @param byrow A logical value indicating whether to arrange the plots by row in the combined plot.
#' Default is `TRUE`.
#'
#' @seealso \code{\link{RunDEtest}}
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' # pancreas_sub <- RunDEtest(
#' #  pancreas_sub,
#' #  group_by = "CellType"
#' # )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   ncol = 2
#' )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05",
#'   ncol = 2
#' )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group_by = "CellType",
#'   x_metric = "avg_log2FC",
#'   ncol = 2
#' )
VolcanoPlot <- function(
    srt,
    group_by = NULL,
    test.use = "wilcox",
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    x_metric = "diff_pct",
    palette = "RdBu",
    palcolor = NULL,
    pt.size = 1,
    pt.alpha = 1,
    cols.highlight = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    nlabel = 5,
    features_label = NULL,
    label.fg = "black",
    label.bg = "white",
    label.bg.r = 0.1,
    label.size = 4,
    aspect.ratio = NULL,
    xlab = x_metric,
    ylab = "-log10(p-adjust)",
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE) {
  if (is.null(group_by)) {
    group_by <- "custom"
  }
  layer <- paste0("DEtest_", group_by)
  if (
    !layer %in% names(srt@tools) ||
      length(grep(pattern = "AllMarkers", names(srt@tools[[layer]]))) == 0
  ) {
    log_message(
      "Cannot find the DEtest result for the group '",
      group_by,
      "'. You may perform RunDEtest first.",
      message_type = "error"
    )
  }
  index <- grep(
    pattern = paste0("AllMarkers_", test.use),
    names(srt@tools[[layer]])
  )[1]
  if (is.na(index)) {
    log_message(
      "Cannot find the 'AllMarkers_", test.use, "' in the DEtest result.",
      message_type = "error"
    )
  }
  de <- names(srt@tools[[layer]])[index]
  de_df <- srt@tools[[layer]][[de]]
  de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
  de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
  de_df[, "DE"] <- FALSE
  de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE

  x_upper <- stats::quantile(
    de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])],
    c(0.99, 1)
  )
  x_lower <- stats::quantile(
    de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])],
    c(0.01, 0)
  )
  x_upper <- ifelse(x_upper[1] > 0, x_upper[1], x_upper[2])
  x_lower <- ifelse(x_lower[1] < 0, x_lower[1], x_lower[2])
  if (x_upper > 0 & x_lower < 0) {
    value_range <- min(abs(c(x_upper, x_lower)), na.rm = TRUE)
    x_upper <- value_range
    x_lower <- -value_range
  }

  de_df[, "border"] <- FALSE
  de_df[de_df[["avg_log2FC"]] > x_upper, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] > x_upper, "avg_log2FC"] <- x_upper
  de_df[de_df[["avg_log2FC"]] < x_lower, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] < x_lower, "avg_log2FC"] <- x_lower

  de_df[, "y"] <- -log10(de_df[, "p_val_adj"])
  if (x_metric == "diff_pct") {
    de_df[, "x"] <- de_df[, "diff_pct"]
    de_df[de_df[, "avg_log2FC"] < 0, "y"] <- -de_df[
      de_df[, "avg_log2FC"] < 0,
      "y"
    ]
    de_df <- de_df[
      order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), ,
      drop = FALSE
    ]
  } else if (x_metric == "avg_log2FC") {
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df[de_df[, "diff_pct"] < 0, "y"] <- -de_df[de_df[, "diff_pct"] < 0, "y"]
    de_df <- de_df[
      order(abs(de_df[, "diff_pct"]), decreasing = FALSE, na.last = FALSE), ,
      drop = FALSE
    ]
  }
  de_df[, "distance"] <- de_df[, "x"]^2 + de_df[, "y"]^2

  plist <- list()
  for (group in levels(de_df[["group1"]])) {
    df <- de_df[de_df[["group1"]] == group, , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }
    x_nudge <- diff(range(df$x)) * 0.05
    df[, "label"] <- FALSE
    if (is.null(features_label)) {
      df[df[["y"]] >= 0, ][
        utils::head(order(df[df[["y"]] >= 0, "distance"], decreasing = TRUE), nlabel),
        "label"
      ] <- TRUE
      df[df[["y"]] < 0, ][
        utils::head(order(df[df[["y"]] < 0, "distance"], decreasing = TRUE), nlabel),
        "label"
      ] <- TRUE
    } else {
      df[df[["gene"]] %in% features_label, "label"] <- TRUE
    }
    jitter <- position_jitter(width = 0.2, height = 0.2, seed = 11)
    color_by <- ifelse(x_metric == "diff_pct", "avg_log2FC", "diff_pct")
    p <- ggplot() +
      geom_point(
        data = df[!df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x, y = y, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_point(
        data = df[!df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x, y = y, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha,
        position = jitter
      ) +
      geom_point(
        data = df[df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x, y = y),
        color = cols.highlight,
        size = sizes.highlight + stroke.highlight,
        alpha = alpha.highlight
      ) +
      geom_point(
        data = df[df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x, y = y),
        color = cols.highlight,
        size = sizes.highlight + stroke.highlight,
        alpha = alpha.highlight,
        position = jitter
      ) +
      geom_point(
        data = df[df[["DE"]] & !df[["border"]], , drop = FALSE],
        aes(x = x, y = y, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha
      ) +
      geom_point(
        data = df[df[["DE"]] & df[["border"]], , drop = FALSE],
        aes(x = x, y = y, color = .data[[color_by]]),
        size = pt.size,
        alpha = pt.alpha,
        position = jitter
      ) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      geom_vline(xintercept = 0, color = "grey", linetype = 2) +
      ggrepel::geom_text_repel(
        data = df[df[["label"]], , drop = FALSE],
        aes(x = x, y = y, label = gene),
        min.segment.length = 0,
        max.overlaps = 100,
        segment.colour = "grey40",
        color = label.fg,
        bg.color = label.bg,
        bg.r = label.bg.r,
        size = label.size,
        force = 20,
        nudge_x = ifelse(df[df[["label"]], "y"] >= 0, -x_nudge, x_nudge)
      ) +
      labs(x = xlab, y = ylab) +
      scale_color_gradientn(
        name = ifelse(x_metric == "diff_pct", "log2FC", "diff_pct"),
        colors = palette_scop(palette = palette, palcolor = palcolor),
        values = rescale(unique(c(
          min(c(df[, color_by], 0), na.rm = TRUE),
          0,
          max(df[, color_by], na.rm = TRUE)
        ))),
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          title.hjust = 0,
          order = 1
        )
      ) +
      scale_y_continuous(labels = abs) +
      facet_wrap(~group1) +
      do.call(theme_use, theme_args) +
      theme(aspect.ratio = aspect.ratio)
    plist[[group]] <- p
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
