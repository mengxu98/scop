#' @title Differential Expression Test Plot
#'
#' @md
#' @inheritParams CellDimPlot
#' @param srt An object of class `Seurat` containing the results of differential expression analysis.
#' @param res A `data.frame` or `data.table` with differential expression results.
#' When `res` is provided, `srt` will be ignored.
#' The data.frame must contain columns: `gene`, `group1` (factor or character),
#' `avg_log2FC`, `p_val_adj`, and optionally `pct.1` and `pct.2` for calculating `diff_pct`.
#' @param test.use A character string specifying the type of statistical test to use.
#' Default is `"wilcox"`.
#' @param plot_type Type of plot to create. Options: `"volcano"`, `"manhattan"`, or `"ring"`.
#' Default is `"volcano"`.
#' @param DE_threshold A character string specifying the threshold for differential expression (used to highlight significant genes in all plot types).
#' Default is `"avg_log2FC > 0 & p_val_adj < 0.05"`.
#' @param x_metric A character string specifying the metric to use for the x-axis (only for volcano plot).
#' Default is `"diff_pct"`.
#' @param y_metric A character string specifying the metric to use for the y-axis (only for Manhattan plot, not used currently).
#' Options: `"p_val"` or `"p_val_adj"`. Default is `"p_val_adj"`.
#' @param x_order A character string specifying how to order genes on x-axis (only for Manhattan plot, not used currently).
#' Options: `"gene"` (alphabetical by gene name) or `"index"` (by data order). Default is `"gene"`.
#' @param palette Color palette name.
#' Available palettes can be found in [thisplot::show_palettes].
#' Default is `"RdBu"`.
#' @param group_palette Palette for cell types (groups) in Manhattan plot.
#' Default is `"Chinese"`.
#' @param group_palcolor Custom colors for cell types (groups) in Manhattan plot.
#' Default is `NULL`.
#' @param pt.size The size of the points.
#' Default is `1`.
#' @param cols.highlight A character string specifying the color for highlighted points.
#' Default is `"black"`.
#' @param sizes.highlight The size of the highlighted points.
#' Default is `1`.
#' @param alpha.highlight The transparency of the highlighted points.
#' Default is `1`.
#' @param stroke.highlight The stroke width for the highlighted points.
#' Default is `0.5`.
#' @param nlabel An integer value specifying the number of labeled points per group.
#' Default is `5`.
#' @param features_label A character vector specifying the feature labels to plot.
#' Default is `NULL`.
#' @param label.fg A character string specifying the color for the labels' foreground.
#' Default is `"black"`.
#' @param label.bg A character string specifying the color for the labels' background.
#' Default is `"white"`.
#' @param label.bg.r The radius of the rounding of the labels' background.
#' Default is `0.1`.
#' @param label.size The size of the labels.
#' Default is `4`.
#' @param aspect.ratio Aspect ratio of the panel.
#' Default is `NULL`.
#' @param xlab A character string specifying the x-axis label.
#' @param ylab A character string specifying the y-axis label.
#' @param theme_use Theme to use for the plot.
#' Default is `"theme_scop"`.
#' @param theme_args A list of additional arguments to pass to the theme function.
#' Default is `list()`.
#' @param combine Whether to combine multiple plots into one.
#' Default is `TRUE`.
#' @param nrow Number of rows for combined plots.
#' Default is `NULL`.
#' @param ncol Number of columns for combined plots.
#' Default is `NULL`.
#' @param byrow Whether to fill plots by row.
#' Default is `TRUE`.
#' @param manhattan.bg Background color for Manhattan plot.
#' Default is `"white"`.
#' @param jitter_width Horizontal jitter range for points in Manhattan plot.
#' Default is `0.5`.
#' @param jitter_height Vertical jitter range for points in Manhattan plot.
#' Default is `0.4`.
#' @param tile_height Height of the cell-type track in ring plot.
#' Default is `0.3`.
#' @param tile_gap Gap between the track and nudged points in ring plot.
#' Default is `0.1`.
#' @param ring_segments Whether to draw segment lines between cell types in ring plot.
#' Default is `TRUE`.
#' @param seed Random seed for jitter in ring plot.
#' Default is `11`.
#'
#' @seealso [RunDEtest], [VolcanoPlot], [DEtestManhattanPlot], [DEtestRingPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   only.pos = FALSE
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "volcano",
#'   ncol = 2
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "manhattan"
#' )
#'
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "ring"
#' )
#'
#' de_results1 <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
#' DEtestPlot(
#'   res = de_results1,
#'   plot_type = "volcano",
#'   ncol = 2
#' )
#'
#' de_results2 <- Seurat::FindMarkers(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   ident.1 = "Ductal",
#'   ident.2 = "Endocrine"
#' )
#' DEtestPlot(
#'   res = de_results2,
#'   plot_type = "volcano"
#' )
#'
#' de_results3 <- Seurat::FindAllMarkers(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' DEtestPlot(
#'   res = de_results3,
#'   plot_type = "volcano",
#'   ncol = 2
#' )
DEtestPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
    plot_type = c("volcano", "manhattan", "ring"),
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    x_metric = "diff_pct",
    y_metric = c("p_val_adj", "p_val"),
    x_order = c("gene", "index"),
    palette = "RdBu",
    palcolor = NULL,
    group_palette = "Chinese",
    group_palcolor = NULL,
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
    xlab = NULL,
    ylab = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    manhattan.bg = "white",
    jitter_width = 0.5,
    jitter_height = 0.4,
    tile_height = 0.3,
    tile_gap = 0.1,
    ring_segments = TRUE,
    seed = 11) {
  plot_type <- match.arg(plot_type)
  y_metric <- match.arg(y_metric)
  x_order <- match.arg(x_order)

  if (plot_type == "volcano") {
    return(
      VolcanoPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        DE_threshold = DE_threshold,
        x_metric = x_metric,
        palette = palette,
        palcolor = palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        aspect.ratio = aspect.ratio,
        xlab = xlab %||% x_metric,
        ylab = ylab %||% "-log10(p-adjust)",
        theme_use = theme_use,
        theme_args = theme_args,
        combine = combine,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      )
    )
  }

  if (plot_type == "manhattan") {
    return(
      DEtestManhattanPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        DE_threshold = DE_threshold,
        group_palette = group_palette,
        group_palcolor = group_palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args,
        manhattan.bg = manhattan.bg,
        jitter_width = jitter_width,
        jitter_height = jitter_height,
        aspect.ratio = aspect.ratio,
        xlab = xlab,
        ylab = ylab
      )
    )
  }
  if (plot_type == "ring") {
    return(
      DEtestRingPlot(
        srt = srt,
        group.by = group.by,
        test.use = test.use,
        res = res,
        DE_threshold = DE_threshold,
        group_palette = group_palette,
        group_palcolor = group_palcolor,
        pt.size = pt.size,
        pt.alpha = pt.alpha,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        alpha.highlight = alpha.highlight,
        stroke.highlight = stroke.highlight,
        nlabel = nlabel,
        features_label = features_label,
        label.fg = label.fg,
        label.bg = label.bg,
        label.bg.r = label.bg.r,
        label.size = label.size,
        palette = palette,
        palcolor = palcolor,
        theme_use = theme_use,
        theme_args = theme_args,
        tile_height = tile_height,
        tile_gap = tile_gap,
        jitter_width = jitter_width,
        ring_segments = ring_segments,
        seed = seed
      )
    )
  }
  invisible(NULL)
}

get_de_data <- function(
    srt,
    group.by,
    test.use,
    DE_threshold,
    res = NULL) {
  if (!is.null(res)) {
    de_df <- as.data.frame(res)
    if (!"gene" %in% colnames(de_df)) {
      if (nrow(de_df) > 0 && !is.null(rownames(de_df))) {
        de_df[, "gene"] <- rownames(de_df)
      } else {
        log_message(
          "Missing required column: {.val gene}. Please provide gene names in a 'gene' column or as row names.",
          message_type = "error"
        )
      }
    }
    required_cols <- c("avg_log2FC", "p_val_adj")
    missing_cols <- setdiff(required_cols, colnames(de_df))
    if (length(missing_cols) > 0) {
      log_message(
        "Missing required columns in data.frame: {.val {missing_cols}}",
        message_type = "error"
      )
    }
    if (!"group1" %in% colnames(de_df)) {
      if ("cluster" %in% colnames(de_df)) {
        de_df[, "group1"] <- de_df[, "cluster"]
      } else {
        de_df[, "group1"] <- "All"
      }
    }
    if (!is.factor(de_df[["group1"]])) {
      de_df[["group1"]] <- factor(
        de_df[["group1"]],
        levels = unique(de_df[["group1"]])
      )
    }
    if ("pct.1" %in% colnames(de_df) && "pct.2" %in% colnames(de_df)) {
      de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
    } else {
      de_df[, "diff_pct"] <- 0
    }
    de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
    de_df[, "DE"] <- FALSE
    de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE
    return(list(de_df = de_df))
  }

  if (is.null(group.by)) {
    group.by <- "custom"
  }
  layer <- paste0("DEtest_", group.by)
  if (
    !layer %in% names(srt@tools) ||
      length(grep(pattern = "AllMarkers", names(srt@tools[[layer]]))) == 0
  ) {
    log_message(
      "Cannot find the DEtest result for the group {.val {group.by}}. Perform {.fn RunDEtest} first",
      message_type = "error"
    )
  }
  index <- grep(
    pattern = paste0("AllMarkers_", test.use),
    names(srt@tools[[layer]])
  )[1]
  if (is.na(index)) {
    log_message(
      "Cannot find the {.val AllMarkers_{test.use}} in the DEtest result",
      message_type = "error"
    )
  }
  de <- names(srt@tools[[layer]])[index]
  de_df <- srt@tools[[layer]][[de]]
  de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
  de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
  de_df[, "DE"] <- FALSE
  de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE
  list(de_df = de_df)
}

clip_log2fc_symmetric <- function(df, fc_col = "avg_log2FC") {
  fc <- df[[fc_col]][is.finite(df[[fc_col]])]
  if (length(fc) == 0) {
    return(list(df = df, fc_lim = c(-1, 1)))
  }
  x_upper <- stats::quantile(fc, c(0.99, 1))
  x_lower <- stats::quantile(fc, c(0.01, 0))
  x_upper <- ifelse(x_upper[1] > 0, x_upper[1], x_upper[2])
  x_lower <- ifelse(x_lower[1] < 0, x_lower[1], x_lower[2])
  if (x_upper > 0 && x_lower < 0) {
    value_range <- min(abs(c(x_upper, x_lower)), na.rm = TRUE)
    x_upper <- value_range
    x_lower <- -value_range
  }
  df[df[[fc_col]] > x_upper, fc_col] <- x_upper
  df[df[[fc_col]] < x_lower, fc_col] <- x_lower
  list(df = df, fc_lim = c(x_lower, x_upper))
}

filter_de_markers <- function(de_df, log2FC_cutoff, pvalue_cutoff) {
  de_df_marker <- de_df[
    abs(de_df[, "avg_log2FC"]) >= log2FC_cutoff &
      de_df[, "p_val"] < pvalue_cutoff, ,
    drop = FALSE
  ]
  if (nrow(de_df_marker) == 0) {
    log_message(
      "No genes pass the threshold. Please adjust the threshold.",
      message_type = "warning"
    )
    return(NULL)
  }
  de_df_marker
}

get_top_markers_for_label <- function(de_df_marker, cluster_levels, nlabel, features_label) {
  top_marker_list <- lapply(cluster_levels, function(x) {
    tmp <- de_df_marker[de_df_marker[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) {
      return(NULL)
    }
    if (is.null(features_label)) {
      top_max <- if (nrow(tmp) > 0) {
        tmp[utils::head(order(tmp[, "avg_log2FC"], decreasing = TRUE), nlabel), , drop = FALSE]
      } else {
        tmp[0, , drop = FALSE]
      }
      top_min <- if (nrow(tmp) > 0) {
        tmp[utils::head(order(tmp[, "avg_log2FC"], decreasing = FALSE), nlabel), , drop = FALSE]
      } else {
        tmp[0, , drop = FALSE]
      }
      rbind(top_max, top_min)
    } else {
      tmp[tmp[, "gene"] %in% features_label, , drop = FALSE]
    }
  })
  do.call(rbind, top_marker_list[!sapply(top_marker_list, is.null)])
}

#' @title DEtest Manhattan Plot
#'
#' @description
#' Draw a Manhattan-style plot of differential expression results by cell type.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [VolcanoPlot], [DEtestRingPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   only.pos = FALSE
#' )
#' DEtestManhattanPlot(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
DEtestManhattanPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    group_palette = "Chinese",
    group_palcolor = NULL,
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
    palette = "RdBu",
    palcolor = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    manhattan.bg = "white",
    jitter_width = 0.5,
    jitter_height = 0.4,
    aspect.ratio = NULL,
    xlab = NULL,
    ylab = NULL) {
  if (is.null(group.by)) {
    group.by <- "custom"
  }
  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res)
  de_df <- data_res$de_df
  clip_res <- clip_log2fc_symmetric(de_df)
  de_df_marker <- clip_res$df
  fc_lim <- clip_res$fc_lim
  de_df_marker[, "type"] <- ifelse(
    de_df_marker[, "avg_log2FC"] >= 0,
    "sigUp",
    "sigDown"
  )
  cluster_levels <- levels(de_df_marker[["group1"]])
  n_m <- nrow(de_df_marker)
  de_df_marker[, "x_num"] <- as.numeric(de_df_marker[, "group1"])
  de_df_marker[, "x_plot"] <- de_df_marker[, "x_num"] + (stats::runif(n_m) - 0.5) * jitter_width
  de_df_marker[, "y_plot"] <- de_df_marker[, "avg_log2FC"] + (stats::runif(n_m) - 0.5) * jitter_height
  top_marker <- get_top_markers_for_label(de_df_marker, cluster_levels, nlabel, features_label)
  back_data_list <- lapply(cluster_levels, function(x) {
    tmp <- de_df_marker[de_df_marker[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) {
      return(NULL)
    }
    data.frame(
      cluster = x,
      min = min(tmp[, "avg_log2FC"], na.rm = TRUE) - 0.2,
      max = max(tmp[, "avg_log2FC"], na.rm = TRUE) + 0.2
    )
  })
  back_data <- do.call(rbind, back_data_list[!sapply(back_data_list, is.null)])
  back_data[, "x_num"] <- match(back_data[, "cluster"], cluster_levels)
  tile_colors <- palette_colors(
    x = cluster_levels,
    palette = group_palette,
    palcolor = group_palcolor
  )
  tile_data <- data.frame(
    group1 = cluster_levels,
    x = seq_along(cluster_levels),
    y = 0
  )
  p1 <- ggplot(de_df_marker, aes(x = x_plot, y = y_plot)) +
    geom_col(
      data = back_data,
      aes(x = x_num, y = min),
      fill = manhattan.bg,
      inherit.aes = FALSE
    ) +
    geom_col(
      data = back_data,
      aes(x = x_num, y = max),
      fill = manhattan.bg,
      inherit.aes = FALSE
    )
  p2 <- p1 +
    geom_point(
      aes(x = x_plot, y = y_plot),
      color = cols.highlight,
      size = sizes.highlight + stroke.highlight,
      alpha = alpha.highlight,
      inherit.aes = FALSE,
      data = de_df_marker
    ) +
    geom_point(aes(color = avg_log2FC), size = pt.size, alpha = pt.alpha) +
    scale_color_gradientn(
      name = "log2FC",
      colors = palette_colors(palette = palette, palcolor = palcolor),
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0,
        order = 1
      )
    )
  xlab_use <- if (!is.null(res) && group.by == "custom") NULL else (xlab %||% group.by)
  p3 <- p2 +
    scale_x_continuous(
      breaks = seq_along(cluster_levels), labels = cluster_levels
    ) +
    scale_y_continuous(n.breaks = 6) +
    labs(
      x = xlab_use,
      y = ylab %||% "Average log2FoldChange"
    ) +
    do.call(theme_use, theme_args) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      aspect.ratio = aspect.ratio
    )
  p4 <- p3 +
    geom_tile(
      aes(x = x, y = y, fill = group1),
      color = "black",
      height = 0.5,
      show.legend = FALSE,
      inherit.aes = FALSE,
      data = tile_data
    ) +
    scale_fill_manual(values = tile_colors) +
    geom_text(
      aes(x = x, y = y, label = group1),
      inherit.aes = FALSE,
      data = tile_data,
      size = 3,
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = top_marker,
      aes(x = x_plot, y = y_plot, label = gene),
      inherit.aes = FALSE,
      min.segment.length = 0,
      max.overlaps = 100,
      segment.colour = "grey40",
      color = label.fg,
      bg.color = label.bg,
      bg.r = label.bg.r,
      size = label.size,
      force = 20
    )
  p4 +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#' @title DEtest Ring Plot
#'
#' @description
#' Draw a circular (ring) plot of differential expression results by cell type.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [VolcanoPlot], [DEtestManhattanPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   only.pos = FALSE
#' )
#' DEtestRingPlot(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
DEtestRingPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    group_palette = "Chinese",
    group_palcolor = NULL,
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
    palette = "RdBu",
    palcolor = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    tile_height = 0.3,
    tile_gap = 0.1,
    jitter_width = 0.5,
    ring_segments = TRUE,
    seed = 11) {
  check_r("geomtextpath", verbose = FALSE)
  if (is.null(group.by)) {
    group.by <- "custom"
  }
  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res)
  de_df <- data_res$de_df
  clip_res <- clip_log2fc_symmetric(de_df)
  de_df_marker <- clip_res$df
  fc_lim <- clip_res$fc_lim
  cluster_levels <- levels(de_df_marker[["group1"]])
  top_marker <- get_top_markers_for_label(
    de_df_marker, cluster_levels, nlabel, features_label
  )
  n_grp <- length(cluster_levels)
  ring_r0 <- 2
  max_abs_fc <- max(abs(de_df_marker[, "avg_log2FC"]), na.rm = TRUE)
  ring_k <- if (max_abs_fc > 0) 1.2 / max_abs_fc else 0.2
  set.seed(seed)
  n_m <- nrow(de_df_marker)
  de_df_marker[, "x_num"] <- as.numeric(de_df_marker[, "group1"])
  de_df_marker[, "x_angle"] <- de_df_marker[, "x_num"] + (stats::runif(n_m) - 0.5) * jitter_width
  de_df_marker[, "y_radius"] <- ring_r0 + ring_k * de_df_marker[, "avg_log2FC"]
  band_lo <- ring_r0 - tile_height / 2
  band_hi <- ring_r0 + tile_height / 2
  in_band <- de_df_marker[, "y_radius"] >= band_lo & de_df_marker[, "y_radius"] <= band_hi
  de_df_marker[in_band & de_df_marker[, "y_radius"] < ring_r0, "y_radius"] <- band_lo - tile_gap
  de_df_marker[in_band & de_df_marker[, "y_radius"] >= ring_r0, "y_radius"] <- band_hi + tile_gap
  top_marker <- merge(
    top_marker,
    de_df_marker[, c("gene", "group1", "x_angle", "y_radius")],
    by = c("gene", "group1"),
    all.x = TRUE
  )
  tile_colors <- palette_colors(
    x = cluster_levels,
    palette = group_palette,
    palcolor = group_palcolor
  )
  tile_data <- data.frame(
    x = seq_len(n_grp),
    y = ring_r0,
    group1 = cluster_levels,
    label = cluster_levels
  )
  npt <- 40
  path_margin <- 0.04
  path_df <- do.call(rbind, lapply(seq_len(n_grp), function(i) {
    x_start <- i - 0.5 + path_margin
    x_end <- i + 0.5 - path_margin
    data.frame(
      x = x_start + (0:(npt - 1)) / (npt - 1) * (x_end - x_start),
      y = ring_r0,
      label = cluster_levels[i],
      group = i
    )
  }))
  p_ring <- ggplot(de_df_marker, aes(x = x_angle, y = y_radius))
  if (isTRUE(ring_segments)) {
    p_ring <- p_ring +
      geom_vline(
        xintercept = seq(0.5, n_grp - 0.5, by = 1),
        color = "grey85",
        linewidth = 0.5
      )
  }
  p_ring <- p_ring +
    geom_hline(
      yintercept = ring_r0,
      linetype = 2,
      color = "grey40",
      linewidth = 0.5
    ) +
    geom_point(
      aes(x = x_angle, y = y_radius),
      color = cols.highlight,
      size = sizes.highlight + stroke.highlight,
      alpha = alpha.highlight,
      inherit.aes = FALSE,
      data = de_df_marker
    ) +
    geom_point(
      aes(color = avg_log2FC),
      size = pt.size,
      alpha = pt.alpha
    ) +
    scale_color_gradientn(
      name = "log2FC",
      colors = palette_colors(palette = palette, palcolor = palcolor),
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0,
        order = 1
      )
    ) +
    scale_x_continuous(
      limits = c(0.5, n_grp + 0.5),
      breaks = seq_len(n_grp), labels = NULL
    ) +
    scale_y_continuous(limits = c(0, NA), n.breaks = 5) +
    coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    geom_tile(
      data = tile_data,
      aes(x = x, y = y, fill = group1),
      color = "black",
      height = tile_height,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = tile_colors) +
    geomtextpath::geom_textpath(
      aes(x = x, y = y, label = label, group = group),
      data = path_df,
      inherit.aes = FALSE,
      size = 3,
      color = "black",
      upright = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = top_marker,
      aes(x = x_angle, y = y_radius, label = gene),
      inherit.aes = FALSE,
      min.segment.length = 0,
      max.overlaps = 100,
      segment.colour = "grey40",
      color = label.fg,
      bg.color = label.bg,
      bg.r = label.bg.r,
      size = label.size,
      force = 5
    ) +
    labs(x = NULL, y = NULL, title = NULL, subtitle = NULL, caption = NULL) +
    do.call(theme_use, theme_args) +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      legend.position = "right",
      legend.background = element_blank()
    )
  p_ring
}

#' @title Volcano Plot
#'
#' @description
#' Generate a volcano plot based on differential expression analysis results.
#'
#' @md
#' @inheritParams DEtestPlot
#'
#' @seealso [DEtestPlot], [RunDEtest], [DEtestManhattanPlot], [DEtestRingPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05",
#'   ncol = 2
#' )
#'
#' VolcanoPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   x_metric = "avg_log2FC",
#'   ncol = 2
#' )
VolcanoPlot <- function(
    srt,
    group.by = NULL,
    test.use = "wilcox",
    res = NULL,
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
  data_res <- get_de_data(srt, group.by, test.use, DE_threshold, res)
  de_df <- data_res$de_df

  clip_res <- clip_log2fc_symmetric(de_df)
  de_df <- clip_res$df
  x_upper <- clip_res$fc_lim[2]
  x_lower <- clip_res$fc_lim[1]
  de_df[, "border"] <- (de_df[["avg_log2FC"]] >= x_upper) | (de_df[["avg_log2FC"]] <= x_lower)

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
        colors = palette_colors(palette = palette, palcolor = palcolor),
        values = scales::rescale(unique(c(
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
      do.call(theme_use, theme_args) +
      theme(aspect.ratio = aspect.ratio)
    if (length(levels(de_df[["group1"]])) > 1) {
      p <- p + facet_wrap(~group1)
    }
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
