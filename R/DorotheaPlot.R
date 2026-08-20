#' @title Plot transcription factor activity
#'
#' @description
#' Visualize TF activity stored by [RunDorothea()]. Comparison plots
#' (`"bar"`, `"lollipop"`, `"volcano"`) test two groups and show the mean
#' activity difference `group1 - group2`. `"heatmap"` uses [GroupHeatmap()]
#' to summarize activity across all groups in `group.by`. `"dim"` compares TF
#' activity with TF expression on embeddings, `"stat"` draws per-group
#' activity distributions with [FeatureStatPlot()], and `"targets"` shows a
#' regulon-target volcano for one TF.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param srt A `Seurat` object containing results from [RunDorothea()].
#' @param group.by Metadata column used to define groups.
#' @param group1,group2 Two group labels to compare. Required for `"bar"`,
#' `"lollipop"`, `"volcano"`, and `"targets"`. Positive logFC means higher
#' TF activity (or target expression) in `group1`. Ignored for `"heatmap"`,
#' `"dim"`, and `"stat"`.
#' @param plot_type Plot type. `"bar"` and `"lollipop"` show signed activity
#' differences, `"volcano"` shows TF-level significance versus logFC,
#' `"heatmap"` draws a [GroupHeatmap()] of TF activity, `"dim"` compares
#' activity and expression on embeddings, `"stat"` draws activity
#' distributions, and `"targets"` shows a regulon-target volcano.
#' @param tool_name Name of the `srt@tools` entry created by [RunDorothea()].
#' @param assay_name Assay used for `"heatmap"`, `"dim"`, and `"stat"`. If
#' `NULL`, the assay stored by [RunDorothea()] is used, or `"dorothea"`.
#' @param features TFs to plot. If `NULL`, `"bar"`/`"lollipop"` use the top
#' `top_n` TFs, `"volcano"` shows all tested TFs, heatmaps and `"stat"` use
#' the `top_n` most variable TFs, `"dim"` uses TFs present in both activity
#' and expression assays, and `"targets"` uses the TF with the largest
#' absolute activity difference.
#' @param top_n Number of TFs to show when `features = NULL`. Ignored for
#' `"volcano"`, which always plots every tested TF. Set `NULL` to show all
#' tested TFs in other comparison plots.
#' @param test.use Statistical test used for each TF or target gene in
#' comparison plots.
#' @param p.adjust.method Method passed to [stats::p.adjust].
#' @param color.by P-value column used for comparison-plot color scales.
#' @param rank.by Metric used to select top TFs when `features = NULL` in
#' comparison plots.
#' @param sort.by Metric used to order TFs in `"bar"` and `"lollipop"` plots.
#' @param p_floor Lower bound used before `-log10()` transformation.
#' @param padjustCutoff Adjusted p-value cutoff used to mark significant TFs
#' in `"volcano"` plots and non-supporting targets in `"targets"` plots.
#' @param palette,palcolor Palette passed to `palette_colors()` for comparison
#' plots.
#' @param heatmap_palette,heatmap_palcolor Palette passed to [GroupHeatmap()]
#' when `plot_type = "heatmap"`.
#' @param group_palette,group_palcolor Group annotation palette passed to
#' [GroupHeatmap()] and [CellDimPlot()].
#' @param exp_method Scaling method passed to [GroupHeatmap()] for
#' `plot_type = "heatmap"`.
#' @param heatmap_args Additional arguments passed to [GroupHeatmap()].
#' @param dim_args Additional arguments passed to [FeatureDimPlot()] and
#' [CellDimPlot()] when `plot_type = "dim"`.
#' @param stat_args Additional arguments passed to [FeatureStatPlot()] when
#' `plot_type = "stat"`.
#' @param compare_expression Whether `"dim"` also plots TF gene expression
#' beside TF activity.
#' @param expression_assay,expression_layer Assay and layer used for TF
#' expression in `"dim"` and for target genes in `"targets"`.
#' @param stat_plot_type Distribution plot type passed to [FeatureStatPlot()]
#' when `plot_type = "stat"`.
#' @param nlabel Number of significant TFs labeled in `"volcano"` plots, or
#' significant target genes labeled in `"targets"` plots.
#' @param reduction Reduction used by `"dim"` plots. If `NULL`, the default
#' reduction of `srt` is used.
#' @param bar_width Width of bars in `"bar"` plots.
#' @param point_size Point size in `"lollipop"`, `"volcano"`, and
#' `"targets"` plots.
#' @param title,xlab,ylab,fill.title Axis, plot, and legend titles.
#' @param flip Whether to draw comparison bar/lollipop plots horizontally.
#' @param cols Optional three-color vector used instead of `palette` for
#' diverging comparison scales. Kept for backward compatibility.
#' @param return_data Whether to return a list with the plot and statistics.
#'
#' @return For `"bar"`, `"lollipop"`, `"volcano"`, `"dim"`, `"stat"`, and
#' `"targets"`, a `ggplot` or patchwork object, or a list with `plot` and
#' `data` when `return_data = TRUE`. For `"heatmap"`, the list returned by
#' [GroupHeatmap()].
#'
#' @seealso [RunDorothea], [GroupHeatmap], [FeatureDimPlot], [FeatureStatPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub, verbose = FALSE)
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   layer = "counts",
#'   species = "Mus_musculus",
#'   method = "ulm",
#'   minsize = 5
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   plot_type = "bar",
#'   top_n = 20
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   plot_type = "lollipop",
#'   top_n = 20
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   plot_type = "volcano"
#' )
#'
#' ht <- DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   top_n = 20
#' )
#' ht$plot
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   features = "Sox9",
#'   plot_type = "dim"
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   features = c("Sox9", "Neurod1", "Pdx1"),
#'   plot_type = "stat",
#'   stat_plot_type = "violin"
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   features = "Sox9",
#'   plot_type = "targets"
#' )
DorotheaPlot <- function(
  srt,
  group.by,
  group1 = NULL,
  group2 = NULL,
  plot_type = c(
    "bar",
    "lollipop",
    "heatmap",
    "volcano",
    "dim",
    "stat",
    "targets"
  ),
  tool_name = "Dorothea",
  assay_name = NULL,
  features = NULL,
  top_n = 30,
  test.use = c("wilcox.test", "t.test"),
  p.adjust.method = "BH",
  color.by = c("p_val", "p_val_adj"),
  rank.by = c("abs_logFC", "p_val", "p_val_adj", "logFC"),
  sort.by = c("logFC", "abs_logFC", "p_val", "p_val_adj"),
  p_floor = .Machine$double.xmin,
  padjustCutoff = 0.05,
  palette = "RdBu",
  palcolor = NULL,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  exp_method = c("zscore", "raw"),
  heatmap_args = list(),
  dim_args = list(),
  stat_args = list(),
  compare_expression = TRUE,
  expression_assay = NULL,
  expression_layer = "data",
  stat_plot_type = c("violin", "box", "dot", "bar", "col"),
  nlabel = 10,
  reduction = NULL,
  bar_width = 0.72,
  point_size = 3.2,
  aspect.ratio = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  fill.title = NULL,
  flip = TRUE,
  cols = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  return_data = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(group.by) != 1L || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (
    !tool_name %in% names(srt@tools) || is.null(srt@tools[[tool_name]]$scores)
  ) {
    log_message(
      "No TF activity scores found in {.code srt@tools[[{tool_name}]]$scores}",
      message_type = "error"
    )
  }
  plot_type <- match.arg(plot_type)
  test.use <- match.arg(test.use)
  color.by <- match.arg(color.by)
  rank.by <- match.arg(rank.by)
  sort.by <- match.arg(sort.by)
  exp_method <- match.arg(exp_method)
  stat_plot_type <- match.arg(stat_plot_type)
  if (!is.null(top_n)) {
    if (length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      log_message(
        "{.arg top_n} must be a positive integer or {.code NULL}",
        message_type = "error"
      )
    }
    top_n <- as.integer(top_n)
  }
  if (!is.null(cols) && length(cols) != 3L) {
    log_message(
      "{.arg cols} must contain three colors: low, midpoint, and high",
      message_type = "error"
    )
  }
  if (!is.list(heatmap_args)) {
    log_message(
      "{.arg heatmap_args} must be a list",
      message_type = "error"
    )
  }
  if (!is.list(dim_args)) {
    log_message(
      "{.arg dim_args} must be a list",
      message_type = "error"
    )
  }
  if (!is.list(stat_args)) {
    log_message(
      "{.arg stat_args} must be a list",
      message_type = "error"
    )
  }
  if (length(nlabel) != 1L || is.na(nlabel) || nlabel < 0) {
    log_message(
      "{.arg nlabel} must be a non-negative integer",
      message_type = "error"
    )
  }
  nlabel <- as.integer(nlabel)

  scores <- as.matrix(srt@tools[[tool_name]]$scores)
  network_info <- srt@tools[[tool_name]]$network_info
  if (is.null(network_info)) {
    network_info <- list(
      source = "dorothea",
      label = "DoRothEA",
      signed = TRUE
    )
  }
  network_label <- network_info$label %||% "DoRothEA"
  assay_name <- assay_name %||%
    srt@tools[[tool_name]]$parameters$assay_name %||%
    "dorothea"
  expression_assay <- dorothea_expression_assay(srt, expression_assay)

  if (identical(plot_type, "heatmap")) {
    ht <- dorothea_plot_heatmap(
      srt = srt,
      scores = scores,
      group.by = group.by,
      assay_name = assay_name,
      features = features,
      top_n = top_n,
      exp_method = exp_method,
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      heatmap_args = heatmap_args,
      title = title,
      network_label = network_label,
      verbose = verbose
    )
    if (isTRUE(return_data)) {
      ht$data <- list(
        features = ht$feature_metadata[["features"]] %||%
          rownames(scores)
      )
    }
    return(ht)
  }

  if (identical(plot_type, "dim")) {
    out <- dorothea_plot_dim(
      srt = srt,
      scores = scores,
      group.by = group.by,
      assay_name = assay_name,
      features = features,
      top_n = top_n,
      compare_expression = compare_expression,
      expression_assay = expression_assay,
      expression_layer = expression_layer,
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      legend.position = legend.position,
      legend.direction = legend.direction,
      dim_args = dim_args,
      title = title,
      network_label = network_label,
      verbose = verbose
    )
    if (isTRUE(return_data)) {
      return(out)
    }
    return(out$plot)
  }

  if (identical(plot_type, "stat")) {
    out <- dorothea_plot_stat(
      srt = srt,
      scores = scores,
      group.by = group.by,
      assay_name = assay_name,
      features = features,
      top_n = top_n,
      stat_plot_type = stat_plot_type,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      legend.position = legend.position,
      legend.direction = legend.direction,
      aspect.ratio = aspect.ratio,
      title = title,
      ylab = ylab,
      stat_args = stat_args,
      network_label = network_label,
      verbose = verbose
    )
    if (isTRUE(return_data)) {
      return(out)
    }
    return(out$plot)
  }

  if (identical(plot_type, "targets")) {
    out <- dorothea_plot_targets(
      srt = srt,
      scores = scores,
      group.by = group.by,
      group1 = group1,
      group2 = group2,
      tool_name = tool_name,
      features = features,
      expression_assay = expression_assay,
      expression_layer = expression_layer,
      test.use = test.use,
      p.adjust.method = p.adjust.method,
      p_floor = p_floor,
      padjustCutoff = padjustCutoff,
      palette = palette,
      palcolor = palcolor,
      point_size = point_size,
      nlabel = nlabel,
      title = title,
      xlab = xlab,
      ylab = ylab,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction,
      network_info = network_info,
      verbose = verbose
    )
    if (isTRUE(return_data)) {
      return(out)
    }
    return(out$plot)
  }

  stat_df <- dorothea_compare_activity(
    srt = srt,
    scores = scores,
    group.by = group.by,
    group1 = group1,
    group2 = group2,
    features = features,
    top_n = if (identical(plot_type, "volcano")) NULL else top_n,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    rank.by = rank.by,
    sort.by = sort.by,
    p_floor = p_floor,
    network_label = network_label,
    verbose = verbose
  )
  if (isTRUE(flip) && plot_type %in% c("bar", "lollipop")) {
    stat_df$TF <- factor(stat_df$TF, levels = rev(levels(stat_df$TF)))
  }

  fill_col <- switch(color.by,
    p_val = "signed_neglog10_p_val",
    p_val_adj = "signed_neglog10_p_val_adj"
  )
  title <- title %||% paste(group1, "vs.", group2)
  if (plot_type %in% c("bar", "lollipop")) {
    ylab <- ylab %||% paste0(group1, " - ", group2)
  }
  fill.title_sig <- fill.title %||%
    ifelse(color.by == "p_val", "-log10(p)", "-log10(padj)")
  theme_use <- apply_plot_theme(theme_use, theme_args)
  fill_colors <- dorothea_diverging_colors(
    palette = palette,
    palcolor = palcolor,
    cols = cols
  )

  p <- switch(plot_type,
    bar = dorothea_plot_bar(
      stat_df = stat_df,
      fill_colors = fill_colors,
      bar_width = bar_width,
      title = title,
      xlab = xlab,
      ylab = ylab,
      fill.title = fill.title %||% ylab,
      flip = flip,
      theme_use = theme_use,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    ),
    lollipop = dorothea_plot_lollipop(
      stat_df = stat_df,
      fill_col = fill_col,
      fill_colors = fill_colors,
      point_size = point_size,
      title = title,
      xlab = xlab,
      ylab = ylab,
      fill.title = fill.title_sig,
      flip = flip,
      theme_use = theme_use,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    ),
    volcano = dorothea_plot_volcano(
      stat_df = stat_df,
      padjustCutoff = padjustCutoff,
      palette = palette,
      palcolor = palcolor,
      point_size = point_size,
      nlabel = nlabel,
      title = title,
      xlab = xlab,
      ylab = ylab,
      theme_use = theme_use,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  )

  if (isTRUE(return_data)) {
    return(list(plot = p, data = stat_df))
  }
  p
}

dorothea_compare_activity <- function(
  srt,
  scores,
  group.by,
  group1,
  group2,
  features,
  top_n,
  test.use,
  p.adjust.method,
  rank.by,
  sort.by,
  p_floor,
  verbose,
  network_label = "DoRothEA"
) {
  if (
    length(group1) != 1L ||
      length(group2) != 1L ||
      any(is.na(c(group1, group2)))
  ) {
    log_message(
      "{.arg group1} and {.arg group2} must each be a single group label",
      message_type = "error"
    )
  }
  groups <- as.character(srt@meta.data[[group.by]])
  names(groups) <- rownames(srt@meta.data)
  cells <- intersect(colnames(scores), names(groups))
  cells <- cells[groups[cells] %in% c(group1, group2)]
  if (length(cells) == 0L) {
    log_message(
      "No cells from {.val {group1}} or {.val {group2}} are shared between {.val {network_label}} scores and metadata",
      message_type = "error"
    )
  }
  cells1 <- cells[groups[cells] == group1]
  cells2 <- cells[groups[cells] == group2]
  if (length(cells1) == 0L || length(cells2) == 0L) {
    log_message(
      "Both {.arg group1} and {.arg group2} must contain cells",
      message_type = "error"
    )
  }

  features_is_null <- is.null(features)
  if (features_is_null) {
    features <- rownames(scores)
  } else {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, rownames(scores))
    if (length(missing_features) > 0L) {
      log_message(
        "Dropping TFs not found in {.val {network_label}} scores: {.val {missing_features}}",
        message_type = "warning",
        verbose = verbose
      )
    }
    features <- intersect(features, rownames(scores))
  }
  if (length(features) == 0L) {
    log_message(
      "No TFs are available for plotting",
      message_type = "error"
    )
  }

  log_message(
    "Compare {.val {network_label}} TF activity: {.val {group1}} vs {.val {group2}}",
    verbose = verbose
  )
  mat1 <- scores[features, cells1, drop = FALSE]
  mat2 <- scores[features, cells2, drop = FALSE]
  mean1 <- Matrix::rowMeans(mat1, na.rm = TRUE)
  mean2 <- Matrix::rowMeans(mat2, na.rm = TRUE)
  p_val <- vapply(
    features,
    function(tf) {
      x <- as.numeric(mat1[tf, ])
      y <- as.numeric(mat2[tf, ])
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      if (length(x) < 1L || length(y) < 1L) {
        return(NA_real_)
      }
      tryCatch(
        switch(test.use,
          wilcox.test = stats::wilcox.test(x, y)$p.value,
          t.test = stats::t.test(x, y)$p.value
        ),
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )
  stat_df <- data.frame(
    TF = features,
    group1 = group1,
    group2 = group2,
    mean1 = as.numeric(mean1[features]),
    mean2 = as.numeric(mean2[features]),
    logFC = as.numeric(mean1[features] - mean2[features]),
    p_val = p_val,
    stringsAsFactors = FALSE
  )
  stat_df$p_val[!is.finite(stat_df$p_val) | stat_df$p_val < 0] <- NA_real_
  stat_df$p_val_adj <- stats::p.adjust(stat_df$p_val, method = p.adjust.method)
  stat_df$neglog10_p_val <- -log10(pmax(stat_df$p_val, p_floor, na.rm = TRUE))
  stat_df$neglog10_p_val_adj <- -log10(pmax(
    stat_df$p_val_adj,
    p_floor,
    na.rm = TRUE
  ))
  stat_df$signed_neglog10_p_val <- sign(stat_df$logFC) * stat_df$neglog10_p_val
  stat_df$signed_neglog10_p_val_adj <- sign(stat_df$logFC) *
    stat_df$neglog10_p_val_adj
  stat_df$neglog10_p_val[!is.finite(stat_df$neglog10_p_val)] <- 0
  stat_df$neglog10_p_val_adj[!is.finite(stat_df$neglog10_p_val_adj)] <- 0
  stat_df$signed_neglog10_p_val[!is.finite(stat_df$signed_neglog10_p_val)] <- 0
  stat_df$signed_neglog10_p_val_adj[
    !is.finite(stat_df$signed_neglog10_p_val_adj)
  ] <- 0

  if (isTRUE(features_is_null) && !is.null(top_n)) {
    top_n <- min(as.integer(top_n), nrow(stat_df))
    if (rank.by == "abs_logFC") {
      rank_order <- order(-abs(stat_df$logFC), stat_df$p_val, na.last = TRUE)
    } else if (rank.by == "logFC") {
      rank_order <- order(-stat_df$logFC, stat_df$p_val, na.last = TRUE)
    } else {
      rank_order <- order(
        stat_df[[rank.by]],
        -abs(stat_df$logFC),
        na.last = TRUE
      )
    }
    stat_df <- stat_df[rank_order[seq_len(top_n)], , drop = FALSE]
  }

  if (sort.by == "abs_logFC") {
    plot_order <- order(-abs(stat_df$logFC), -stat_df$logFC, na.last = TRUE)
  } else if (sort.by == "logFC") {
    plot_order <- order(stat_df$logFC, decreasing = TRUE)
  } else {
    plot_order <- order(stat_df[[sort.by]], -abs(stat_df$logFC), na.last = TRUE)
  }
  stat_df <- stat_df[plot_order, , drop = FALSE]
  stat_df$TF <- factor(stat_df$TF, levels = stat_df$TF)
  stat_df
}

dorothea_diverging_colors <- function(palette, palcolor, cols = NULL) {
  if (!is.null(cols)) {
    return(unname(cols))
  }
  palette_colors(
    palette = palette,
    palcolor = palcolor
  )
}

dorothea_fill_scale <- function(fill_colors, name) {
  ggplot2::scale_fill_gradientn(
    colours = fill_colors,
    name = name,
    na.value = "grey80",
    guide = ggplot2::guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      title.hjust = 0
    )
  )
}

dorothea_color_scale <- function(fill_colors, name) {
  ggplot2::scale_color_gradientn(
    colours = fill_colors,
    name = name,
    na.value = "grey80",
    guide = ggplot2::guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      title.hjust = 0
    )
  )
}

dorothea_comparison_theme <- function(
  theme_use,
  flip,
  aspect.ratio,
  legend.position,
  legend.direction
) {
  theme_use +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction,
      panel.grid.major = ggplot2::element_line(colour = "grey80", linetype = 2),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1)
    )
}

dorothea_plot_bar <- function(
  stat_df,
  fill_colors,
  bar_width,
  title,
  xlab,
  ylab,
  fill.title,
  flip,
  theme_use,
  aspect.ratio,
  legend.position,
  legend.direction
) {
  p <- ggplot2::ggplot(
    stat_df,
    ggplot2::aes(x = .data[["TF"]], y = .data[["logFC"]], fill = .data[["logFC"]])
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.25) +
    dorothea_fill_scale(fill_colors, fill.title) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    dorothea_comparison_theme(
      theme_use = theme_use,
      flip = flip,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }
  p
}

dorothea_plot_lollipop <- function(
  stat_df,
  fill_col,
  fill_colors,
  point_size,
  title,
  xlab,
  ylab,
  fill.title,
  flip,
  theme_use,
  aspect.ratio,
  legend.position,
  legend.direction
) {
  p <- ggplot2::ggplot(
    stat_df,
    ggplot2::aes(x = .data[["TF"]], y = .data[["logFC"]])
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data[["TF"]],
        y = 0,
        yend = .data[["logFC"]]
      ),
      color = "grey35",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = .data[[fill_col]],
        size = abs(.data[["logFC"]])
      ),
      shape = 21,
      color = "black",
      stroke = 0.35
    ) +
    ggplot2::scale_size_continuous(
      name = "|logFC|",
      range = c(2, max(point_size, 2) + 2)
    ) +
    dorothea_fill_scale(fill_colors, fill.title) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, override.aes = list(fill = "grey30")),
      fill = ggplot2::guide_colorbar(
        order = 2,
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0
      )
    ) +
    dorothea_comparison_theme(
      theme_use = theme_use,
      flip = flip,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }
  p
}

dorothea_volcano_y <- function(stat_df) {
  y <- as.numeric(stat_df$neglog10_p_val_adj)
  y[!is.finite(y)] <- 0
  finite <- y[is.finite(stat_df$p_val_adj) & stat_df$p_val_adj > 0]
  if (length(finite) == 0L) {
    cap <- max(y, na.rm = TRUE)
    if (!is.finite(cap) || cap <= 0) {
      cap <- 10
    }
    return(pmin(y, cap))
  }
  cap <- 10 * ceiling(max(finite) / 10)
  if (!is.finite(cap) || cap <= 0) {
    cap <- max(finite)
  }
  pmin(y, cap)
}

dorothea_volcano_direction <- function(stat_df, cutoff) {
  direction <- rep("NS", nrow(stat_df))
  sig <- is.finite(stat_df$p_val_adj) & stat_df$p_val_adj <= cutoff
  direction[sig & is.finite(stat_df$logFC) & stat_df$logFC > 0] <- "Up"
  direction[sig & is.finite(stat_df$logFC) & stat_df$logFC < 0] <- "Down"
  factor(direction, levels = c("Down", "NS", "Up"))
}

dorothea_volcano_colors <- function(palette, palcolor) {
  div <- unname(dorothea_diverging_colors(palette, palcolor))
  if (length(div) < 2L) {
    div <- c("#2166AC", "#B2182B")
  }
  c(Down = div[[1L]], NS = "grey75", Up = div[[length(div)]])
}

dorothea_volcano_label_df <- function(df, nlabel, y_col, sig) {
  n_sig <- sum(sig, na.rm = TRUE)
  nlabel <- min(as.integer(nlabel), n_sig, nrow(df))
  if (!is.finite(nlabel) || nlabel <= 0L) {
    return(df[integer(), , drop = FALSE])
  }
  rank_score <- abs(df$logFC) * df[[y_col]]
  rank_score[!is.finite(rank_score) | !sig] <- -Inf
  df[order(-rank_score), , drop = FALSE][seq_len(nlabel), , drop = FALSE]
}

dorothea_volcano_repel <- function(label_df, label_col) {
  ggrepel::geom_text_repel(
    data = label_df,
    ggplot2::aes(
      x = .data[["logFC"]],
      y = .data[["y_plot"]],
      label = .data[[label_col]]
    ),
    inherit.aes = FALSE,
    min.segment.length = 0.2,
    max.overlaps = Inf,
    segment.colour = "grey40",
    size = 3.2,
    seed = 42,
    show.legend = FALSE
  )
}

dorothea_plot_volcano <- function(
  stat_df,
  padjustCutoff,
  palette,
  palcolor,
  point_size,
  nlabel,
  title,
  xlab,
  ylab,
  theme_use,
  aspect.ratio,
  legend.position,
  legend.direction
) {
  cutoff <- padjustCutoff %||% 0.05
  stat_df$y_plot <- dorothea_volcano_y(stat_df)
  stat_df$direction <- dorothea_volcano_direction(stat_df, cutoff)
  cols <- dorothea_volcano_colors(palette, palcolor)
  pt_size <- if (nrow(stat_df) >= 50L) min(point_size, 1.8) else point_size
  label_df <- dorothea_volcano_label_df(
    df = stat_df,
    nlabel = nlabel,
    y_col = "y_plot",
    sig = stat_df$direction != "NS"
  )
  p <- ggplot2::ggplot(
    stat_df,
    ggplot2::aes(
      x = .data[["logFC"]],
      y = .data[["y_plot"]],
      color = .data[["direction"]]
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_hline(
      yintercept = -log10(cutoff),
      color = "grey75",
      linewidth = 0.35,
      linetype = 2
    ) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = cols,
      breaks = names(cols),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = xlab %||% paste0(unique(stat_df$group1), " - ", unique(stat_df$group2)),
      y = ylab %||% "-log10(adjusted p-value)",
      color = NULL
    ) +
    theme_use +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (nrow(label_df) > 0L) {
    p <- p + dorothea_volcano_repel(label_df, "TF")
  }
  p
}

dorothea_plot_heatmap <- function(
  srt,
  scores,
  group.by,
  assay_name,
  features,
  top_n,
  exp_method,
  heatmap_palette,
  heatmap_palcolor,
  group_palette,
  group_palcolor,
  heatmap_args,
  title,
  verbose,
  network_label = "DoRothEA"
) {
  features <- dorothea_resolve_features(
    scores = scores,
    features = features,
    top_n = top_n,
    verbose = verbose
  )

  if (!assay_name %in% SeuratObject::Assays(srt)) {
    srt <- dorothea_attach_assay(
      srt = srt,
      scores = scores,
      assay_name = assay_name
    )
  }

  log_message(
    "Draw {.val {network_label}} TF activity heatmap for {.val {length(features)}} TFs",
    verbose = verbose
  )
  args <- list(
    srt = srt,
    features = features,
    group.by = group.by,
    assay = assay_name,
    layer = "data",
    aggregate_fun = base::mean,
    exp_method = exp_method,
    exp_legend_title = if (identical(exp_method, "zscore")) {
      "Activity z"
    } else {
      "Mean activity"
    },
    lib_normalize = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_title = title %||% paste(network_label, "TF activity"),
    heatmap_palette = heatmap_palette,
    heatmap_palcolor = heatmap_palcolor,
    group_palette = group_palette,
    group_palcolor = group_palcolor
  )
  args[names(heatmap_args)] <- heatmap_args
  do.call(GroupHeatmap, args)
}

dorothea_resolve_features <- function(scores, features, top_n, verbose) {
  if (is.null(features)) {
    tf_var <- apply(scores, 1, stats::var, na.rm = TRUE)
    tf_var <- tf_var[is.finite(tf_var)]
    if (length(tf_var) == 0L) {
      log_message(
        "No TFs with finite activity variance are available for plotting",
        message_type = "error"
      )
    }
    n_keep <- if (is.null(top_n)) length(tf_var) else min(top_n, length(tf_var))
    features <- names(sort(tf_var, decreasing = TRUE))[seq_len(n_keep)]
  } else {
    features <- unique(as.character(features))
    missing_features <- setdiff(features, rownames(scores))
    if (length(missing_features) > 0L) {
      log_message(
        "Dropping TFs not found in DoRothEA scores: {.val {missing_features}}",
        message_type = "warning",
        verbose = verbose
      )
    }
    features <- intersect(features, rownames(scores))
  }
  if (length(features) == 0L) {
    log_message(
      "No TFs are available for plotting",
      message_type = "error"
    )
  }
  features
}

dorothea_expression_assay <- function(srt, expression_assay = NULL) {
  assays <- SeuratObject::Assays(srt)
  if (!is.null(expression_assay)) {
    if (!expression_assay %in% assays) {
      log_message(
        "{.arg expression_assay} {.val {expression_assay}} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    return(expression_assay)
  }
  if ("RNA" %in% assays) {
    return("RNA")
  }
  SeuratObject::DefaultAssay(srt)
}

dorothea_match_feature <- function(feature, rownames_vec) {
  if (feature %in% rownames_vec) {
    return(feature)
  }
  hit <- rownames_vec[tolower(rownames_vec) == tolower(feature)]
  if (length(hit) == 1L) {
    return(hit[[1L]])
  }
  NA_character_
}

dorothea_ensure_assay <- function(srt, scores, assay_name) {
  if (!assay_name %in% SeuratObject::Assays(srt)) {
    srt <- dorothea_attach_assay(
      srt = srt,
      scores = scores,
      assay_name = assay_name
    )
  }
  srt
}

dorothea_plot_dim <- function(
  srt,
  scores,
  group.by,
  assay_name,
  features,
  top_n,
  compare_expression,
  expression_assay,
  expression_layer,
  reduction,
  palette,
  palcolor,
  group_palette,
  group_palcolor,
  theme_use,
  theme_args,
  legend.position,
  legend.direction,
  dim_args,
  title,
  verbose,
  network_label = "DoRothEA"
) {
  expr_names <- tryCatch(
    rownames(srt[[expression_assay]]),
    error = function(e) character()
  )
  if (is.null(features) && isTRUE(compare_expression) && length(expr_names) > 0L) {
    shared <- rownames(scores)[
      tolower(rownames(scores)) %in% tolower(expr_names)
    ]
    if (length(shared) > 0L) {
      features <- dorothea_resolve_features(
        scores = scores[shared, , drop = FALSE],
        features = NULL,
        top_n = top_n,
        verbose = verbose
      )
    }
  }
  features <- dorothea_resolve_features(
    scores = scores,
    features = features,
    top_n = if (is.null(features)) top_n else NULL,
    verbose = verbose
  )
  srt <- dorothea_ensure_assay(srt, scores, assay_name)

  log_message(
    "Draw {.val {network_label}} embedding plots for {.val {length(features)}} TFs",
    verbose = verbose
  )

  cell_formals <- names(formals(CellDimPlot))
  feat_formals <- names(formals(FeatureDimPlot))
  cell_args <- list(
    srt = srt,
    group.by = group.by,
    reduction = reduction,
    palette = group_palette,
    palcolor = group_palcolor,
    label = TRUE,
    title = "Cell types",
    theme_use = theme_use,
    theme_args = theme_args,
    legend.position = legend.position,
    legend.direction = legend.direction
  )
  extra_cell <- dim_args[intersect(names(dim_args), cell_formals)]
  extra_cell <- extra_cell[setdiff(names(extra_cell), c("srt", "group.by"))]
  cell_args[names(extra_cell)] <- extra_cell
  p_cell <- do.call(CellDimPlot, cell_args)

  tf_plots <- lapply(features, function(tf) {
    act_args <- list(
      srt = srt,
      features = tf,
      assay = assay_name,
      layer = "data",
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      bg_cutoff = -Inf,
      title = paste(tf, "activity"),
      legend.title = "Activity",
      theme_use = theme_use,
      theme_args = theme_args,
      legend.position = legend.position,
      legend.direction = legend.direction,
      combine = TRUE
    )
    extra_feat <- dim_args[intersect(names(dim_args), feat_formals)]
    extra_feat <- extra_feat[setdiff(
      names(extra_feat),
      c("srt", "features", "assay", "layer", "title", "legend.title")
    )]
    act_args[names(extra_feat)] <- extra_feat
    p_act <- do.call(FeatureDimPlot, act_args)

    expr_feature <- dorothea_match_feature(tf, expr_names)
    if (!isTRUE(compare_expression) || is.na(expr_feature)) {
      if (isTRUE(compare_expression) && is.na(expr_feature)) {
        log_message(
          "No matching expression feature for TF {.val {tf}} in assay {.val {expression_assay}}",
          message_type = "warning",
          verbose = verbose
        )
      }
      return(p_act)
    }
    exp_args <- act_args
    exp_args$features <- expr_feature
    exp_args$assay <- expression_assay
    exp_args$layer <- expression_layer
    exp_args$palette <- dim_args$expression_palette %||% "Spectral"
    exp_args$palcolor <- dim_args$expression_palcolor
    exp_args$bg_cutoff <- dim_args$expression_bg_cutoff %||% 0
    exp_args$title <- paste(expr_feature, "expression")
    exp_args$legend.title <- "Expression"
    p_exp <- do.call(FeatureDimPlot, exp_args)
    patchwork::wrap_plots(p_act, p_exp, ncol = 2)
  })

  if (length(features) == 1L && isTRUE(compare_expression)) {
    expr_feature <- dorothea_match_feature(features[[1L]], expr_names)
    if (!is.na(expr_feature)) {
      p <- patchwork::wrap_plots(p_cell, tf_plots[[1L]], ncol = 2, widths = c(1, 2))
      if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
      }
      return(list(plot = p, data = features))
    }
  }
  p <- patchwork::wrap_plots(
    c(list(p_cell), tf_plots),
    ncol = 1
  )
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(title = title)
  }
  list(plot = p, data = features)
}

dorothea_plot_stat <- function(
  srt,
  scores,
  group.by,
  assay_name,
  features,
  top_n,
  stat_plot_type,
  group_palette,
  group_palcolor,
  theme_use,
  theme_args,
  legend.position,
  legend.direction,
  aspect.ratio,
  title,
  ylab,
  stat_args,
  verbose,
  network_label = "DoRothEA"
) {
  features <- dorothea_resolve_features(
    scores = scores,
    features = features,
    top_n = top_n,
    verbose = verbose
  )
  srt <- dorothea_ensure_assay(srt, scores, assay_name)
  log_message(
    "Draw {.val {network_label}} activity distributions for {.val {length(features)}} TFs",
    verbose = verbose
  )
  args <- list(
    srt = srt,
    stat.by = features,
    group.by = group.by,
    assay = assay_name,
    layer = "data",
    plot_type = stat_plot_type,
    palette = group_palette,
    palcolor = group_palcolor,
    ylab = ylab %||% "TF activity",
    title = title,
    theme_use = theme_use,
    theme_args = theme_args,
    legend.position = legend.position,
    legend.direction = legend.direction,
    aspect.ratio = aspect.ratio
  )
  extra <- stat_args[intersect(names(stat_args), names(formals(FeatureStatPlot)))]
  extra <- extra[setdiff(names(extra), c("srt", "stat.by", "assay", "layer"))]
  args[names(extra)] <- extra
  list(plot = do.call(FeatureStatPlot, args), data = features)
}

dorothea_plot_targets <- function(
  srt,
  scores,
  group.by,
  group1,
  group2,
  tool_name,
  features,
  expression_assay,
  expression_layer,
  test.use,
  p.adjust.method,
  p_floor,
  padjustCutoff,
  palette,
  palcolor,
  point_size,
  nlabel,
  title,
  xlab,
  ylab,
  theme_use,
  theme_args,
  aspect.ratio,
  legend.position,
  legend.direction,
  network_info,
  verbose
) {
  regulons <- srt@tools[[tool_name]]$regulons
  network_label <- network_info$label %||% "DoRothEA"
  signed_network <- !identical(network_info$signed, FALSE)
  if (is.null(regulons) || !"tf" %in% colnames(regulons)) {
    log_message(
      "No regulon network found in {.code srt@tools[[{tool_name}]]$regulons}. Re-run {.fn RunDorothea}.",
      message_type = "error"
    )
  }
  if (is.null(features)) {
    stat_df <- dorothea_compare_activity(
      srt = srt,
      scores = scores,
      group.by = group.by,
      group1 = group1,
      group2 = group2,
      features = NULL,
      top_n = 1L,
      test.use = test.use,
      p.adjust.method = p.adjust.method,
      rank.by = "abs_logFC",
      sort.by = "abs_logFC",
      p_floor = p_floor,
      network_label = network_label,
      verbose = FALSE
    )
    features <- as.character(stat_df$TF[[1L]])
  }
  if (length(features) != 1L) {
    log_message(
      "{.arg plot_type = 'targets'} accepts a single TF in {.arg features}",
      message_type = "error"
    )
  }
  tf <- features[[1L]]
  tf_use <- dorothea_match_feature(tf, unique(as.character(regulons$tf)))
  if (is.na(tf_use)) {
    log_message(
      "TF {.val {tf}} is not present in the stored {.val {network_label}} regulons",
      message_type = "error"
    )
  }
  net <- regulons[as.character(regulons$tf) == tf_use, , drop = FALSE]
  expr <- GetAssayData5(
    srt,
    assay = expression_assay,
    layer = expression_layer
  )
  if (nrow(expr) == 0L && !identical(expression_layer, "counts")) {
    expr <- GetAssayData5(
      srt,
      assay = expression_assay,
      layer = "counts"
    )
  }
  expr <- as.matrix(expr)
  target_map <- vapply(
    as.character(net$target),
    function(gene) dorothea_match_feature(gene, rownames(expr)),
    character(1)
  )
  keep <- !is.na(target_map)
  if (!any(keep)) {
    log_message(
      "None of the {.val {tf_use}} targets are present in assay {.val {expression_assay}}",
      message_type = "error"
    )
  }
  if (any(!keep)) {
    log_message(
      "Dropping {.val {sum(!keep)}} {.val {tf_use}} targets missing from assay {.val {expression_assay}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  net <- net[keep, , drop = FALSE]
  net$target_expr <- unname(target_map[keep])
  target_mat <- expr[unique(net$target_expr), , drop = FALSE]
  log_message(
    "Draw {.val {network_label}} regulon-target volcano for {.val {tf_use}} ({.val {nrow(net)}} targets)",
    verbose = verbose
  )
  target_stat <- dorothea_compare_activity(
    srt = srt,
    scores = target_mat,
    group.by = group.by,
    group1 = group1,
    group2 = group2,
    features = rownames(target_mat),
    top_n = NULL,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    rank.by = "abs_logFC",
    sort.by = "abs_logFC",
    p_floor = p_floor,
    network_label = network_label,
    verbose = FALSE
  )
  plot_df <- merge(
    net,
    target_stat,
    by.x = "target_expr",
    by.y = "TF",
    all.x = TRUE,
    sort = FALSE
  )
  plot_df$mor <- as.numeric(plot_df$mor)
  if (signed_network) {
    plot_df$support <- ifelse(
      is.finite(plot_df$logFC) & is.finite(plot_df$mor) & plot_df$mor != 0,
      ifelse(sign(plot_df$mor) * sign(plot_df$logFC) > 0, "Support", "Oppose"),
      "NS"
    )
    support_levels <- c("Oppose", "NS", "Support")
    size_title <- "|MoR|"
  } else {
    plot_df$support <- ifelse(
      is.finite(plot_df$logFC),
      ifelse(plot_df$logFC > 0, "Higher group1", "Higher group2"),
      "NS"
    )
    support_levels <- c("Higher group2", "NS", "Higher group1")
    size_title <- "Importance"
  }
  cutoff <- padjustCutoff %||% 0.05
  plot_df$support[
    !is.finite(plot_df$p_val_adj) | plot_df$p_val_adj > cutoff
  ] <- "NS"
  cols <- palette_colors(
    support_levels,
    palette = palette,
    palcolor = palcolor
  )
  if (is.null(names(cols)) || !all(support_levels %in% names(cols))) {
    cols <- unname(cols)[seq_len(min(3L, length(cols)))]
    if (length(cols) < 3L) {
      cols <- rep(cols, length.out = 3L)
    }
    names(cols) <- support_levels
  } else {
    cols <- cols[support_levels]
  }
  theme_use <- apply_plot_theme(theme_use, theme_args)
  plot_df$y_plot <- dorothea_volcano_y(plot_df)
  label_df <- dorothea_volcano_label_df(
    df = plot_df,
    nlabel = nlabel,
    y_col = "y_plot",
    sig = plot_df$support != "NS"
  )
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["logFC"]],
      y = .data[["y_plot"]],
      color = .data[["support"]],
      size = abs(.data[["mor"]])
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_hline(
      yintercept = -log10(cutoff),
      color = "grey75",
      linewidth = 0.35,
      linetype = 2
    ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE, name = NULL) +
    ggplot2::scale_size_continuous(name = size_title, range = c(1.6, max(point_size, 1.6) + 2)) +
    ggplot2::labs(
      title = title %||% paste0(network_label, ": ", tf_use, " targets (", group1, " vs. ", group2, ")"),
      x = xlab %||% paste0(group1, " - ", group2),
      y = ylab %||% "-log10(adjusted p-value)"
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, override.aes = list(color = "grey30")),
      color = ggplot2::guide_legend(order = 2, override.aes = list(size = 3))
    ) +
    theme_use +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (nrow(label_df) > 0L) {
    p <- p + dorothea_volcano_repel(label_df, "target_expr")
  }
  list(plot = p, data = plot_df)
}
