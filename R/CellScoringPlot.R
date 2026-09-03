#' @title Plot CellScoring results from a Seurat object
#'
#' @inheritParams CellDimPlot
#'
#' @param srt A `Seurat` object.
#' @param method Cell-scoring method to retrieve and label. Accepts the same
#'   values as the `method` argument of [CellScoring()].
#' @param features Score names to plot. For metadata input these are metadata
#'   columns; for `scores` input these are matrix row or column names. A named
#'   character vector uses its names as display labels. When both `features`
#'   and `scores` are `NULL`, the most recent available result for `method`
#'   recorded by [CellScoring()] is used. Older AUCell objects fall back to
#'   numeric metadata.
#' @param scores Optional score matrix or data frame, with features by cells or
#'   cells by features. When `NULL`, scores are read from `srt` metadata or the
#'   matching [CellScoring()] record.
#' @param group.by Metadata column used for cell groups. `NULL` uses active
#'   Seurat identities.
#' @param thresholds Optional score thresholds. Supply one value per feature,
#'   a named vector keyed by display or source feature name, or `"auto"` to
#'   infer thresholds for AUCell scores. When `NULL`, `add_bar = TRUE`, and
#'   `method = "AUCell"`, automatic thresholds are used.
#' @param add_box Add the score boxplot.
#' @param add_bar Add the threshold-proportion bar chart. Defaults to `TRUE`.
#'   Without explicit `thresholds`, thresholds are inferred automatically only
#'   for AUCell. Other methods require explicit thresholds or `add_bar = FALSE`.
#' @param bar.text.size Text size for counts and percentages inside threshold
#'   bars.
#' @param bar.text.color Text color for counts and percentages inside threshold
#'   bars.
#' @param box.text.size Text size for the mean-score label in boxplots.
#' @param box.text.color Text color for the mean-score label in boxplots.
#' @param ncol Number of physical columns. `NULL` uses one column per group or
#'   feature slot when `nrow` is also `NULL`, or the minimum required by `nrow`.
#' @param nrow Number of physical rows. Rows occur in UMAP/statistic pairs and
#'   must therefore be even. `NULL` chooses the minimum required value.
#' @param row.heights Relative heights of the UMAP and statistic rows within
#'   each row pair. Use `c(0.47, 0.53)` for the compact benchmark layout.
#' @param point.size Point size in rasterized UMAPs. `NULL` lets
#'   [CellDimPlot()] choose a readable size from the number of plotted cells.
#' @param point.fraction Fraction of cells shown as jittered points in each
#'   boxplot. Boxplot statistics always use all cells.
#' @param point.alpha Alpha for jittered points.
#' @param boxplot.y.range If `TRUE`, use the reference boxplot y-axis range
#'   of 0--0.5 when it contains all scores; otherwise use the score range.
#' @param group.palette,group.palcolor Group color palette arguments passed to
#'   [thisplot::palette_colors()].
#' @param score.palette,score.palcolor Continuous score palette name and optional
#'   custom colors, resolved with [thisplot::palette_colors()]. The default
#'   `"Spectral"` matches other continuous feature plots in scop.
#' @param theme_use Theme used by [CellDimPlot()] for the main group UMAP.
#' @param auc_theme_use Theme used by [CellDimPlot()] for AUC UMAPs.
#' @param highlight.label Whether to label highlighted groups on score UMAPs.
#'   The highlighted group is inferred from each feature's source name, falling
#'   back to its display label. Defaults to `FALSE`; outlines are still drawn.
#' @param highlight.mark.type Highlight outline geometry. `"ellipse"` draws
#'   a confidence ellipse; `"hull"` follows the selected cells more
#'   closely with [ggforce::geom_mark_hull()].
#' @param highlight.mark.color Color of highlighted group outlines. `NULL`
#'   inherits the matching group colors, consistent with [CellDimPlot()].
#' @param highlight.label.color Highlight label color; defaults to grey.
#' @param highlight.mean.color Mean-score line color. `NULL` inherits
#'   `highlight.label.color`. The mean-score label text is styled with
#'   `box.text.size` and `box.text.color`.
#' @param highlight.color Compatibility shortcut that sets all highlight
#'   colors when their newer arguments are not supplied.
#' @param ellipse.level Confidence level for highlighted group ellipses.
#' @param show.score.legend Which score UMAPs show a continuous legend: `TRUE`
#'   shows the first, `FALSE` shows none, or supply feature names.
#' @param legend.position Position of collected score and threshold legends.
#'   Legends are horizontal at `"top"` or `"bottom"`, and vertical at
#'   `"left"` or `"right"`.
#' @param group.title Optional group legend title.
#' @param threshold.colors Colors for cells above and below each threshold.
#' @param ... Additional arguments forwarded to [CellDimPlot()].
#' @param combine Return a patchwork when `TRUE`; otherwise return named plot
#'   components before layout spacers are inserted.
#' @param seed Seed used for automatic threshold inference and jitter-point
#'   sampling.
#' @param verbose Whether to report plotting progress.
#'
#' @return A `patchwork` object, or a named list of `ggplot`/`patchwork`
#'   components when `combine = FALSE`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub, verbose = FALSE)
#' reduction <- DefaultReduction(pancreas_sub)
#' SeuratObject::Key(pancreas_sub[[reduction]]) <- "UMAP_"
#' genesets <- list(
#'   Alpha = c("Gcg", "Ttr", "Mafb", "Arx", "Slc7a2"),
#'   Beta = c("Ins1", "Ins2", "Iapp", "Pdx1"),
#'   Ductal = c("Krt19", "Krt8", "Krt18", "Sox9", "Krt7")
#' )
#' pancreas_sub <- CellScoring(
#'   pancreas_sub,
#'   features = genesets,
#'   method = "AUCell",
#'   classification = FALSE
#' )
#' CellScoringPlot(
#'   pancreas_sub,
#'   method = "AUCell",
#'   group.by = "SubCellType"
#' )
CellScoringPlot <- function(
  srt,
  method = "AUCell",
  features = NULL,
  scores = NULL,
  group.by = NULL,
  reduction = NULL,
  dims = c(1, 2),
  cells = NULL,
  thresholds = NULL,
  add_box = TRUE,
  add_bar = TRUE,
  bar.text.size = 2,
  bar.text.color = "white",
  box.text.size = 2,
  box.text.color = "grey50",
  ncol = NULL,
  nrow = NULL,
  row.heights = c(1, 1),
  point.size = NULL,
  raster = NULL,
  raster.dpi = 512,
  point.fraction = 0.1,
  point.alpha = 0.2,
  boxplot.y.range = TRUE,
  group.palette = "Chinese",
  group.palcolor = NULL,
  score.palette = "Spectral",
  score.palcolor = NULL,
  theme_use = "theme_scop",
  auc_theme_use = "theme_scop",
  highlight.label = FALSE,
  highlight.mark.type = c("ellipse", "hull"),
  highlight.mark.color = NULL,
  highlight.label.color = "grey40",
  highlight.mean.color = NULL,
  highlight.color = NULL,
  ellipse.level = 0.98,
  show.score.legend = TRUE,
  legend.position = "bottom",
  group.title = NULL,
  threshold.colors = c(
    "Above AUC threshold" = "black",
    "Below AUC threshold" = "grey50"
  ),
  ...,
  combine = TRUE,
  seed = 42,
  verbose = TRUE
) {
  if (length(method) != 1L || is.na(method)) {
    log_message(
      "{.arg method} must name one CellScoring method",
      message_type = "error"
    )
  }
  method <- normalize_gene_set_scoring_method(method, arg_name = "method")
  score_label <- if (identical(method, "AUCell")) "AUC" else method
  highlight.mark.type <- match.arg(highlight.mark.type)
  legend.position <- match.arg(
    legend.position,
    choices = c("bottom", "top", "left", "right")
  )
  mark_color_missing <- missing(highlight.mark.color)
  label_color_missing <- missing(highlight.label.color)
  mean_color_missing <- missing(highlight.mean.color)
  if (!is.null(highlight.color)) {
    if (mark_color_missing) {
      highlight.mark.color <- highlight.color
    }
    if (label_color_missing) {
      highlight.label.color <- highlight.color
    }
    if (mean_color_missing) highlight.mean.color <- highlight.color
  }
  if (is.null(highlight.mean.color)) {
    highlight.mean.color <- highlight.label.color
  }
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(dims) != 2L || anyNA(dims) || any(dims < 1)) {
    log_message(
      "{.arg dims} must contain two positive dimension indices",
      message_type = "error"
    )
  }
  scoreplot_validate_layout(ncol, "ncol")
  scoreplot_validate_layout(nrow, "nrow")
  if (!is.null(nrow) && nrow %% 2L != 0L) {
    log_message(
      "{.arg nrow} must be even because rows occur in UMAP/statistic pairs",
      message_type = "error"
    )
  }
  if (
    !is.numeric(row.heights) ||
      length(row.heights) != 2L ||
      anyNA(row.heights) ||
      any(!is.finite(row.heights)) ||
      any(row.heights <= 0)
  ) {
    log_message(
      "{.arg row.heights} must contain two positive finite values",
      message_type = "error"
    )
  }
  if (
    !is.numeric(point.fraction) ||
      length(point.fraction) != 1L ||
      is.na(point.fraction) ||
      point.fraction < 0 ||
      point.fraction > 1
  ) {
    log_message(
      "{.arg point.fraction} must be between 0 and 1",
      message_type = "error"
    )
  }
  if (
    !is.logical(add_box) ||
      length(add_box) != 1L ||
      is.na(add_box) ||
      !is.logical(add_bar) ||
      length(add_bar) != 1L ||
      is.na(add_bar)
  ) {
    log_message(
      "{.arg add_box} and {.arg add_bar} must be single logical values",
      message_type = "error"
    )
  }
  if (
    !is.numeric(bar.text.size) ||
      length(bar.text.size) != 1L ||
      is.na(bar.text.size) ||
      !is.finite(bar.text.size) ||
      bar.text.size <= 0
  ) {
    log_message(
      "{.arg bar.text.size} must be one positive finite number",
      message_type = "error"
    )
  }
  if (
    !is.character(bar.text.color) ||
      length(bar.text.color) != 1L ||
      is.na(bar.text.color) ||
      !nzchar(bar.text.color)
  ) {
    log_message(
      "{.arg bar.text.color} must be one non-empty color string",
      message_type = "error"
    )
  }
  if (
    !is.numeric(box.text.size) ||
      length(box.text.size) != 1L ||
      is.na(box.text.size) ||
      !is.finite(box.text.size) ||
      box.text.size <= 0
  ) {
    log_message(
      "{.arg box.text.size} must be one positive finite number",
      message_type = "error"
    )
  }
  if (
    !is.character(box.text.color) ||
      length(box.text.color) != 1L ||
      is.na(box.text.color) ||
      !nzchar(box.text.color)
  ) {
    log_message(
      "{.arg box.text.color} must be one non-empty color string",
      message_type = "error"
    )
  }
  if (
    !is.logical(highlight.label) ||
      length(highlight.label) != 1L ||
      is.na(highlight.label)
  ) {
    log_message(
      "{.arg highlight.label} must be a single logical value",
      message_type = "error"
    )
  }
  if (
    !is.null(raster) &&
      (!is.logical(raster) || length(raster) != 1L || is.na(raster))
  ) {
    log_message(
      "{.arg raster} must be TRUE, FALSE, or NULL",
      message_type = "error"
    )
  }

  all_cells <- colnames(srt)
  if (is.null(cells)) {
    cells <- all_cells
  } else {
    cells <- unique(as.character(cells))
    missing_cells <- setdiff(cells, all_cells)
    if (length(missing_cells)) {
      log_message(
        "{.arg cells} contains cells absent from {.arg srt}: {.val {missing_cells}}",
        message_type = "error"
      )
    }
  }
  if (!length(cells)) {
    log_message("No cells are available for plotting", message_type = "error")
  }

  if (is.null(scores) && is.null(features)) {
    stored_input <- scoreplot_input(srt, method = method)
    if (!is.null(stored_input)) {
      scores <- stored_input$scores
      features <- stored_input$features
      if (isTRUE(verbose)) {
        log_message(
          "Using the latest {method} result recorded by {.fn CellScoring}",
          message_type = "info"
        )
      }
    } else if (!identical(method, "AUCell")) {
      log_message(
        "No available {.fn CellScoring} result was recorded for method {.val {method}}; supply {.arg features} or {.arg scores}",
        message_type = "error"
      )
    }
  }
  score_input <- scoreplot_scores(srt, scores, features, cells)
  scores <- score_input$scores
  feature_labels <- score_input$labels
  feature_sources <- score_input$sources
  thresholds <- scoreplot_thresholds(
    thresholds,
    feature_labels,
    feature_sources,
    scores = scores,
    add_bar = add_bar,
    method = method,
    seed = seed,
    verbose = verbose
  )
  groups <- scoreplot_groups(srt, group.by, cells)

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  }
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "Reduction {.val {reduction}} is absent from {.arg srt}",
      message_type = "error"
    )
  }
  reduction_embedding <- SeuratObject::Embeddings(srt[[reduction]])
  if (max(dims) > ncol(reduction_embedding)) {
    log_message(
      "{.arg dims} exceeds the dimensions available in reduction {.val {reduction}}",
      message_type = "error"
    )
  }
  group_palcolor <- scoreplot_named_colors(group.palcolor, levels(groups))
  group_colors <- thisplot::palette_colors(
    levels(groups),
    palette = group.palette,
    palcolor = group_palcolor,
    type = "discrete",
    matched = TRUE
  )
  score_colors <- thisplot::palette_colors(
    n = 1000,
    palette = score.palette,
    palcolor = score.palcolor,
    type = "continuous",
    matched = FALSE
  )
  highlights <- scoreplot_highlights(
    feature_labels,
    feature_sources,
    levels(groups)
  )
  feature_colors <- vapply(
    highlights,
    function(highlight) {
      if (length(highlight)) {
        unname(group_colors[[highlight]])
      } else {
        "black"
      }
    },
    character(1)
  )
  legend_features <- scoreplot_legend_features(
    show.score.legend,
    feature_labels,
    feature_sources
  )

  group_umap <- scoreplot_group_umap(
    srt = srt,
    group.by = group.by,
    groups = groups,
    reduction = reduction,
    dims = dims,
    cells = cells,
    group.palette = group.palette,
    group.palcolor = group_palcolor,
    point.size = point.size,
    raster = raster,
    raster.dpi = raster.dpi,
    title = if (is.null(group.title)) "Celltype" else group.title,
    theme_use = theme_use,
    verbose = verbose,
    ...
  )
  group_legend <- scoreplot_group_legend(
    groups,
    group_colors,
    if (is.null(group.title)) "Celltype" else group.title,
    scoreplot_group_labels(groups, feature_labels, feature_sources)
  )

  umap_plots <- vector("list", length(feature_labels))
  stat_plots <- vector("list", length(feature_labels))
  for (i in seq_along(feature_labels)) {
    if (isTRUE(verbose)) {
      log_message(
        "Plotting {method} feature {.val {feature_labels[[i]]}}",
        message_type = "info"
      )
    }
    umap_plots[[i]] <- scoreplot_umap(
      srt = srt,
      group.by = group.by,
      reduction = reduction,
      dims = dims,
      cells = cells,
      score = scores[, i],
      group = groups,
      feature = feature_labels[[i]],
      title.color = feature_colors[[i]],
      score.colors = score_colors,
      score.label = score_label,
      highlight.group = highlights[[i]],
      group.colors = group_colors,
      highlight.label = highlight.label,
      highlight.mark.type = highlight.mark.type,
      highlight.mark.color = highlight.mark.color,
      highlight.label.color = highlight.label.color,
      ellipse.level = ellipse.level,
      show.legend = feature_labels[[i]] %in% legend_features,
      legend.position = legend.position,
      point.size = point.size,
      raster = raster,
      raster.dpi = raster.dpi,
      theme_use = auc_theme_use,
      verbose = verbose,
      ...
    )
    stat_plots[[i]] <- scoreplot_stat(
      score = scores[, i],
      group = groups,
      group.colors = group_colors,
      threshold = if (is.null(thresholds)) NULL else thresholds[[i]],
      add_box = add_box,
      add_bar = add_bar,
      highlight.group = highlights[[i]],
      highlight.mean.color = highlight.mean.color,
      threshold.colors = threshold.colors,
      bar.text.size = bar.text.size,
      bar.text.color = bar.text.color,
      box.text.size = box.text.size,
      box.text.color = box.text.color,
      point.fraction = point.fraction,
      point.alpha = point.alpha,
      boxplot.y.range = boxplot.y.range,
      seed = seed + i,
      legend.position = legend.position,
      score.label = score_label
    )
  }
  names(umap_plots) <- paste0("umap_", feature_labels)
  names(stat_plots) <- paste0("stats_", feature_labels)

  components <- c(
    list(group_umap = group_umap, group_legend = group_legend),
    umap_plots,
    stat_plots
  )
  if (!isTRUE(combine)) {
    return(components)
  }

  slot_count <- length(feature_labels) + 1L
  if (is.null(ncol)) {
    ncol <- if (is.null(nrow)) {
      slot_count
    } else {
      ceiling(slot_count / (nrow / 2L))
    }
  }
  required_pairs <- ceiling(slot_count / ncol)
  if (is.null(nrow)) {
    nrow <- as.integer(required_pairs * 2L)
  }
  if (nrow / 2L * ncol < slot_count) {
    log_message(
      "The requested layout has too few paired slots for {slot_count} group/feature slots",
      message_type = "error"
    )
  }

  tops <- lapply(
    c(list(group_umap), unname(umap_plots)),
    scoreplot_panel_title
  )
  bottoms <- c(list(group_legend), unname(stat_plots))
  physical <- list()
  pair_rows <- nrow / 2L
  for (block in seq_len(pair_rows)) {
    index <- ((block - 1L) * ncol + 1L):(block * ncol)
    top_row <- lapply(index, function(i) {
      if (i <= length(tops)) tops[[i]] else patchwork::plot_spacer()
    })
    bottom_row <- lapply(index, function(i) {
      if (i <= length(bottoms)) bottoms[[i]] else patchwork::plot_spacer()
    })
    if (block < pair_rows) {
      bottom_row <- lapply(bottom_row, function(p) {
        p &
          ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank()
          )
      })
    }
    physical <- c(physical, top_row, bottom_row)
  }
  patchwork::wrap_plots(physical, ncol = ncol, nrow = nrow, byrow = TRUE) +
    patchwork::plot_layout(
      guides = "collect",
      heights = rep(row.heights, pair_rows)
    ) &
    scoreplot_legend_theme(legend.position) &
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "pt"))
}

scoreplot_group_umap <- function(
  srt,
  group.by,
  groups,
  reduction,
  dims,
  cells,
  group.palette,
  group.palcolor,
  point.size,
  raster,
  raster.dpi,
  title,
  theme_use,
  verbose,
  ...
) {
  dots <- list(...)
  if (
    length(dots) &&
      (is.null(names(dots)) || any(!nzchar(names(dots))))
  ) {
    log_message(
      "Arguments passed via {.arg ...} must be named",
      message_type = "error"
    )
  }
  controlled <- c(
    "srt",
    "group.by",
    "reduction",
    "dims",
    "cells",
    "combine",
    "legend.position",
    "pt.size",
    "raster",
    "raster.dpi"
  )
  conflicts <- intersect(names(dots), controlled)
  if (length(conflicts)) {
    log_message(
      "Arguments passed via {.arg ...} cannot override CellScoringPlot structural arguments: {.val {conflicts}}",
      message_type = "error"
    )
  }
  plot_srt <- srt
  group_column <- group.by
  if (is.null(group_column)) {
    group_column <- ".scop_scoring_group"
    while (group_column %in% colnames(plot_srt@meta.data)) {
      group_column <- paste0(group_column, "_")
    }
    values <- stats::setNames(
      rep(NA_character_, ncol(plot_srt)),
      colnames(plot_srt)
    )
    values[cells] <- as.character(groups)
    plot_srt[[group_column]] <- factor(
      values[colnames(plot_srt)],
      levels = levels(groups)
    )
  }
  defaults <- list(
    srt = plot_srt,
    group.by = group_column,
    reduction = reduction,
    dims = dims,
    cells = cells,
    show_stat = FALSE,
    label = FALSE,
    pt.size = point.size,
    pt.alpha = 1,
    palette = group.palette,
    palcolor = group.palcolor,
    raster = raster,
    raster.dpi = rep(raster.dpi, length.out = 2L),
    aspect.ratio = 1,
    title = title,
    legend.position = "none",
    theme_use = theme_use,
    theme_args = list(base_size = scoreplot_base_size()),
    combine = TRUE,
    force = TRUE,
    verbose = verbose
  )
  plot <- do.call(CellDimPlot, utils::modifyList(defaults, dots))
  plot &
    ggplot2::guides(
      color = "none",
      fill = "none",
      shape = "none",
      linetype = "none",
      linewidth = "none",
      size = "none"
    )
}

scoreplot_umap <- function(
  srt,
  group.by,
  reduction,
  dims,
  cells,
  score,
  group,
  feature,
  title.color,
  score.colors,
  score.label,
  highlight.group,
  group.colors,
  highlight.label,
  highlight.mark.type,
  highlight.mark.color,
  highlight.label.color,
  ellipse.level,
  show.legend,
  legend.position,
  point.size,
  raster,
  raster.dpi,
  theme_use,
  verbose,
  ...
) {
  dots <- list(...)
  if (
    length(dots) &&
      (is.null(names(dots)) || any(!nzchar(names(dots))))
  ) {
    log_message(
      "Arguments passed via {.arg ...} must be named",
      message_type = "error"
    )
  }
  controlled <- c(
    "srt",
    "group.by",
    "reduction",
    "dims",
    "cells",
    "combine",
    "legend.position",
    "pt.size",
    "raster",
    "raster.dpi"
  )
  conflicts <- intersect(names(dots), controlled)
  if (length(conflicts)) {
    log_message(
      "Arguments passed via {.arg ...} cannot override CellScoringPlot structural arguments: {.val {conflicts}}",
      message_type = "error"
    )
  }
  plot_srt <- srt
  group_column <- group.by
  if (is.null(group_column)) {
    group_column <- ".scop_scoring_group"
    while (group_column %in% colnames(plot_srt@meta.data)) {
      group_column <- paste0(group_column, "_")
    }
    values <- stats::setNames(
      rep(NA_character_, ncol(plot_srt)),
      colnames(plot_srt)
    )
    values[cells] <- as.character(group)
    plot_srt[[group_column]] <- factor(
      values[colnames(plot_srt)],
      levels = levels(group)
    )
  }
  defaults <- list(
    srt = plot_srt,
    group.by = group_column,
    reduction = reduction,
    dims = dims,
    cells = cells,
    show_stat = FALSE,
    label = FALSE,
    pt.size = point.size,
    pt.alpha = 1,
    palette = "Chinese",
    raster = raster,
    raster.dpi = rep(raster.dpi, length.out = 2L),
    aspect.ratio = 1,
    title = feature,
    legend.position = if (show.legend) legend.position else "none",
    legend.direction = if (scoreplot_horizontal_legend(legend.position)) {
      "horizontal"
    } else {
      "vertical"
    },
    theme_use = theme_use,
    theme_args = list(base_size = scoreplot_base_size()),
    combine = TRUE,
    force = TRUE,
    verbose = verbose
  )
  plot <- do.call(CellDimPlot, utils::modifyList(defaults, dots))

  score_by_cell <- stats::setNames(as.numeric(score), cells)
  plot$data$group <- plot$data$group.by
  plot$data$score <- unname(score_by_cell[rownames(plot$data)])
  plot$data$Dim1 <- plot$data$x
  plot$data$Dim2 <- plot$data$y
  plot$data$group.by <- plot$data$score
  plot$data <- plot$data[
    order(plot$data$score, na.last = TRUE),
    ,
    drop = FALSE
  ]
  for (i in seq_along(plot$layers)) {
    layer_data <- plot$layers[[i]]$data
    if (
      is.data.frame(layer_data) &&
        "group.by" %in% colnames(layer_data) &&
        nrow(layer_data)
    ) {
      layer_data$group.by <- unname(score_by_cell[rownames(layer_data)])
      layer_data <- layer_data[
        order(layer_data$group.by, na.last = TRUE),
        ,
        drop = FALSE
      ]
      plot$layers[[i]]$data <- layer_data
    }
  }
  color_scale <- vapply(
    plot$scales$scales,
    function(scale) identical(scale$aesthetics, "colour"),
    logical(1)
  )
  plot$scales$scales <- plot$scales$scales[!color_scale]

  finite <- is.finite(plot$data$score)
  limits <- range(plot$data$score[finite], na.rm = TRUE)
  if (!length(limits) || any(!is.finite(limits))) {
    limits <- c(0, 1)
  }
  if (diff(limits) == 0) {
    limits <- limits + c(-0.5, 0.5) * max(abs(limits), 1) * 1e-6
  }

  plot <- plot +
    ggplot2::scale_color_gradientn(
      name = if (show.legend) score.label else NULL,
      colors = score.colors,
      limits = limits,
      na.value = "grey90"
    ) +
    ggplot2::labs(color = if (show.legend) score.label else NULL) +
    ggplot2::guides(
      color = if (show.legend) {
        scoreplot_colorbar(legend.position)
      } else {
        "none"
      }
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = title.color
      ),
      legend.position = if (show.legend) legend.position else "none"
    ) +
    scoreplot_legend_theme(if (show.legend) legend.position else "none")

  highlighted <- plot$data[
    as.character(plot$data$group) %in% highlight.group,
    ,
    drop = FALSE
  ]
  for (highlight in highlight.group) {
    subset <- highlighted[
      as.character(highlighted$group) == highlight,
      ,
      drop = FALSE
    ]
    if (
      nrow(subset) < 3L ||
        length(unique(subset$Dim1)) <= 1L ||
        length(unique(subset$Dim2)) <= 1L
    ) {
      next
    }
    mark_color <- scoreplot_highlight_color(
      highlight.mark.color,
      highlight,
      group.colors[[highlight]]
    )
    if (identical(highlight.mark.type, "hull")) {
      plot <- plot +
        ggforce::geom_mark_hull(
          data = subset,
          ggplot2::aes(x = Dim1, y = Dim2),
          inherit.aes = FALSE,
          concavity = 3.5,
          expand = grid::unit(0.0276, "in"),
          radius = grid::unit(0.0276, "in"),
          color = mark_color,
          fill = NA,
          linewidth = 0.55,
          linetype = "3313",
          show.legend = FALSE
        )
    } else {
      plot <- plot +
        ggplot2::stat_ellipse(
          data = subset,
          ggplot2::aes(x = Dim1, y = Dim2),
          inherit.aes = FALSE,
          type = "t",
          level = ellipse.level,
          color = mark_color,
          linewidth = 0.6,
          linetype = "3313"
        )
    }
    if (isTRUE(highlight.label)) {
      label_color <- scoreplot_highlight_color(
        highlight.label.color,
        highlight,
        "grey40"
      )
      covariance <- stats::cov(
        cbind(subset$Dim1, subset$Dim2),
        use = "complete.obs"
      )
      scale_factor <- sqrt(stats::qf(ellipse.level, 2, nrow(subset) - 1L) * 2)
      label_x <- mean(subset$Dim1, na.rm = TRUE) +
        scale_factor * sqrt(covariance[1, 1])
      label_y <- mean(subset$Dim2, na.rm = TRUE)
      label_hjust <- if (label_x >= mean(range(plot$data$Dim1, na.rm = TRUE))) {
        1
      } else {
        0
      }
      plot <- plot +
        ggplot2::annotate(
          "text",
          x = label_x,
          y = label_y,
          label = highlight,
          color = label_color,
          size = 3.5,
          hjust = label_hjust,
          vjust = 0.5
        )
    }
  }
  plot
}

scoreplot_highlight_color <- function(value, group, fallback) {
  if (is.null(value)) {
    return(unname(fallback))
  }
  if (!is.null(names(value)) && group %in% names(value)) {
    return(unname(value[[group]]))
  }
  unname(value[[1]])
}

scoreplot_named_colors <- function(colors, labels) {
  if (
    is.null(colors) || is.null(names(colors)) || !any(nzchar(names(colors)))
  ) {
    return(colors)
  }
  if (all(labels %in% names(colors))) {
    return(unname(colors[labels]))
  }
  colors
}

scoreplot_stat <- function(
  score,
  group,
  group.colors,
  threshold,
  add_box,
  add_bar,
  highlight.group,
  highlight.mean.color,
  threshold.colors,
  bar.text.size,
  bar.text.color,
  box.text.size,
  box.text.color,
  point.fraction,
  point.alpha,
  boxplot.y.range,
  seed,
  legend.position,
  score.label
) {
  data <- data.frame(group = group, score = as.numeric(score))
  data <- data[is.finite(data$score) & !is.na(data$group), , drop = FALSE]
  box <- NULL
  if (isTRUE(add_box)) {
    box <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = group, y = score, fill = group)
    ) +
      ggplot2::geom_boxplot(width = 0.9, outlier.shape = NA) +
      ggplot2::scale_fill_manual(
        values = group.colors,
        drop = FALSE,
        guide = "none"
      ) +
      ggplot2::labs(x = NULL, y = score.label) +
      ggplot2::theme_bw(base_size = scoreplot_base_size()) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-5, 5, -5, 5)
      )
    if (point.fraction > 0 && nrow(data)) {
      set.seed(seed)
      count <- max(1L, floor(nrow(data) * point.fraction))
      point_data <- data[
        sample.int(nrow(data), min(count, nrow(data))),
        ,
        drop = FALSE
      ]
      box <- box +
        ggplot2::geom_jitter(
          data = point_data,
          width = 0.3,
          height = 0,
          size = 0.2,
          alpha = point.alpha
        )
    }
    score_range <- range(data$score, na.rm = TRUE)
    if (
      isTRUE(boxplot.y.range) &&
        score_range[[1]] >= 0 &&
        score_range[[2]] <= 0.5
    ) {
      box <- box + ggplot2::coord_cartesian(ylim = c(0, 0.5))
    }
    if (length(highlight.group)) {
      mean_score <- mean(
        data$score[as.character(data$group) %in% highlight.group],
        na.rm = TRUE
      )
      if (is.finite(mean_score)) {
        mean_color <- scoreplot_highlight_color(
          highlight.mean.color,
          highlight.group[[1]],
          "grey40"
        )
        label_margin <- grid::unit(5, "pt")
        mean_label <- grid::textGrob(
          sprintf("Mean %s: %.3f", score.label, mean_score),
          x = label_margin,
          y = grid::unit(1, "npc") - label_margin,
          just = c("left", "top"),
          gp = grid::gpar(
            col = box.text.color,
            fontsize = box.text.size * 72 / 25.4
          )
        )
        box <- box +
          ggplot2::geom_hline(
            yintercept = mean_score,
            linetype = "dashed",
            color = mean_color
          ) +
          ggplot2::annotation_custom(
            mean_label,
            xmin = -Inf,
            xmax = Inf,
            ymin = -Inf,
            ymax = Inf
          )
      }
    }
  }

  bar <- NULL
  if (isTRUE(add_bar)) {
    status_levels <- c(
      paste("Above", score.label, "threshold"),
      paste("Below", score.label, "threshold")
    )
    status_colors <- stats::setNames(
      rep(unname(threshold.colors), length.out = 2L),
      status_levels
    )
    total <- table(data$group)
    above <- table(factor(
      data$group[data$score >= threshold],
      levels = levels(data$group)
    ))
    threshold_data <- data.frame(
      group = rep(levels(data$group), each = 2L),
      status = factor(
        rep(
          status_levels,
          times = length(levels(data$group))
        ),
        levels = status_levels
      ),
      count = as.vector(rbind(as.numeric(above), as.numeric(total - above)))
    )
    threshold_data$total <- rep(as.numeric(total), each = 2L)
    threshold_data$proportion <- threshold_data$count /
      pmax(threshold_data$total, 1)
    threshold_data$label <- ifelse(
      threshold_data$status == status_levels[[1]],
      sprintf(
        "%d/%d\n(%.1f%%)",
        threshold_data$count,
        threshold_data$total,
        100 * threshold_data$count / pmax(threshold_data$total, 1)
      ),
      ""
    )
    bar <- ggplot2::ggplot(
      threshold_data,
      ggplot2::aes(x = group, y = proportion, fill = status)
    ) +
      ggplot2::geom_col(width = 0.9, position = "stack") +
      ggplot2::geom_text(
        data = threshold_data[
          threshold_data$status == status_levels[[1]],
          ,
          drop = FALSE
        ],
        ggplot2::aes(y = 0.5, label = label),
        angle = 90,
        lineheight = 0.8,
        size = bar.text.size,
        color = bar.text.color,
        vjust = 0.5,
        hjust = 0.5
      ) +
      ggplot2::scale_fill_manual(
        values = status_colors,
        drop = FALSE,
        name = "Status"
      ) +
      ggplot2::labs(x = "Celltype", y = "Proportion") +
      ggplot2::theme_bw(base_size = scoreplot_base_size()) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-5, 5, 10, 5)
      ) +
      scoreplot_legend_theme(legend.position)
  }
  if (isTRUE(add_box) && isTRUE(add_bar)) {
    return(box / bar + patchwork::plot_layout(heights = c(0.6, 0.4)))
  }
  if (isTRUE(add_box)) {
    return(box)
  }
  if (isTRUE(add_bar)) {
    return(bar)
  }
  patchwork::plot_spacer()
}

scoreplot_horizontal_legend <- function(legend.position) {
  legend.position %in% c("bottom", "top")
}

scoreplot_panel_title <- function(plot) {
  if (!inherits(plot, "ggplot") || is.null(plot$labels$title)) {
    return(plot)
  }
  title <- plot$labels$title
  title_theme <- plot$theme$plot.title
  plot$labels$title <- NULL
  plot$coordinates$clip <- "off"
  plot +
    ggplot2::annotation_custom(
      grid::textGrob(
        title,
        x = 0,
        y = grid::unit(1, "npc") + grid::unit(2, "pt"),
        just = c("left", "bottom"),
        gp = grid::gpar(
          col = title_theme$colour %||% "black",
          fontsize = title_theme$size %||% scoreplot_base_size()
        )
      ),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}

scoreplot_colorbar <- function(legend.position) {
  horizontal <- scoreplot_horizontal_legend(legend.position)
  ggplot2::guide_colorbar(
    title.position = if (horizontal) "left" else "top",
    title.hjust = 0.5,
    title.vjust = if (horizontal) 1 else 0.5,
    theme = ggplot2::theme(
      legend.text = ggplot2::element_text(
        margin = ggplot2::margin(t = if (horizontal) 1 else 0, unit = "pt")
      )
    ),
    barwidth = if (horizontal) grid::unit(42, "pt") else grid::unit(5, "pt"),
    barheight = if (horizontal) grid::unit(4, "pt") else grid::unit(42, "pt")
  )
}

scoreplot_legend_theme <- function(legend.position) {
  horizontal <- scoreplot_horizontal_legend(legend.position)
  ggplot2::theme(
    legend.position = legend.position,
    legend.direction = if (horizontal) "horizontal" else "vertical",
    legend.box = if (horizontal) "horizontal" else "vertical",
    legend.key.size = grid::unit(7, "pt"),
    legend.text = ggplot2::element_text(size = 6),
    legend.title = ggplot2::element_text(size = 7),
    legend.spacing.x = grid::unit(3, "pt"),
    legend.spacing.y = grid::unit(2, "pt")
  )
}

scoreplot_group_legend <- function(
  group,
  colors,
  title,
  labels
) {
  levels_group <- levels(group)
  data <- data.frame(
    x = 0.3,
    y = seq_along(levels_group),
    group = factor(levels_group, levels = levels_group),
    label = labels
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = group), size = 2.5, shape = 15) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      x = 0.5,
      hjust = 0,
      vjust = 0.5,
      color = "black",
      size = 2.1
    ) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE, guide = "none") +
    ggplot2::coord_cartesian(
      xlim = c(0.2, 2.5),
      ylim = c(0.5, length(levels_group) + 1.2),
      clip = "off"
    ) +
    ggplot2::theme_void(base_size = scoreplot_base_size()) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(-5, 5, 10, 5)
    )
  if (!is.null(title)) {
    plot <- plot +
      ggplot2::annotate(
        "text",
        x = 0.4,
        y = length(levels_group) + 0.8,
        label = title,
        hjust = 0,
        vjust = 0,
        size = 2.45
      )
  }
  plot
}

scoreplot_input <- function(srt, method = "AUCell") {
  provenance <- srt@misc[["CellScoring"]]
  records <- if (is.list(provenance)) provenance$records else NULL
  if (!is.list(records) || !length(records)) {
    return(NULL)
  }
  for (i in rev(seq_along(records))) {
    record <- records[[i]]
    if (!identical(record$method, method)) {
      next
    }
    features <- record$features
    if (
      !is.character(features) ||
        !length(features) ||
        is.null(names(features)) ||
        any(!nzchar(names(features)))
    ) {
      next
    }
    sources <- unname(features)
    if (all(sources %in% colnames(srt@meta.data))) {
      return(list(scores = NULL, features = features))
    }
    assay <- record$assay
    if (
      is.character(assay) &&
        length(assay) == 1L &&
        nzchar(assay) &&
        assay %in% names(srt@assays)
    ) {
      scores <- GetAssayData5(srt, assay = assay, layer = "counts")
      assay_features <- record$assay_features
      if (!is.character(assay_features) || !length(assay_features)) {
        assay_features <- features
      }
      if (all(unname(assay_features) %in% rownames(scores))) {
        return(list(scores = scores, features = assay_features))
      }
    }
  }
  NULL
}

scoreplot_scores <- function(srt, score_input, features, cells) {
  if (is.null(score_input)) {
    metadata <- srt@meta.data[cells, , drop = FALSE]
    if (is.null(features)) {
      numeric_columns <- names(metadata)[vapply(
        metadata,
        is.numeric,
        logical(1)
      )]
      if (!length(numeric_columns)) {
        log_message(
          "No numeric score columns were found in Seurat metadata",
          message_type = "error"
        )
      }
      features <- numeric_columns
    }
    feature_map <- scoreplot_feature_map(features)
    missing <- setdiff(feature_map$sources, colnames(metadata))
    if (length(missing)) {
      log_message(
        "Score columns are absent from Seurat metadata: {.val {missing}}",
        message_type = "error"
      )
    }
    non_numeric <- feature_map$sources[
      !vapply(metadata[feature_map$sources], is.numeric, logical(1))
    ]
    if (length(non_numeric)) {
      log_message(
        "Score columns must be numeric: {.val {non_numeric}}",
        message_type = "error"
      )
    }
    scores <- as.matrix(metadata[, feature_map$sources, drop = FALSE])
  } else {
    score_input <- as.matrix(score_input)
    storage.mode(score_input) <- "double"
    if (is.null(rownames(score_input)) || is.null(colnames(score_input))) {
      log_message(
        "{.arg scores} must have both row and column names",
        message_type = "error"
      )
    }
    cells_in_rows <- all(cells %in% rownames(score_input))
    cells_in_columns <- all(cells %in% colnames(score_input))
    if (cells_in_columns && !cells_in_rows) {
      source_names <- rownames(score_input)
      if (is.null(features)) {
        features <- source_names
      }
      feature_map <- scoreplot_feature_map(features)
      missing <- setdiff(feature_map$sources, source_names)
      if (length(missing)) {
        log_message(
          "Features are absent from {.arg scores}: {.val {missing}}",
          message_type = "error"
        )
      }
      scores <- t(score_input[feature_map$sources, cells, drop = FALSE])
    } else if (cells_in_rows && !cells_in_columns) {
      source_names <- colnames(score_input)
      if (is.null(features)) {
        features <- source_names
      }
      feature_map <- scoreplot_feature_map(features)
      missing <- setdiff(feature_map$sources, source_names)
      if (length(missing)) {
        log_message(
          "Features are absent from {.arg scores}: {.val {missing}}",
          message_type = "error"
        )
      }
      scores <- score_input[cells, feature_map$sources, drop = FALSE]
    } else if (cells_in_columns && cells_in_rows) {
      log_message(
        "Cannot infer {.arg scores} orientation because cell names occur on both axes",
        message_type = "error"
      )
    } else {
      log_message(
        "Neither axis of {.arg scores} contains all requested cell names",
        message_type = "error"
      )
    }
  }
  colnames(scores) <- feature_map$labels
  rownames(scores) <- cells
  list(
    scores = scores,
    labels = feature_map$labels,
    sources = feature_map$sources
  )
}

scoreplot_feature_map <- function(features) {
  labels <- names(features)
  features <- as.character(features)
  if (!length(features) || anyNA(features) || any(!nzchar(features))) {
    log_message(
      "{.arg features} must contain at least one non-empty feature name",
      message_type = "error"
    )
  }
  if (is.null(labels)) {
    labels <- features
  }
  labels[!nzchar(labels)] <- features[!nzchar(labels)]
  if (anyDuplicated(labels)) {
    log_message("Feature display labels must be unique", message_type = "error")
  }
  list(labels = labels, sources = unname(features))
}

scoreplot_groups <- function(srt, group.by, cells) {
  if (is.null(group.by)) {
    values <- SeuratObject::Idents(srt)[cells]
  } else {
    if (length(group.by) != 1L || !group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} must name one Seurat metadata column",
        message_type = "error"
      )
    }
    values <- srt@meta.data[cells, group.by]
  }
  if (anyNA(values)) {
    log_message("Cell groups contain missing values", message_type = "error")
  }
  if (is.factor(values)) {
    factor(values, levels = levels(droplevels(values)))
  } else {
    factor(values, levels = unique(as.character(values)))
  }
}

scoreplot_thresholds <- function(
  thresholds,
  labels,
  sources,
  scores,
  add_bar,
  method,
  seed,
  verbose
) {
  automatic <- identical(thresholds, "auto") ||
    (is.null(thresholds) && isTRUE(add_bar))
  if (automatic) {
    if (!isTRUE(add_bar)) {
      return(NULL)
    }
    if (!identical(method, "AUCell")) {
      log_message(
        "automatic thresholds are only available for AUCell; supply explicit {.arg thresholds} for method {.val {method}}",
        message_type = "error"
      )
    }
    return(scoreplot_auto_thresholds(scores, labels, seed, verbose))
  }
  if (is.null(thresholds)) {
    return(NULL)
  }
  if (
    !is.numeric(thresholds) || anyNA(thresholds) || any(!is.finite(thresholds))
  ) {
    log_message(
      "{.arg thresholds} must be {.val \"auto\"} or contain finite numeric values",
      message_type = "error"
    )
  }
  if (is.null(names(thresholds))) {
    if (length(thresholds) == 1L) {
      thresholds <- rep(thresholds, length(labels))
    }
    if (length(thresholds) != length(labels)) {
      log_message(
        "Unnamed {.arg thresholds} must have length 1 or one value per feature",
        message_type = "error"
      )
    }
    return(stats::setNames(as.numeric(thresholds), labels))
  }
  result <- vapply(
    seq_along(labels),
    function(i) {
      key <- if (labels[[i]] %in% names(thresholds)) {
        labels[[i]]
      } else {
        sources[[i]]
      }
      if (!key %in% names(thresholds)) {
        log_message(
          "No threshold was supplied for feature {.val {labels[[i]]}}",
          message_type = "error"
        )
      }
      thresholds[[key]]
    },
    numeric(1)
  )
  stats::setNames(result, labels)
}

scoreplot_auto_thresholds <- function(scores, labels, seed, verbose) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    log_message(
      "Package {.pkg AUCell} is required to infer automatic thresholds",
      message_type = "error"
    )
  }
  threshold_results <- AUCell::aucellResults(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(AUC = t(scores))
    )
  )
  set.seed(seed)
  thresholds <- AUCell::getThresholdSelected(
    AUCell::AUCell_exploreThresholds(
      threshold_results,
      nCores = 1,
      plotHist = FALSE,
      assignCells = TRUE,
      verbose = verbose
    )
  )
  if (
    !setequal(names(thresholds), labels) ||
      anyNA(thresholds) ||
      any(!is.finite(thresholds))
  ) {
    log_message(
      "AUCell did not return one finite automatic threshold per feature",
      message_type = "error"
    )
  }
  stats::setNames(as.numeric(thresholds[labels]), labels)
}

scoreplot_highlights <- function(labels, sources, groups) {
  result <- lapply(seq_along(labels), function(i) {
    if (sources[[i]] %in% groups) {
      sources[[i]]
    } else if (labels[[i]] %in% groups) {
      labels[[i]]
    } else {
      character()
    }
  })
  stats::setNames(result, labels)
}

scoreplot_group_labels <- function(group, labels, sources) {
  groups <- levels(group)
  display <- stats::setNames(groups, groups)
  highlights <- scoreplot_highlights(labels, sources, groups)
  for (i in seq_along(highlights)) {
    if (length(highlights[[i]])) {
      display[[highlights[[i]]]] <- labels[[i]]
    }
  }
  counts <- tabulate(as.integer(group), nbins = length(groups))
  stats::setNames(
    sprintf(
      "%s (n = %s)",
      unname(display[groups]),
      formatC(counts, format = "d", big.mark = ",")
    ),
    groups
  )
}

scoreplot_legend_features <- function(value, labels, sources) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    return(if (value) labels[[1]] else character())
  }
  value <- as.character(value)
  matched <- labels[labels %in% value | sources %in% value]
  unknown <- setdiff(value, c(labels, sources))
  if (length(unknown)) {
    log_message(
      "{.arg show.score.legend} contains unknown features: {.val {unknown}}",
      message_type = "error"
    )
  }
  matched
}

scoreplot_validate_layout <- function(value, name) {
  if (is.null(value)) {
    return(invisible(NULL))
  }
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      value < 1 ||
      value != as.integer(value)
  ) {
    log_message(
      "{.arg {name}} must be a positive integer or NULL",
      message_type = "error"
    )
  }
  invisible(NULL)
}

scoreplot_base_size <- function() {
  7
}
