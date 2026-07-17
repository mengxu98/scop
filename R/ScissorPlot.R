#' @title Plot Scissor results
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams FeatureHeatmap
#' @inheritParams thisplot::StatPlot
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object after [RunScissor].
#' @param plot_type Plot type. `"umap"` shows embedding panels, `"heatmap"`
#' shows a `FeatureHeatmap`, and statistical views such as `"bar"` and
#' `"upset"` are drawn with [thisplot::StatPlot].
#' @param prefix Prefix used by [RunScissor].
#' @param group.by Optional metadata column shown together with Scissor status.
#' For `"heatmap"`, the default is the Scissor status column. For statistical
#' plots, it is passed to [thisplot::StatPlot], except that `"upset"` uses it
#' to split Scissor status distributions by group.
#' @param nfeatures Number of features to show when `plot_type = "heatmap"` and
#' `features = NULL`.
#' @param feature_method Method used to rank heatmap features when
#' `features = NULL`. `"variance"` ranks genes by variance in selected cells,
#' `"status_diff"` ranks by the largest mean-expression difference between
#' Scissor status groups, `"coef_cor"` ranks by absolute correlation with
#' Scissor coefficients, and `"input_order"` keeps the [RunScissor] input order.
#' @param tool_name Name of the `srt@tools` entry created by [RunScissor].
#' @param status Scissor status levels included in the heatmap.
#' @param include.background Whether to include background cells in heatmap and
#' background status in upset plots.
#' @param upset_top_n Maximum number of `group.by` levels to show in
#' `plot_type = "upset"`. The most frequent levels are kept. `NULL` keeps all.
#' @param combine Whether to combine UMAP panels or StatPlot panels.
#' @param nrow,ncol,byrow Layout parameters passed to `patchwork::wrap_plots()`
#' or [thisplot::StatPlot].
#' @param theme_use,theme_args Theme used by [thisplot::StatPlot].
#' @param ... Additional arguments passed to [CellDimPlot], [FeatureDimPlot],
#' [FeatureHeatmap], or [thisplot::StatPlot], depending on `plot_type`.
#'
#' @return A ggplot/patchwork object for embedding and statistical plots, or a
#' list returned by [FeatureHeatmap] for `plot_type = "heatmap"`.
#' @export
#'
#' @seealso [ScissorPlot]
#'
#' @examples
#' data(panc8_sub)
#' data(islet_bulk)
#' panc8_sub <- RunScissor(
#'   panc8_sub,
#'   bulk_dataset = islet_bulk,
#'   condition.by = "condition",
#'   positive = "bfa",
#'   family = "binomial",
#'   features = head(intersect(
#'     rownames(panc8_sub),
#'     rownames(SummarizedExperiment::assay(islet_bulk, "counts"))
#'   ), 1000),
#'   alpha = 0.2,
#'   cutoff = 0.5
#' )
#' panc8_sub <- standard_scop(panc8_sub, verbose = FALSE)
#'
#' ScissorPlot(
#'   panc8_sub,
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   group.by = "celltype",
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#'
#' ht <- ScissorPlot(
#'   panc8_sub,
#'   plot_type = "heatmap",
#'   group.by = "celltype"
#' )
#' ht$plot
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "bar",
#'   group.by = "celltype"
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "upset"
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "rose",
#'   label = TRUE
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "ring",
#'   label = TRUE
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "pie",
#'   label = TRUE
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   plot_type = "dot",
#'   label = TRUE
#' )
ScissorPlot <- function(
  srt,
  plot_type = c(
    "umap",
    "heatmap",
    "bar",
    "upset",
    "rose",
    "ring",
    "pie",
    "dot"
  ),
  reduction = NULL,
  prefix = "Scissor",
  group.by = NULL,
  split.by = NULL,
  features = NULL,
  nfeatures = 50,
  feature_method = c("variance", "status_diff", "coef_cor", "input_order"),
  tool_name = "Scissor",
  status = c("Scissor+", "Scissor-"),
  include.background = TRUE,
  upset_top_n = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  max_cells = 100,
  cell_order = NULL,
  exp_method = "zscore",
  stat_type = c("percent", "count"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "RdBu",
  group_palette = "Chinese",
  group_palcolor = NULL,
  cell_annotation = NULL,
  cell_annotation_palette = "Chinese",
  cell_annotation_palcolor = NULL,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
) {
  context <- scissor_plot_context(srt = srt, prefix = prefix)
  srt <- context$srt
  coef_col <- context$coef_col
  status_col <- context$status_col
  status_colors <- context$status_colors
  plot_type <- match.arg(plot_type)
  feature_method <- match.arg(feature_method)

  if (identical(plot_type, "heatmap")) {
    return(scissor_heatmap_plot(
      srt = srt,
      features = features,
      nfeatures = nfeatures,
      feature_method = feature_method,
      prefix = prefix,
      tool_name = tool_name,
      status = status,
      include.background = include.background,
      cells = cells,
      group.by = group.by,
      split.by = split.by,
      layer = layer,
      assay = assay,
      max_cells = max_cells,
      cell_order = cell_order,
      exp_method = exp_method,
      heatmap_palette = heatmap_palette,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      cell_annotation = cell_annotation,
      cell_annotation_palette = cell_annotation_palette,
      cell_annotation_palcolor = cell_annotation_palcolor,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      coef_col = coef_col,
      status_col = status_col,
      status_colors = status_colors,
      verbose = verbose,
      ...
    ))
  }

  if (!identical(plot_type, "umap")) {
    stat_type <- match.arg(stat_type)
    return(scissor_stat_plot(
      srt = srt,
      plot_type = plot_type,
      status_col = status_col,
      status_colors = status_colors,
      group.by = group.by,
      split.by = split.by,
      include.background = include.background,
      upset_top_n = upset_top_n,
      stat_type = stat_type,
      palette = palette,
      palcolor = palcolor %||% status_colors,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }

  scissor_umap_plot(
    srt = srt,
    reduction = reduction,
    coef_col = coef_col,
    status_col = status_col,
    status_colors = status_colors,
    group.by = group.by,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palette = palette,
    palcolor = palcolor,
    ...
  )
}

scissor_plot_context <- function(srt, prefix) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat}",
      message_type = "error"
    )
  }

  coef_col <- paste0(prefix, "_coef")
  status_col <- paste0(prefix, "_status")
  missing_cols <- setdiff(c(coef_col, status_col), colnames(srt@meta.data))
  if (length(missing_cols) > 0L) {
    log_message(
      "Missing Scissor results: {.val {missing_cols}}, please run {.fn RunScissor} first",
      message_type = "error"
    )
  }
  status_colors <- c(
    "Scissor+" = "#D70440",
    "Scissor-" = "#2C7BB6",
    "Background" = "#D9D9D9"
  )
  srt@meta.data[[status_col]] <- factor(
    as.character(srt@meta.data[[status_col]]),
    levels = names(status_colors)
  )
  list(
    srt = srt,
    coef_col = coef_col,
    status_col = status_col,
    status_colors = status_colors
  )
}

scissor_umap_plot <- function(
  srt,
  reduction,
  coef_col,
  status_col,
  status_colors,
  group.by,
  combine,
  nrow,
  ncol,
  byrow,
  pt.size,
  pt.alpha,
  palette,
  palcolor,
  ...
) {
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "Column {.val {group.by}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  }

  plist <- list()
  plist[["Status"]] <- CellDimPlot(
    srt,
    group.by = status_col,
    reduction = reduction,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    palcolor = status_colors,
    combine = FALSE,
    ...
  )[[1]] +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = "Scissor"),
      fill = ggplot2::guide_legend(title = "Scissor")
    )

  plist[["Coefficient"]] <- FeatureDimPlot(
    srt,
    features = coef_col,
    reduction = reduction,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    combine = FALSE,
    ...
  )[[1]]
  plist[["Coefficient"]] <- scissor_drop_aes_scales(
    plist[["Coefficient"]],
    c("colour", "color")
  ) +
    ggplot2::scale_colour_gradient2(
      low = status_colors[["Scissor-"]],
      mid = "#F7F7F7",
      high = status_colors[["Scissor+"]],
      midpoint = 0
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(title = "Coefficient")
    )

  if (!is.null(group.by)) {
    plist[["Group"]] <- CellDimPlot(
      srt,
      group.by = group.by,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      palette = palette,
      palcolor = palcolor,
      combine = FALSE,
      ...
    )[[1]]
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1L) {
      if (is.null(nrow) && is.null(ncol)) {
        nrow <- 1L
        ncol <- length(plist)
      }
      return(patchwork::wrap_plots(
        plotlist = plist,
        nrow = nrow,
        ncol = ncol,
        byrow = byrow
      ))
    }
    return(plist[[1]])
  }
  plist
}

scissor_stat_plot <- function(
  srt,
  plot_type,
  status_col,
  status_colors,
  group.by,
  split.by,
  include.background,
  upset_top_n,
  stat_type,
  palette,
  palcolor,
  combine,
  nrow,
  ncol,
  byrow,
  theme_use,
  theme_args,
  ...
) {
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "Column {.val {group.by}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.null(split.by) && !split.by %in% colnames(srt@meta.data)) {
    log_message(
      "Column {.val {split.by}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (
    !is.null(upset_top_n) &&
      (length(upset_top_n) != 1L ||
        !is.numeric(upset_top_n) ||
        is.na(upset_top_n) ||
        !is.finite(upset_top_n) ||
        upset_top_n < 1L)
  ) {
    log_message(
      "{.arg upset_top_n} must be a positive number or {.val NULL}",
      message_type = "error"
    )
  }
  if (identical(plot_type, "upset")) {
    upset_data <- scissor_upset_data(
      meta.data = srt@meta.data,
      status_col = status_col,
      group.by = group.by,
      split.by = split.by,
      status_colors = status_colors,
      include.background = include.background,
      upset_top_n = upset_top_n
    )
    return(thisplot::StatPlot(
      meta.data = upset_data$meta.data,
      stat.by = upset_data$stat.by,
      stat_level = rep("yes", length(upset_data$stat.by)),
      split.by = upset_data$split.by,
      plot_type = "upset",
      stat_type = "count",
      palette = palette,
      palcolor = palcolor,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ))
  }
  thisplot::StatPlot(
    meta.data = srt@meta.data,
    stat.by = status_col,
    group.by = group.by,
    split.by = split.by,
    plot_type = plot_type,
    stat_type = stat_type,
    palette = palette,
    palcolor = palcolor,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    theme_use = theme_use,
    theme_args = theme_args,
    ...
  )
}

scissor_upset_data <- function(
  meta.data,
  status_col,
  group.by,
  split.by,
  status_colors,
  include.background,
  upset_top_n
) {
  meta_use <- meta.data
  status_columns <- c(
    "Scissor+" = "Scissor_positive",
    "Scissor-" = "Scissor_negative",
    "Background" = "Scissor_background"
  )
  status_keep <- c("Scissor+", "Scissor-")
  if (isTRUE(include.background)) {
    status_keep <- c(status_keep, "Background")
  }
  for (status in names(status_columns)) {
    meta_use[[status_columns[[status]]]] <- ifelse(
      meta_use[[status_col]] == status,
      "yes",
      NA
    )
  }
  stat.by <- unname(status_columns[intersect(
    status_keep,
    names(status_colors)
  )])

  if (!is.null(group.by)) {
    group_values <- as.character(meta_use[[group.by]])
    if (!is.null(upset_top_n)) {
      group_counts <- sort(table(group_values), decreasing = TRUE)
      keep_groups <- names(utils::head(group_counts, upset_top_n))
      meta_use <- meta_use[group_values %in% keep_groups, , drop = FALSE]
      group_values <- as.character(meta_use[[group.by]])
    }
    split_col <- paste0("Scissor_", group.by, "_distribution")
    if (!is.null(split.by)) {
      split_values <- as.character(meta_use[[split.by]])
      meta_use[[split_col]] <- factor(
        paste(group_values, split_values, sep = " | "),
        levels = unique(paste(group_values, split_values, sep = " | "))
      )
    } else {
      meta_use[[split_col]] <- factor(
        group_values,
        levels = unique(group_values)
      )
    }
    split.by <- split_col
  }

  list(meta.data = meta_use, stat.by = stat.by, split.by = split.by)
}

scissor_heatmap_plot <- function(
  srt,
  features,
  nfeatures,
  feature_method,
  prefix,
  tool_name,
  status,
  include.background,
  cells,
  group.by,
  split.by,
  layer,
  assay,
  max_cells,
  cell_order,
  exp_method,
  heatmap_palette,
  group_palette,
  group_palcolor,
  cell_annotation,
  cell_annotation_palette,
  cell_annotation_palcolor,
  show_row_names,
  show_column_names,
  cluster_rows,
  cluster_columns,
  coef_col,
  status_col,
  status_colors,
  verbose,
  ...
) {
  status <- match.arg(status, choices = names(status_colors), several.ok = TRUE)
  if (isTRUE(include.background)) {
    status <- union(status, "Background")
  }

  cells_status <- rownames(srt@meta.data)[
    srt@meta.data[[status_col]] %in% status
  ]
  cells <- cells %||% cells_status
  cells <- intersect(cells, cells_status)
  if (length(cells) == 0L) {
    log_message(
      "No cells match requested Scissor status levels",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  data_layer <- suppressWarnings(
    tryCatch(
      GetAssayData5(srt, assay = assay, layer = "data"),
      error = function(e) NULL
    )
  )
  if (is.null(data_layer) || nrow(data_layer) == 0L || ncol(data_layer) == 0L) {
    srt <- Seurat::NormalizeData(srt, assay = assay, verbose = FALSE)
  }
  if (is.null(features)) {
    features <- scissor_select_heatmap_features(
      srt = srt,
      cells = cells,
      nfeatures = nfeatures,
      feature_method = feature_method,
      tool_name = tool_name,
      assay = assay,
      layer = layer,
      coef_col = coef_col,
      status_col = status_col
    )
  }
  if (length(features) == 0L) {
    log_message(
      "No heatmap features are available",
      message_type = "error"
    )
  }

  group.by <- group.by %||% status_col
  if (!all(group.by %in% colnames(srt@meta.data))) {
    log_message(
      "Some {.arg group.by} columns are not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.null(split.by) && !split.by %in% colnames(srt@meta.data)) {
    log_message(
      "Column {.val {split.by}} not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!identical(group.by, status_col)) {
    annotation_info <- scissor_add_status_cell_annotation(
      cell_annotation = cell_annotation,
      cell_annotation_palette = cell_annotation_palette,
      cell_annotation_palcolor = cell_annotation_palcolor,
      status_col = status_col,
      status_colors = status_colors
    )
    cell_annotation <- annotation_info$cell_annotation
    cell_annotation_palette <- annotation_info$cell_annotation_palette
    cell_annotation_palcolor <- annotation_info$cell_annotation_palcolor
  }
  if (is.null(group_palcolor) && identical(group.by, status_col)) {
    group_palcolor <- status_colors
  }
  if (is.null(cell_order)) {
    coef <- srt@meta.data[cells, coef_col]
    cell_order <- names(sort(stats::setNames(coef, cells), decreasing = TRUE))
  }

  FeatureHeatmap(
    srt = srt,
    features = features,
    cells = cells,
    group.by = group.by,
    split.by = split.by,
    max_cells = max_cells,
    cell_order = cell_order,
    layer = layer,
    assay = assay,
    exp_method = exp_method,
    heatmap_palette = heatmap_palette,
    group_palette = group_palette,
    group_palcolor = group_palcolor,
    cell_annotation = cell_annotation,
    cell_annotation_palette = cell_annotation_palette,
    cell_annotation_palcolor = cell_annotation_palcolor,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    verbose = verbose,
    ...
  )
}

scissor_drop_aes_scales <- function(plot, aesthetics) {
  if (!inherits(plot, c("gg", "ggplot")) || is.null(plot$scales)) {
    return(plot)
  }
  plot$scales$scales <- Filter(
    function(scale) {
      !any(tolower(scale$aesthetics) %in% aesthetics)
    },
    plot$scales$scales
  )
  plot
}

scissor_add_status_cell_annotation <- function(
  cell_annotation,
  cell_annotation_palette,
  cell_annotation_palcolor,
  status_col,
  status_colors
) {
  annotation_original <- cell_annotation
  palette_original <- scissor_normalize_cell_annotation_palette(
    palette = cell_annotation_palette,
    n = length(annotation_original)
  )
  palcolor_original <- scissor_normalize_cell_annotation_palcolor(
    palcolor = cell_annotation_palcolor,
    n = length(annotation_original)
  )
  cell_annotation <- unique(c(status_col, annotation_original))
  cell_annotation_palette <- rep("Chinese", length(cell_annotation))
  cell_annotation_palcolor <- vector("list", length(cell_annotation))
  for (i in seq_along(cell_annotation)) {
    if (identical(cell_annotation[[i]], status_col)) {
      cell_annotation_palcolor[i] <- list(status_colors)
    } else {
      cell_annotation_palette[[i]] <- palette_original[[
        match(cell_annotation[[i]], annotation_original)
      ]]
      cell_annotation_palcolor[i] <- list(
        palcolor_original[[match(cell_annotation[[i]], annotation_original)]]
      )
    }
  }
  list(
    cell_annotation = cell_annotation,
    cell_annotation_palette = cell_annotation_palette,
    cell_annotation_palcolor = cell_annotation_palcolor
  )
}

scissor_normalize_cell_annotation_palette <- function(palette, n) {
  if (n == 0L) {
    return(character())
  }
  if (length(palette) == 0L) {
    return(rep("Chinese", n))
  }
  if (length(palette) == 1L) {
    return(rep(palette, n))
  }
  if (length(palette) != n) {
    log_message(
      "{.arg cell_annotation_palette} must be length 1 or the same length as {.arg cell_annotation}",
      message_type = "error"
    )
  }
  palette
}

scissor_normalize_cell_annotation_palcolor <- function(palcolor, n) {
  if (n == 0L) {
    return(list())
  }
  if (is.null(palcolor) || length(palcolor) == 0L) {
    return(rep(list(NULL), n))
  }
  if (!is.list(palcolor)) {
    return(rep(list(palcolor), n))
  }
  if (length(palcolor) == 1L) {
    return(rep(palcolor, n))
  }
  if (length(palcolor) != n) {
    log_message(
      "{.arg cell_annotation_palcolor} must be length 1 or the same length as {.arg cell_annotation}",
      message_type = "error"
    )
  }
  palcolor
}

scissor_select_heatmap_features <- function(
  srt,
  cells,
  nfeatures,
  feature_method,
  tool_name,
  assay,
  layer,
  coef_col,
  status_col
) {
  tool <- srt@tools[[tool_name]]
  features <- tool$input_summary$features
  if (length(features) == 0L) {
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
  }
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  if (length(features) == 0L) {
    features <- rownames(mat)
  }
  features <- intersect(features, rownames(mat))
  if (length(features) == 0L) {
    return(features)
  }
  if (identical(feature_method, "input_order")) {
    return(utils::head(features, nfeatures))
  }

  mat_select <- mat[features, cells, drop = FALSE]
  if (identical(feature_method, "status_diff")) {
    scores <- scissor_feature_status_diff(
      mat = mat_select,
      status = srt@meta.data[cells, status_col]
    )
  } else if (identical(feature_method, "coef_cor")) {
    scores <- scissor_feature_coef_cor(
      mat = mat_select,
      coef = srt@meta.data[cells, coef_col]
    )
  } else {
    scores <- scissor_feature_variance(mat_select)
  }
  scores[!is.finite(scores)] <- 0
  utils::head(names(sort(scores, decreasing = TRUE)), nfeatures)
}

scissor_feature_variance <- function(mat) {
  Matrix::rowMeans(mat^2) - Matrix::rowMeans(mat)^2
}

scissor_feature_group_range <- function(means) {
  check_r("matrixStats", verbose = FALSE)
  out <- matrixStats::rowMaxs(means, na.rm = TRUE) -
    matrixStats::rowMins(means, na.rm = TRUE)
  names(out) <- rownames(means)
  out
}

scissor_feature_status_diff <- function(mat, status) {
  status <- droplevels(factor(status))
  levels_use <- levels(status)
  levels_use <- levels_use[!is.na(levels_use)]
  if (length(levels_use) < 2L) {
    return(scissor_feature_variance(mat))
  }
  means <- vapply(
    levels_use,
    function(level) Matrix::rowMeans(mat[, status == level, drop = FALSE]),
    numeric(nrow(mat))
  )
  rownames(means) <- rownames(mat)
  scissor_feature_group_range(means)
}

scissor_feature_coef_cor <- function(mat, coef) {
  coef <- as.numeric(coef)
  keep <- is.finite(coef)
  if (sum(keep) < 3L || stats::var(coef[keep]) == 0) {
    return(scissor_feature_variance(mat))
  }
  mat <- mat[, keep, drop = FALSE]
  coef <- coef[keep] - mean(coef[keep])
  numerator <- as.numeric(mat %*% coef)
  denom <- sqrt(scissor_feature_variance(mat) * ncol(mat) * sum(coef^2))
  stats::setNames(abs(numerator / denom), rownames(mat))
}
