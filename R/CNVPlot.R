#' @title Plot CNV results
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams SpatialSpotPlot
#' @param plot_type Plot type: `"heatmap"`, `"dim"`, `"spatial"`, `"bar"`, or `"tree"`.
#' Heatmaps show cells in rows and features in columns. Use `cluster_columns`
#' for cells or `cluster_rows` for features; both default to `FALSE`.
#' @param tool_name Name of the stored CNV results.
#' @param method CNV method. `NULL` uses the active method.
#' @param group.by A metadata column identifying cell groups.
#' @param value Metadata column for `"dim"`, `"spatial"`, or `"bar"`.
#' @param cluster_tree_by Build the tree from cells or CNV clusters.
#' @param palette,palcolor Colors for discrete groups.
#' @param heatmap_palette,heatmap_palcolor Colors for heatmap scores.
#' Defaults to muted blue, white, and red.
#' @param title,subtitle Plot title and ggplot subtitle.
#' @param heatmap_limits Color limits. `NULL` uses `c(0.75, 1.25)` for inferCNV
#' and `c(-1, 1)` otherwise. `"auto"` uses the 95th percentile of absolute
#' deviations from the method's midpoint (1 or 0), with the method limits as
#' fallback. A numeric pair sets the limits and must span the midpoint.
#' @param ... Additional parameters to other functions.
#'
#' @return A `ComplexHeatmap::Heatmap` for `"heatmap"`, an `hclust` for
#' `"tree"`, or a `ggplot`, `patchwork`, or list of `ggplot` objects otherwise.
#'
#' @seealso [RunCNV()]
#' @export
#' @inherit RunCNV examples source
CNVPlot <- function(
  object,
  plot_type = c("heatmap", "dim", "spatial", "bar", "tree"),
  tool_name = "CNV",
  method = NULL,
  group.by = NULL,
  reduction = NULL,
  value = NULL,
  cluster_tree_by = c("cluster", "cell"),
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  heatmap_limits = NULL,
  srt = NULL,
  ...
) {
  check_r("withr", verbose = FALSE)
  withr::local_preserve_seed()
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  plot_type <- match.arg(plot_type)
  cluster_tree_by <- match.arg(cluster_tree_by)
  bundle <- cnv_get_method_bundle(srt, tool_name = tool_name, method = method)
  method <- bundle$method %||% method %||% srt@tools[[tool_name]]$active_method

  switch(plot_type,
    heatmap = {
      if (missing(heatmap_palette) && is.null(heatmap_palcolor)) {
        heatmap_palcolor <- c("#285E89", "#F9F9F7", "#AB403C")
      }
      cnv_plot_heatmap(
        srt = srt,
        bundle = bundle,
        method = method,
        group.by = group.by,
        palette = palette,
        palcolor = palcolor,
        heatmap_palette = heatmap_palette,
        heatmap_palcolor = heatmap_palcolor,
        heatmap_limits = heatmap_limits,
        title = title,
        ...
      )
    },
    dim = cnv_plot_dim(
      srt = srt,
      bundle = bundle,
      method = method,
      group.by = group.by,
      reduction = reduction,
      value = value,
      palette = palette,
      palcolor = palcolor,
      title = title,
      subtitle = subtitle,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ),
    spatial = cnv_plot_spatial(
      srt = srt,
      bundle = bundle,
      method = method,
      group.by = group.by,
      value = value,
      palette = palette,
      palcolor = palcolor,
      title = title,
      subtitle = subtitle,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ),
    bar = cnv_plot_bar(
      srt = srt,
      bundle = bundle,
      method = method,
      group.by = group.by,
      value = value,
      palette = palette,
      palcolor = palcolor,
      title = title,
      subtitle = subtitle,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    ),
    tree = cnv_plot_tree(
      bundle = bundle,
      by = cluster_tree_by,
      ...
    )
  )
}

cnv_get_method_bundle <- function(srt, tool_name, method = NULL) {
  tool <- srt@tools[[tool_name]]
  if (is.null(tool) || is.null(tool$methods)) {
    log_message(
      "Cannot find CNV results in {.code srt@tools[[{tool_name}]]}",
      message_type = "error"
    )
  }
  method <- method %||% tool$active_method
  if (is.null(method) || !method %in% names(tool$methods)) {
    log_message(
      "{.arg method} must be one of stored CNV methods: {.val {names(tool$methods)}}",
      message_type = "error"
    )
  }
  tool$methods[[method]]
}

cnv_plot_heatmap <- function(
  srt,
  bundle,
  method,
  group.by = NULL,
  max_cells = 2000,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  palette = "Chinese",
  palcolor = NULL,
  title = NULL,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  heatmap_limits = NULL,
  ...
) {
  check_r("ComplexHeatmap", verbose = FALSE)

  check_r("circlize", verbose = FALSE)
  mat <- cnv_require_matrix(bundle)
  cells <- intersect(colnames(mat), colnames(srt))
  if (length(cells) == 0L) {
    log_message(
      "CNV matrix contains no cells present in {.arg srt}",
      message_type = "error"
    )
  }
  if (!is.null(max_cells) && length(cells) > max_cells) {
    set.seed(11)
    cells <- sample(cells, max_cells)
  }
  mat <- mat[, cells, drop = FALSE]
  bin_info <- bundle$bin_info
  if (identical(method, "fastCNV") &&
      (is.null(bin_info) || nrow(bin_info) != nrow(mat) ||
        !any(!is.na(bin_info$chr) & nzchar(bin_info$chr)))) {
    windows <- rownames(mat)
    if (all(grepl("^(?:chr)?(?:[1-9]|1[0-9]|2[0-3]|X|Y)[.][pq][0-9]+$", windows,
                  ignore.case = TRUE, perl = TRUE))) {
      bin_info <- data.frame(chr = sub("[.].*$", "", windows), start = seq_along(windows))
    }
  }
  if (!is.null(bin_info) && nrow(bin_info) == nrow(mat)) {
    row_order <- cnv_bin_order(bin_info)
    mat <- mat[row_order, , drop = FALSE]
    bin_info <- bin_info[row_order, , drop = FALSE]
  }
  if (is.null(bin_info) || nrow(bin_info) != nrow(mat)) {
    row_split <- NULL
  } else {
    row_split <- bin_info$chr
    row_split[is.na(row_split) | !nzchar(as.character(row_split))] <- "unknown"
    row_split <- factor(
      as.character(row_split),
      levels = unique(as.character(row_split))
    )
  }

  group.by <- group.by %||% cnv_default_group_col(srt, method = method)
  group_values <- group_colors <- NULL
  if (!is.null(group.by)) {
    cnv_validate_metadata(srt, group.by, "group.by")
    group_values <- srt@meta.data[cells, group.by, drop = TRUE]
    if (!is.factor(group_values)) {
      group_values <- factor(
        as.character(group_values),
        levels = unique(as.character(group_values))
      )
    }
    group_order <- order(group_values, na.last = TRUE)
    mat <- mat[, group_order, drop = FALSE]
    group_values <- group_values[group_order]
    group_colors <- list(group = cnv_discrete_colors(group_values, palette = palette, palcolor = palcolor))
  }

  is_infercnv <- identical(method, "infercnv")
  neutral <- if (is_infercnv) 1 else 0
  default_limits <- if (is_infercnv) c(0.75, 1.25) else c(-1, 1)
  limits <- cnv_heatmap_limits(mat, heatmap_limits, default_limits, neutral)
  colors <- palette_colors(
    n = 101L,
    palette = heatmap_palette,
    palcolor = heatmap_palcolor
  )
  legend_at <- c(limits[1L], neutral, limits[2L])
  n_side <- (length(colors) + 1L) %/% 2L
  breaks <- c(
    seq(limits[1L], neutral, length.out = n_side),
    seq(neutral, limits[2L], length.out = n_side)[-1L]
  )
  ticks <- format(legend_at, digits = 4L, trim = TRUE)
  if (anyDuplicated(ticks)) {
    ticks <- format(legend_at, digits = 17L, trim = TRUE)
  }

  chr <- if (!is.null(row_split) && any(row_split != "unknown")) row_split
  Heatmap <- get_namespace_fun("ComplexHeatmap", "Heatmap")
  HeatmapAnnotation <- get_namespace_fun("ComplexHeatmap", "HeatmapAnnotation")
  legend_style <- list(
    legend_gp = heatmap_border_gp(TRUE, "black"),
    border = heatmap_legend_border(TRUE, "black")
  )
  args <- list(
    matrix = t(mat),
    name = paste0(method, " CNV"),
    col = circlize::colorRamp2(breaks, colors),
    column_title = if (!identical(title, "")) title %||% paste0(method, " CNV"),
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_title = rep("", max(1L, length(unique(group_values)))),
    row_split = group_values,
    column_split = chr,
    row_gap = grid::unit(0.5, "mm"),
    column_gap = grid::unit(0.12, "mm"),
    cluster_rows = cluster_columns,
    cluster_columns = cluster_rows,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    show_row_names = show_column_names,
    show_column_names = show_row_names,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    use_raster = use_raster,
    border = FALSE,
    left_annotation = if (!is.null(group_values)) HeatmapAnnotation(
      group = group_values,
      col = group_colors,
      which = "row",
      show_annotation_name = FALSE,
      simple_anno_size = grid::unit(1.5, "mm"),
      annotation_legend_param = list(group = legend_style)
    ),
    top_annotation = if (!is.null(chr)) HeatmapAnnotation(
      chromosome = cnv_chromosome_strip(chr, method),
      show_annotation_name = FALSE
    ),
    heatmap_legend_param = c(
      list(
        title = switch(method,
          infercnv = "Relative expression",
          fastCNV = "Window score",
          copykat = "CNA score",
          scevan = "Relative CNA",
          numbat = "Signed CNV evidence",
          "CNV"
        ),
        at = legend_at,
        labels = ticks
      ),
      legend_style
    )
  )
  dots <- list(...)
  if (is.list(dots$heatmap_legend_param)) {
    args$heatmap_legend_param <- utils::modifyList(args$heatmap_legend_param, dots$heatmap_legend_param)
    if (!is.null(dots$heatmap_legend_param$at) && is.null(dots$heatmap_legend_param$labels)) {
      args$heatmap_legend_param$labels <- as.character(dots$heatmap_legend_param$at)
    }
    dots$heatmap_legend_param <- NULL
  }
  do.call(Heatmap, c(args[setdiff(names(args), names(dots))], dots))
}

cnv_heatmap_limits <- function(mat, limits, defaults, neutral) {
  if (is.null(limits)) {
    return(defaults)
  }
  if (identical(limits, "auto")) {
    deviations <- abs(mat[is.finite(mat)] - neutral)
    if (!length(deviations)) {
      return(defaults)
    }
    span <- stats::quantile(deviations, 0.95, names = FALSE)
    limits <- neutral + c(-span, span)
    return(if (limits[1L] < neutral && limits[2L] > neutral) limits else defaults)
  }
  if (!is.numeric(limits) || is.complex(limits) || length(limits) != 2L ||
      any(!is.finite(limits)) || limits[1L] >= neutral || limits[2L] <= neutral) {
    log_message(
      "{.arg heatmap_limits} must be {.val auto} or two finite values that strictly span the neutral value ({neutral})",
      message_type = "error"
    )
  }
  limits
}

cnv_chromosome_strip <- function(chromosomes, method) {
  labels <- sub("^chr", "", levels(chromosomes), ignore.case = TRUE)
  if (method %in% c("fastCNV", "copykat")) {
    labels[labels == "23"] <- "X"
  }
  n <- length(labels)
  if (n > 12L) {
    keep <- unique(c(seq.int(1L, n, by = 3L), n))
    if (length(keep) > 1L && diff(tail(keep, 2L)) == 1L) {
      keep <- keep[-(length(keep) - 1L)]
    }
    labels[-keep] <- ""
  }
  widths <- as.numeric(table(chromosomes)) / length(chromosomes)
  labels[widths < ifelse(nchar(labels) > 1L, 0.032, 0.020)] <- ""
  names(labels) <- levels(chromosomes)
  fills <- stats::setNames(rep_len(c("#D5DCE3", "#ECEFF2"), n), names(labels))
  anno_block <- get_namespace_fun("ComplexHeatmap", "anno_block")
  callback <- if ("panel_fun" %in% names(formals(anno_block))) "panel_fun" else "graphics"
  args <- list(height = grid::unit(2.6, "mm"))
  args[[callback]] <- function(index, levels) {
    chromosome <- as.character(chromosomes[index[1L]])
    grid::grid.rect(gp = grid::gpar(fill = fills[[chromosome]], col = NA))
    grid::grid.text(labels[[chromosome]], gp = grid::gpar(fontsize = 6, col = "#4B5869"))
  }
  do.call(anno_block, args)
}

cnv_plot_dim <- function(
  srt,
  bundle,
  method,
  group.by = NULL,
  reduction = NULL,
  value = NULL,
  palette = "Chinese",
  palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  value <- value %||% cnv_default_value_col(srt, method = method)
  cnv_validate_metadata(srt, value, "value")
  value_data <- srt@meta.data[[value]]
  if (is.numeric(value_data) || is.integer(value_data)) {
    FeatureDimPlot(
      object = srt,
      features = value,
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      title = title %||% value,
      subtitle = subtitle,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    )
  } else {
    palcolor <- cnv_discrete_colors(
      value_data,
      palette = palette,
      palcolor = palcolor
    )
    CellDimPlot(
      object = srt,
      group.by = value,
      reduction = reduction,
      palette = palette,
      palcolor = palcolor,
      title = title %||% value,
      subtitle = subtitle,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    )
  }
}

cnv_plot_spatial <- function(
  srt,
  bundle,
  method,
  group.by = NULL,
  value = NULL,
  palette = "Chinese",
  palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  value <- value %||% group.by %||% cnv_default_value_col(srt, method = method)
  cnv_validate_metadata(srt, value, "value")
  value_data <- srt@meta.data[[value]]
  if (!is.numeric(value_data) && !is.integer(value_data)) {
    palcolor <- cnv_discrete_colors(
      value_data,
      palette = palette,
      palcolor = palcolor
    )
  }
  p <- SpatialSpotPlot(
    object = srt,
    group.by = value,
    palette = palette,
    palcolor = palcolor,
    theme_use = theme_use,
    theme_args = theme_args,
    ...
  )
  if (inherits(p, "ggplot")) {
    p <- p + ggplot2::labs(title = title %||% value, subtitle = subtitle)
  }
  p
}

cnv_plot_bar <- function(
  srt,
  bundle,
  method,
  group.by = NULL,
  value = NULL,
  palette = "Chinese",
  palcolor = NULL,
  position = c("fill", "stack", "dodge"),
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  position <- match.arg(position)
  value <- value %||%
    cnv_col(srt, method, "prediction") %||%
    cnv_col(srt, method, "cluster")
  if (is.null(value) || !value %in% colnames(srt@meta.data)) {
    log_message(
      "No CNV prediction or cluster metadata column is available for bar plotting",
      message_type = "error"
    )
  }
  group.by <- group.by %||% "All.cells"
  df <- srt@meta.data
  if (identical(group.by, "All.cells")) {
    df[[group.by]] <- "All.cells"
  } else {
    cnv_validate_metadata(srt, group.by, "group.by")
  }
  df <- df[, c(group.by, value), drop = FALSE]
  colnames(df) <- c("group", "value")
  df$value <- as.character(df$value)
  df$value[is.na(df$value) | !nzchar(df$value)] <- "NA"
  df$group <- as.character(df$group)
  value_levels <- cnv_order_cnv_levels(unique(df$value))
  group_levels <- unique(df$group)
  df$value <- factor(df$value, levels = value_levels)
  df$group <- factor(df$group, levels = group_levels)
  summary_df <- as.data.frame(
    table(group = df$group, value = df$value),
    stringsAsFactors = FALSE
  )
  summary_df$n <- as.numeric(summary_df$Freq)
  summary_df$Freq <- NULL
  totals <- tapply(summary_df$n, summary_df$group, sum)
  summary_df$fraction <- summary_df$n / totals[as.character(summary_df$group)]
  summary_df$y <- if (identical(position, "fill")) {
    summary_df$fraction
  } else {
    summary_df$n
  }
  colors <- cnv_discrete_colors(
    summary_df$value,
    palette = palette,
    palcolor = palcolor
  )
  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = .data$group, y = .data$y, fill = .data$value)
  ) +
    ggplot2::geom_col(
      position = if (identical(position, "fill")) "stack" else position,
      width = 0.72,
      ...
    ) +
    ggplot2::scale_fill_manual(values = colors, name = value) +
    ggplot2::labs(
      x = group.by,
      y = if (identical(position, "fill")) "Fraction" else "Cells",
      title = title %||% "CNV composition",
      subtitle = subtitle
    ) +
    ggplot2::scale_y_continuous(
      labels = if (identical(position, "fill")) {
        function(x) paste0(round(x * 100), "%")
      } else {
        ggplot2::waiver()
      },
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    apply_plot_theme(
      theme_use = theme_use,
      theme_args = theme_args,
      fallback = ggplot2::theme
    )
}

cnv_plot_tree <- function(
  bundle,
  by = c("cluster", "cell"),
  distance = c("euclidean", "correlation"),
  hclust_method = "ward.D2",
  ...
) {
  by <- match.arg(by)
  distance <- match.arg(distance)
  mat <- cnv_require_matrix(bundle)
  if (ncol(mat) < 2L) {
    log_message(
      "CNV tree requires at least two cells",
      message_type = "error"
    )
  }
  profile <- t(mat)
  if (identical(by, "cluster")) {
    cell_info <- bundle$cell_info
    cluster <- cell_info$cluster
    names(cluster) <- cell_info$cell
    cluster <- cluster[rownames(profile)]
    if (all(is.na(cluster) | !nzchar(as.character(cluster)))) {
      log_message(
        "CNV cluster information is unavailable; use {.arg cluster_tree_by = 'cell'}",
        message_type = "error"
      )
    }
    profile <- stats::aggregate(
      profile,
      by = list(cluster = as.character(cluster)),
      FUN = mean,
      na.rm = TRUE
    )
    rownames(profile) <- profile$cluster
    profile$cluster <- NULL
    profile <- as.matrix(profile)
  }
  profile[!is.finite(profile)] <- 0
  dist_obj <- if (identical(distance, "correlation")) {
    stats::as.dist(1 - stats::cor(t(profile), use = "pairwise.complete.obs"))
  } else {
    stats::dist(profile, method = "euclidean")
  }
  stats::hclust(dist_obj, method = hclust_method, ...)
}

cnv_require_matrix <- function(bundle) {
  mat <- bundle$cnv_matrix
  if (is.null(mat)) {
    log_message(
      "CNV matrix was not stored for this result. Re-run {.fn RunCNV} with {.arg store_matrix = TRUE}",
      message_type = "error"
    )
  }
  mat <- as.matrix(mat)
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    log_message(
      "Stored CNV matrix is empty",
      message_type = "error"
    )
  }
  mat
}

cnv_default_group_col <- function(srt, method) {
  for (candidate in c(
    cnv_col(srt, method, "prediction"),
    cnv_col(srt, method, "cluster")
  )) {
    if (is.null(candidate)) {
      next
    }
    values <- srt@meta.data[[candidate]]
    if (any(!is.na(values) & nzchar(as.character(values)))) {
      return(candidate)
    }
  }
  NULL
}

cnv_default_value_col <- function(srt, method) {
  pred <- cnv_col(srt, method, "prediction")
  if (!is.null(pred)) {
    vals <- srt@meta.data[[pred]]
    if (any(!is.na(vals) & nzchar(as.character(vals)))) {
      return(pred)
    }
  }
  score <- cnv_col(srt, method, "score")
  if (!is.null(score)) {
    return(score)
  }
  log_message(
    "No CNV metadata columns were found for method {.val {method}}",
    message_type = "error"
  )
}

cnv_col <- function(srt, method, suffix) {
  candidates <- c(paste0("CNV_", method, "_", suffix), paste0("CNV_", suffix))
  hit <- candidates[candidates %in% colnames(srt@meta.data)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

cnv_bin_order <- function(bin_info) {
  chr <- as.character(bin_info$chr %||% rep(NA_character_, nrow(bin_info)))
  chr_clean <- sub("^chr", "", chr, ignore.case = TRUE)
  chr_rank <- match(
    chr_clean,
    c(as.character(seq_len(22)), "X", "Y", "M", "MT")
  )
  missing_chr <- is.na(chr_rank)
  if (any(missing_chr)) {
    base_rank <- if (any(!missing_chr)) max(chr_rank[!missing_chr]) else 1000L
    chr_rank[missing_chr] <- base_rank + seq_len(sum(missing_chr))
  }
  chr_rank[!is.finite(chr_rank)] <- seq_len(length(chr_rank))[
    !is.finite(chr_rank)
  ]
  start <- suppressWarnings(as.numeric(bin_info$start))
  start[!is.finite(start)] <- seq_len(length(start))[!is.finite(start)]
  order(chr_rank, start, seq_along(chr_rank), na.last = TRUE)
}

cnv_discrete_colors <- function(values, palette = "Chinese", palcolor = NULL) {
  values <- cnv_order_cnv_levels(unique(as.character(values)))
  values <- values[!is.na(values)]
  fallback <- palette_colors(values, palette = palette, palcolor = palcolor)
  defaults <- c(
    aneuploid = "#C44E52",
    tumor = "#C44E52",
    malignant = "#C44E52",
    altered = "#C44E52",
    gain = "#E76F51",
    amplification = "#E76F51",
    amp = "#E76F51",
    diploid = "#4C72B0",
    normal = "#4C72B0",
    neutral = "#BDBDBD",
    loss = "#3B6FB6",
    deletion = "#3B6FB6",
    del = "#3B6FB6",
    `NA` = "#BDBDBD"
  )
  out <- fallback
  keys <- tolower(values)
  matched <- keys %in% names(defaults)
  out[matched] <- defaults[keys[matched]]
  stats::setNames(unname(out), values)
}

cnv_order_cnv_levels <- function(values) {
  values <- as.character(values)
  priority <- c(
    "aneuploid",
    "tumor",
    "malignant",
    "altered",
    "gain",
    "amplification",
    "amp",
    "loss",
    "deletion",
    "del",
    "diploid",
    "normal",
    "neutral",
    "NA"
  )
  keys <- tolower(values)
  rank <- match(keys, tolower(priority))
  rank[is.na(rank)] <- length(priority) + seq_len(sum(is.na(rank)))
  values[order(rank, values)]
}
