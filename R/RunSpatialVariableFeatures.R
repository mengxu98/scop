#' @title Run spatial variable feature detection
#'
#' @description
#' Score genes by spot-level spatial autocorrelation. The native `"moran"` and
#' `"geary"` methods use a lightweight coordinate KNN graph. `"SPARKX"` and
#' `"nnSVG"` use optional external backends when their packages are installed.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used for expression values.
#' @param features Features to score. If `NULL`, current variable features are
#' used; if no variable features are present, all assay features are used.
#' @param method Spatial variable feature detection method.
#' @param coordinate_space Coordinate system used for distance-sensitive
#' analysis. The default is raw, unscaled acquisition coordinates. Use
#' `"legacy_display"` explicitly to reproduce the display-scaled coordinates
#' used before scop 0.9.0. Distance thresholds and weights use the selected
#' coordinate units; `k` is a unitless neighbor count.
#' @param k Number of nearest spatial neighbors per spot.
#' @param nfeatures Number of top spatial features stored in
#' `srt@misc[["SpatialVariableFeatures"]]`.
#' @param min_spots Minimum number of spots with non-zero expression required
#' for a feature to be tested.
#' @param nperm Number of label permutations used for empirical p values. The
#' default `0` skips p-value calculation.
#' @param set_variable_features Whether to set the top spatial features as
#' variable features for `assay`.
#' @param store_results Whether to store the full result in `srt@tools`.
#' @param seed Random seed used for permutation tests.
#' @param ... Additional arguments passed to external backends.
#'
#' @return A `Seurat` object with spatial variable feature results stored in
#' `srt@tools[["SpatialVariableFeatures"]]` and top feature names stored in
#' `srt@misc[["SpatialVariableFeatures"]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- Seurat::NormalizeData(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   verbose = FALSE
#' )
#' spatial <- Seurat::FindVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 100,
#'   verbose = FALSE
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   features = Seurat::VariableFeatures(spatial, assay = "Spatial")[1:2]
#' )
#'
#' spatial <- RunSpatialVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 50
#' )
#' SpatialVariableFeaturePlot(spatial, plot_type = "combined", nfeatures = 2)
RunSpatialVariableFeatures <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  method = c("moran", "geary", "SPARKX", "nnSVG"),
  image = NULL,
  coord.cols = c("x", "y"),
  k = 6,
  nfeatures = 2000,
  min_spots = 5,
  nperm = 0,
  set_variable_features = TRUE,
  store_results = TRUE,
  verbose = TRUE,
  seed = 11,
  coordinate_space = c("raw", "legacy_display"),
  ...
) {
  coordinate_space <- match.arg(coordinate_space)
  log_message(
    "Running spatial variable feature detection",
    message_type = "running",
    verbose = verbose
  )
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    log_message(
      "{.arg k} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.numeric(nfeatures) || length(nfeatures) != 1L || is.na(nfeatures) || nfeatures < 1) {
    log_message(
      "{.arg nfeatures} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.numeric(min_spots) || length(min_spots) != 1L || is.na(min_spots) || min_spots < 1) {
    log_message(
      "{.arg min_spots} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.numeric(nperm) || length(nperm) != 1L || is.na(nperm) || nperm < 0) {
    log_message(
      "{.arg nperm} must be a non-negative number",
      message_type = "error"
    )
  }
  k <- as.integer(k)
  nfeatures <- as.integer(nfeatures)
  min_spots <- as.integer(min_spots)
  nperm <- as.integer(nperm)
  set.seed(seed)

  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  spots <- intersect(colnames(srt), rownames(coords))
  if (length(spots) == 0L) {
    log_message(
      "No spatial coordinates match spots in {.arg srt}",
      message_type = "error"
    )
  }
  coords <- coords[spots, , drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  coords <- coords[keep_coords, , drop = FALSE]
  spots <- rownames(coords)
  if (length(spots) < 3L) {
    log_message(
      "At least three spots with finite coordinates are required",
      message_type = "error"
    )
  }
  if (k >= length(spots)) {
    k <- length(spots) - 1L
  }

  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
    if (length(features) == 0L) {
      features <- rownames(expr)
    }
  }
  features <- unique(features)
  features <- intersect(features, rownames(expr))
  if (length(features) == 0L) {
    log_message(
      "No requested {.arg features} are present in assay {.val {assay}}",
      message_type = "error"
    )
  }

  expr <- expr[features, spots, drop = FALSE]
  expressed_spots <- Matrix::rowSums(expr > 0)
  keep_features <- expressed_spots >= min_spots
  if (!any(keep_features)) {
    log_message(
      "No features are expressed in at least {.val {min_spots}} spots",
      message_type = "error"
    )
  }
  expr <- expr[keep_features, , drop = FALSE]
  expressed_spots <- expressed_spots[keep_features]

  expr_mat <- as.matrix(expr)
  edges <- NULL
  if (method %in% c("moran", "geary")) {
    edges <- spatial_variable_knn_edges(coords, k = k)
    result <- spatial_variable_run_knn(
      expr = expr_mat,
      edges = edges,
      method = method,
      nperm = nperm
    )
  } else if (identical(method, "SPARKX")) {
    result <- spatial_variable_run_sparkx(
      expr = expr_mat,
      coords = coords,
      ...
    )
  } else {
    result <- spatial_variable_run_nnsvg(
      expr = expr_mat,
      coords = coords,
      assay = assay,
      ...
    )
  }
  result <- spatial_variable_finalize_result(
    result = result,
    expr = expr_mat,
    expressed_spots = expressed_spots,
    method = method
  )
  top_features <- utils::head(result$feature[is.finite(result$score)], nfeatures)
  srt@misc[["SpatialVariableFeatures"]] <- top_features
  if (isTRUE(set_variable_features) && length(top_features) > 0L) {
    SeuratObject::VariableFeatures(srt, assay = assay) <- top_features
  }
  if (isTRUE(store_results)) {
    score_lookup <- stats::setNames(result$score, result$feature)
    srt@tools[["SpatialVariableFeatures"]] <- list(
      result = result,
      coords = coords,
      edges = edges,
      summary = list(
        n_features = nrow(result),
        top_features = scop_spatial_feature_summary(top_features, scores = score_lookup)
      ),
      parameters = list(
        assay = assay,
        layer = layer,
        method = method,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        k = k,
        nfeatures = nfeatures,
        min_spots = min_spots,
        nperm = nperm,
        seed = seed,
        set_variable_features = set_variable_features
      )
    )
    srt@tools[["SpatialVariableFeatures"]] <- spatial_result_build(
      bundle = srt@tools[["SpatialVariableFeatures"]],
      method = "SpatialVariableFeatures",
      result_type = "feature_pattern",
      provenance = list(
        producer = "RunSpatialVariableFeatures",
        backend_id = if (identical(method, "moran")) "core" else method
      )
    )
  }
  log_message(
    "Stored {.val {length(top_features)}} spatial variable features",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spatial_variable_run_knn <- function(expr, edges, method, nperm = 0) {
  scores <- spatial_variable_score_matrix(
    expr = expr,
    edges = edges,
    method = method,
    nperm = nperm
  )
  data.frame(
    feature = rownames(expr),
    statistic = scores$statistic,
    score = scores$score,
    p_value = scores$p_value,
    q_value = if (all(is.na(scores$p_value))) {
      NA_real_
    } else {
      stats::p.adjust(scores$p_value, method = "BH")
    },
    stringsAsFactors = FALSE
  )
}

spatial_variable_run_sparkx <- function(expr, coords, ...) {
  sparkx <- spatial_variable_get_fun("SPARK", "sparkx")
  out <- sparkx(
    count_in = expr,
    locus_in = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    ...
  )
  res <- if (is.data.frame(out)) {
    out
  } else if (is.list(out) && is.data.frame(out$res_mtest)) {
    out$res_mtest
  } else if (is.list(out) && is.data.frame(out$result)) {
    out$result
  } else {
    log_message(
      "{.pkg SPARK} returned an unsupported SPARK-X result format",
      message_type = "error"
    )
  }
  res <- as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  features <- spatial_variable_result_features(res, rownames(expr))
  p_value <- spatial_variable_pick_numeric(res, c("combinedPval", "combined_pvalue", "p_value", "pvalue", "pval"))
  q_value <- spatial_variable_pick_numeric(res, c("adjustedPval", "adjusted_pvalue", "q_value", "qvalue", "padj", "fdr"))
  statistic <- spatial_variable_pick_numeric(res, c("statistic", "score", "combinedStat", "combined_stat"))
  data.frame(
    feature = features,
    statistic = statistic,
    p_value = p_value,
    q_value = q_value,
    stringsAsFactors = FALSE
  )
}

spatial_variable_run_nnsvg <- function(expr, coords, assay, ...) {
  spatial_variable_require_package("SpatialExperiment")
  spatial_variable_require_package("SummarizedExperiment")
  spatial_variable_require_package("S4Vectors")
  spatial_variable_require_package("nnSVG")
  spe <- spatial_variable_make_spe(expr = expr, coords = coords, assay = assay)
  extra_args <- list(...)
  args <- c(list(spe), extra_args)
  if (!"assay_name" %in% names(extra_args)) {
    args$assay_name <- "counts"
  }
  out <- do.call(spatial_variable_get_fun("nnSVG", "nnSVG"), args)
  res <- spatial_variable_row_data(out)
  features <- spatial_variable_result_features(res, rownames(expr))
  statistic <- spatial_variable_pick_numeric(res, c("LR_stat", "LR.stat", "statistic", "score", "prop_sv", "prop.sv"))
  p_value <- spatial_variable_pick_numeric(res, c("pval", "p_value", "p.value"))
  q_value <- spatial_variable_pick_numeric(res, c("padj", "q_value", "q.value", "fdr"))
  rank <- spatial_variable_pick_numeric(res, c("rank", "Rank"))
  data.frame(
    feature = features,
    statistic = statistic,
    p_value = p_value,
    q_value = q_value,
    rank = rank,
    stringsAsFactors = FALSE
  )
}

spatial_variable_make_spe <- function(expr, coords, assay = NULL) {
  spatial_variable_require_package("SpatialExperiment")
  SpatialExperiment <- get_namespace_fun("SpatialExperiment", "SpatialExperiment")
  spe <- SpatialExperiment(
    assays = list(counts = expr),
    spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE])
  )
  if (!is.null(assay)) {
    SummarizedExperiment::rowData(spe)$feature <- rownames(expr)
  }
  spe
}

spatial_variable_row_data <- function(x) {
  as.data.frame(SummarizedExperiment::rowData(x), stringsAsFactors = FALSE, check.names = FALSE)
}

spatial_variable_finalize_result <- function(result, expr, expressed_spots, method) {
  result <- as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"feature" %in% colnames(result)) {
    log_message(
      "Spatial variable feature result is missing a {.field feature} column",
      message_type = "error"
    )
  }
  result$feature <- as.character(result$feature)
  result <- result[result$feature %in% rownames(expr), , drop = FALSE]
  result <- result[!duplicated(result$feature), , drop = FALSE]
  if (nrow(result) == 0L) {
    log_message(
      "Spatial variable feature backend returned no tested features",
      message_type = "error"
    )
  }
  for (col in c("statistic", "score", "p_value", "q_value", "rank")) {
    if (!col %in% colnames(result)) {
      result[[col]] <- NA_real_
    }
    result[[col]] <- suppressWarnings(as.numeric(result[[col]]))
  }
  missing_q <- is.na(result$q_value) & is.finite(result$p_value)
  if (any(missing_q)) {
    result$q_value[missing_q] <- stats::p.adjust(result$p_value[missing_q], method = "BH")
  }
  if (all(!is.finite(result$score))) {
    result$score <- spatial_variable_score_from_significance(
      q_value = result$q_value,
      p_value = result$p_value,
      statistic = result$statistic,
      rank = result$rank
    )
  }
  result$mean <- rowMeans(expr[result$feature, , drop = FALSE])
  result$variance <- fast_row_vars(expr[result$feature, , drop = FALSE])
  result$n_spots <- as.integer(expressed_spots[result$feature])
  result$method <- method
  score_order <- ifelse(is.finite(result$score), result$score, -Inf)
  p_order <- ifelse(is.finite(result$p_value), result$p_value, Inf)
  q_order <- ifelse(is.finite(result$q_value), result$q_value, Inf)
  result <- result[order(-score_order, p_order, q_order), , drop = FALSE]
  result$rank <- seq_len(nrow(result))
  rownames(result) <- NULL
  spatial_variable_reorder_cols(
    result,
    c("feature", "rank", "method", "statistic", "score", "p_value", "q_value", "mean", "variance", "n_spots")
  )
}

spatial_variable_score_from_significance <- function(q_value, p_value, statistic, rank) {
  sig <- q_value
  sig[!is.finite(sig)] <- p_value[!is.finite(sig)]
  if (any(is.finite(sig))) {
    positive <- sig[is.finite(sig) & sig > 0]
    floor_value <- if (length(positive) > 0L) min(positive) * 0.1 else .Machine$double.xmin
    sig[is.finite(sig) & sig <= 0] <- floor_value
    return(-log10(sig))
  }
  if (any(is.finite(rank))) {
    return(-rank)
  }
  statistic
}

spatial_variable_result_features <- function(df, fallback) {
  feature_col <- intersect(c("feature", "features", "gene", "genes", "gene_id", "geneid"), colnames(df))
  if (length(feature_col) > 0L) {
    return(as.character(df[[feature_col[[1L]]]]))
  }
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df) && !all(grepl("^[0-9]+$", rn))) {
    return(as.character(rn))
  }
  utils::head(fallback, nrow(df))
}

spatial_variable_pick_numeric <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[hit[[1L]]]]))
}

spatial_variable_reorder_cols <- function(df, first_cols) {
  first_cols <- intersect(first_cols, colnames(df))
  df[, c(first_cols, setdiff(colnames(df), first_cols)), drop = FALSE]
}

spatial_variable_require_package <- function(pkg) {
  repo <- if (identical(pkg, "SPARK")) "xzhoulab/SPARK" else pkg
  status <- tryCatch(check_r(repo, verbose = FALSE), error = function(e) FALSE)
  if (!isTRUE(unname(unlist(status))[1])) {
    log_message(
      "Please install required package before running this function: {.val {pkg}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spatial_variable_get_fun <- function(pkg, fun) {
  spatial_variable_require_package(pkg)
  tryCatch(
    get_namespace_fun(pkg, fun),
    error = function(e) {
      log_message(
        "{.pkg {pkg}} does not export required function {.fn {fun}}",
        message_type = "error"
      )
    }
  )
}

#' @title Plot spatial variable feature results
#'
#' @description
#' Visualize normalized results produced by [RunSpatialVariableFeatures()]. The
#' summary view shows feature ranks and, when finite p- or q-values are stored,
#' significance. Results without permutation statistics are shown without a
#' significance size mapping or legend. The surface view reuses
#' [SpatialSpotPlot()] to draw spatial expression for selected features.
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param plot_type Plot type: `"summary"`, `"surface"`, or `"combined"`.
#' @param features Features to plot. If `NULL`, top features from the stored
#' spatial variable feature result are used.
#' @param nfeatures Number of top features used when `features = NULL`.
#' @param score_col Result column used for the summary x-axis.
#'
#' @return A `ggplot` or `patchwork` object.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- Seurat::NormalizeData(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   verbose = FALSE
#' )
#' spatial <- RunSpatialVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 10,
#'   verbose = FALSE
#' )
#' SpatialVariableFeaturePlot(spatial, plot_type = "summary")
SpatialVariableFeaturePlot <- function(
  srt,
  plot_type = c("summary", "surface", "combined"),
  features = NULL,
  nfeatures = 10,
  score_col = "score",
  assay = NULL,
  layer = NULL,
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  pt.size = NULL,
  pt.alpha = 0.9,
  stroke = 0.1,
  palette = "Spectral",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  stored <- spatial_variable_get_stored_result(srt)
  result <- stored$result
  features <- spatial_variable_plot_features(result, features = features, nfeatures = nfeatures)
  if (identical(plot_type, "summary")) {
    return(spatial_variable_summary_plot(
      result = result,
      features = features,
      score_col = score_col,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
  if (is.null(layer)) {
    layer <- stored$parameters$layer %||% "data"
  }
  if (identical(plot_type, "surface")) {
    return(spatial_variable_surface_plot(
      srt = srt,
      features = features,
      assay = assay,
      layer = layer,
      image = image,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      coord.cols = coord.cols,
      flip.y = flip.y,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    ))
  }
  spatial_variable_require_package("patchwork")
  summary_plot <- spatial_variable_summary_plot(
    result = result,
    features = features,
    score_col = score_col,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args
  )
  surface_plot <- spatial_variable_surface_plot(
    srt = srt,
    features = features,
    assay = assay,
    layer = layer,
    image = image,
    overlay_image = overlay_image,
    image.alpha = image.alpha,
    coord.cols = coord.cols,
    flip.y = flip.y,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    stroke = stroke,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args,
    combine = TRUE,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow
  )
  wrap_plots <- spatial_variable_get_fun("patchwork", "wrap_plots")
  wrap_plots(summary_plot, surface_plot, ncol = 1)
}

spatial_variable_get_stored_result <- function(srt) {
  stored <- srt@tools[["SpatialVariableFeatures"]]
  if (is.null(stored) || !is.list(stored) || is.null(stored$result)) {
    log_message(
      "No spatial variable feature result is stored in {.code srt@tools[['SpatialVariableFeatures']]}",
      message_type = "error"
    )
  }
  if (!is.data.frame(stored$result) || nrow(stored$result) == 0L) {
    log_message(
      "Stored spatial variable feature result is empty",
      message_type = "error"
    )
  }
  stored
}

spatial_variable_plot_features <- function(result, features = NULL, nfeatures = 10) {
  if (is.null(features)) {
    features <- utils::head(result$feature, nfeatures)
  }
  features <- unique(as.character(features))
  features <- features[!is.na(features) & nzchar(features)]
  if (length(features) == 0L) {
    log_message("No {.arg features} are available for plotting", message_type = "error")
  }
  missing <- setdiff(features, result$feature)
  if (length(missing) > 0L) {
    log_message(
      "Requested feature(s) are not present in stored spatial variable feature results: {.val {missing}}",
      message_type = "error"
    )
  }
  features
}

spatial_variable_summary_plot <- function(
  result,
  features,
  score_col,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args
) {
  if (!score_col %in% colnames(result)) {
    log_message(
      "{.arg score_col} {.val {score_col}} is not present in the stored result",
      message_type = "error"
    )
  }
  df <- result[result$feature %in% features, , drop = FALSE]
  df$feature <- factor(df$feature, levels = rev(features))
  df[[".score"]] <- suppressWarnings(as.numeric(df[[score_col]]))
  has_significance <- ("q_value" %in% colnames(df) && any(is.finite(df$q_value))) ||
    ("p_value" %in% colnames(df) && any(is.finite(df$p_value)))
  df[[".significance"]] <- if ("q_value" %in% colnames(df) && any(is.finite(df$q_value))) {
    -log10(pmax(df$q_value, .Machine$double.xmin))
  } else if ("p_value" %in% colnames(df) && any(is.finite(df$p_value))) {
    -log10(pmax(df$p_value, .Machine$double.xmin))
  } else {
    NA_real_
  }
  cols <- palette_colors(as.character(features), palette = palette, palcolor = palcolor)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[".score"]], y = .data[["feature"]], color = .data[["feature"]])) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[[".score"]], yend = .data[["feature"]]),
      linewidth = 0.35,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::labs(
      x = score_col,
      y = NULL,
      color = "Feature"
    ) +
    spatial_variable_plot_theme(theme_use = theme_use, theme_args = theme_args) +
    ggplot2::theme(legend.position = legend.position)
  if (isTRUE(has_significance)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(size = .data[[".significance"]]),
        na.rm = TRUE
      ) +
      ggplot2::labs(size = "-log10(q/p)")
  } else {
    p <- p + ggplot2::geom_point(na.rm = TRUE)
  }
  p
}

spatial_variable_surface_plot <- function(
  srt,
  features,
  assay,
  layer,
  image,
  overlay_image,
  image.alpha,
  coord.cols,
  flip.y,
  pt.size,
  pt.alpha,
  stroke,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args,
  combine,
  nrow,
  ncol,
  byrow
) {
  if (is.null(theme_use)) {
    theme_use <- ggplot2::theme_minimal
  }
  SpatialSpotPlot(
    srt = srt,
    features = features,
    assay = assay,
    layer = layer,
    image = image,
    overlay_image = overlay_image,
    image.alpha = image.alpha,
    coord.cols = coord.cols,
    flip.y = flip.y,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    stroke = stroke,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow
  )
}

spatial_variable_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(ggplot2::theme_minimal())
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_fun <- if (is.character(theme_use)) {
    get(theme_use, mode = "function", inherits = TRUE)
  } else {
    theme_use
  }
  do.call(theme_fun, theme_args)
}

spatial_variable_knn_edges <- function(coords, k = 6) {
  k <- min(as.integer(k), nrow(coords) - 1L)
  graph <- spatial_graph_compute(
    coords = coords,
    method = "knn",
    k = k,
    directed = TRUE,
    weight = "binary"
  )
  graph$edges[, c("from", "to"), drop = FALSE]
}

spatial_variable_score_matrix <- function(
  expr,
  edges,
  method = c("moran", "geary"),
  nperm = 0
) {
  method <- match.arg(method)
  statistic <- rep(NA_real_, nrow(expr))
  p_value <- rep(NA_real_, nrow(expr))
  for (i in seq_len(nrow(expr))) {
    x <- as.numeric(expr[i, ])
    statistic[[i]] <- spatial_variable_score_vector(
      x = x,
      edges = edges,
      method = method
    )
    if (nperm > 0L && is.finite(statistic[[i]])) {
      perm_scores <- replicate(
        nperm,
        spatial_variable_score_vector(
          x = sample(x),
          edges = edges,
          method = method
        )
      )
      perm_scores <- perm_scores[is.finite(perm_scores)]
      if (length(perm_scores) > 0L) {
        if (identical(method, "moran")) {
          p_value[[i]] <- (sum(perm_scores >= statistic[[i]]) + 1) /
            (length(perm_scores) + 1)
        } else {
          p_value[[i]] <- (sum(perm_scores <= statistic[[i]]) + 1) /
            (length(perm_scores) + 1)
        }
      }
    }
  }
  score <- if (identical(method, "moran")) {
    statistic
  } else {
    1 - statistic
  }
  list(
    statistic = statistic,
    score = score,
    p_value = p_value
  )
}

spatial_variable_score_vector <- function(x, edges, method = c("moran", "geary")) {
  method <- match.arg(method)
  keep <- is.finite(x)
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  idx <- seq_along(x)
  keep_edge <- keep[edges$from] & keep[edges$to]
  if (!any(keep_edge)) {
    return(NA_real_)
  }
  edges <- edges[keep_edge, , drop = FALSE]
  idx_keep <- idx[keep]
  x_keep <- x[keep]
  x_centered <- x - mean(x_keep)
  denom <- sum((x_keep - mean(x_keep))^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  n <- length(idx_keep)
  w <- nrow(edges)
  if (identical(method, "moran")) {
    n / w * sum(x_centered[edges$from] * x_centered[edges$to]) / denom
  } else {
    (n - 1) / (2 * w) *
      sum((x[edges$from] - x[edges$to])^2) / denom
  }
}
