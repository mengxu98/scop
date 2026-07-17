#' @title Run STdeconvolve reference-free spatial deconvolution
#'
#' @description
#' Estimate spot-level topic proportions from a spatial `Seurat` object using
#' the optional `STdeconvolve` package.
#'
#' @md
#' @inheritParams RunRCTD
#' @param layer Assay layer used as STdeconvolve input.
#' @param k Number of topics. If `NULL`, models are fit over `k_candidates` and
#' `STdeconvolve::optimalModel()` is used to choose a model.
#' @param k_candidates Candidate topic numbers used when `k = NULL`.
#' @param opt Optimal model selector passed to `STdeconvolve::optimalModel()`.
#' @param clean_counts Whether to call `STdeconvolve::cleanCounts()`.
#' @param clean_counts_params Additional parameters passed to
#' `STdeconvolve::cleanCounts()`.
#' @param restrict_corpus Whether to call `STdeconvolve::restrictCorpus()`.
#' @param restrict_corpus_params Additional parameters passed to
#' `STdeconvolve::restrictCorpus()`.
#' @param fit_lda_params Additional parameters passed to
#' `STdeconvolve::fitLDA()`.
#' @param get_beta_theta_params Additional parameters passed to
#' `STdeconvolve::getBetaTheta()`.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param round_counts Whether to round non-integer counts before model fitting.
#'
#' @return A `Seurat` object with topic proportions in metadata and detailed
#' results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' topic_weights <- data.frame(
#'   STdeconvolve_prop_topic_1 = seq(0.75, 0.20, length.out = ncol(spatial)),
#'   STdeconvolve_prop_topic_2 = seq(0.20, 0.70, length.out = ncol(spatial)),
#'   STdeconvolve_prop_topic_3 = 0.10,
#'   row.names = colnames(spatial)
#' )
#' topic_weights <- topic_weights / rowSums(topic_weights)
#' spatial <- Seurat::AddMetaData(spatial, topic_weights)
#' spatial$STdeconvolve_dominant_type <- sub(
#'   "^STdeconvolve_prop_",
#'   "",
#'   colnames(topic_weights)[max.col(topic_weights)]
#' )
#' spatial$STdeconvolve_max_prop <- apply(topic_weights, 1, max)
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "STdeconvolve_dominant_type",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' spatial <- RunSTdeconvolve(
#'   spatial,
#'   assay = "Spatial",
#'   features = rownames(spatial)[1:300],
#'   k = 3,
#'   verbose = FALSE
#' )
#' STdeconvolvePlot(
#'   spatial,
#'   topics = 1:2,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'   STdeconvolvePlot(
#'     spatial,
#'     plot_type = "pie",
#'     overlay_image = FALSE,
#'     coord.cols = c("x", "y")
#'   )
RunSTdeconvolve <- function(
  srt,
  assay = NULL,
  layer = "counts",
  features = NULL,
  k = NULL,
  k_candidates = 2:9,
  opt = "min",
  clean_counts = TRUE,
  clean_counts_params = list(),
  restrict_corpus = TRUE,
  restrict_corpus_params = list(),
  fit_lda_params = list(),
  get_beta_theta_params = list(),
  prefix = "STdeconvolve",
  tool_name = "STdeconvolve",
  store_results = TRUE,
  round_counts = TRUE,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  stdeconvolve_assert_scalar_string(prefix, "prefix")
  stdeconvolve_assert_scalar_string(tool_name, "tool_name")
  stdeconvolve_validate_param_list(clean_counts_params, "clean_counts_params")
  stdeconvolve_validate_param_list(restrict_corpus_params, "restrict_corpus_params")
  stdeconvolve_validate_param_list(fit_lda_params, "fit_lda_params")
  stdeconvolve_validate_param_list(get_beta_theta_params, "get_beta_theta_params")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  features_use <- features %||% rownames(srt[[assay]])
  features_use <- intersect(features_use, rownames(srt[[assay]]))
  if (length(features_use) == 0L) {
    log_message(
      "No features are available for {.fn RunSTdeconvolve}",
      message_type = "error"
    )
  }

  counts <- rctd_get_count_matrix(
    srt,
    assay = assay,
    layer = layer,
    features = features_use,
    data_label = "Spatial",
    round_counts = round_counts,
    verbose = verbose
  )
  counts <- stdeconvolve_filter_counts(counts)
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    log_message(
      "No non-zero features or spots remain for {.fn RunSTdeconvolve}",
      message_type = "error"
    )
  }

  extra_fit_params <- list(...)
  if (length(extra_fit_params) > 0L) {
    stdeconvolve_validate_param_list(extra_fit_params, "...")
    fit_lda_params <- c(fit_lda_params, extra_fit_params)
  }
  topic_k <- stdeconvolve_resolve_k(k = k, k_candidates = k_candidates)

  log_message(
    "Run {.pkg STdeconvolve} with {.val {nrow(counts)}} features and {.val {ncol(counts)}} spatial spots",
    verbose = verbose
  )
  backend <- stdeconvolve_run_backend(
    counts = counts,
    k = topic_k,
    opt = opt,
    clean_counts = clean_counts,
    clean_counts_params = clean_counts_params,
    restrict_corpus = restrict_corpus,
    restrict_corpus_params = restrict_corpus_params,
    fit_lda_params = fit_lda_params,
    get_beta_theta_params = get_beta_theta_params
  )

  theta <- stdeconvolve_orient_theta(backend$theta, spot_ids = colnames(counts))
  colnames(theta) <- make.unique(make.names(colnames(theta)), sep = "_")
  weight_summary <- scop_spatial_finalize_weights(
    weights = theta,
    all_spots = colnames(srt)
  )
  # Retain all original spots in the stored result.  Spots removed after
  # corpus restriction have unknown topic weights and are represented by NA.
  theta <- weight_summary$full_weights
  srt <- scop_spatial_add_deconv_metadata(
    srt,
    weights = theta,
    prefix = prefix,
    metadata = weight_summary
  )

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      theta = theta,
      beta = backend$beta,
      corpus = backend$corpus,
      model = backend$model,
      models = backend$models,
      selected_k = backend$selected_k,
      features = rownames(counts),
      summary = scop_spatial_weight_summary(theta),
      parameters = list(
        assay = assay,
        layer = layer,
        k = k,
        k_candidates = k_candidates,
        opt = opt,
        clean_counts = clean_counts,
        clean_counts_params = clean_counts_params,
        restrict_corpus = restrict_corpus,
        restrict_corpus_params = restrict_corpus_params,
        fit_lda_params = fit_lda_params,
        get_beta_theta_params = get_beta_theta_params,
        prefix = prefix,
        tool_name = tool_name,
        round_counts = round_counts
      )
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "STdeconvolve",
      result_type = "deconvolution",
      provenance = list(producer = "RunSTdeconvolve", backend_id = "stdeconvolve")
    )
  }

  log_message(
    "{.pkg STdeconvolve} topic proportions stored in metadata columns with prefix {.val {prefix}_prop_}",
    verbose = verbose
  )
  srt
}

#' @title Plot STdeconvolve topic proportions
#'
#' @description
#' Plot schema-v1 topic proportions without rerunning the optional backend.
#' Multi-topic point maps use one shared proportion scale and reserve title
#' space above the automatic layout.
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param topics Topic names, topic numbers, or metadata columns to plot. If
#' `NULL`, all topics in the stored result are used.
#' @param tool_name Exact schema-v1 result key written by
#' `RunSTdeconvolve()`.
#' @param prefix Metadata prefix used by `RunSTdeconvolve()`. If `NULL`, the
#' prefix is read from the stored result.
#' @param combine Whether to combine point plots. If `FALSE`, a named list of
#' plots is returned.
#' @param nrow,ncol,byrow Point-plot layout controls. When both `nrow` and
#' `ncol` are `NULL`, a near-square layout with at most three columns is used.
#' @param ... Additional arguments passed to `SpatialSpotPlot()`.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#' @export
#'
#' @examples
#' check_r("JEFworks-Lab/STdeconvolve", verbose = FALSE)
#' data(visium_human_pancreas_sub)
#' spatial <- RunSTdeconvolve(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   features = rownames(visium_human_pancreas_sub)[1:300],
#'   k = 3,
#'   prefix = "STFull",
#'   tool_name = "STdeconvolveFull",
#'   verbose = FALSE
#' )
#' STdeconvolvePlot(
#'   spatial,
#'   tool_name = "STdeconvolveFull",
#'   topics = 1:2,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
STdeconvolvePlot <- function(
  srt,
  tool_name = "STdeconvolve",
  topics = NULL,
  prefix = NULL,
  plot_type = c("point", "pie"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  stdeconvolve_assert_scalar_string(tool_name, "tool_name")
  plot_type <- match.arg(plot_type)
  stored <- GetSpatialResult(srt, tool_name = tool_name)
  if (!identical(stored$method, "STdeconvolve")) {
    log_message(
      "Stored result {.val {tool_name}} was not produced by {.fn RunSTdeconvolve}",
      message_type = "error"
    )
  }
  prefix <- prefix %||% stored$parameters$prefix
  stdeconvolve_assert_scalar_string(prefix, "prefix")
  theta <- stdeconvolve_plot_theta(
    theta = stored$theta,
    spot_ids = colnames(srt),
    tool_name = tool_name
  )
  topic_names <- stdeconvolve_resolve_plot_topics(
    topics = topics,
    topic_names = colnames(theta),
    prefix = prefix
  )
  values <- theta[, topic_names, drop = FALSE]
  colnames(values) <- paste0(prefix, "_prop_", make.names(topic_names))
  if (identical(plot_type, "pie")) {
    return(SpatialSpotPlot(srt, values = values, plot_type = "pie", ...))
  }
  plots <- SpatialSpotPlot(
    srt,
    values = values,
    plot_type = "point",
    combine = FALSE,
    ...
  )
  if (isFALSE(combine)) {
    return(plots)
  }
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- min(3L, ceiling(sqrt(length(plots))))
  }
  finite_values <- values[is.finite(values)]
  common_limit <- if (length(finite_values) == 0L) 1 else max(finite_values)
  common_limit <- max(common_limit, .Machine$double.eps)
  legend_title <- list(...)$legend.title %||% "Proportion"
  plots <- Map(
    function(plot, topic) {
      plot <- stdeconvolve_set_point_scale(
        plot = plot,
        limits = c(0, common_limit),
        title = legend_title
      )
      plot +
        ggplot2::labs(title = topic) +
        ggplot2::theme(plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 4)))
    },
    plots,
    topic_names
  )
  patchwork::wrap_plots(
    plots,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(tool_name, " topic proportions"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 8))
      )
    )
}

stdeconvolve_set_point_scale <- function(plot, limits, title) {
  scale_index <- which(vapply(
    plot$scales$scales,
    function(scale) any(scale$aesthetics %in% c("colour", "color")),
    logical(1)
  ))
  if (length(scale_index) != 1L) {
    log_message(
      "Unable to identify the continuous topic-proportion color scale",
      message_type = "error"
    )
  }
  plot$scales$scales[[scale_index]]$limits <- limits
  plot$scales$scales[[scale_index]]$name <- title
  plot$labels$colour <- title
  plot
}

stdeconvolve_plot_theta <- function(theta, spot_ids, tool_name) {
  if (is.null(theta) || (!is.matrix(theta) && !inherits(theta, "Matrix"))) {
    log_message(
      "Stored result {.val {tool_name}} does not contain a matrix-like {.val theta}",
      message_type = "error"
    )
  }
  theta <- as.matrix(theta)
  if (nrow(theta) == 0L || ncol(theta) == 0L) {
    log_message("Stored {.val theta} is empty", message_type = "error")
  }
  if (
    is.null(rownames(theta)) || anyNA(rownames(theta)) ||
      any(!nzchar(rownames(theta))) || anyDuplicated(rownames(theta))
  ) {
    log_message(
      "Stored {.val theta} must have unique, non-missing spatial spot names",
      message_type = "error"
    )
  }
  if (
    is.null(colnames(theta)) || anyNA(colnames(theta)) ||
      any(!nzchar(colnames(theta))) || anyDuplicated(colnames(theta))
  ) {
    log_message(
      "Stored {.val theta} must have unique, non-missing topic names",
      message_type = "error"
    )
  }
  missing_spots <- setdiff(spot_ids, rownames(theta))
  extra_spots <- setdiff(rownames(theta), spot_ids)
  if (length(missing_spots) > 0L || length(extra_spots) > 0L) {
    log_message(
      "Stored {.val theta} spot identities are stale or incomplete: {.val {length(missing_spots)}} missing and {.val {length(extra_spots)}} unknown",
      message_type = "error"
    )
  }
  theta[spot_ids, , drop = FALSE]
}

stdeconvolve_resolve_plot_topics <- function(topics, topic_names, prefix) {
  if (is.null(topics)) {
    return(topic_names)
  }
  topics_chr <- as.character(topics)
  out <- vapply(topics_chr, function(topic) {
    if (topic %in% topic_names) {
      return(topic)
    }
    if (grepl("^[0-9]+$", topic)) {
      idx <- as.integer(topic)
      if (idx >= 1L && idx <= length(topic_names)) {
        return(topic_names[[idx]])
      }
    }
    prefix_value <- paste0(prefix, "_prop_", make.names(topic_names))
    hit <- match(topic, prefix_value)
    if (!is.na(hit)) {
      return(topic_names[[hit]])
    }
    NA_character_
  }, character(1))
  if (anyNA(out)) {
    log_message(
      "Requested {.arg topics} were not found: {.val {topics_chr[is.na(out)]}}",
      message_type = "error"
    )
  }
  unname(out)
}

stdeconvolve_run_backend <- function(
  counts,
  k,
  opt,
  clean_counts,
  clean_counts_params,
  restrict_corpus,
  restrict_corpus_params,
  fit_lda_params,
  get_beta_theta_params
) {
  check_r("JEFworks-Lab/STdeconvolve", verbose = FALSE)
  clean_counts_fun <- get_namespace_fun("STdeconvolve", "cleanCounts")
  restrict_fun <- get_namespace_fun("STdeconvolve", "restrictCorpus")
  fit_lda <- get_namespace_fun("STdeconvolve", "fitLDA")
  optimal_model <- get_namespace_fun("STdeconvolve", "optimalModel")
  get_beta_theta <- get_namespace_fun("STdeconvolve", "getBetaTheta")

  corpus <- counts
  if (isTRUE(clean_counts)) {
    corpus <- do.call(clean_counts_fun, c(list(counts = corpus), clean_counts_params))
  }
  if (isTRUE(restrict_corpus)) {
    corpus <- do.call(restrict_fun, c(list(counts = corpus), restrict_corpus_params))
  }
  corpus_mat <- stdeconvolve_extract_corpus(corpus)
  # `restrictCorpus()` can remove every retained gene from individual spots.
  # Drop those spots before fitLDA(), which requires a non-zero row per spot
  # after the matrix is transposed below.
  corpus_mat <- corpus_mat[, Matrix::colSums(corpus_mat) > 0, drop = FALSE]
  if (nrow(corpus_mat) == 0L || ncol(corpus_mat) == 0L) {
    log_message(
      "{.pkg STdeconvolve} corpus is empty after preprocessing",
      message_type = "error"
    )
  }

  lda_input <- t(as.matrix(corpus_mat))
  models <- do.call(fit_lda, c(list(counts = lda_input, Ks = k), fit_lda_params))
  model <- do.call(optimal_model, list(models = models, opt = opt))
  result <- do.call(get_beta_theta, c(list(lda = model), get_beta_theta_params))
  stdeconvolve_extract_result(
    result = result,
    model = model,
    models = models,
    corpus = corpus_mat,
    requested_k = k
  )
}

stdeconvolve_extract_result <- function(result, model, models, corpus, requested_k) {
  theta <- result$theta %||% result$Theta %||% result$deconProp
  beta <- result$beta %||% result$Beta %||% result$geneExpr
  if (is.null(theta)) {
    log_message(
      "{.pkg STdeconvolve} did not return topic proportion matrix {.val theta}",
      message_type = "error"
    )
  }
  selected_k <- ncol(as.matrix(theta))
  if (length(requested_k) == 1L) {
    selected_k <- requested_k
  }
  list(
    theta = theta,
    beta = beta,
    corpus = corpus,
    model = model,
    models = models,
    selected_k = selected_k
  )
}

stdeconvolve_extract_corpus <- function(corpus) {
  if (is.list(corpus) && !is.null(corpus$corpus)) {
    corpus <- corpus$corpus
  }
  if (is.data.frame(corpus)) {
    corpus <- as.matrix(corpus)
  }
  if (!is.matrix(corpus) && !inherits(corpus, "Matrix")) {
    log_message(
      "{.pkg STdeconvolve} corpus must be a matrix-like object",
      message_type = "error"
    )
  }
  corpus
}

stdeconvolve_filter_counts <- function(counts) {
  counts <- methods::as(counts, "dgCMatrix")
  keep_features <- Matrix::rowSums(counts) > 0
  keep_spots <- Matrix::colSums(counts) > 0
  counts[keep_features, keep_spots, drop = FALSE]
}

stdeconvolve_orient_theta <- function(theta, spot_ids) {
  theta <- as.matrix(theta)
  rn_match <- if (is.null(rownames(theta))) 0L else sum(rownames(theta) %in% spot_ids)
  cn_match <- if (is.null(colnames(theta))) 0L else sum(colnames(theta) %in% spot_ids)
  if (cn_match > rn_match) {
    theta <- t(theta)
  }
  if (is.null(rownames(theta))) {
    log_message(
      "{.pkg STdeconvolve} topic proportions must contain spatial spot names",
      message_type = "error"
    )
  }
  spots_use <- spot_ids[spot_ids %in% rownames(theta)]
  if (length(spots_use) == 0L) {
    log_message(
      "{.pkg STdeconvolve} topic proportions could not be matched to spatial spot names",
      message_type = "error"
    )
  }
  theta <- theta[spots_use, , drop = FALSE]
  if (is.null(colnames(theta)) || any(!nzchar(colnames(theta)))) {
    colnames(theta) <- paste0("topic_", seq_len(ncol(theta)))
  }
  colnames(theta) <- make.unique(as.character(colnames(theta)), sep = "_")
  theta
}

stdeconvolve_resolve_k <- function(k = NULL, k_candidates = 2:9) {
  k_use <- k %||% k_candidates
  if (!is.numeric(k_use) || length(k_use) == 0L || any(is.na(k_use) | k_use < 2)) {
    log_message(
      "{.arg k} or {.arg k_candidates} must contain topic numbers >= 2",
      message_type = "error"
    )
  }
  unique(as.integer(k_use))
}

stdeconvolve_validate_param_list <- function(x, arg_name) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg_name}} must be a list",
      message_type = "error"
    )
  }
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      "{.arg {arg_name}} must contain named arguments only",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

stdeconvolve_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
}
