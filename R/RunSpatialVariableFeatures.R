#' @title Run spatial variable feature detection
#'
#' @description
#' Score genes by spot-level spatial autocorrelation using a lightweight
#' coordinate KNN graph. Moran's I is ranked high-to-low, while Geary's C is
#' converted to `1 - C` for the stored ranking score.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used for expression values.
#' @param features Features to score. If `NULL`, current variable features are
#' used; if no variable features are present, all assay features are used.
#' @param method Spatial autocorrelation statistic.
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
#'
#' @return A `Seurat` object with spatial variable feature results stored in
#' `srt@tools[["SpatialVariableFeatures"]]` and top feature names stored in
#' `srt@misc[["SpatialVariableFeatures"]]`.
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
#' spatial <- RunSpatialVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 50
#' )
#' SpatialSpotPlot(
#'   spatial,
#'   features = spatial@misc[["SpatialVariableFeatures"]][1:2]
#' )
RunSpatialVariableFeatures <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  method = c("moran", "geary"),
  image = NULL,
  coord.cols = c("x", "y"),
  k = 6,
  nfeatures = 2000,
  min_spots = 5,
  nperm = 0,
  set_variable_features = TRUE,
  store_results = TRUE,
  verbose = TRUE,
  seed = 11
) {
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

  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
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

  edges <- spatial_variable_knn_edges(coords, k = k)
  expr_mat <- as.matrix(expr)
  scores <- spatial_variable_score_matrix(
    expr = expr_mat,
    edges = edges,
    method = method,
    nperm = nperm
  )
  result <- data.frame(
    feature = rownames(expr_mat),
    statistic = scores$statistic,
    score = scores$score,
    p_value = scores$p_value,
    q_value = if (all(is.na(scores$p_value))) {
      NA_real_
    } else {
      stats::p.adjust(scores$p_value, method = "BH")
    },
    mean = rowMeans(expr_mat),
    variance = fast_row_vars(expr_mat),
    n_spots = as.integer(expressed_spots[rownames(expr_mat)]),
    stringsAsFactors = FALSE
  )
  result <- result[order(result$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  top_features <- utils::head(result$feature[is.finite(result$score)], nfeatures)
  srt@misc[["SpatialVariableFeatures"]] <- top_features
  if (isTRUE(set_variable_features) && length(top_features) > 0L) {
    SeuratObject::VariableFeatures(srt, assay = assay) <- top_features
  }
  if (isTRUE(store_results)) {
    srt@tools[["SpatialVariableFeatures"]] <- list(
      result = result,
      coords = coords,
      edges = edges,
      parameters = list(
        assay = assay,
        layer = layer,
        method = method,
        image = image,
        coord.cols = coord.cols,
        k = k,
        nfeatures = nfeatures,
        min_spots = min_spots,
        nperm = nperm,
        set_variable_features = set_variable_features
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

spatial_variable_knn_edges <- function(coords, k = 6) {
  coord_mat <- as.matrix(coords[, c("x", "y"), drop = FALSE])
  dmat <- as.matrix(stats::dist(coord_mat))
  diag(dmat) <- Inf
  nn <- t(apply(dmat, 1L, function(x) utils::head(order(x), k)))
  if (k == 1L) {
    nn <- matrix(nn, ncol = 1L)
  }
  data.frame(
    from = rep(seq_len(nrow(coords)), each = k),
    to = as.vector(t(nn)),
    stringsAsFactors = FALSE
  )
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
