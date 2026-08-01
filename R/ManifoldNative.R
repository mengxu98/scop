native_manifold_metric_id <- function(metric) {
  metric <- tolower(metric)
  if (identical(metric, "angular")) {
    metric <- "cosine"
  }
  metric_id <- match(metric, c("euclidean", "manhattan", "cosine", "hamming"))
  if (is.na(metric_id)) {
    log_message(
      "Native manifold backend does not support distance {.val {metric}}",
      message_type = "error"
    )
  }
  metric_id
}

native_manifold_prepare_data <- function(data, apply_pca = TRUE) {
  data <- as.matrix(data)
  storage.mode(data) <- "double"
  if (any(!is.finite(data))) {
    log_message(
      "Native manifold backends require finite input values",
      message_type = "error"
    )
  }
  if (isTRUE(apply_pca) && ncol(data) > 100L) {
    rank <- min(100L, ncol(data), nrow(data) - 1L)
    data <- stats::prcomp(
      data,
      center = TRUE,
      scale. = FALSE,
      rank. = rank
    )$x[, seq_len(rank), drop = FALSE]
  }
  data
}

native_pacmap_prepare_data <- function(data, apply_pca = TRUE) {
  data <- as.matrix(data)
  storage.mode(data) <- "double"
  if (any(!is.finite(data))) {
    log_message(
      "Native PaCMAP requires finite input values",
      message_type = "error"
    )
  }
  if (isTRUE(apply_pca) && ncol(data) > 100L) {
    rank <- min(100L, ncol(data), nrow(data) - 1L)
    return(stats::prcomp(
      data,
      center = TRUE,
      scale. = FALSE,
      rank. = rank
    )$x[, seq_len(rank), drop = FALSE])
  }
  data <- data - min(data)
  maximum <- max(data)
  if (is.finite(maximum) && maximum > 0) {
    data <- data / maximum
  }
  sweep(data, 2L, colMeans(data), "-")
}

native_manifold_initialization <- function(
  data,
  n_components,
  init = c("pca", "random"),
  seed = 11L
) {
  init <- match.arg(init)
  set.seed(seed)
  if (identical(init, "random")) {
    return(matrix(
      stats::rnorm(nrow(data) * n_components, sd = 1e-4),
      nrow = nrow(data),
      ncol = n_components
    ))
  }
  rank <- min(n_components, ncol(data), nrow(data) - 1L)
  initial <- stats::prcomp(
    data,
    center = TRUE,
    scale. = FALSE,
    rank. = rank
  )$x[, seq_len(rank), drop = FALSE]
  if (rank < n_components) {
    initial <- cbind(
      initial,
      matrix(0, nrow(initial), n_components - rank)
    )
  }
  initial <- scale(initial, center = TRUE, scale = FALSE)
  initial[!is.finite(initial)] <- 0
  as.matrix(initial) * 0.01
}
