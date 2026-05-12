run_cpp_knn <- function(
  reference,
  query = NULL,
  k,
  metric = c("euclidean", "cosine"),
  exclude_self = FALSE,
  n_threads = 0L
) {
  metric <- match.arg(metric)
  reference <- as.matrix(reference)
  query <- as.matrix(query %||% reference)
  storage.mode(reference) <- "double"
  storage.mode(query) <- "double"
  if (ncol(reference) != ncol(query)) {
    log_message(
      "{.arg reference} and {.arg query} must have the same number of columns",
      message_type = "error"
    )
  }

  result <- knn_topk_cpp(
    reference = reference,
    query = query,
    k = as.integer(k),
    metric = metric,
    exclude_self = isTRUE(exclude_self),
    n_threads = as.integer(n_threads)
  )
  rownames(result[["idx"]]) <- rownames(query)
  rownames(result[["dist"]]) <- rownames(query)
  result
}

run_knn_topk <- function(
  reference,
  query = NULL,
  k,
  metric = c("euclidean", "cosine"),
  backend = "cpp",
  exclude_self = FALSE,
  n_threads = 0L
) {
  metric <- match.arg(metric)
  if (!identical(backend, "cpp")) {
    log_message(
      "{.arg backend} must be {.val cpp}",
      message_type = "error"
    )
  }
  run_cpp_knn(
    reference = reference,
    query = query,
    k = k,
    metric = metric,
    exclude_self = exclude_self,
    n_threads = n_threads
  )
}
