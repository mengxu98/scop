run_sparse_topk_by_column <- function(x, k, decreasing = TRUE) {
  if (!run_sparse_topk_by_column_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  }
  sparse_topk_by_column(
    mat = x,
    k = as.integer(k),
    decreasing = isTRUE(decreasing)
  )
}

run_sparse_topk_by_column_available <- function() {
  exists("sparse_topk_by_column", mode = "function") &&
    isTRUE(is.loaded("_scop_sparse_topk_by_column"))
}

run_dense_topk_by_column <- function(x, k, decreasing = FALSE) {
  if (!run_dense_topk_by_column_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  dense_topk_by_column(
    mat = x,
    k = as.integer(k),
    decreasing = isTRUE(decreasing)
  )
}

run_dense_topk_by_column_available <- function() {
  exists("dense_topk_by_column", mode = "function") &&
    isTRUE(is.loaded("_scop_dense_topk_by_column"))
}
