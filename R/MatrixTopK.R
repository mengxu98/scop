run_sparse_topk_by_column <- function(x, k, decreasing = TRUE) {
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  }
  sparse_topk_by_column(
    mat = x,
    k = as.integer(k),
    decreasing = isTRUE(decreasing)
  )
}

run_dense_topk_by_column <- function(x, k, decreasing = FALSE) {
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
