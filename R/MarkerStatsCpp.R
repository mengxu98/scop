run_sparse_wilcox <- function(x, n_group1, min.expression = 0) {
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  }
  p_val <- wilcox_rank_sum_sparse(
    mat = x,
    n_group1 = as.integer(n_group1),
    min_expression = as.numeric(min.expression)
  )
  names(p_val) <- rownames(x)
  p_val
}

run_sparse_wilcox_all_cells <- function(x, n_group1) {
  if (!inherits(x, "dgCMatrix")) {
    x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  }
  p_val <- wilcox_rank_sum_sparse_all_cells(
    mat = x,
    n_group1 = as.integer(n_group1)
  )
  names(p_val) <- rownames(x)
  p_val
}
