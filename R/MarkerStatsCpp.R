run_sparse_wilcox <- function(x, n_group1, min.expression = 0) {
  if (!run_sparse_wilcox_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
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

run_sparse_wilcox_available <- function() {
  exists("wilcox_rank_sum_sparse", mode = "function") &&
    isTRUE(is.loaded("_scop_wilcox_rank_sum_sparse"))
}

run_sparse_wilcox_all_cells <- function(x, n_group1) {
  if (!run_sparse_wilcox_all_cells_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
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

run_sparse_wilcox_all_cells_available <- function() {
  exists("wilcox_rank_sum_sparse_all_cells", mode = "function") &&
    isTRUE(is.loaded("_scop_wilcox_rank_sum_sparse_all_cells"))
}
