test_that("C++ GRNBoost2 sparse path matches dense path", {
  expr_dense <- matrix(
    c(
      0, 1, 0, 2, 0, 3,
      2, 0, 1, 0, 3, 0,
      0, 0, 2, 2, 1, 1,
      3, 2, 0, 0, 1, 0,
      1, 0, 3, 1, 0, 2,
      0, 3, 1, 0, 2, 1
    ),
    nrow = 6,
    byrow = TRUE
  )
  colnames(expr_dense) <- paste0("Gene", seq_len(ncol(expr_dense)))
  rownames(expr_dense) <- paste0("Cell", seq_len(nrow(expr_dense)))
  expr_sparse <- Matrix::Matrix(expr_dense, sparse = TRUE)
  regulators <- colnames(expr_dense)[1:3]
  targets <- colnames(expr_dense)[4:6]

  dense <- RunGRNBoost2(
    expr_dense,
    genes_in = "columns",
    regulators = regulators,
    targets = targets,
    backend = "cpp",
    n_rounds = 12,
    max_depth = 2,
    max_features = 0.8,
    subsample = 0.8,
    early_stop_window_length = 4,
    seed = 7,
    verbose = FALSE
  )
  auto_dense <- RunGRNBoost2(
    expr_sparse,
    genes_in = "columns",
    regulators = regulators,
    targets = targets,
    backend = "cpp",
    n_rounds = 12,
    max_depth = 2,
    max_features = 0.8,
    subsample = 0.8,
    early_stop_window_length = 4,
    seed = 7,
    verbose = FALSE
  )
  old_options <- options(scop.grnboost_dense_max_entries = 0)
  on.exit(options(old_options), add = TRUE)
  sparse <- RunGRNBoost2(
    expr_sparse,
    genes_in = "columns",
    regulators = regulators,
    targets = targets,
    backend = "cpp",
    n_rounds = 12,
    max_depth = 2,
    max_features = 0.8,
    subsample = 0.8,
    early_stop_window_length = 4,
    seed = 7,
    verbose = FALSE
  )

  dense <- dense[order(dense$TF, dense$target), , drop = FALSE]
  auto_dense <- auto_dense[order(auto_dense$TF, auto_dense$target), , drop = FALSE]
  sparse <- sparse[order(sparse$TF, sparse$target), , drop = FALSE]
  rownames(dense) <- NULL
  rownames(auto_dense) <- NULL
  rownames(sparse) <- NULL

  expect_equal(auto_dense, dense, tolerance = 1e-12)
  expect_equal(sparse, dense, tolerance = 1e-12)
})

test_that("C++ GRNBoost2 preserves sparse zero ordering on larger matrices", {
  old_options <- options(scop.grnboost_dense_max_entries = 0)
  on.exit(options(old_options), add = TRUE)
  set.seed(42)
  expr_dense <- matrix(stats::rpois(300L * 150L, lambda = 0.3), nrow = 300L)
  expr_dense[expr_dense < 1L] <- 0
  rownames(expr_dense) <- paste0("Cell", seq_len(nrow(expr_dense)))
  colnames(expr_dense) <- paste0("Gene", seq_len(ncol(expr_dense)))
  regulators <- colnames(expr_dense)[1:50]
  targets <- colnames(expr_dense)[51:150]

  dense <- RunGRNBoost2(
    expr_dense, genes_in = "columns", regulators = regulators, targets = targets,
    backend = "cpp", n_rounds = 30L, max_depth = 3L, max_features = 0.1,
    subsample = 0.9, early_stop_window_length = 25L, seed = 1234L,
    verbose = FALSE
  )
  sparse <- RunGRNBoost2(
    Matrix::Matrix(expr_dense, sparse = TRUE), genes_in = "columns",
    regulators = regulators, targets = targets, backend = "cpp", n_rounds = 30L,
    max_depth = 3L, max_features = 0.1, subsample = 0.9,
    early_stop_window_length = 25L, seed = 1234L, verbose = FALSE
  )
  dense <- dense[order(dense$TF, dense$target), , drop = FALSE]
  sparse <- sparse[order(sparse$TF, sparse$target), , drop = FALSE]
  rownames(dense) <- NULL
  rownames(sparse) <- NULL

  expect_equal(sparse, dense, tolerance = 1e-12)
})
