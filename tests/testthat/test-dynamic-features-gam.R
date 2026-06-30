test_that("dynamic_features_gam preserves failed feature slots", {
  if (identical(Sys.info()[["sysname"]], "Darwin") && getRversion() >= "4.6.0") {
    skip("mgcv can segfault while loading on macOS ARM with R >= 4.6")
  }
  skip_if_not_installed("mgcv")

  cells <- paste0("cell", seq_len(30))
  t_ordered <- stats::setNames(seq(0, 1, length.out = length(cells)), cells)
  y_ordered <- rbind(
    good = 1 + sin(t_ordered * pi),
    partial = 1 + cos(t_ordered * pi),
    bad = rep(NA_real_, length(cells))
  )
  y_ordered["partial", c(3, 12)] <- NA_real_
  y_ordered["partial", 20] <- Inf
  colnames(y_ordered) <- cells
  family <- c(good = "gaussian", partial = "gaussian", bad = "gaussian")
  y_libsize <- stats::setNames(rep(1, length(cells)), cells)

  out <- dynamic_features_gam(
    y_ordered = y_ordered,
    t_ordered = t_ordered,
    features = rownames(y_ordered),
    gene = rownames(y_ordered),
    meta = character(),
    family = family,
    layer = "data",
    y_libsize = y_libsize,
    padjust_method = "fdr",
    cores = 1,
    verbose = FALSE
  )

  expect_equal(
    colnames(out$fitted_matrix),
    c("pseudotime", "good", "partial", "bad")
  )
  expect_equal(rownames(out$DynamicFeatures), c("good", "partial", "bad"))
  expect_false(all(is.na(out$fitted_matrix[, "partial"])))
  expect_false(is.na(out$DynamicFeatures["partial", "pvalue"]))
  expect_true(all(is.na(out$fitted_matrix[, "bad"])))
  expect_true(is.na(out$DynamicFeatures["bad", "pvalue"]))
})

test_that("dynamic feature sparse row unique counts match dense apply", {
  sparse <- Matrix::Matrix(
    c(
      0, 1, 1, 2,
      0, 0, 0, 0,
      3, 0, 3, 4
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(sparse) <- paste0("gene", seq_len(nrow(sparse)))
  dense_counts <- apply(as.matrix(sparse), 1, function(row) length(unique(row)))

  expect_equal(dynamic_row_unique_counts(sparse), dense_counts)
  expect_equal(dynamic_row_unique_counts(as.matrix(sparse)), dense_counts)
})

test_that("dynamic raw matrix matches legacy construction", {
  t_ordered <- stats::setNames(c(0.2, 0.5, 0.9), paste0("cell", 1:3))
  y_ordered <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("geneA", "geneB"), names(t_ordered))
  )
  legacy <- as.matrix(
    cbind(
      data.frame(pseudotime = t_ordered),
      Matrix::t(y_ordered)
    )
  )

  expect_equal(dynamic_raw_matrix(y_ordered, t_ordered), legacy)
})
