test_that("dynamic_features_gam preserves failed feature slots", {
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
