test_that("Seurat dynamic fitting preserves storage", {
  set.seed(42)
  x <- matrix(rpois(6 * 60, 5), 6, dimnames = list(paste0("g", 1:6), paste0("c", 1:60)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(x, sparse = TRUE))
  object$time <- seq(0, 1, length.out = 60)
  out <- RunDynamicFeatures(object, "time",
    features = rownames(object),
    fit_method = "pretsa", verbose = FALSE
  )
  reference <- thisutils::fit_trends(log1p(x), object$time, method = "pretsa")
  table <- reference$statistics
  names(table)[names(table) == "n_above_min"] <- "exp_ncells"
  expect_equal(out@tools$DynamicFeatures_time$DynamicFeatures, table)
  expect_equal(out@tools$DynamicFeatures_time$fitted_matrix[, -1], t(reference$fitted))
  expect_equal(out@tools$DynamicFeatures_time$raw_matrix[, -1], t(x))
})
