test_that("SCENIC edge correlation matches stats::cor without a full pair matrix", {
  expr <- cbind(
    TF1 = c(1, 2, 3, 4, 5),
    TF2 = c(5, 4, 3, 2, 1),
    TargetA = c(2, 4, 6, 8, 10),
    TargetB = c(1, 1, 1, 1, 1),
    Unused = c(4, 3, 2, 1, 0)
  )
  adjacency <- data.frame(
    TF = c("TF1", "TF2", "TF1", "missing"),
    target = c("TargetA", "TargetA", "TargetB", "TargetA"),
    importance = c(0.9, 0.8, 0.7, 0.6),
    stringsAsFactors = FALSE
  )

  actual <- scenic_add_correlation(adjacency, expr, rho_threshold = 0.5)
  expected_rho <- c(1, -1, 0)
  expect_equal(nrow(actual), 3L)
  expect_equal(actual$rho, expected_rho, tolerance = 1e-12)
  expect_equal(actual$regulation, c(1L, -1L, 0L))
})

test_that("SCENIC edge correlation treats non-finite values as legacy zero rho", {
  expr <- cbind(
    TF1 = c(1, 2, NA, 4),
    TargetA = c(2, 4, 6, 8)
  )
  adjacency <- data.frame(TF = "TF1", target = "TargetA")
  actual <- scenic_add_correlation(adjacency, expr)
  expect_equal(actual$rho, 0)
  expect_equal(actual$regulation, 0L)
})

test_that("SCENICPlus edge correlation matches the legacy pair matrix", {
  expr <- rbind(
    TF1 = c(1, 2, 3, 4, 5),
    TF2 = c(5, 4, 3, 2, 1),
    TargetA = c(2, 4, 6, 8, 10),
    TargetB = c(1, 1, 1, 1, 1)
  )
  adjacency <- data.frame(
    TF = c("TF1", "TF2", "TF1", "missing"),
    target = c("TargetA", "TargetA", "TargetB", "TargetA"),
    importance = c(0.9, 0.8, 0.7, 0.6),
    stringsAsFactors = FALSE
  )
  legacy <- function(adjacency, expr) {
    corr <- stats::cor(t(expr), method = "pearson", use = "pairwise.complete.obs")
    rho <- vapply(seq_len(nrow(adjacency)), function(i) {
      tf <- adjacency[["TF"]][[i]]
      target <- adjacency[["target"]][[i]]
      if (!tf %in% rownames(corr) || !target %in% colnames(corr)) return(NA_real_)
      as.numeric(corr[tf, target])
    }, numeric(1))
    rho
  }
  actual <- scenicplus_add_tfg_cor(adjacency, expr)
  expect_equal(
    actual$rho,
    suppressWarnings(legacy(adjacency, expr)),
    tolerance = 1e-12
  )

  expr[[1L, 3L]] <- NA_real_
  actual_nonfinite <- suppressWarnings(scenicplus_add_tfg_cor(adjacency, expr))
  expect_equal(
    actual_nonfinite$rho,
    suppressWarnings(legacy(adjacency, expr)),
    tolerance = 1e-12
  )
})
