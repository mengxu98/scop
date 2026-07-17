test_that("scomm_assign_labels matches thresholded per-cell selection", {
  probabilities <- matrix(
    c(0.8, 0.6, 0.2, 0.4, 0.7, 0.7, 0.3, 0.1, 0.9),
    nrow = 3,
    dimnames = list(NULL, c("A", "B", "C"))
  )
  thresholds <- c(A = 0.5, B = 0.65, C = 0.95)
  legacy <- apply(probabilities, 1, function(x) {
    valid <- x >= thresholds[names(x)]
    if (any(valid)) names(which.max(x[valid]))[[1]] else "unclassified"
  })

  expect_identical(scomm_assign_labels(probabilities, thresholds), legacy)
})

test_that("scomm_assign_labels retains the legacy missing-value path", {
  probabilities <- matrix(
    c(0.8, 0.4, 0.1, 0.2, NA_real_, 0.7),
    nrow = 2,
    dimnames = list(NULL, c("A", "B", "C"))
  )
  thresholds <- c(A = 0.5, B = 0.5, C = 0.5)

  expect_identical(
    scomm_assign_labels(probabilities, thresholds),
    c("A", "C")
  )
})

test_that("add_prediction_meta computes maximum probabilities by row", {
  counts <- Matrix::Matrix(
    matrix(1:6, nrow = 2,
      dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  probabilities <- data.frame(
    A = c(0.1, NA_real_, 0.6),
    B = c(0.7, NA_real_, 0.4),
    check.names = FALSE,
    row.names = colnames(srt)
  )
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    .package = "scop"
  )

  out <- add_prediction_meta(srt, c("B", "unclassified", "A"), probabilities)

  expect_equal(
    unname(out$srt[[out$probability_col, drop = TRUE]]),
    c(0.7, -Inf, 0.6),
    tolerance = 1e-12
  )
})
