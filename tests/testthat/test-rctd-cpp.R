test_that("RCTD sparse quality helper matches R sparse summaries", {
  st_dense <- matrix(
    c(
      1, 0, 2,
      0, 0, 0,
      3, 0, 0,
      0, 5, 0,
      0, 0, 0
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("Spot", 1:3))
  )
  ref_dense <- matrix(
    c(
      0, 4, 0, 1,
      2, 0, 0, 0,
      0, 0, 0, 0,
      1, 1, 0, 0,
      0, 0, 7, 0
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:4))
  )
  st <- methods::as(Matrix::Matrix(st_dense, sparse = TRUE), "dgCMatrix")
  ref <- methods::as(Matrix::Matrix(ref_dense, sparse = TRUE), "dgCMatrix")

  out <- rctd_sparse_quality_cpp(st, ref)
  keep <- Matrix::rowSums(st) > 0 & Matrix::rowSums(ref) > 0

  expect_equal(out$keep_features, unname(keep))
  expect_equal(
    out$st_numi,
    unname(Matrix::colSums(st[keep, , drop = FALSE]))
  )
  expect_equal(
    out$ref_numi,
    unname(Matrix::colSums(ref[keep, , drop = FALSE]))
  )
  expect_equal(rctd_nonzero_shared_features(st, ref), rownames(st)[keep])
})

test_that("RCTD weight normalization helper matches R behavior", {
  weights <- matrix(
    c(
      1, 2, 3,
      NA, -1, Inf,
      0, 0, 0,
      4, 0, 4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Spot", 1:4), paste0("Type", 1:3))
  )
  expected <- weights
  expected[!is.finite(expected) | expected < 0] <- 0
  totals <- rowSums(expected)
  keep <- is.finite(totals) & totals > 0
  expected[keep, ] <- sweep(expected[keep, , drop = FALSE], 1, totals[keep], "/")
  expected[!keep, ] <- 0

  out <- rctd_normalize_weights(weights)
  expect_equal(out, expected)
  expect_equal(rownames(out), rownames(weights))
  expect_equal(colnames(out), colnames(weights))
})

test_that("RCTD metadata helper matches metadata fill and dominant summaries", {
  weights <- matrix(
    c(
      0.25, 0.75,
      0.00, 0.00,
      0.60, 0.40
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Spot1", "Spot3", "Spot4"), c("Alpha/Beta", "Delta"))
  )
  all_spots <- paste0("Spot", 1:5)
  out <- rctd_metadata_cpp(weights, all_spots)

  expected_weights <- matrix(
    NA_real_,
    nrow = length(all_spots),
    ncol = ncol(weights),
    dimnames = list(all_spots, colnames(weights))
  )
  expected_weights[rownames(weights), ] <- weights
  expected_dominant <- c("Delta", NA, NA, "Alpha/Beta", NA)
  expected_max <- c(0.75, NA, 0, 0.60, NA)

  expect_equal(out$weights, expected_weights)
  expect_equal(out$dominant, expected_dominant)
  expect_equal(out$max_prop, expected_max)
})

test_that("RCTD finalize helper matches separate normalization and metadata helpers", {
  weights <- matrix(
    c(
      1, 3, 0,
      NA, -2, Inf,
      0, 0, 0,
      4, 1, 5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Spot", c(1, 3, 4, 6)), paste0("Type", 1:3))
  )
  all_spots <- paste0("Spot", 1:6)

  expected_weights <- rctd_normalize_weights(weights)
  expected_meta <- rctd_metadata_cpp(expected_weights, all_spots)
  out <- rctd_finalize_weights_cpp(weights, all_spots)

  expect_equal(out$weights, expected_weights)
  expect_equal(out$full_weights, expected_meta$weights)
  expect_equal(out$dominant, expected_meta$dominant)
  expect_equal(out$max_prop, expected_meta$max_prop)
})
