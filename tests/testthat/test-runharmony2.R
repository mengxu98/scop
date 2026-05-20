test_that("Harmony result extraction supports fields and callable methods", {
  zcorr <- matrix(seq_len(6), nrow = 2)
  clusters <- matrix(seq_len(8), nrow = 2)

  old_object <- list(Z_corr = zcorr, R = clusters)
  expect_equal(
    scop:::extract_harmony_component(
      old_object,
      field = "Z_corr",
      method = "getZcorr",
      label = "corrected embeddings"
    ),
    zcorr
  )

  new_object <- new.env(parent = emptyenv())
  new_object$getZcorr <- function() zcorr
  new_object$getR <- function() clusters
  expect_equal(
    scop:::extract_harmony_component(
      new_object,
      field = "Z_corr",
      method = "getZcorr",
      label = "corrected embeddings"
    ),
    zcorr
  )
  expect_equal(
    scop:::extract_harmony_component(
      new_object,
      field = "R",
      method = "getR",
      label = "soft cluster assignments"
    ),
    clusters
  )
})
