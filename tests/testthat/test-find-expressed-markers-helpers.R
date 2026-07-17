test_that("conserved-marker sign filtering matches legacy all() semantics", {
  logfc <- rbind(
    Positive = c(0.2, 0.4, 0.1),
    Negative = c(-0.2, -0.4, -0.1),
    Mixed = c(-0.1, 0.3, 0.2),
    Missing = c(NA_real_, 0.3, 0.2),
    Missing_with_false = c(NA_real_, -0.3, 0.2)
  )
  legacy <- function(x, only.pos) {
    positive <- apply(x > 0, 1L, all)
    if (only.pos) return(positive)
    apply(x < 0, 1L, all) | positive
  }

  expect_identical(
    find_conserved_marker_sign_keep(logfc, only.pos = TRUE),
    legacy(logfc, TRUE)
  )
  expect_identical(
    find_conserved_marker_sign_keep(logfc, only.pos = FALSE),
    legacy(logfc, FALSE)
  )
})

test_that("conserved-marker maximum p-values match legacy apply", {
  pvalues <- rbind(
    GeneA = c(0.01, 0.04, 0.02),
    GeneB = c(NA_real_, 0.03, 0.02),
    GeneC = c(0.4, 0.3, 0.5)
  )

  expect_identical(
    find_conserved_marker_max_pval(pvalues),
    apply(pvalues, 1L, max)
  )
})
