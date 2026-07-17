test_that("Palantir branch selection preserves finite-probability handling", {
  dat <- data.frame(
    BranchA = c(0.4, NA_real_, Inf, -Inf, 0.5),
    BranchB = c(0.4, 0.9, -Inf, NA_real_, 0.5),
    BranchC = c(0.1, 0.2, NA_real_, NA_real_, 0.4),
    row.names = paste0("Cell", seq_len(5))
  )
  branches <- colnames(dat)
  legacy <- function(x) {
    probs <- as.matrix(x[, branches, drop = FALSE])
    probs[is.na(probs)] <- -Inf
    selected <- branches[max.col(probs, ties.method = "first")]
    has_prob <- apply(probs, 1L, function(values) any(is.finite(values)))
    selected[!has_prob] <- NA_character_
    factor(selected, levels = branches)
  }

  expect_identical(
    palantir_branch_selection(dat, branches),
    legacy(dat)
  )
})
