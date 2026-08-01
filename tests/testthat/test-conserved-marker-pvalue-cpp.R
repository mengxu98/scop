test_that("conserved-marker p-value C++ methods agree with R", {
  combine <- getFromNamespace("find_conserved_marker_combine_pval", "scop")
  set.seed(52)
  pvalues <- matrix(runif(500 * 12), nrow = 500)
  rownames(pvalues) <- paste0("g", seq_len(nrow(pvalues)))
  pvalues[1, 1] <- -1
  pvalues[2, 2] <- 2
  methods <- c(
    "maximump",
    "minimump",
    "wilkinsonp",
    "meanp",
    "sump",
    "votep"
  )

  for (method in methods) {
    reference <- combine(pvalues, method = method, backend = "r")
    native <- combine(pvalues, method = method, backend = "cpp")
    expect_identical(native, reference)
  }
})

test_that("conserved-marker p-value C++ handles insufficient valid values", {
  combine <- getFromNamespace("find_conserved_marker_combine_pval", "scop")
  pvalues <- matrix(c(0.1, NA, 2, -1, 0.2, 0.3), nrow = 2, byrow = TRUE)
  rownames(pvalues) <- c("a", "b")

  expect_true(is.na(combine(pvalues, "maximump", "cpp")[[1]]))
  expect_true(is.na(combine(pvalues, "meanp", "cpp")[[2]]))
})
