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
    expect_equal(native, reference, tolerance = 1e-12)
  }
})

test_that("conserved-marker p-value C++ handles insufficient valid values", {
  combine <- getFromNamespace("find_conserved_marker_combine_pval", "scop")
  pvalues <- matrix(c(0.1, NA, 2, -1, 0.2, 0.3), nrow = 2, byrow = TRUE)
  rownames(pvalues) <- c("a", "b")

  expect_true(is.na(combine(pvalues, "maximump", "cpp")[[1]]))
  expect_true(is.na(combine(pvalues, "meanp", "cpp")[[2]]))
})

test_that("conserved-marker sum p-values agree with an independent high-precision reference", {
  combine <- getFromNamespace("find_conserved_marker_combine_pval", "scop")

  # Irwin-Hall reference values computed by evaluating the alternating
  # polynomial with 250-digit decimal arithmetic, not the C++ recurrence.
  pvalues <- outer(c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 1), rep(1, 100))
  actual <- unname(combine(pvalues, "sump", "cpp"))
  expected <- c(
    0,
    6.243339283753968e-13,
    0.00025065623009830237,
    0.5,
    0.99974934376990170,
    0.99999999999937567,
    1
  )
  expect_true(all(is.finite(actual) & actual >= 0 & actual <= 1))
  expect_true(all(diff(actual) > 0))
  expect_identical(actual[c(1L, 4L, 7L)], c(0, 0.5, 1))
  expect_equal(actual[2:3] / expected[2:3], c(1, 1), tolerance = 1e-12)
  expect_equal(actual[5:6], expected[5:6], tolerance = 1e-15)
  expect_equal(actual + rev(actual), rep(1, length(actual)), tolerance = 1e-15)

  pvalues_50 <- outer(c(0.466, 0.534, 0.6, 0.7), rep(1, 50))
  actual_50 <- unname(combine(pvalues_50, "sump"))
  expected_50 <- c(
    0.20301575157035084, 0.79698424842964916,
    0.99299508839848781, 0.99999973115347334
  )
  expect_equal(actual_50, expected_50, tolerance = 1e-15)
  r_50 <- unname(combine(pvalues_50, "sump", "r"))
  expect_equal(r_50, expected_50, tolerance = 1e-15)

  r_100 <- unname(combine(pvalues, "sump", "r"))
  expect_identical(r_100[c(1L, 4L, 7L)], c(0, 0.5, 1))
  expect_equal(r_100[2:3] / expected[2:3], c(1, 1), tolerance = 1e-12)
  expect_equal(r_100[5:6], expected[5:6], tolerance = 1e-15)

  tiny <- matrix(0.001, nrow = 1L, ncol = 100L)
  for (backend in c("cpp", "r")) {
    expect_equal(
      unname(combine(tiny, "sump", backend)) / 1.071510288125467e-258,
      1,
      tolerance = 1e-12
    )
  }
})

test_that("conserved-marker sum p-values preserve valid-value filtering", {
  combine <- getFromNamespace("find_conserved_marker_combine_pval", "scop")
  pvalues <- rbind(
    valid = c(0.125, 0.25, NA, NaN, Inf, -Inf, -1, 2),
    insufficient = c(0.2, NA, NA, NaN, Inf, -Inf, -1, 2),
    absent = c(NA, NA, NA, NaN, Inf, -Inf, -1, 2)
  )
  # For two uniforms and x <= 1 the exact CDF is x^2 / 2.
  for (backend in c("cpp", "r")) {
    expect_no_warning(result <- combine(pvalues, "sump", backend))
    expect_equal(
      result,
      c(valid = 0.375^2 / 2, insufficient = NA_real_, absent = NA_real_),
      tolerance = 1e-15
    )
  }
})
