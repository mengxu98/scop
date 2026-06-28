make_cibersort_mock <- function(
  n_genes = 10,
  n_cell_types = 3,
  n_samples = 5,
  seed = 1
) {
  set.seed(seed)
  signature <- matrix(
    stats::runif(n_genes * n_cell_types, 50, 500),
    nrow = n_genes,
    ncol = n_cell_types,
    dimnames = list(
      sprintf("Gene%02d", seq_len(n_genes)),
      sprintf("Cell%02d", seq_len(n_cell_types))
    )
  )
  fractions <- matrix(
    stats::runif(n_cell_types * n_samples),
    nrow = n_cell_types,
    ncol = n_samples,
    dimnames = list(colnames(signature), sprintf("Sample%02d", seq_len(n_samples)))
  )
  fractions <- sweep(fractions, 2, colSums(fractions), "/")
  mixture <- signature %*% fractions
  rownames(mixture) <- rownames(signature)
  colnames(mixture) <- colnames(fractions)
  list(signature = signature, mixture = mixture, fractions = t(fractions))
}

test_that("RunCIBERSORT cpp backend returns a valid bundle", {
  dat <- make_cibersort_mock()
  out <- RunCIBERSORT(
    count_matrix = dat$mixture,
    sig_matrix = dat$signature,
    backend = "cpp",
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  )

  expect_equal(out$status, "success")
  expect_equal(out$parameters$backend, "cpp")
  expect_equal(out$details$engine, "scop::cibersort_cpp")
  expect_false(out$details$package_backend)
  prop <- out$details$proportion_matrix
  expect_equal(rownames(prop), colnames(dat$mixture))
  expect_equal(colnames(prop), colnames(dat$signature))
  expect_true(all(is.finite(prop)))
  expect_true(all(prop >= 0))
  expect_equal(unname(rowSums(prop)), rep(1, nrow(prop)), tolerance = 1e-8)
  expect_true(all(c("P-value", "Correlation", "RMSE") %in% colnames(out$details$statistics)))
  expect_true(nrow(out$results) > 0)
})

test_that("RunCIBERSORT cpp backend supports permutations and quantile normalization", {
  dat <- make_cibersort_mock(seed = 2)
  out <- RunCIBERSORT(
    count_matrix = dat$mixture,
    sig_matrix = dat$signature,
    backend = "cpp",
    perm = 10,
    QN = TRUE,
    n_threads = 2,
    seed = 99,
    verbose = FALSE
  )

  expect_equal(out$status, "success")
  stats <- out$details$statistics
  expect_true(all(is.finite(stats[["P-value"]])))
  expect_true(all(stats[["P-value"]] >= 0 & stats[["P-value"]] <= 1))
  expect_true(all(is.finite(stats[["Correlation"]])))
  expect_true(all(is.finite(stats[["RMSE"]])))
})

test_that("RunCIBERSORT cpp backend stores SummarizedExperiment metadata", {
  skip_if_not_installed("SummarizedExperiment")
  dat <- make_cibersort_mock(seed = 3)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = dat$mixture),
    colData = S4Vectors::DataFrame(
      condition = rep(c("A", "B"), length.out = ncol(dat$mixture))
    )
  )
  out <- RunCIBERSORT(
    object = se,
    sig_matrix = dat$signature,
    backend = "cpp",
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  )

  store <- S4Vectors::metadata(out)[["Deconvolution"]]
  expect_equal(store$input$backend, "cpp")
  expect_equal(store$status$status, "success")
  expect_true(nrow(store$results) > 0)
})

test_that("RunBayesPrism accepts sample.by from RunDeconvolution dispatch", {
  expect_true("sample.by" %in% names(formals(scop:::RunBayesPrism)))
})

test_that("RunCIBERSORT cpp backend matches the R backend on small deterministic data", {
  skip_if_not_installed("CIBERSORT")
  dat <- make_cibersort_mock(seed = 4)
  old <- suppressWarnings(RunCIBERSORT(
    count_matrix = dat$mixture,
    sig_matrix = dat$signature,
    backend = "r",
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  ))
  new <- RunCIBERSORT(
    count_matrix = dat$mixture,
    sig_matrix = dat$signature,
    backend = "cpp",
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  )

  skip_if_not(identical(old$status, "success"))
  old_fraction <- old$details$proportion_matrix
  new_fraction <- new$details$proportion_matrix[
    rownames(old_fraction),
    colnames(old_fraction),
    drop = FALSE
  ]
  expect_gt(stats::cor(as.vector(old_fraction), as.vector(new_fraction)), 0.995)
  expect_lt(mean(abs(old_fraction - new_fraction)), 0.01)
  expect_gt(mean(max.col(old_fraction) == max.col(new_fraction)), 0.95)
})
