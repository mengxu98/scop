make_stdeconvolve_seurat <- function() {
  counts <- matrix(
    c(
      10, 8, 1, 0,
      0, 2, 9, 8,
      6, 0, 1, 0,
      1, 7, 2, 5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$col <- c(1, 2, 1, 2)
  srt$row <- c(1, 1, 2, 2)
  srt
}

with_mock_stdeconvolve <- function(code) {
  clean_fun <- function(counts, ...) counts
  restrict_fun <- function(counts, ...) list(corpus = counts)
  fit_fun <- function(pixels, Ks, ...) {
    expect_equal(Ks, 2L)
    expect_equal(rownames(pixels), paste0("Spot", 1:4))
    list(k2 = list(pixels = pixels, Ks = Ks))
  }
  optimal_fun <- function(models, opt = "min", ...) {
    expect_identical(opt, "min")
    models[[1L]]
  }
  beta_theta_fun <- function(deconvolved, ...) {
    list(
      theta = matrix(
        c(
          0.80, 0.20,
          0.35, 0.65,
          0.10, 0.90,
          0.55, 0.45
        ),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(paste0("Spot", 1:4), c("topic_1", "topic_2"))
      ),
      beta = matrix(
        seq_len(8) / 8,
        nrow = 4,
        dimnames = list(paste0("Gene", 1:4), c("topic_1", "topic_2"))
      )
    )
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "STdeconvolve")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "STdeconvolve")
      switch(name,
        cleanCounts = clean_fun,
        restrictCorpus = restrict_fun,
        fitLDA = fit_fun,
        optimalModel = optimal_fun,
        getBetaTheta = beta_theta_fun,
        stop("unexpected function")
      )
    }
  )
  force(code)
}

test_that("RunSTdeconvolve writes topic proportions and tool results", {
  srt <- make_stdeconvolve_seurat()
  with_mock_stdeconvolve({
    out <- RunSTdeconvolve(srt, k = 2, verbose = FALSE)
  })

  expect_equal(unname(out$STdeconvolve_prop_topic_1), c(0.80, 0.35, 0.10, 0.55))
  expect_equal(unname(out$STdeconvolve_prop_topic_2), c(0.20, 0.65, 0.90, 0.45))
  expect_equal(unname(out$STdeconvolve_dominant_type), c("topic_1", "topic_2", "topic_2", "topic_1"))
  expect_equal(unname(out$STdeconvolve_max_prop), c(0.80, 0.65, 0.90, 0.55))
  expect_true("STdeconvolve" %in% names(out@tools))
  expect_equal(out@tools$STdeconvolve$selected_k, 2L)
  expect_equal(colnames(out@tools$STdeconvolve$theta), c("topic_1", "topic_2"))
})

test_that("RunSTdeconvolve validates inputs before backend work", {
  srt <- make_stdeconvolve_seurat()
  expect_error(
    RunSTdeconvolve(matrix(1, nrow = 2), verbose = FALSE),
    "Seurat"
  )
  with_mock_stdeconvolve({
    expect_error(
      RunSTdeconvolve(srt, features = "AbsentGene", verbose = FALSE),
      "No features"
    )
    expect_error(
      RunSTdeconvolve(srt, k = 1, verbose = FALSE),
      "topic numbers"
    )
    expect_error(
      RunSTdeconvolve(srt, fit_lda_params = list(1), verbose = FALSE),
      "named arguments"
    )
  })
})

test_that("STdeconvolvePlot uses SCOP spatial plotting", {
  srt <- make_stdeconvolve_seurat()
  srt$STdeconvolve_prop_topic_1 <- c(0.8, 0.4, 0.1, 0.6)
  srt$STdeconvolve_prop_topic_2 <- c(0.2, 0.6, 0.9, 0.4)
  srt$STdeconvolve_dominant_type <- c("topic_1", "topic_2", "topic_2", "topic_1")
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      if (identical(packages, "scatterpie")) {
        testthat::skip_if_not_installed("scatterpie")
      }
      invisible(TRUE)
    }
  )

  p1 <- STdeconvolvePlot(
    srt,
    topics = 1,
    overlay_image = FALSE
  )
  p2 <- STdeconvolvePlot(
    srt,
    plot_type = "pie",
    overlay_image = FALSE
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
