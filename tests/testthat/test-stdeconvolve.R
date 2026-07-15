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
  fit_fun <- function(counts, Ks, ...) {
    expect_equal(Ks, 2L)
    expect_equal(rownames(counts), paste0("Spot", 1:4))
    list(k2 = list(counts = counts, Ks = Ks))
  }
  optimal_fun <- function(models, opt = "min", ...) {
    expect_identical(opt, "min")
    models[[1L]]
  }
  beta_theta_fun <- function(lda, ...) {
    expect_named(lda, c("counts", "Ks"))
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
    .package = "scop",
    check_r = function(packages, ...) {
      expect_identical(packages, "JEFworks-Lab/STdeconvolve")
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
  expect_named(out@tools$STdeconvolve$summary, c("n_spots", "n_types", "dominant_counts", "max_prop"))
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
  with_mock_stdeconvolve({
    srt <- RunSTdeconvolve(srt, k = 2, verbose = FALSE)
  })
  testthat::local_mocked_bindings(
    .package = "scop",
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

test_that("STdeconvolvePlot reads custom tool names and stored prefixes", {
  srt <- make_stdeconvolve_seurat()
  with_mock_stdeconvolve({
    out <- RunSTdeconvolve(
      srt,
      k = 2,
      prefix = "STFull",
      tool_name = "STdeconvolveFull",
      verbose = FALSE
    )
  })

  p <- STdeconvolvePlot(
    out,
    tool_name = "STdeconvolveFull",
    topics = "STFull_prop_topic_2",
    overlay_image = FALSE
  )
  expect_s3_class(p, "ggplot")
  expect_identical(GetSpatialResult(out, tool_name = "STdeconvolveFull")$parameters$prefix, "STFull")
})

test_that("STdeconvolvePlot uses readable automatic layouts and supports lists", {
  srt <- make_stdeconvolve_seurat()
  theta <- matrix(
    seq_len(ncol(srt) * 7),
    nrow = ncol(srt),
    dimnames = list(colnames(srt), paste0("topic_", 1:7))
  )
  theta <- theta / rowSums(theta)
  srt@tools$STSeven <- spatial_result_build(
    bundle = list(
      theta = theta,
      parameters = list(prefix = "STSeven"),
      summary = scop_spatial_weight_summary(theta)
    ),
    method = "STdeconvolve",
    result_type = "deconvolution",
    provenance = list(producer = "RunSTdeconvolve", backend_id = "stdeconvolve")
  )

  plots <- STdeconvolvePlot(
    srt,
    tool_name = "STSeven",
    combine = FALSE,
    overlay_image = FALSE
  )
  combined <- STdeconvolvePlot(
    srt,
    tool_name = "STSeven",
    overlay_image = FALSE
  )
  expect_length(plots, 7L)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_s3_class(combined, "patchwork")
  combined_plots <- combined$patches$plots
  expect_length(combined_plots, 6L)
  limits <- lapply(combined_plots, function(plot) {
    color_scale <- Filter(
      function(scale) any(scale$aesthetics %in% c("colour", "color")),
      plot$scales$scales
    )
    color_scale[[1L]]$limits
  })
  expect_true(all(vapply(limits, identical, logical(1), limits[[1L]])))
  expect_identical(combined$patches$annotation$title, "STSeven topic proportions")
})

test_that("STdeconvolvePlot rejects missing, stale, and malformed stored results", {
  srt <- make_stdeconvolve_seurat()
  expect_error(STdeconvolvePlot(srt, tool_name = "missing"), "No stored spatial result")

  theta <- matrix(
    0.5,
    nrow = ncol(srt),
    ncol = 2,
    dimnames = list(colnames(srt), c("topic_1", "topic_2"))
  )
  make_bundle <- function(value) spatial_result_build(
    bundle = list(theta = value, parameters = list(prefix = "Bad")),
    method = "STdeconvolve",
    result_type = "deconvolution",
    provenance = list(producer = "RunSTdeconvolve", backend_id = "stdeconvolve")
  )

  srt@tools$StaleST <- make_bundle(theta[-1, , drop = FALSE])
  expect_error(
    STdeconvolvePlot(srt, tool_name = "StaleST", overlay_image = FALSE),
    "stale or incomplete"
  )
  duplicated_theta <- theta
  rownames(duplicated_theta)[2] <- rownames(duplicated_theta)[1]
  srt@tools$DuplicatedST <- make_bundle(duplicated_theta)
  expect_error(
    STdeconvolvePlot(srt, tool_name = "DuplicatedST", overlay_image = FALSE),
    "unique, non-missing"
  )
  srt@tools$EmptyST <- make_bundle(theta[0, , drop = FALSE])
  expect_error(
    STdeconvolvePlot(srt, tool_name = "EmptyST", overlay_image = FALSE),
    "empty"
  )
})
