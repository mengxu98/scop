make_mistyr_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 2, 1, 3,
      0, 3, 0, 4, 1,
      1, 1, 1, 1, 1,
      2, 0, 2, 0, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:5))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$x <- c(0, 1, 2, 0, 1)
  srt$y <- c(0, 0, 0, 1, 1)
  srt
}

with_mock_mistyr <- function(code) {
  fake_create <- function(expr) {
    expect_true(is.data.frame(expr))
    expect_equal(colnames(expr), c("Gene1", "Gene2", "Gene3"))
    list(intraview = list(data = expr), misty.uniqueid = "mock")
  }
  fake_juxta <- function(current.views, positions, neighbor.thr, cached, verbose) {
    expect_equal(neighbor.thr, 2)
    expect_false(cached)
    current.views$juxtaview <- list(data = positions)
    current.views
  }
  fake_para <- function(current.views, positions, l, zoi, family, approx, nn, cached, verbose) {
    expect_equal(l, 3)
    expect_equal(zoi, 0)
    expect_identical(family, "gaussian")
    expect_equal(approx, 1)
    expect_null(nn)
    current.views$paraview.3 <- list(data = positions)
    current.views
  }
  fake_run <- function(
    views,
    results.folder,
    seed,
    target.subset,
    bypass.intra,
    cv.folds,
    cached,
    append,
    num.trees = NULL,
    ...
  ) {
    expect_true(all(c("intraview", "juxtaview", "paraview.3") %in% names(views)))
    expect_equal(seed, 7)
    expect_equal(target.subset, "Gene1")
    expect_false(bypass.intra)
    expect_equal(cv.folds, 3)
    expect_false(cached)
    expect_false(append)
    expect_equal(num.trees, 25)
    results.folder
  }
  fake_collect <- function(folders) {
    list(
      improvements = data.frame(target = c("Gene1", "Gene2"), measure = "gain.R2", value = c(1, 2)),
      contributions = data.frame(target = "Gene1", view = "paraview.3", value = 0.4),
      importances.aggregated = list(intra = data.frame(Predictor = "Gene2", Gene1 = 1))
    )
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "mistyR")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "mistyR")
      switch(name,
        create_initial_view = fake_create,
        add_juxtaview = fake_juxta,
        add_paraview = fake_para,
        run_misty = fake_run,
        collect_results = fake_collect,
        stop("unexpected mistyR function: ", name)
      )
    }
  )
  force(code)
}

test_that("RunMistyR builds views, collects results, and stores schema", {
  srt <- make_mistyr_seurat()
  with_mock_mistyr({
    out <- RunMistyR(
      srt,
      assay = "RNA",
      layer = "data",
      features = paste0("Gene", 1:3),
      coord.cols = c("x", "y"),
      views = c("juxta", "para"),
      para_l = 3,
      juxta_neighbor_thr = 2,
      results_folder = tempfile("mistyr-test-"),
      seed = 7,
      target_subset = "Gene1",
      cv_folds = 3,
      store_views = TRUE,
      num.trees = 25,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("MistyR" %in% names(out@tools))
  bundle <- out@tools$MistyR
  expect_named(bundle, c("results", "summary", "results_folder", "features", "cells", "parameters", "views"))
  expect_equal(bundle$features, paste0("Gene", 1:3))
  expect_equal(bundle$summary$n_targets, 2)
  expect_equal(bundle$summary$n_improvement_records, 2)
  expect_equal(bundle$summary$n_contribution_records, 1)
  expect_equal(bundle$summary$importance_views, "intra")
  expect_equal(bundle$parameters$views, c("juxta", "para"))
  expect_equal(bundle$parameters$backend_args$num.trees, 25)
})

test_that("RunMistyR validates features and views", {
  srt <- make_mistyr_seurat()
  expect_error(
    RunMistyR(
      srt,
      assay = "RNA",
      features = "missing",
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "features"
  )
  expect_error(
    RunMistyR(
      srt,
      assay = "RNA",
      features = "Gene1",
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "At least two features"
  )
  expect_error(
    RunMistyR(
      srt,
      assay = "RNA",
      views = "bad",
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "views"
  )
})

test_that("RunMistyR can skip result storage", {
  srt <- make_mistyr_seurat()
  with_mock_mistyr({
    out <- RunMistyR(
      srt,
      assay = "RNA",
      layer = "data",
      features = paste0("Gene", 1:3),
      coord.cols = c("x", "y"),
      views = c("juxta", "para"),
      para_l = 3,
      juxta_neighbor_thr = 2,
      seed = 7,
      target_subset = "Gene1",
      cv_folds = 3,
      store_results = FALSE,
      num.trees = 25,
      verbose = FALSE
    )
  })

  expect_false("MistyR" %in% names(out@tools))
})
