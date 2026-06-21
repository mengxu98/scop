make_spotsweeper_seurat <- function(include_mito = TRUE) {
  genes <- if (isTRUE(include_mito)) {
    c("Gene1", "Gene2", "MT-ND1", "Gene4")
  } else {
    c("Gene1", "Gene2", "Gene3", "Gene4")
  }
  counts <- matrix(
    c(
      10, 1, 8, 7, 5, 6,
      0, 6, 2, 5, 4, 1,
      1, 2, 3, 1, 2, 3,
      5, 0, 0, 2, 0, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(genes, paste0("Spot", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(0, 1, 2, 0, 1, 2)
  srt$y <- c(0, 0, 0, 1, 1, 1)
  srt$sample <- rep(c("S1", "S2"), each = 3)
  srt
}

skip_if_no_spotsweeper_infra <- function() {
  testthat::skip_if_not_installed("SpatialExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
}

with_mock_spotsweeper <- function(code) {
  fake_local <- function(
    spe,
    metric,
    direction,
    n_neighbors,
    samples,
    log,
    cutoff,
    workers
  ) {
    expect_true(metric %in% colnames(as.data.frame(SummarizedExperiment::colData(spe))))
    expect_true(direction %in% c("lower", "higher", "both"))
    expect_lte(n_neighbors, 2)
    expect_identical(workers, 1L)
    cdata <- as.data.frame(SummarizedExperiment::colData(spe))
    outlier <- rep(FALSE, nrow(cdata))
    names(outlier) <- colnames(spe)
    if (identical(metric, "nCount_RNA")) {
      outlier["Spot1"] <- TRUE
    }
    if (identical(metric, "percent.mito")) {
      outlier["Spot2"] <- TRUE
    }
    cdata[[paste0(metric, "_outliers")]] <- outlier
    cdata[[paste0(metric, "_z")]] <- ifelse(outlier, cutoff + 1, 0)
    SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(cdata)
    spe
  }
  fake_artifact <- function(
    spe,
    mito_percent,
    mito_sum,
    samples,
    n_order,
    shape,
    log,
    name
  ) {
    cdata <- as.data.frame(SummarizedExperiment::colData(spe))
    expect_true(mito_percent %in% colnames(cdata))
    expect_true(mito_sum %in% colnames(cdata))
    expect_equal(length(unique(cdata[[samples]])), 1)
    expect_identical(shape, "hexagonal")
    cdata$artifact <- colnames(spe) %in% c("Spot4", "Spot6")
    cdata$k18 <- seq_len(nrow(cdata))
    SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(cdata)
    spe
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_true(all(c(
        "SpotSweeper", "SpatialExperiment",
        "SummarizedExperiment", "S4Vectors"
      ) %in% packages))
      invisible(TRUE)
    },
    spot_sweeper_get_fun = function(fun) {
      switch(fun,
        localOutliers = fake_local,
        findArtifacts = fake_artifact,
        stop("unexpected SpotSweeper function")
      )
    }
  )
  force(code)
}

test_that("RunSpotSweeper validates Seurat input before backend work", {
  expect_error(
    RunSpotSweeper(matrix(1, nrow = 2), verbose = FALSE),
    "Seurat"
  )
})

test_that("RunSpotSweeper stores local outlier and artifact metadata", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  with_mock_spotsweeper({
    out <- RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      n_order = 2,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true(all(c(
    "SpotSweeper_QC",
    "SpotSweeper_local_outlier_qc",
    "SpotSweeper_artifact_qc",
    "SpotSweeper_nCount_RNA_outlier",
    "SpotSweeper_nFeature_RNA_outlier",
    "SpotSweeper_percent.mito_outlier",
    "SpotSweeper_nCount_RNA_z",
    "SpotSweeper_artifact"
  ) %in% colnames(out@meta.data)))
  expect_equal(as.character(out$SpotSweeper_QC[1:4]), c("Fail", "Fail", "Pass", "Fail"))
  expect_true(out$SpotSweeper_nCount_RNA_outlier[1])
  expect_true(out$SpotSweeper_percent.mito_outlier[2])
  expect_true(out$SpotSweeper_artifact[4])
  expect_true("SpotSweeper" %in% names(out@tools))
  expect_equal(out@tools$SpotSweeper$metrics, c("nCount_RNA", "nFeature_RNA", "percent.mito"))
  expect_equal(out@tools$SpotSweeper$directions[["percent.mito"]], "higher")
  expect_true(all(c("local_outliers", "artifacts", "coords", "parameters") %in% names(out@tools$SpotSweeper)))
})

test_that("RunSpotSweeper skips artifact detection when mitochondrial signal is unavailable", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat(include_mito = FALSE)
  with_mock_spotsweeper({
    out <- RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      n_order = 2,
      verbose = FALSE
    )
  })

  expect_true("SpotSweeper_artifact" %in% colnames(out@meta.data))
  expect_false(any(out$SpotSweeper_artifact))
  expect_true(all(out$SpotSweeper_artifact_qc == "Pass"))
  expect_true("skip_reason" %in% colnames(out@tools$SpotSweeper$artifacts))
  expect_true(all(!is.na(out@tools$SpotSweeper$artifacts$skip_reason)))
})

test_that("RunSpotSweeper supports return_filtered and store_results = FALSE", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  with_mock_spotsweeper({
    out <- RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      run_artifact = FALSE,
      return_filtered = TRUE,
      store_results = FALSE,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_equal(colnames(out), c("Spot3", "Spot4", "Spot5", "Spot6"))
  expect_false("SpotSweeper" %in% names(out@tools))
  expect_true(all(out$SpotSweeper_QC == "Pass"))
})

test_that("RunSpotSweeper output is directly plottable with SpatialSpotPlot", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  with_mock_spotsweeper({
    out <- RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      run_artifact = FALSE,
      verbose = FALSE
    )
  })

  p <- SpatialSpotPlot(
    out,
    group.by = "SpotSweeper_QC",
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = "theme_blank"
  )
  expect_s3_class(p, "ggplot")
})
