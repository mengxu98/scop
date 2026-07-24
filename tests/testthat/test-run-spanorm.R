make_spanorm_seurat <- function() {
  counts <- matrix(
    c(
      5, 4, 0, 0, 0,
      0, 0, 4, 5, 4,
      1, 1, 1, 1, 1,
      3, 0, 3, 0, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:5))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$col <- c(1, 2, 3, 4, 5)
  srt$row <- c(1, 1, 2, 2, 3)
  srt
}

with_mock_spanorm <- function(code) {
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      if (identical(package, "SpatialExperiment")) {
        return(getExportedValue(package, name))
      }
      expect_identical(package, "SpaNorm")
      expect_identical(name, "SpaNorm")
      function(spe, sample.p = NULL, ...) {
        expect_equal(sample.p, 0.5)
        counts <- SummarizedExperiment::assay(spe, "counts")
        SummarizedExperiment::assay(spe, "logcounts") <- log1p(counts) + 1
        spe
      }
    }
  )
  force(code)
}

test_that("RunSpaNorm creates a new assay and stores normalized schema", {
  skip_if_not_installed("SpatialExperiment")
  srt <- make_spanorm_seurat()
  original_counts <- GetAssayData5(srt, assay = "RNA", layer = "counts")

  with_mock_spanorm({
    out <- RunSpaNorm(
      srt,
      assay = "RNA",
      layer = "counts",
      coord.cols = c("col", "row"),
      new_assay = "SpaNormMock",
      tool_name = "SpaNormMock",
      store_spe = TRUE,
      sample.p = 0.5,
      verbose = FALSE
    )
  })

  expect_true("SpaNormMock" %in% SeuratObject::Assays(out))
  normalized <- GetAssayData5(out, assay = "SpaNormMock", layer = "data")
  expect_equal(as.matrix(normalized), as.matrix(log1p(original_counts) + 1))
  expect_equal(
    as.matrix(GetAssayData5(out, assay = "RNA", layer = "counts")),
    as.matrix(original_counts)
  )
  expect_true("SpaNormMock" %in% names(out@tools))
  expect_equal(out@tools$SpaNormMock$parameters$assay, "RNA")
  expect_equal(out@tools$SpaNormMock$parameters$layer, "counts")
  expect_null(out@tools$SpaNormMock$parameters$image)
  expect_equal(out@tools$SpaNormMock$parameters$coord.cols, c("col", "row"))
  expect_identical(out@tools$SpaNormMock$source$coordinate_space, "raw")
  expect_equal(out@tools$SpaNormMock$parameters$new_assay, "SpaNormMock")
  expect_equal(out@tools$SpaNormMock$parameters$backend_args$sample.p, 0.5)
  expect_equal(out@tools$SpaNormMock$features, rownames(original_counts))
  expect_equal(out@tools$SpaNormMock$cells, colnames(original_counts))
  expect_s4_class(out@tools$SpaNormMock$spe, "SpatialExperiment")
})

test_that("RunSpaNorm can skip tool storage", {
  skip_if_not_installed("SpatialExperiment")
  srt <- make_spanorm_seurat()
  with_mock_spanorm({
    out <- RunSpaNorm(
      srt,
      assay = "RNA",
      layer = "counts",
      coord.cols = c("col", "row"),
      store_results = FALSE,
      sample.p = 0.5,
      verbose = FALSE
    )
  })

  expect_true("SpaNorm" %in% SeuratObject::Assays(out))
  expect_null(out@tools$SpaNorm)
})

test_that("RunSpaNorm validates coordinates and backend output clearly", {
  skip_if_not_installed("SpatialExperiment")
  srt <- make_spanorm_seurat()
  expect_error(
    RunSpaNorm(
      srt,
      assay = "RNA",
      coord.cols = c("missing_x", "missing_y"),
      verbose = FALSE
    ),
    "Spatial coordinates"
  )

  testthat::local_mocked_bindings(
    check_r = function(packages, ...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      if (identical(package, "SpatialExperiment")) {
        return(getExportedValue(package, name))
      }
      function(spe, ...) spe
    }
  )
  expect_error(
    RunSpaNorm(
      srt,
      assay = "RNA",
      coord.cols = c("col", "row"),
      verbose = FALSE
    ),
    "logcounts"
  )
})

test_that("spatial wrapper exports remain available", {
  expect_true(is.function(RunSpaNorm))
  expect_true(is.function(RunSmoothClust))
})
