make_spatial_deconvolution_plot_object <- function() {
  counts <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- SeuratObject::CreateSeuratObject(counts)
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  srt
}

add_spatial_deconvolution_result <- function(srt, key) {
  proportions <- matrix(
    c(
      0.8, 0.2,
      0.6, 0.4,
      0.3, 0.7,
      0.1, 0.9
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(colnames(srt), c("Alpha", "Beta"))
  )
  srt@tools[[key]] <- list(
    proportions = proportions,
    cells = colnames(srt),
    coordinate_contract_version = 2L,
    parameters = list(
      coordinate_space = "raw",
      coordinate_contract_version = 2L
    )
  )
  srt
}

test_that("SpatialDeconvolutionPlot consumes each compatible stored family", {
  keys <- c("RCTD", "CARD", "SPOTlight", "SpatialDWLS")
  for (key in keys) {
    srt <- add_spatial_deconvolution_result(
      make_spatial_deconvolution_plot_object(),
      key = paste0(key, "Custom")
    )
    plots <- SpatialDeconvolutionPlot(
      srt,
      tool_name = paste0(key, "Custom"),
      combine = FALSE,
      overlay_image = FALSE
    )
    expect_named(plots, c("Alpha", "Beta"))
    expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
    expect_true(all(vapply(
      plots,
      function(plot) identical(plot$scales$scales[[1L]]$limits, c(0, 1)),
      logical(1)
    )))
    expect_s3_class(SpatialDeconvolutionPlot(
      srt,
      tool_name = paste0(key, "Custom"),
      plot_type = "dominant",
      overlay_image = FALSE
    ), "ggplot")
  }
})

test_that("SpatialDeconvolutionPlot requires an explicit result and supports layouts", {
  srt <- add_spatial_deconvolution_result(
    make_spatial_deconvolution_plot_object(),
    "CARD"
  )
  expect_error(
    SpatialDeconvolutionPlot(srt),
    "tool_name.*explicit"
  )
  expect_s3_class(SpatialDeconvolutionPlot(
    srt,
    tool_name = "CARD",
    overlay_image = FALSE,
    nrow = 1,
    byrow = FALSE
  ), "patchwork")
  expect_s3_class(SpatialDeconvolutionPlot(
    srt,
    tool_name = "CARD",
    cell_types = "Alpha",
    overlay_image = FALSE
  ), "ggplot")
  expect_error(
    SpatialDeconvolutionPlot(srt, tool_name = "CARD", cell_types = "Missing"),
    "Unknown.*cell_types"
  )
  srt <- add_spatial_deconvolution_result(srt, "RCTD")
  expect_s3_class(SpatialDeconvolutionPlot(
    srt,
    tool_name = "RCTD",
    cell_types = "Beta",
    overlay_image = FALSE
  ), "ggplot")
})

test_that("SpatialDeconvolutionPlot rejects empty, partial, stale, and wrong results", {
  srt <- make_spatial_deconvolution_plot_object()
  srt <- add_spatial_deconvolution_result(srt, "CustomCARD")

  bad <- srt
  bad@tools$CustomCARD$proportions <- NULL
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "does not contain")

  bad <- srt
  bad@tools$CustomCARD$proportions <- bad@tools$CustomCARD$proportions[-1, , drop = FALSE]
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "stale or incomplete")

  bad <- srt
  rownames(bad@tools$CustomCARD$proportions)[[1L]] <- "UnknownSpot"
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "stale or incomplete")

  bad <- srt
  bad@tools$CustomCARD$proportions <- bad@tools$CustomCARD$proportions[, FALSE, drop = FALSE]
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "empty")

  bad <- srt
  bad@tools$CustomCARD$proportions[[1L]] <- 1.2
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "between zero and one")

  bad <- srt
  bad@tools$CustomCARD$proportions <- apply(
    bad@tools$CustomCARD$proportions,
    2,
    as.character
  )
  expect_error(SpatialDeconvolutionPlot(bad, "CustomCARD"), "must be numeric")

  expect_error(SpatialDeconvolutionPlot(srt, "MissingKey"), "not a spatial deconvolution result")
})

test_that("SpatialDeconvolutionPlot rejects old coordinate-dependent results", {
  srt <- add_spatial_deconvolution_result(
    make_spatial_deconvolution_plot_object(),
    "OldRCTD"
  )
  srt@tools$OldRCTD$coordinate_contract_version <- NULL
  srt@tools$OldRCTD$parameters$coordinate_contract_version <- NULL
  srt@tools$OldRCTD$backend_api <- "mock"
  expect_error(
    SpatialDeconvolutionPlot(srt, "OldRCTD"),
    "rerun.*OldRCTD"
  )
})

test_that("SpatialDeconvolutionPlot accepts legacy coordinate-independent SpatialDWLS results", {
  srt <- add_spatial_deconvolution_result(
    make_spatial_deconvolution_plot_object(),
    "OldSpatialDWLS"
  )
  srt@tools$OldSpatialDWLS$coordinate_contract_version <- NULL
  srt@tools$OldSpatialDWLS$parameters$coordinate_contract_version <- NULL
  srt@tools$OldSpatialDWLS$parameters$coordinate_dependent <- FALSE
  srt@tools$OldSpatialDWLS$coords <- data.frame(
    x = srt$col,
    y = srt$row,
    row.names = colnames(srt)
  )
  expect_s3_class(
    SpatialDeconvolutionPlot(
      srt,
      tool_name = "OldSpatialDWLS",
      overlay_image = FALSE
    ),
    "patchwork"
  )
})

test_that("RunRCTD writes custom-key plot-ready proportions", {
  srt <- make_spatial_deconvolution_plot_object()
  reference_counts <- matrix(
    c(8, 6, 1, 1, 1, 1, 7, 9, 5, 4, 2, 1),
    nrow = 3,
    dimnames = list(rownames(srt), paste0("Cell", 1:4))
  )
  reference_counts <- methods::as(
    Matrix::Matrix(reference_counts, sparse = TRUE),
    "dgCMatrix"
  )
  reference <- SeuratObject::CreateSeuratObject(reference_counts)
  reference$celltype <- c("Alpha", "Alpha", "Beta", "Beta")
  testthat::local_mocked_bindings(
    rctd_run_spacexr = function(st_counts, ...) {
      list(
        weights = matrix(
          rep(c(0.7, 0.3), each = ncol(st_counts)),
          nrow = ncol(st_counts),
          dimnames = list(colnames(st_counts), c("Alpha", "Beta"))
        ),
        metadata = list(),
        api = "mock",
        object = list()
      )
    },
    .package = "scop"
  )
  out <- RunRCTD(
    srt,
    reference,
    reference_label = "celltype",
    min_cells = 1,
    tool_name = "RCTDCustom",
    coord.cols = c("col", "row"),
    verbose = FALSE
  )
  expect_identical(out@tools$RCTDCustom$parameters$coordinate_space, "raw")
  expect_false("RCTD" %in% names(out@tools))
  expect_identical(rownames(out@tools$RCTDCustom$proportions), colnames(out))
  expect_s3_class(SpatialDeconvolutionPlot(
    out,
    tool_name = "RCTDCustom",
    cell_types = "Alpha",
    overlay_image = FALSE
  ), "ggplot")
})
