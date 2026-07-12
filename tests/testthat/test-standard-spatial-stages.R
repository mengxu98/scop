make_standard_spatial_stage_object <- function() {
  counts <- matrix(
    c(3, 1, 0, 2, 0, 4, 1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  srt
}

test_that("standard spatial workflow records completed and skipped stages truthfully", {
  original <- getFromNamespace("standard_spatial_scop", "scop")
  testthat::local_mocked_bindings(
    .package = "scop",
    standard_scop = function(srt, ...) srt,
    RunSpotQC = function(srt, ...) {
      srt$SpotQC <- "Pass"
      srt
    },
    RunSpatialVariableFeatures = function(srt, ...) {
      srt@tools$SpatialVariableFeatures <- list(result = data.frame(feature = rownames(srt)))
      srt
    }
  )

  out <- suppressWarnings(original(
    make_standard_spatial_stage_object(),
    assay = "RNA",
    do_spot_qc = TRUE,
    do_spatial_variable_features = TRUE,
    do_spatial_cluster = FALSE,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = NULL,
    verbose = FALSE
  ))
  workflow <- out@tools$standard_spatial_scop
  expect_identical(workflow$status, "partial")
  expect_identical(
    workflow$stages$status,
    c("completed", "completed", "skipped", "skipped")
  )
  expect_identical(
    workflow$stages$actual_method[workflow$stages$stage == "quality_control"],
    "RunSpotQC"
  )
  expect_identical(
    workflow$stages$result_tool_key[workflow$stages$stage == "spatial_variable_features"],
    "SpatialVariableFeatures"
  )
  deconv <- workflow$stages[workflow$stages$stage == "deconvolution", , drop = FALSE]
  expect_true(deconv$requested)
  expect_identical(deconv$status, "skipped")
  expect_match(deconv$reason, "unavailable")
})

test_that("standard spatial workflow exposes failed stage diagnostics on errors", {
  original <- getFromNamespace("standard_spatial_scop", "scop")
  testthat::local_mocked_bindings(
    .package = "scop",
    standard_scop = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(...) stop("synthetic SVF failure")
  )
  condition <- tryCatch(
    original(
      make_standard_spatial_stage_object(),
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = TRUE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      verbose = FALSE
    ),
    error = identity
  )
  expect_s3_class(condition, "error")
  stages <- attr(condition, "standard_spatial_stages")
  expect_s3_class(stages, "data.frame")
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]
  expect_identical(svf$status, "failed")
  expect_identical(svf$actual_method, "RunSpatialVariableFeatures")
  expect_match(svf$reason, "synthetic SVF failure")
})

test_that("standard spatial workflow marks deconvolution preflight failures", {
  original <- getFromNamespace("standard_spatial_scop", "scop")
  testthat::local_mocked_bindings(
    .package = "scop",
    standard_scop = function(srt, ...) srt
  )
  srt <- make_standard_spatial_stage_object()
  condition <- tryCatch(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = TRUE,
      deconvolution_method = "RCTD",
      reference = srt,
      reference_label = NULL,
      verbose = FALSE
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]
  expect_identical(deconv$status, "failed")
  expect_identical(deconv$actual_method, "RunRCTD")
  expect_match(deconv$reason, "reference_label")
})
