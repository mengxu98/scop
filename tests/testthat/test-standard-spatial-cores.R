make_standard_spatial_cores_object <- function() {
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

test_that("public spatial workflow forwards cores to its helper", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    run_standard_spatial_workflow = function(...) {
      captured$args <- list(...)
      list(workflow = "spatial")
    },
    .package = "scop"
  )

  out <- RunStandardWorkflow(
    srt = NULL,
    workflow = "spatial",
    cores = 3L,
    verbose = FALSE
  )

  expect_identical(out, list(workflow = "spatial"))
  expect_identical(captured$args$cores, 3L)
})

test_that("spatial helper forwards cores to nested single-cell preprocessing", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) {
      captured$args <- list(...)
      srt
    },
    .package = "scop"
  )

  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  out <- original(
    make_standard_spatial_cores_object(),
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = FALSE,
    do_deconvolution = FALSE,
    cores = 3L,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_identical(captured$args$cores, 3L)
})

test_that("spatial workflow rejects invalid cores before running stages", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  expect_error(
    original(
      make_standard_spatial_cores_object(),
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      cores = 0L,
      verbose = FALSE
    ),
    "cores.*positive integer"
  )
})
