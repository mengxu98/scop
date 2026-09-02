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
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpotQC = function(srt, ...) {
      srt$SpotQC <- "Pass"
      srt
    },
    RunSpatialVariableFeatures = function(srt, set_variable_features, ...) {
      expect_false(set_variable_features)
      srt@tools$SpatialVariableFeatures <- list(result = data.frame(feature = rownames(srt)))
      srt
    }
  )

  input <- make_standard_spatial_stage_object()
  SeuratObject::VariableFeatures(input) <- rownames(input)[1:2]
  out <- suppressWarnings(original(
    input,
    assay = "RNA",
    do_spot_qc = TRUE,
    do_spatial_variable_features = TRUE,
    do_spatial_cluster = FALSE,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = NULL,
    verbose = FALSE
  ))
  workflow <- out@tools$run_standard_spatial_workflow
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
  svf <- workflow$stages[workflow$stages$stage == "spatial_variable_features", , drop = FALSE]
  expect_identical(svf$variable_features_before, 2L)
  expect_identical(svf$variable_features_after, 2L)
  expect_false(svf$set_variable_features)
  expect_identical(SeuratObject::VariableFeatures(out), rownames(input)[1:2])
  deconv <- workflow$stages[workflow$stages$stage == "deconvolution", , drop = FALSE]
  expect_true(deconv$requested)
  expect_identical(deconv$status, "skipped")
  expect_match(deconv$reason, "unavailable")
})

test_that("standard spatial workflow exposes failed stage diagnostics on errors", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
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

test_that("BayesSpace planning failures carry clustering stage diagnostics", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
      srt
    },
    RunBayesSpace = function(...) {
      producer_calls <<- c(producer_calls, "RunBayesSpace")
      stop("BayesSpace producer was reached")
    }
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  cases <- list(
    malformed_cluster = list(
      params = list(cluster_colname = c("one", "two"), init_colname = NULL),
      pattern = "bayesspace_params\\$cluster_colname"
    ),
    malformed_init = list(
      params = list(cluster_colname = "CustomCluster", init_colname = ""),
      pattern = "bayesspace_params\\$init_colname"
    ),
    preprocessing_collision = list(
      params = list(cluster_colname = "Standardclusters", init_colname = NULL),
      pattern = "metadata outputs collide"
    )
  )

  for (case in cases) {
    condition <- tryCatch(
      original(
        srt,
        assay = "RNA",
        do_spot_qc = FALSE,
        do_spatial_variable_features = FALSE,
        do_spatial_cluster = TRUE,
        spatial_q = 2,
        bayesspace_params = case$params,
        do_deconvolution = FALSE,
        verbose = FALSE
      ),
      error = identity
    )
    stages <- attr(condition, "standard_spatial_stages")
    clustering <- stages[
      stages$stage == "spatial_clustering",
      ,
      drop = FALSE
    ]

    expect_s3_class(condition, "error")
    expect_match(conditionMessage(condition), case$pattern)
    expect_s3_class(stages, "data.frame")
    expect_true(clustering$requested)
    expect_identical(clustering$status, "failed")
    expect_identical(clustering$actual_method, "RunBayesSpace")
    expect_false(isTRUE(clustering$result_stored))
    expect_match(clustering$reason, case$pattern)
    expect_length(producer_calls, 0L)
    expect_identical(srt, caller_before)
  }
})

test_that("BayesSpace q inference failures retain clustering stage diagnostics", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  preprocessing_called <- FALSE
  bayesspace_called <- FALSE
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) {
      preprocessing_called <<- TRUE
      srt
    },
    RunBayesSpace = function(...) {
      bayesspace_called <<- TRUE
      stop("BayesSpace producer was reached")
    }
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  condition <- tryCatch(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = NULL,
      bayesspace_params = list(
        cluster_colname = "CustomCluster",
        init_colname = NULL
      ),
      do_deconvolution = FALSE,
      do_normalization = FALSE,
      do_HVF_finding = FALSE,
      do_scaling = FALSE,
      verbose = FALSE
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  clustering <- stages[
    stages$stage == "spatial_clustering",
    ,
    drop = FALSE
  ]

  expect_s3_class(condition, "error")
  expect_identical(clustering$status, "failed")
  expect_identical(clustering$actual_method, "RunBayesSpace")
  expect_match(clustering$reason, "Unable to infer.*spatial_q")
  expect_true(preprocessing_called)
  expect_false(bayesspace_called)
  expect_identical(srt, caller_before)
})

test_that("standard spatial workflow marks deconvolution preflight failures", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
      srt
    },
    RunRCTD = function(...) {
      producer_calls <<- c(producer_calls, "RunRCTD")
      stop("deconvolution backend was reached")
    }
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
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
  expect_length(producer_calls, 0L)
  expect_identical(srt, caller_before)
})

test_that("deconvolution runtime setup failures retain stage diagnostics", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  backend_called <- FALSE
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(...) {
      backend_called <<- TRUE
      stop("deconvolution backend was reached")
    },
    standard_spatial_clear_outputs = function(...) {
      stop("synthetic deconvolution cleanup failure")
    }
  )
  srt <- make_standard_spatial_stage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  caller_before <- srt
  condition <- tryCatch(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = TRUE,
      reference = srt,
      reference_label = "label",
      do_normalization = FALSE,
      do_HVF_finding = FALSE,
      do_scaling = FALSE,
      verbose = FALSE
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(deconv$status, "failed")
  expect_identical(deconv$actual_method, "RunRCTD")
  expect_match(deconv$reason, "synthetic deconvolution cleanup failure")
  expect_false(backend_called)
  expect_identical(srt, caller_before)
})
