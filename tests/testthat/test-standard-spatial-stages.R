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

test_that("top-level spatial stage controls require logical scalars", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunSpotQC = function(...) producer_calls <<- c(producer_calls, "RunSpotQC"),
    RunStandardWorkflow = function(...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
    },
    RunSpatialVariableFeatures = function(...) {
      producer_calls <<- c(producer_calls, "RunSpatialVariableFeatures")
    },
    RunBayesSpace = function(...) producer_calls <<- c(producer_calls, "RunBayesSpace"),
    RunRCTD = function(...) producer_calls <<- c(producer_calls, "RunRCTD")
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  cases <- list(
    do_spot_qc = "yes",
    do_spatial_variable_features = 1,
    do_spatial_cluster = c(TRUE, FALSE),
    do_deconvolution = NA
  )

  for (arg in names(cases)) {
    args <- list(
      srt = srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      verbose = FALSE
    )
    args[[arg]] <- cases[[arg]]
    expect_error(
      do.call(original, args),
      paste0(arg, ".*TRUE or FALSE")
    )
    expect_length(producer_calls, 0L)
    expect_identical(srt, caller_before)
  }
})

test_that("requested stage parameter lists fail through their stage wrapper", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunSpotQC = function(...) producer_calls <<- c(producer_calls, "RunSpotQC"),
    RunStandardWorkflow = function(...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
    },
    RunSpatialVariableFeatures = function(...) {
      producer_calls <<- c(producer_calls, "RunSpatialVariableFeatures")
    },
    RunBayesSpace = function(...) producer_calls <<- c(producer_calls, "RunBayesSpace")
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  cases <- list(
    quality_control = list(
      flag = "do_spot_qc",
      params = "spot_qc_params",
      actual_method = "RunSpotQC"
    ),
    spatial_variable_features = list(
      flag = "do_spatial_variable_features",
      params = "spatial_variable_features_params",
      actual_method = "RunSpatialVariableFeatures"
    ),
    spatial_clustering = list(
      flag = "do_spatial_cluster",
      params = "bayesspace_params",
      actual_method = "RunBayesSpace"
    )
  )

  for (stage in names(cases)) {
    case <- cases[[stage]]
    args <- list(
      srt = srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      verbose = FALSE
    )
    args[[case$flag]] <- TRUE
    args[[case$params]] <- "invalid"
    condition <- tryCatch(do.call(original, args), error = identity)
    stages <- attr(condition, "standard_spatial_stages")
    failed <- stages[stages$stage == stage, , drop = FALSE]

    expect_s3_class(condition, "error")
    expect_match(conditionMessage(condition), case$params)
    expect_s3_class(stages, "data.frame")
    expect_identical(failed$status, "failed")
    expect_identical(failed$actual_method, case$actual_method)
    expect_false(isTRUE(failed$result_stored))
    expect_length(producer_calls, 0L)
    expect_identical(srt, caller_before)
  }
})

test_that("requested stage methods fail through their stage wrapper", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
    },
    RunBayesSpace = function(...) producer_calls <<- c(producer_calls, "RunBayesSpace"),
    RunRCTD = function(...) producer_calls <<- c(producer_calls, "RunRCTD"),
    RunSPOTlight = function(...) producer_calls <<- c(producer_calls, "RunSPOTlight"),
    RunCell2location = function(...) producer_calls <<- c(producer_calls, "RunCell2location")
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  cases <- list(
    spatial_clustering = list(
      flag = "do_spatial_cluster",
      value = TRUE,
      method_arg = "spatial_cluster_method",
      actual_method = "RunBayesSpace"
    ),
    deconvolution = list(
      flag = "do_deconvolution",
      value = TRUE,
      method_arg = "deconvolution_method",
      actual_method = "deconvolution_method validation"
    ),
    automatic_deconvolution = list(
      stage = "deconvolution",
      flag = "do_deconvolution",
      value = NULL,
      method_arg = "deconvolution_method",
      actual_method = "deconvolution_method validation"
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    stage <- if (is.null(case$stage)) case_name else case$stage
    args <- list(
      srt = srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      verbose = FALSE
    )
    args[case$flag] <- list(case$value)
    args[[case$method_arg]] <- "invalid"
    condition <- tryCatch(do.call(original, args), error = identity)
    stages <- attr(condition, "standard_spatial_stages")
    failed <- stages[stages$stage == stage, , drop = FALSE]

    expect_s3_class(condition, "error")
    expect_match(conditionMessage(condition), "should be")
    expect_s3_class(stages, "data.frame")
    expect_true(failed$requested)
    expect_identical(failed$status, "failed")
    expect_identical(failed$requested_method, "invalid")
    expect_identical(failed$actual_method, case$actual_method)
    expect_length(producer_calls, 0L)
    expect_identical(srt, caller_before)
  }
})

test_that("parameters for unrequested spatial stages are ignored", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  preprocessing_calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) {
      preprocessing_calls <<- preprocessing_calls + 1L
      srt
    }
  )
  out <- original(
    make_standard_spatial_stage_object(),
    assay = "RNA",
    do_spot_qc = FALSE,
    spot_qc_params = "ignored",
    do_spatial_variable_features = FALSE,
    spatial_variable_features_params = "ignored",
    do_spatial_cluster = FALSE,
    spatial_cluster_method = "ignored",
    bayesspace_params = "ignored",
    do_deconvolution = FALSE,
    deconvolution_method = "ignored",
    deconvolution_params = "ignored",
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )

  expect_identical(preprocessing_calls, 1L)
  expect_identical(
    out@tools$run_standard_spatial_workflow$status,
    "completed"
  )
  expect_identical(
    out@tools$run_standard_spatial_workflow$stages$requested_method,
    c("RunSpotQC", "RunSpatialVariableFeatures", "BayesSpace", "RCTD")
  )
})

test_that("empty and NULL spatial prefixes normalize to the empty string", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  plan_targets <- getFromNamespace(
    "standard_spatial_preprocessing_metadata_targets",
    "scop"
  )
  expected_targets <- c(
    "pcaclusters",
    "clusters",
    "pca_SNN_res.0.6",
    "ident"
  )
  expect_identical(
    plan_targets("", "pca", "LogNormalize"),
    expected_targets
  )
  expect_identical(
    plan_targets(NULL, "pca", "LogNormalize"),
    expected_targets
  )
  for (invalid_prefix in list(1, c("one", "two"), NA_character_)) {
    expect_error(
      plan_targets(invalid_prefix, "pca", "LogNormalize"),
      "prefix.*single character string"
    )
  }

  nested_prefixes <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, prefix, ...) {
      nested_prefixes <<- c(nested_prefixes, prefix)
      srt
    }
  )
  for (prefix_value in list("", NULL)) {
    out <- original(
      make_standard_spatial_stage_object(),
      prefix = prefix_value,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      do_normalization = FALSE,
      do_HVF_finding = FALSE,
      do_scaling = FALSE,
      verbose = FALSE
    )
    expect_identical(
      out@tools$run_standard_spatial_workflow$parameters$prefix,
      ""
    )
  }
  expect_identical(nested_prefixes, c("", ""))
})

test_that("duplicate SpotQC rules do not create a planning collision", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  duplicate_rules <- rep("log10_nCount:lower:3", 2L)
  captured <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpotQC = function(srt, outlier_threshold, outlier_n, ...) {
      captured <<- list(
        outlier_threshold = outlier_threshold,
        outlier_n = outlier_n
      )
      srt$SpotQC <- rep("Pass", ncol(srt))
      srt
    }
  )

  out <- original(
    make_standard_spatial_stage_object(),
    assay = "RNA",
    do_spot_qc = TRUE,
    spot_qc_params = list(
      qc_metrics = "outlier",
      outlier_threshold = duplicate_rules,
      outlier_n = 2
    ),
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = FALSE,
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )
  qc <- out@tools$run_standard_spatial_workflow$stages
  qc <- qc[qc$stage == "quality_control", , drop = FALSE]

  expect_identical(captured$outlier_threshold, duplicate_rules)
  expect_identical(captured$outlier_n, 2)
  expect_identical(qc$status, "completed")
  expect_identical(unique(out$SpotQC), "Pass")
})

test_that("distinct SpotQC rules cannot normalize to the same output", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunSpotQC = function(...) {
      producer_calls <<- c(producer_calls, "RunSpotQC")
    },
    RunStandardWorkflow = function(...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
    }
  )
  srt <- make_standard_spatial_stage_object()
  srt[["foo-bar"]] <- seq_len(ncol(srt))
  srt[["foo.bar"]] <- rev(seq_len(ncol(srt)))
  caller_before <- srt

  condition <- tryCatch(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = TRUE,
      spot_qc_params = list(
        qc_metrics = "outlier",
        outlier_threshold = c("foo-bar:lower:1", "foo.bar:lower:1")
      ),
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      verbose = FALSE
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  qc <- stages[stages$stage == "quality_control", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "metadata outputs collide")
  expect_s3_class(stages, "data.frame")
  expect_identical(qc$status, "failed")
  expect_identical(qc$actual_method, "RunSpotQC")
  expect_false(isTRUE(qc$result_stored))
  expect_length(producer_calls, 0L)
  expect_identical(srt, caller_before)
})

test_that("effective SVF storage controls require logical scalars", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_called <- FALSE
  testthat::local_mocked_bindings(
    .package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(...) {
      producer_called <<- TRUE
      stop("SVF producer was reached")
    }
  )
  srt <- make_standard_spatial_stage_object()
  caller_before <- srt
  cases <- list(
    set_variable_features = "yes",
    store_results = c(TRUE, FALSE)
  )

  for (arg in names(cases)) {
    params <- list()
    params[[arg]] <- cases[[arg]]
    condition <- tryCatch(
      original(
        srt,
        assay = "RNA",
        do_spot_qc = FALSE,
        do_spatial_variable_features = TRUE,
        spatial_variable_features_params = params,
        do_spatial_cluster = FALSE,
        do_deconvolution = FALSE,
        do_normalization = FALSE,
        do_HVF_finding = FALSE,
        do_scaling = FALSE,
        verbose = FALSE
      ),
      error = identity
    )
    stages <- attr(condition, "standard_spatial_stages")
    svf <- stages[
      stages$stage == "spatial_variable_features",
      ,
      drop = FALSE
    ]

    expect_s3_class(condition, "error")
    expect_match(conditionMessage(condition), paste0("spatial_variable_features_params\\$", arg))
    expect_identical(svf$status, "failed")
    expect_identical(svf$actual_method, "RunSpatialVariableFeatures")
    expect_false(producer_called)
    expect_identical(srt, caller_before)
  }
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

test_that("SpotQC planning collisions carry quality-control diagnostics", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  producer_calls <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    RunSpotQC = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunSpotQC")
      srt
    },
    RunStandardWorkflow = function(srt, ...) {
      producer_calls <<- c(producer_calls, "RunStandardWorkflow")
      srt
    }
  )
  srt <- make_standard_spatial_stage_object()
  counts <- SeuratObject::LayerData(srt, assay = "RNA", layer = "counts")
  srt[["pcaclusters"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  caller_before <- srt
  condition <- tryCatch(
    original(
      srt,
      prefix = "nCount_",
      assay = "pcaclusters",
      do_spot_qc = TRUE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      do_deconvolution = FALSE,
      linear_reduction = "pca",
      verbose = FALSE
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  qc <- stages[stages$stage == "quality_control", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "metadata outputs collide")
  expect_s3_class(stages, "data.frame")
  expect_true(qc$requested)
  expect_identical(qc$status, "failed")
  expect_identical(qc$actual_method, "RunSpotQC")
  expect_false(isTRUE(qc$result_stored))
  expect_match(qc$reason, "metadata outputs collide")
  expect_length(producer_calls, 0L)
  expect_identical(srt, caller_before)
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
    non_character_cluster = list(
      params = list(cluster_colname = 1, init_colname = NULL),
      pattern = "bayesspace_params\\$cluster_colname"
    ),
    malformed_init = list(
      params = list(cluster_colname = "CustomCluster", init_colname = ""),
      pattern = "bayesspace_params\\$init_colname"
    ),
    non_character_init = list(
      params = list(cluster_colname = "CustomCluster", init_colname = 1),
      pattern = "bayesspace_params\\$init_colname"
    ),
    preprocessing_collision = list(
      params = list(cluster_colname = "Standardclusters", init_colname = NULL),
      pattern = "metadata outputs collide"
    ),
    retained_ident_cluster = list(
      params = list(cluster_colname = "ident", init_colname = NULL),
      pattern = "metadata outputs collide"
    ),
    retained_ident_init = list(
      params = list(cluster_colname = "CustomCluster", init_colname = "ident"),
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
