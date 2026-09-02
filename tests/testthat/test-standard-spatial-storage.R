make_standard_spatial_storage_object <- function() {
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

make_standard_spatial_rerun_object <- function() {
  counts <- matrix(
    c(
      5, 1, 0, 2, 3, 0, 4, 1, 2,
      1, 3, 2, 0, 5, 2, 1, 4, 0,
      3, 1, 2, 5, 0, 2, 1, 4, 3,
      0, 2, 5, 1, 3, 4, 2, 0, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:4), paste0("spot", 1:9))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- rep(0:2, 3)
  srt$row <- rep(0:2, each = 3)
  srt
}

run_storage_workflow <- function(
  srt,
  do_spatial_variable_features = FALSE,
  spatial_variable_features_params = list(),
  do_deconvolution = FALSE,
  deconvolution_method = "RCTD",
  reference = NULL,
  reference_label = NULL,
  deconvolution_params = list()
) {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  original(
    srt,
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = do_spatial_variable_features,
    spatial_variable_features_params = spatial_variable_features_params,
    do_spatial_cluster = FALSE,
    do_deconvolution = do_deconvolution,
    deconvolution_method = deconvolution_method,
    reference = reference,
    reference_label = reference_label,
    deconvolution_params = deconvolution_params,
    verbose = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE
  )
}

add_mock_deconv_outputs <- function(
  srt,
  prefix,
  tool_name = NULL,
  store_results = FALSE
) {
  srt[[paste0(prefix, "_prop_typeA")]] <- rep(0.5, ncol(srt))
  srt[[paste0(prefix, "_dominant_type")]] <- rep("typeA", ncol(srt))
  srt[[paste0(prefix, "_max_prop")]] <- rep(0.5, ncol(srt))
  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      proportions = matrix(
        0.5,
        nrow = ncol(srt),
        ncol = 1,
        dimnames = list(colnames(srt), "typeA")
      )
    )
  }
  srt
}

make_valid_empty_svf_tool <- function() {
  list(
    result = data.frame(
      feature = c("gene1", "gene2"),
      rank = 1:2,
      method = rep("moran", 2),
      statistic = c(NA_real_, NA_real_),
      score = c(NA_real_, NaN),
      p_value = c(NA_real_, NA_real_),
      q_value = c(NA_real_, NA_real_),
      mean = c(1, 2),
      variance = c(0, 0),
      n_spots = c(4L, 4L)
    ),
    summary = list(
      n_features = 2L,
      top_features = character(),
      top_feature_summary = data.frame(
        feature = character(),
        rank = integer(),
        score = numeric()
      )
    ),
    parameters = list(set_variable_features = TRUE)
  )
}

test_that("stored spatial variable feature results are reported", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_false(set_variable_features)
      expect_true(store_results)
      expect_null(srt@tools[["SpatialVariableFeatures"]])
      srt@tools[["SpatialVariableFeatures"]] <- list(
        result = data.frame(feature = "gene1")
      )
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    make_standard_spatial_storage_object(),
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(store_results = TRUE)
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_true(isTRUE(svf$result_stored))
  expect_identical(svf$result_tool_key, "SpatialVariableFeatures")
  expect_match(svf$result_location, "tools")
})

test_that("storage-off SVF ignores and preserves an older tool result", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_false(set_variable_features)
      expect_false(store_results)
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  stale_tool <- list(result = data.frame(feature = "stale"))
  srt@tools[["SpatialVariableFeatures"]] <- stale_tool

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = FALSE,
      set_variable_features = FALSE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_false(isTRUE(svf$result_stored))
  expect_true(is.na(svf$result_tool_key))
  expect_true(is.na(svf$result_location))
  expect_identical(out@tools[["SpatialVariableFeatures"]], stale_tool)
  expect_identical(out@tools$run_standard_spatial_workflow$status, "completed")
})

test_that("variable-feature-only output is certified without a tool", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_true(set_variable_features)
      expect_false(store_results)
      top_features <- c("gene1", "gene2")
      SeuratObject::VariableFeatures(srt, assay = assay) <- top_features
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  stale_tool <- list(result = data.frame(feature = "stale"))
  srt@tools[["SpatialVariableFeatures"]] <- stale_tool

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = FALSE,
      set_variable_features = TRUE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_true(isTRUE(svf$result_stored))
  expect_true(is.na(svf$result_tool_key))
  expect_match(svf$result_location, "VariableFeatures")
  expect_identical(
    SeuratObject::VariableFeatures(out, assay = "RNA"),
    c("gene1", "gene2")
  )
  expect_identical(out@tools[["SpatialVariableFeatures"]], stale_tool)
})

test_that("SVF bookkeeping follows the effective assay override", {
  srt <- make_standard_spatial_storage_object()
  alt_counts <- SeuratObject::LayerData(srt, assay = "RNA", layer = "counts")
  srt[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt_counts)
  SeuratObject::VariableFeatures(srt, assay = "RNA") <- "gene1"
  SeuratObject::VariableFeatures(srt, assay = "ALT") <- c("gene2", "gene3")

  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      features,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_identical(assay, "ALT")
      expect_identical(features, c("gene2", "gene3"))
      expect_length(
        suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
        0L
      )
      expect_identical(
        standard_spatial_variable_features(srt, assay = "RNA"),
        "gene1"
      )
      expect_true(set_variable_features)
      expect_false(store_results)
      SeuratObject::VariableFeatures(srt, assay = assay) <- "gene3"
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      assay = "ALT",
      store_results = FALSE,
      set_variable_features = TRUE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_identical(svf$variable_features_before, 2L)
  expect_identical(svf$variable_features_after, 1L)
  expect_match(svf$result_location, 'VariableFeatures\\("ALT"\\)')
  expect_identical(
    standard_spatial_variable_features(out, assay = "RNA"),
    "gene1"
  )
  expect_identical(
    standard_spatial_variable_features(out, assay = "ALT"),
    "gene3"
  )
})

test_that("NULL SVF assay resolves to the producer default assay", {
  srt <- make_standard_spatial_storage_object()
  alt_counts <- SeuratObject::LayerData(srt, assay = "RNA", layer = "counts")
  srt[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt_counts)
  SeuratObject::DefaultAssay(srt) <- "ALT"
  SeuratObject::VariableFeatures(srt, assay = "ALT") <- "gene2"

  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(srt, assay, features, ...) {
      expect_identical(assay, "ALT")
      expect_identical(features, "gene2")
      SeuratObject::VariableFeatures(srt, assay = assay) <- "gene3"
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      assay = NULL,
      store_results = FALSE,
      set_variable_features = TRUE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$variable_features_before, 1L)
  expect_identical(svf$variable_features_after, 1L)
  expect_match(svf$result_location, 'VariableFeatures\\("ALT"\\)')
})

test_that("invalid effective SVF assays fail through the stage wrapper", {
  producer_calls <- 0L
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(...) {
      producer_calls <<- producer_calls + 1L
      stop("SVF producer was reached")
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  SeuratObject::VariableFeatures(srt, assay = "RNA") <- c("gene1", "gene2")
  srt@tools[["SpatialVariableFeatures"]] <- list(result = "stale")
  caller_before <- srt
  invalid_assays <- list(
    missing = list(value = "missing_assay", reason = "not present"),
    malformed = list(
      value = c("RNA", "another_assay"),
      reason = "single non-empty\\s+string"
    )
  )

  for (case_name in names(invalid_assays)) {
    case <- invalid_assays[[case_name]]
    condition <- tryCatch(
      run_storage_workflow(
        srt,
        do_spatial_variable_features = TRUE,
        spatial_variable_features_params = list(assay = case$value)
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
    expect_s3_class(stages, "data.frame")
    expect_identical(svf$status, "failed", info = case_name)
    expect_identical(
      svf$actual_method,
      "RunSpatialVariableFeatures",
      info = case_name
    )
    expect_false(isTRUE(svf$result_stored), info = case_name)
    expect_true(is.na(svf$result_location), info = case_name)
    expect_match(svf$reason, case$reason, info = case_name)
    expect_identical(srt, caller_before, info = case_name)
  }

  expect_identical(producer_calls, 0L)
})

test_that("SVF output replaces stale selection and preserves HVF metadata", {
  srt <- make_standard_spatial_storage_object()
  assay_object <- srt[["RNA"]]
  skip_if_not(inherits(assay_object, "StdAssay"))
  variable <- c(TRUE, FALSE, TRUE)
  rank <- c(1L, NA_integer_, 2L)
  names(variable) <- names(rank) <- rownames(assay_object)
  assay_object[["vf_mock_counts_variable"]] <- variable
  assay_object[["vf_mock_counts_rank"]] <- rank
  srt[["RNA"]] <- assay_object
  SeuratObject::VariableFeatures(srt, assay = "RNA") <- c("gene1", "gene2")
  hvf_before <- srt[["RNA"]][[]][
    ,
    c("vf_mock_counts_variable", "vf_mock_counts_rank"),
    drop = FALSE
  ]

  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      features,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_identical(features, c("gene1", "gene2"))
      expect_true(set_variable_features)
      expect_false(store_results)
      expect_false(any(startsWith(colnames(srt[[assay]][[]]), "vf_")))
      expect_length(
        suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
        0L
      )
      SeuratObject::VariableFeatures(srt, assay = assay) <- c("gene2", "gene3")
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = FALSE,
      set_variable_features = TRUE
    )
  )
  hvf_after <- out[["RNA"]][[]][
    ,
    c("vf_mock_counts_variable", "vf_mock_counts_rank"),
    drop = FALSE
  ]
  output_features <- suppressWarnings(
    SeuratObject::VariableFeatures(out, assay = "RNA")
  )
  output_features <- output_features[!is.na(output_features)]

  expect_identical(output_features, c("gene2", "gene3"))
  expect_identical(hvf_after, hvf_before)
})

test_that("stale variable features do not satisfy a quiet producer", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      features,
      store_results,
      ...
    ) {
      expect_false(store_results)
      expect_null(srt@tools[["SpatialVariableFeatures"]])
      expect_identical(features, "gene1")
      expect_length(
        suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
        0L
      )
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  SeuratObject::VariableFeatures(srt, assay = "RNA") <- "gene1"

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_spatial_variable_features = TRUE,
      spatial_variable_features_params = list(
        store_results = FALSE,
        set_variable_features = TRUE
      )
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(svf$status, "failed")
  expect_false(isTRUE(svf$result_stored))
  expect_match(svf$reason, "VariableFeatures")
  expect_identical(
    suppressWarnings(SeuratObject::VariableFeatures(srt, assay = "RNA"))[[1L]],
    "gene1"
  )
})

test_that("native spatial variable feature workflow can be rerun unchanged", {
  skip_if_not_installed("BiocNeighbors")
  public_workflow <- RunStandardWorkflow
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    .package = "scop"
  )
  args <- list(
    workflow = "spatial",
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      layer = "counts",
      coord.cols = c("col", "row"),
      k = 2,
      nfeatures = 3,
      min_spots = 1,
      nperm = 0,
      backend = "r",
      store_results = TRUE,
      set_variable_features = FALSE
    ),
    do_spatial_cluster = FALSE,
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )

  first <- do.call(
    public_workflow,
    c(list(srt = make_standard_spatial_rerun_object()), args)
  )
  first_tool <- first@tools[["SpatialVariableFeatures"]]
  second <- do.call(public_workflow, c(list(srt = first), args))
  second_tool <- second@tools[["SpatialVariableFeatures"]]

  expect_identical(first_tool, second_tool)
  expect_identical(
    second@tools$run_standard_spatial_workflow$status,
    "completed"
  )
})

test_that("native variable-feature-only output is verified without a tool", {
  skip_if_not_installed("BiocNeighbors")
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    make_standard_spatial_rerun_object(),
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      layer = "counts",
      coord.cols = c("col", "row"),
      k = 2,
      nfeatures = 3,
      min_spots = 1,
      nperm = 0,
      backend = "r",
      store_results = FALSE,
      set_variable_features = TRUE
    ),
    do_spatial_cluster = FALSE,
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_true(isTRUE(svf$result_stored))
  expect_match(svf$result_location, "VariableFeatures")
  expect_null(out@tools[["SpatialVariableFeatures"]])
  expect_length(
    suppressWarnings(SeuratObject::VariableFeatures(out, assay = "RNA")),
    3L
  )
})

test_that("NULL SVF controls resolve to producer defaults", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      set_variable_features = TRUE,
      store_results = TRUE,
      ...
    ) {
      expect_true(set_variable_features)
      expect_true(store_results)
      SeuratObject::VariableFeatures(srt, assay = assay) <- "gene1"
      srt@tools[["SpatialVariableFeatures"]] <- list(
        summary = list(top_features = "gene1")
      )
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    make_standard_spatial_storage_object(),
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = NULL,
      set_variable_features = NULL
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_true(isTRUE(svf$result_stored))
  expect_identical(svf$result_tool_key, "SpatialVariableFeatures")
  expect_true(svf$set_variable_features)
  expect_match(svf$result_location, "VariableFeatures")
})

test_that("fresh non-finite SVF scores certify a valid empty selection", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(
      srt,
      assay,
      set_variable_features,
      store_results,
      ...
    ) {
      expect_true(set_variable_features)
      expect_true(store_results)
      expect_null(srt@tools[["SpatialVariableFeatures"]])
      expect_length(
        suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
        0L
      )
      srt@tools[["SpatialVariableFeatures"]] <- make_valid_empty_svf_tool()
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    make_standard_spatial_storage_object(),
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = TRUE,
      set_variable_features = TRUE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(svf$status, "completed")
  expect_true(isTRUE(svf$result_stored))
  expect_identical(svf$result_tool_key, "SpatialVariableFeatures")
  expect_match(svf$result_location, "SpatialVariableFeatures")
  expect_identical(svf$variable_features_after, 0L)
  expect_identical(out@tools$run_standard_spatial_workflow$status, "completed")
})

test_that("restored HVF statistics do not revive a valid empty SVF selection", {
  srt <- make_standard_spatial_storage_object()
  assay_object <- srt[["RNA"]]
  skip_if_not(inherits(assay_object, "StdAssay"))
  variable <- c(TRUE, FALSE, TRUE)
  rank <- c(1L, NA_integer_, 2L)
  names(variable) <- names(rank) <- rownames(assay_object)
  assay_object[["vf_mock_counts_variable"]] <- variable
  assay_object[["vf_mock_counts_rank"]] <- rank
  srt[["RNA"]] <- assay_object
  SeuratObject::VariableFeatures(srt, assay = "RNA") <- c("gene1", "gene3")
  hvf_before <- srt[["RNA"]][[]][
    ,
    c("vf_mock_counts_variable", "vf_mock_counts_rank"),
    drop = FALSE
  ]

  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(srt, assay, ...) {
      expect_false(any(startsWith(colnames(srt[[assay]][[]]), "vf_")))
      srt@tools[["SpatialVariableFeatures"]] <- make_valid_empty_svf_tool()
      srt
    },
    .package = "scop"
  )

  out <- run_storage_workflow(
    srt,
    do_spatial_variable_features = TRUE,
    spatial_variable_features_params = list(
      store_results = TRUE,
      set_variable_features = TRUE
    )
  )
  feature_metadata <- out[["RNA"]][[]]
  hvf_after <- feature_metadata[
    ,
    c("vf_mock_counts_variable", "vf_mock_counts_rank"),
    drop = FALSE
  ]
  raw_active_features <- suppressWarnings(
    SeuratObject::VariableFeatures(out, assay = "RNA")
  )
  active_features <- standard_spatial_variable_features(out, assay = "RNA")
  historical_hvfs <- suppressWarnings(
    SeuratObject::VariableFeatures(
      out[["RNA"]],
      method = "mock",
      layer = "counts"
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_identical(hvf_after, hvf_before)
  expect_identical(active_features, character())
  expect_false(any(raw_active_features %in% rownames(out[["RNA"]]), na.rm = TRUE))
  expect_identical(historical_hvfs, c("gene1", "gene3"))
  expect_true("var.features" %in% colnames(feature_metadata))
  expect_false(any(feature_metadata$var.features, na.rm = TRUE))
  expect_true(all(is.na(feature_metadata$var.features.rank)))
  expect_identical(svf$variable_features_after, 0L)
})

test_that("malformed SVF tools cannot certify an empty selection", {
  valid_tool <- make_valid_empty_svf_tool()
  malformed_tools <- list(
    list_top_features = within(valid_tool, {
      summary$top_features <- list()
    }),
    logical_scores = within(valid_tool, {
      result$score <- c(NA, NA)
    }),
    finite_score = within(valid_tool, {
      result$score <- c(NA_real_, 1)
    }),
    incomplete_result = within(valid_tool, {
      result$method <- NULL
    }),
    inconsistent_summary = within(valid_tool, {
      summary$n_features <- 1L
    })
  )

  for (case_name in names(malformed_tools)) {
    candidate <- malformed_tools[[case_name]]
    testthat::local_mocked_bindings(
      RunStandardWorkflow = function(srt, ...) srt,
      RunSpatialVariableFeatures = function(srt, assay, ...) {
        expect_null(srt@tools[["SpatialVariableFeatures"]])
        expect_length(
          suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
          0L
        )
        srt@tools[["SpatialVariableFeatures"]] <- candidate
        srt
      },
      .package = "scop"
    )

    expect_error(
      run_storage_workflow(
        make_standard_spatial_storage_object(),
        do_spatial_variable_features = TRUE,
        spatial_variable_features_params = list(
          store_results = TRUE,
          set_variable_features = TRUE
        )
      ),
      "completed without the expected",
      info = case_name
    )
  }
})

test_that("a stale empty SVF tool cannot certify an empty selection", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(srt, assay, ...) {
      expect_null(srt@tools[["SpatialVariableFeatures"]])
      expect_length(
        suppressWarnings(SeuratObject::VariableFeatures(srt, assay = assay)),
        0L
      )
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  stale_tool <- list(
    result = data.frame(feature = "gene1", score = NA_real_),
    summary = list(top_features = character())
  )
  srt@tools[["SpatialVariableFeatures"]] <- stale_tool

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_spatial_variable_features = TRUE,
      spatial_variable_features_params = list(
        store_results = TRUE,
        set_variable_features = TRUE
      )
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(svf$status, "failed")
  expect_false(isTRUE(svf$result_stored))
  expect_match(svf$reason, "SpatialVariableFeatures")
  expect_identical(srt@tools[["SpatialVariableFeatures"]], stale_tool)
})

test_that("a stale SVF tool does not satisfy a quiet stored producer", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpatialVariableFeatures = function(srt, ...) {
      expect_null(srt@tools[["SpatialVariableFeatures"]])
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt@tools[["SpatialVariableFeatures"]] <- list(result = "stale")

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_spatial_variable_features = TRUE,
      spatial_variable_features_params = list(store_results = TRUE)
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  svf <- stages[stages$stage == "spatial_variable_features", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(svf$status, "failed")
  expect_false(isTRUE(svf$result_stored))
  expect_match(svf$reason, "SpatialVariableFeatures")
  expect_identical(srt@tools[["SpatialVariableFeatures"]]$result, "stale")
})

test_that("metadata-only deconvolution records actual output columns", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, tool_name, store_results, ...) {
      expect_identical(prefix, "Custom")
      expect_identical(tool_name, "CustomTool")
      expect_false(store_results)
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = "Custom",
      tool_name = "CustomTool",
      store_results = FALSE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_identical(deconv$status, "completed")
  expect_true(isTRUE(deconv$result_stored))
  expect_true(is.na(deconv$result_tool_key))
  expect_match(deconv$result_metadata_key, "Custom_prop_typeA")
  expect_match(deconv$result_location, "Custom_dominant_type")
})

test_that("non-cell2location methods preserve same-prefix abundance annotations", {
  producer <- function(srt, prefix, store_results, ...) {
    expect_false(store_results)
    expect_true(paste0(prefix, "_abundance_annotation") %in%
      colnames(srt@meta.data))
    expect_false(any(startsWith(
      colnames(srt@meta.data),
      paste0(prefix, "_prop_")
    )))
    expect_false(
      paste0(prefix, "_dominant_type") %in% colnames(srt@meta.data)
    )
    expect_false(paste0(prefix, "_max_prop") %in% colnames(srt@meta.data))
    add_mock_deconv_outputs(srt, prefix = prefix)
  }
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = producer,
    RunSPOTlight = producer,
    .package = "scop"
  )

  for (method in c("RCTD", "SPOTlight")) {
    srt <- make_standard_spatial_storage_object()
    srt$label <- factor(rep("typeA", ncol(srt)))
    srt$Custom_abundance_annotation <- seq_len(ncol(srt))
    srt <- add_mock_deconv_outputs(srt, prefix = "Custom")

    out <- run_storage_workflow(
      srt,
      do_deconvolution = TRUE,
      deconvolution_method = method,
      reference = srt,
      reference_label = "label",
      deconvolution_params = list(
        prefix = "Custom",
        store_results = FALSE
      )
    )
    stages <- out@tools$run_standard_spatial_workflow$stages
    deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

    expect_identical(deconv$status, "completed", info = method)
    expect_identical(
      unname(out$Custom_abundance_annotation),
      seq_len(ncol(out)),
      info = method
    )
  }
})

test_that("Cell2location clears and replaces stale abundance metadata", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunCell2location = function(srt, prefix, store_results, ...) {
      expect_false(store_results)
      expect_false(any(startsWith(
        colnames(srt@meta.data),
        paste0(prefix, "_abundance_")
      )))
      expect_false(any(startsWith(
        colnames(srt@meta.data),
        paste0(prefix, "_prop_")
      )))
      srt[[paste0(prefix, "_abundance_typeA")]] <- rep(3, ncol(srt))
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  srt$Custom_abundance_stale <- 99
  srt <- add_mock_deconv_outputs(srt, prefix = "Custom")

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "Cell2location",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = "Custom",
      store_results = FALSE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_identical(deconv$status, "completed")
  expect_false("Custom_abundance_stale" %in% colnames(out@meta.data))
  expect_true("Custom_abundance_typeA" %in% colnames(out@meta.data))
  expect_match(deconv$result_metadata_key, "Custom_abundance_typeA")
})

test_that("Cell2location cannot complete without fresh abundance metadata", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunCell2location = function(srt, prefix, ...) {
      expect_false(any(startsWith(
        colnames(srt@meta.data),
        paste0(prefix, "_abundance_")
      )))
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  srt$Custom_abundance_stale <- 99
  caller_before <- srt

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_deconvolution = TRUE,
      deconvolution_method = "Cell2location",
      reference = srt,
      reference_label = "label",
      deconvolution_params = list(
        prefix = "Custom",
        store_results = FALSE
      )
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(deconv$status, "failed")
  expect_match(deconv$reason, "Custom_abundance")
  expect_identical(srt, caller_before)
})

test_that("stored deconvolution records a custom tool and metadata", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, tool_name, store_results, ...) {
      expect_true(store_results)
      add_mock_deconv_outputs(
        srt,
        prefix = prefix,
        tool_name = tool_name,
        store_results = store_results
      )
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = "Custom",
      tool_name = "CustomTool",
      store_results = TRUE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_identical(deconv$status, "completed")
  expect_true(isTRUE(deconv$result_stored))
  expect_identical(deconv$result_tool_key, "CustomTool")
  expect_match(deconv$result_metadata_key, "Custom_prop_typeA")
  expect_match(deconv$result_location, "CustomTool")
})

test_that("deconvolution rejects workflow-owned result keys before execution", {
  calls <- character()
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) {
      calls <<- c(calls, "RunStandardWorkflow")
      srt
    },
    RunRCTD = function(...) {
      calls <<- c(calls, "RunRCTD")
      stop("deconvolution backend was reached")
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  workflow_owned_keys <- c(
    "run_standard_spatial_workflow",
    "SpotQC",
    "SpatialVariableFeatures",
    "BayesSpace"
  )

  for (tool_name in workflow_owned_keys) {
    expect_error(
      run_storage_workflow(
        srt,
        do_deconvolution = TRUE,
        deconvolution_method = "RCTD",
        reference = srt,
        reference_label = "label",
        deconvolution_params = list(tool_name = tool_name)
      ),
      "reserved",
      fixed = TRUE
    )
  }

  expect_length(calls, 0L)
})

test_that("metadata-only deconvolution allows workflow-owned tool names", {
  called <- FALSE
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, tool_name, store_results, ...) {
      called <<- TRUE
      expect_identical(tool_name, "BayesSpace")
      expect_false(store_results)
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  bayesspace_tool <- list(result = "keep")
  srt@tools[["BayesSpace"]] <- bayesspace_tool

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = "Custom",
      tool_name = "BayesSpace",
      store_results = FALSE
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_true(called)
  expect_identical(deconv$status, "completed")
  expect_true(is.na(deconv$result_tool_key))
  expect_match(deconv$result_metadata_key, "Custom_prop_typeA")
  expect_identical(out@tools[["BayesSpace"]], bayesspace_tool)
})

test_that("planned spatial metadata collisions fail before any producer", {
  calls <- character()
  testthat::local_mocked_bindings(
    RunSpotQC = function(srt, ...) {
      calls <<- c(calls, "RunSpotQC")
      srt
    },
    RunStandardWorkflow = function(srt, ...) {
      calls <<- c(calls, "RunStandardWorkflow")
      srt
    },
    RunBayesSpace = function(srt, ...) {
      calls <<- c(calls, "RunBayesSpace")
      srt
    },
    RunRCTD = function(srt, ...) {
      calls <<- c(calls, "RunRCTD")
      srt
    },
    RunCell2location = function(srt, ...) {
      calls <<- c(calls, "RunCell2location")
      srt
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  expect_error(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = TRUE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = 2,
      bayesspace_params = list(cluster_colname = "SpotQC"),
      do_deconvolution = FALSE,
      verbose = FALSE
    ),
    "metadata outputs collide"
  )
  expect_length(calls, 0L)

  expect_error(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = 2,
      bayesspace_params = list(
        cluster_colname = "SharedDomain",
        init_colname = "SharedDomain"
      ),
      do_deconvolution = FALSE,
      verbose = FALSE
    ),
    "metadata outputs collide"
  )
  expect_length(calls, 0L)

  expect_error(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = 2,
      bayesspace_params = list(
        cluster_colname = "Custom_prop_typeA"
      ),
      do_deconvolution = TRUE,
      deconvolution_method = "RCTD",
      reference = srt,
      reference_label = "label",
      deconvolution_params = list(prefix = "Custom"),
      verbose = FALSE
    ),
    "deconvolution metadata.*collides"
  )
  expect_length(calls, 0L)

  expect_error(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = TRUE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = 2,
      bayesspace_params = list(
        cluster_colname = "CustomCluster",
        init_colname = "SpotQC"
      ),
      do_deconvolution = FALSE,
      verbose = FALSE
    ),
    "metadata outputs collide"
  )
  expect_length(calls, 0L)

  expect_error(
    original(
      srt,
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE,
      spatial_q = 2,
      bayesspace_params = list(
        cluster_colname = "CustomCluster",
        init_colname = "Custom_abundance_typeA"
      ),
      do_deconvolution = TRUE,
      deconvolution_method = "Cell2location",
      reference = srt,
      reference_label = "label",
      deconvolution_params = list(
        prefix = "Custom",
        store_results = FALSE
      ),
      verbose = FALSE
    ),
    "deconvolution metadata.*collides"
  )
  expect_length(calls, 0L)
})

test_that("RCTD preserves an abundance-named BayesSpace init output", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunBayesSpace = function(srt, cluster_colname, init_colname, ...) {
      srt[[cluster_colname]] <- rep("domain1", ncol(srt))
      srt[[init_colname]] <- rep("initial1", ncol(srt))
      srt@tools[["BayesSpace"]] <- list(result = "fresh")
      srt
    },
    RunRCTD = function(srt, prefix, store_results, ...) {
      expect_false(store_results)
      expect_true("Custom_abundance_annotation" %in% colnames(srt@meta.data))
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  out <- original(
    srt,
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "CustomCluster",
      init_colname = "Custom_abundance_annotation"
    ),
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = "Custom",
      store_results = FALSE
    ),
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages

  expect_true(all(stages$status[stages$requested] == "completed"))
  expect_identical(unique(out$Custom_abundance_annotation), "initial1")
})

test_that("unrequested stage names do not create metadata collisions", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpotQC = function(srt, ...) {
      srt$SpotQC <- rep("Pass", ncol(srt))
      srt
    },
    RunBayesSpace = function(srt, cluster_colname, init_colname, ...) {
      expect_identical(cluster_colname, "SpotQC")
      expect_identical(init_colname, "BayesSpace_init")
      srt[[cluster_colname]] <- rep("domain1", ncol(srt))
      srt[[init_colname]] <- rep("initial1", ncol(srt))
      srt@tools[["BayesSpace"]] <- list(result = "fresh")
      srt
    },
    .package = "scop"
  )
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    make_standard_spatial_storage_object(),
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "SpotQC",
      init_colname = NULL
    ),
    do_deconvolution = FALSE,
    verbose = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  domain <- stages[stages$stage == "spatial_clustering", , drop = FALSE]

  expect_identical(domain$status, "completed")
  expect_match(domain$result_metadata_key, "SpotQC")
  expect_match(domain$result_metadata_key, "BayesSpace_init")
  expect_equal(unique(out$SpotQC), "domain1")

  qc_only <- original(
    make_standard_spatial_storage_object(),
    assay = "RNA",
    do_spot_qc = TRUE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = FALSE,
    bayesspace_params = list(init_colname = "SpotQC"),
    do_deconvolution = FALSE,
    verbose = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE
  )
  qc_stages <- qc_only@tools$run_standard_spatial_workflow$stages
  qc <- qc_stages[qc_stages$stage == "quality_control", , drop = FALSE]

  expect_identical(qc$status, "completed")
  expect_equal(unique(qc_only$SpotQC), "Pass")
})

test_that("SPOTlight dispatch satisfies the normalized storage contract", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSPOTlight = function(
      srt,
      reference,
      reference_label,
      prefix,
      tool_name,
      store_results,
      ...
    ) {
      expect_s4_class(reference, "Seurat")
      expect_identical(reference_label, "label")
      expect_identical(prefix, "SPOTlight")
      expect_identical(tool_name, "SPOTlight")
      expect_true(store_results)
      add_mock_deconv_outputs(
        srt,
        prefix = prefix,
        tool_name = tool_name,
        store_results = store_results
      )
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "SPOTlight",
    reference = srt,
    reference_label = "label"
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_identical(deconv$status, "completed")
  expect_identical(deconv$result_tool_key, "SPOTlight")
  expect_match(deconv$result_metadata_key, "SPOTlight_prop_typeA")
})

test_that("stale deconvolution metadata does not satisfy a quiet producer", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, store_results, ...) {
      expect_false(store_results)
      expect_false(any(startsWith(
        colnames(srt@meta.data),
        paste0(prefix, "_prop_")
      )))
      expect_false(
        paste0(prefix, "_dominant_type") %in% colnames(srt@meta.data)
      )
      expect_false(paste0(prefix, "_max_prop") %in% colnames(srt@meta.data))
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))
  srt <- add_mock_deconv_outputs(srt, prefix = "Custom")

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_deconvolution = TRUE,
      deconvolution_method = "RCTD",
      reference = srt,
      reference_label = "label",
      deconvolution_params = list(
        prefix = "Custom",
        store_results = FALSE
      )
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(deconv$status, "failed")
  expect_false(isTRUE(deconv$result_stored))
  expect_match(deconv$reason, "Custom_prop")
  expect_true("Custom_prop_typeA" %in% colnames(srt@meta.data))
  expect_true("Custom_dominant_type" %in% colnames(srt@meta.data))
  expect_true("Custom_max_prop" %in% colnames(srt@meta.data))
})

test_that("NULL deconvolution controls resolve to producer defaults", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(
      srt,
      prefix = "RCTD",
      tool_name = "RCTD",
      store_results = TRUE,
      ...
    ) {
      expect_identical(prefix, "RCTD")
      expect_identical(tool_name, "RCTD")
      expect_true(store_results)
      add_mock_deconv_outputs(
        srt,
        prefix = prefix,
        tool_name = tool_name,
        store_results = store_results
      )
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  out <- run_storage_workflow(
    srt,
    do_deconvolution = TRUE,
    deconvolution_method = "RCTD",
    reference = srt,
    reference_label = "label",
    deconvolution_params = list(
      prefix = NULL,
      tool_name = NULL,
      store_results = NULL
    )
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_identical(deconv$status, "completed")
  expect_identical(deconv$result_tool_key, "RCTD")
  expect_match(deconv$result_metadata_key, "RCTD_prop_typeA")
})

test_that("partial postconditions retain stored metadata diagnostics", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, store_results, ...) {
      expect_true(store_results)
      add_mock_deconv_outputs(srt, prefix = prefix)
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_deconvolution = TRUE,
      reference = srt,
      reference_label = "label"
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(deconv$status, "failed")
  expect_true(isTRUE(deconv$result_stored))
  expect_true(is.na(deconv$result_tool_key))
  expect_match(deconv$result_metadata_key, "RCTD_prop_typeA")
  expect_match(deconv$result_location, "RCTD_dominant_type")
  expect_match(deconv$reason, "tools")
})

test_that("incomplete metadata remains visible when a tool was stored", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunRCTD = function(srt, prefix, tool_name, ...) {
      srt[[paste0(prefix, "_dominant_type")]] <- rep("typeA", ncol(srt))
      srt@tools[[tool_name]] <- list(proportions = matrix(0.5, ncol(srt), 1))
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$label <- factor(rep("typeA", ncol(srt)))

  condition <- tryCatch(
    run_storage_workflow(
      srt,
      do_deconvolution = TRUE,
      reference = srt,
      reference_label = "label"
    ),
    error = identity
  )
  stages <- attr(condition, "standard_spatial_stages")
  deconv <- stages[stages$stage == "deconvolution", , drop = FALSE]

  expect_s3_class(condition, "error")
  expect_identical(deconv$status, "failed")
  expect_true(isTRUE(deconv$result_stored))
  expect_identical(deconv$result_tool_key, "RCTD")
  expect_match(deconv$result_metadata_key, "RCTD_dominant_type")
  expect_match(deconv$reason, "RCTD_prop")
})

test_that("a stale BayesSpace init is removed when fresh init is absent", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunBayesSpace = function(srt, cluster_colname, init_colname, ...) {
      expect_false(cluster_colname %in% colnames(srt@meta.data))
      expect_false(init_colname %in% colnames(srt@meta.data))
      expect_null(srt@tools[["BayesSpace"]])
      srt[[cluster_colname]] <- rep("domain1", ncol(srt))
      srt@tools[["BayesSpace"]] <- list(result = "fresh")
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$CustomCluster <- "stale cluster"
  srt$CustomInit <- "stale init"
  srt@tools[["BayesSpace"]] <- list(result = "stale")
  caller_before <- srt
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    srt,
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "CustomCluster",
      init_colname = "CustomInit"
    ),
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  domain <- stages[stages$stage == "spatial_clustering", , drop = FALSE]

  expect_identical(domain$status, "completed")
  expect_true(isTRUE(domain$result_stored))
  expect_match(domain$result_metadata_key, "CustomCluster")
  expect_false(grepl("CustomInit", domain$result_metadata_key, fixed = TRUE))
  expect_false("CustomInit" %in% colnames(out@meta.data))
  expect_identical(srt, caller_before)
})

test_that("quality control and BayesSpace probes replace stale outputs", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(srt, ...) srt,
    RunSpotQC = function(srt, ...) {
      expect_false("SpotQC" %in% colnames(srt@meta.data))
      srt$SpotQC <- factor(rep("Pass", ncol(srt)), levels = c("Pass", "Fail"))
      srt
    },
    RunBayesSpace = function(
      srt,
      cluster_colname = "BayesSpace_cluster",
      init_colname = "BayesSpace_init",
      ...
    ) {
      expect_false(cluster_colname %in% colnames(srt@meta.data))
      expect_false(init_colname %in% colnames(srt@meta.data))
      expect_null(srt@tools[["BayesSpace"]])
      srt[[cluster_colname]] <- factor(
        rep(c("A", "B"), length.out = ncol(srt))
      )
      srt[[init_colname]] <- factor(
        rep(c("B", "A"), length.out = ncol(srt))
      )
      srt@tools[["BayesSpace"]] <- list(
        colData = data.frame(row.names = colnames(srt))
      )
      srt
    },
    .package = "scop"
  )
  srt <- make_standard_spatial_storage_object()
  srt$SpotQC <- "stale"
  srt$CustomCluster <- "stale"
  srt$CustomInit <- "stale"
  srt@tools[["BayesSpace"]] <- list(result = "stale")
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")

  out <- original(
    srt,
    assay = "RNA",
    do_spot_qc = TRUE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_q = 2,
    bayesspace_params = list(
      cluster_colname = "CustomCluster",
      init_colname = "CustomInit"
    ),
    do_deconvolution = FALSE,
    do_normalization = FALSE,
    do_HVF_finding = FALSE,
    do_scaling = FALSE,
    verbose = FALSE
  )
  stages <- out@tools$run_standard_spatial_workflow$stages
  qc <- stages[stages$stage == "quality_control", , drop = FALSE]
  domain <- stages[stages$stage == "spatial_clustering", , drop = FALSE]

  expect_identical(qc$status, "completed")
  expect_true(isTRUE(qc$result_stored))
  expect_identical(qc$result_metadata_key, "SpotQC")
  expect_identical(domain$status, "completed")
  expect_true(isTRUE(domain$result_stored))
  expect_identical(domain$result_tool_key, "BayesSpace")
  expect_match(domain$result_metadata_key, "CustomCluster")
  expect_match(domain$result_metadata_key, "CustomInit")
})
