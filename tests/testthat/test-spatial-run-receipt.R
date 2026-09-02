make_receipt_spatial_object <- function() {
  counts <- matrix(
    c(
      3, 1, 0, 2,
      0, 4, 1, 0,
      2, 1, 3, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt <- SeuratObject::AddMetaData(srt, metadata = c(0, 2, 1, 3), col.name = "col")
  SeuratObject::AddMetaData(srt, metadata = c(0, 0, 1, 1), col.name = "row")
}

make_receipt_multi_image_object <- function() {
  srt <- make_receipt_spatial_object()
  assay <- SeuratObject::DefaultAssay(srt)
  slice1 <- data.frame(
    x = c(0, 2),
    y = c(0, 0),
    row.names = c("Spot1", "Spot2")
  )
  slice2 <- data.frame(
    x = c(1, 3),
    y = c(1, 1),
    row.names = c("Spot3", "Spot4")
  )
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    slice1,
    type = "centroids",
    assay = assay,
    key = "receipt1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    slice2,
    type = "centroids",
    assay = assay,
    key = "receipt2_"
  )
  srt
}

make_receipt_svf_multi_image_object <- function() {
  srt <- make_receipt_spatial_object()
  alt_counts <- SeuratObject::LayerData(srt, assay = "RNA", layer = "counts")
  rownames(alt_counts) <- paste0("AltGene", seq_len(nrow(alt_counts)))
  srt[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt_counts)
  slice1 <- data.frame(
    x = c(0, 2, 1),
    y = c(0, 0, 1),
    row.names = c("Spot1", "Spot2", "Spot3")
  )
  slice2 <- data.frame(
    x = c(2, 1, 3),
    y = c(0, 1, 1),
    row.names = c("Spot2", "Spot3", "Spot4")
  )
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    slice1,
    type = "centroids",
    assay = "RNA",
    key = "svf1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    slice2,
    type = "centroids",
    assay = "RNA",
    key = "svf2_"
  )
  srt
}

receipt_plain <- function(messages) {
  cli::ansi_strip(paste(messages, collapse = "\n"))
}

expect_receipt_plot_hint <- function(messages, expected_call) {
  plain <- receipt_plain(messages)
  expect_true(grepl(
    paste0("Plot returned object `", expected_call, "`"),
    plain,
    fixed = TRUE
  ))
  expect_true(grepl("<returned_object>", expected_call, fixed = TRUE))
  expect_false(grepl("srt", expected_call, fixed = TRUE))
  invisible(plain)
}

expect_receipt_inspect_hint <- function(messages, expected_call) {
  plain <- receipt_plain(messages)
  expect_true(grepl(
    paste0("Inspect returned object `", expected_call, "`"),
    plain,
    fixed = TRUE
  ))
  expect_true(grepl("<returned_object>", expected_call, fixed = TRUE))
  expect_false(grepl("srt", expected_call, fixed = TRUE))
  invisible(plain)
}

test_that("spatial_run_receipt renders exact inline labels and an explicit result placeholder", {
  plot_call <- 'SpatialSpotPlot(<returned_object>, group.by = "BANKSY_cluster")'
  out <- testthat::capture_messages(
    spatial_run_receipt(
      done = "{.pkg BANKSY} completed for {.val 10} spots",
      scope = "assay {.val Spatial}, layer {.val counts}",
      saved = "metadata column {.var BANKSY_cluster} and returned object tool bundle {.var BANKSY}",
      plot = plot_call,
      verbose = TRUE,
      .envir = environment()
    )
  )
  plain <- receipt_plain(out)
  expect_match(plain, "BANKSY completed for \"10\" spots")
  expect_match(plain, "Scope assay \"Spatial\", layer \"counts\"")
  expect_match(plain, "Saved metadata column `BANKSY_cluster`")
  expect_match(plain, "returned object tool bundle `BANKSY`")
  expect_true(grepl(
    paste0("Plot returned object `", plot_call, "`"),
    plain,
    fixed = TRUE
  ))
  expect_false(grepl("srt", plain, fixed = TRUE))
})

test_that("spatial_run_receipt renders partial output without a warning condition", {
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)
  out <- testthat::expect_warning(
    testthat::capture_messages(
      spatial_run_receipt(
        done = "SpotSweeper completed partially",
        scope = "assay {.val Spatial}",
        inspect = "tool bundle SpotSweeper",
        status = "partial",
        verbose = TRUE,
        .envir = environment()
      )
    ),
    NA
  )
  plain <- receipt_plain(out)
  expect_match(plain, "SpotSweeper completed partially")
  expect_true(grepl(
    "Inspect returned object `tool bundle SpotSweeper`",
    plain,
    fixed = TRUE
  ))
})

test_that("spatial_run_receipt maps partial output to non-success display", {
  message_types <- character()
  testthat::local_mocked_bindings(
    log_message = function(..., message_type = "info") {
      message_types <<- c(message_types, message_type)
      invisible(NULL)
    },
    .package = "scop"
  )
  spatial_run_receipt(
    done = "SpotSweeper completed partially",
    status = "partial",
    verbose = TRUE
  )
  expect_identical(message_types, "running")
})

test_that("spatial_run_receipt is silent when verbose is FALSE", {
  expect_silent(
    spatial_run_receipt(
      done = "{.pkg BANKSY} completed",
      scope = "assay {.val Spatial}",
      saved = "returned object tool bundle {.var BANKSY}",
      plot = "SpatialSpotPlot(<returned_object>)",
      verbose = FALSE,
      .envir = environment()
    )
  )
})

test_that("spatial_run_receipt reuses thisutils verbose NULL semantics", {
  old_options <- options(log_message.verbose = FALSE)
  on.exit(options(old_options), add = TRUE)
  expect_silent(
    spatial_run_receipt(
      done = "{.pkg BANKSY} completed",
      plot = "SpatialSpotPlot(<returned_object>)",
      verbose = NULL,
      .envir = environment()
    )
  )
})

test_that("spatial_run_receipt uses Replaced and prefers Plot over Inspect", {
  plot_call <- 'SpatialNetworkPlot(<returned_object>, graph.name = "knn_k6")'
  out <- testthat::capture_messages(
    spatial_run_receipt(
      done = "{.pkg SpatialNetwork} completed",
      saved = "graph {.val knn_k6}",
      plot = plot_call,
      inspect = "tool bundle SpatialNetwork",
      replaced = TRUE,
      verbose = TRUE,
      .envir = environment()
    )
  )
  plain <- receipt_plain(out)
  expect_match(plain, "Replaced graph")
  expect_true(grepl(
    paste0("Plot returned object `", plot_call, "`"),
    plain,
    fixed = TRUE
  ))
  expect_false(grepl("Inspect", plain, fixed = TRUE))
})

test_that("spatial_run_receipt_quote safely round-trips custom keys", {
  key <- 'domain "quoted"\\path'
  quoted <- spatial_run_receipt_quote(key, "key")
  expect_identical(eval(parse(text = quoted)), key)
  expect_error(spatial_run_receipt_quote(character(), "key"), "key")
})

test_that("spatial_run_receipt validates its visible interface", {
  expect_error(spatial_run_receipt(done = character()), "done")
  expect_error(spatial_run_receipt(done = "ok", scope = 1), "scope")
  expect_error(spatial_run_receipt(done = "ok", plot = c("a", "b")), "plot")
  expect_error(spatial_run_receipt(done = "ok", inspect = c("a", "b")), "inspect")
  expect_error(spatial_run_receipt(done = "ok", replaced = NA), "replaced")
  expect_error(spatial_run_receipt(done = "ok", .envir = NULL), "envir")
})

test_that("RunSpotQC emits an exact receipt with an explicit result placeholder", {
  srt <- make_receipt_spatial_object()
  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = TRUE
    )
  ))
  plain <- receipt_plain(messages)
  expect_true(grepl(
    "Spot QC completed: 4 evaluated, 4 Pass, 0 Fail",
    plain,
    fixed = TRUE
  ))
  expect_true(grepl('Scope assay "RNA", layer "counts"', plain, fixed = TRUE))
  expect_true(grepl("Saved metadata column `SpotQC`", plain, fixed = TRUE))
  expect_identical(as.integer(table(out$SpotQC)), c(4L, 0L))
  expect_receipt_plot_hint(
    messages,
    'SpatialSpotPlot(<returned_object>, group.by = "SpotQC")'
  )

  filtered <- NULL
  filtered_messages <- suppressWarnings(testthat::capture_messages(
    filtered <- RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      return_filtered = TRUE,
      verbose = TRUE
    )
  ))
  filtered_plain <- receipt_plain(filtered_messages)
  expect_true(grepl("4 returned", filtered_plain, fixed = TRUE))
  expect_false(grepl("Plot", filtered_plain, fixed = TRUE))
  expect_equal(ncol(filtered), 4)
})

test_that("RunSpotQC offers an honest Inspect hint without plottable coordinates", {
  counts <- matrix(
    c(3, 1, 0, 2, 0, 4, 1, 0),
    nrow = 2,
    dimnames = list(paste0("Gene", 1:2), paste0("Spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = TRUE
    )
  ))
  inspect_call <- paste0(
    "table(<returned_object>[[\"SpotQC\", drop = TRUE]], ",
    "useNA = \"ifany\")"
  )
  plain <- expect_receipt_inspect_hint(messages, inspect_call)
  expect_false(grepl("Plot returned object", plain, fixed = TRUE))
  expect_true("SpotQC" %in% colnames(out[[]]))
  expect_identical(
    eval(parse(text = sub("<returned_object>", "out", inspect_call, fixed = TRUE))),
    table(out[["SpotQC", drop = TRUE]], useNA = "ifany")
  )
})

test_that("RunSpotQC counts and labels only cells present in the selected assay", {
  full_counts <- matrix(
    1,
    nrow = 2,
    ncol = 6,
    dimnames = list(paste0("Gene", 1:2), paste0("Spot", 1:6))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(full_counts))
  assay_counts <- full_counts[, 1:4, drop = FALSE]
  assay_counts[, "Spot1"] <- 0
  srt[["SUBSET"]] <- suppressWarnings(
    SeuratObject::CreateAssay5Object(counts = assay_counts)
  )

  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "SUBSET",
      qc_metrics = c("umi", "gene", "mito"),
      UMI_threshold = 1,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = TRUE
    )
  ))
  plain <- receipt_plain(messages)

  expect_true(grepl(
    "Spot QC completed: 4 evaluated, 3 Pass, 1 Fail",
    plain,
    fixed = TRUE
  ))
  expect_identical(
    as.character(out$SpotQC[1:4]),
    c("Fail", rep("Pass", 3))
  )
  expect_true(all(is.na(out$SpotQC[5:6])))
  expect_identical(
    as.integer(table(out$SpotQC, useNA = "always")),
    c(3L, 1L, 2L)
  )

  filtered <- suppressWarnings(RunSpotQC(
    srt,
    assay = "SUBSET",
    qc_metrics = c("umi", "gene", "mito"),
    UMI_threshold = 1,
    gene_threshold = 0,
    mito_threshold = 100,
    return_filtered = TRUE,
    verbose = FALSE
  ))
  expect_identical(colnames(filtered), paste0("Spot", 2:4))
})

test_that("RunSpotQC aligns metadata outlier metrics to noncontiguous assay cells", {
  full_counts <- matrix(
    1,
    nrow = 2,
    ncol = 6,
    dimnames = list(paste0("Gene", 1:2), paste0("Spot", 1:6))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(full_counts))
  evaluated_spots <- c("Spot1", "Spot3", "Spot4", "Spot6")
  srt[["SUBSET"]] <- suppressWarnings(
    SeuratObject::CreateAssay5Object(
      counts = full_counts[, evaluated_spots, drop = FALSE]
    )
  )
  srt$custom_metric <- c(0, 100, 0.5, 1, 2, 1.5)

  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "SUBSET",
      qc_metrics = "outlier",
      outlier_threshold = "custom_metric:higher:3",
      outlier_n = 1,
      verbose = TRUE
    )
  ))

  expect_true(grepl(
    "Spot QC completed: 4 evaluated, 4 Pass, 0 Fail",
    receipt_plain(messages),
    fixed = TRUE
  ))
  expect_identical(as.character(out$SpotQC[evaluated_spots]), rep("Pass", 4))
  expect_true(all(is.na(out$SpotQC[c("Spot2", "Spot5")])))
  outlier_col <- make.names("spot_custom_metric:higher:3")
  expect_identical(
    unname(out[[outlier_col, drop = TRUE]][evaluated_spots]),
    rep(FALSE, 4)
  )
  expect_true(all(is.na(out[[outlier_col, drop = TRUE]][c("Spot2", "Spot5")])))
})

test_that("RunSpotQC preserves tri-state outlier rules and aggregates uncertainty", {
  full_counts <- matrix(
    1,
    nrow = 2,
    ncol = 8,
    dimnames = list(paste0("Gene", 1:2), paste0("Spot", 1:8))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(full_counts))
  evaluated_spots <- c("Spot1", "Spot3", "Spot4", "Spot5", "Spot7", "Spot8")
  assay_counts <- full_counts[, evaluated_spots, drop = FALSE]
  assay_counts[, "Spot3"] <- 0
  srt[["SUBSET"]] <- suppressWarnings(
    SeuratObject::CreateAssay5Object(
      counts = assay_counts
    )
  )
  metric_a <- setNames(rep(0, ncol(srt)), colnames(srt))
  metric_b <- metric_a
  metric_a[evaluated_spots] <- c(100, 200, 1, NA, 3, 4)
  metric_b[evaluated_spots] <- c(100, NA, Inf, NA, 2, 3)
  srt$metric_a <- metric_a
  srt$metric_b <- metric_b

  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "SUBSET",
      qc_metrics = c("outlier", "umi"),
      outlier_threshold = c(
        "metric_a:higher:3",
        "metric_b:higher:3"
      ),
      outlier_n = 2,
      UMI_threshold = 1,
      verbose = TRUE
    )
  ))

  expect_identical(
    unname(as.character(out$SpotQC[evaluated_spots])),
    c("Fail", "Fail", "Pass", NA, "Pass", "Pass")
  )
  expect_identical(
    unname(as.character(out$spot_outlier_qc[evaluated_spots])),
    c("Fail", NA, "Pass", NA, "Pass", "Pass")
  )
  expect_true(all(is.na(out$SpotQC[c("Spot2", "Spot6")])))
  raw_a <- make.names("spot_metric_a:higher:3")
  raw_b <- make.names("spot_metric_b:higher:3")
  expect_identical(
    unname(out[[raw_a, drop = TRUE]][evaluated_spots]),
    c(TRUE, TRUE, FALSE, NA, FALSE, FALSE)
  )
  expect_identical(
    unname(out[[raw_b, drop = TRUE]][evaluated_spots]),
    c(TRUE, NA, NA, NA, FALSE, FALSE)
  )
  expect_true(all(is.na(out[[raw_a, drop = TRUE]][c("Spot2", "Spot6")])))
  expect_true(grepl(
    paste0(
      "Spot QC completed: 5 evaluated, 3 Pass, 2 Fail; 1\\s+",
      "NotEvaluated"
    ),
    receipt_plain(messages)
  ))

  filtered <- suppressWarnings(RunSpotQC(
    srt,
    assay = "SUBSET",
    qc_metrics = c("outlier", "umi"),
    outlier_threshold = c(
      "metric_a:higher:3",
      "metric_b:higher:3"
    ),
    outlier_n = 2,
    UMI_threshold = 1,
    return_filtered = TRUE,
    verbose = FALSE
  ))
  expect_identical(
    colnames(filtered),
    c("Spot4", "Spot7", "Spot8")
  )
})

test_that("RunSpotQC validates outlier aggregation and metadata metrics", {
  srt <- make_receipt_spatial_object()
  srt$custom_metric <- c(1, 2, 100, 3)

  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "outlier",
      outlier_threshold = NA_character_,
      verbose = FALSE
    ),
    "outlier_threshold.*non-empty character vector"
  )
  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "outlier",
      outlier_threshold = "custom_metric:higher:3",
      outlier_n = NA_real_,
      verbose = FALSE
    ),
    "outlier_n.*finite positive integer"
  )
  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "outlier",
      outlier_threshold = "custom_metric:higher:3",
      outlier_n = 2,
      verbose = FALSE
    ),
    "cannot exceed.*outlier_threshold"
  )

  srt$factor_metric <- factor(c("1", "2", "100", "3"))
  expect_error(
    suppressWarnings(RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "outlier",
      outlier_threshold = "factor_metric:higher:3",
      verbose = FALSE
    )),
    "factor_metric.*numeric vector"
  )
  expect_error(
    spot_qc_validate_metric(1, metric = "short", n_expected = 2),
    "one\\s+value per selected-assay spot"
  )
  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "umi",
      UMI_threshold = NA_real_,
      verbose = FALSE
    ),
    "UMI_threshold.*finite.*number"
  )
  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "mito",
      mito_threshold = Inf,
      verbose = FALSE
    ),
    "mito_threshold.*finite.*number"
  )
  expect_error(
    RunSpotQC(
      srt,
      assay = "RNA",
      return_filtered = NA,
      verbose = FALSE
    ),
    "return_filtered.*TRUE or FALSE"
  )
})

test_that("RunSpotQC gives a clear error when filtering would return no spots", {
  srt <- make_receipt_spatial_object()
  expect_error(
    suppressWarnings(RunSpotQC(
      srt,
      assay = "RNA",
      qc_metrics = "umi",
      UMI_threshold = 1e9,
      return_filtered = TRUE,
      verbose = FALSE
    )),
    "No spots passed QC.*return_filtered = FALSE"
  )
})

test_that("RunSpotQC is silent on request and failure has no receipt", {
  srt <- make_receipt_spatial_object()
  quiet <- NULL
  quiet_messages <- suppressWarnings(testthat::capture_messages(
    quiet <- RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = FALSE
    )
  ))
  expect_length(quiet_messages, 0L)
  expect_s4_class(quiet, "Seurat")

  failure_messages <- testthat::capture_messages(
    expect_error(RunSpotQC(list(), verbose = TRUE), "Seurat")
  )
  expect_false(grepl("Spot QC completed", receipt_plain(failure_messages), fixed = TRUE))
})

test_that("quiet RunSpotQC skips receipt image and coordinate probes", {
  srt <- make_receipt_spatial_object()
  image_probes <- coordinate_probes <- 0L
  testthat::local_mocked_bindings(
    Images = function(...) {
      image_probes <<- image_probes + 1L
      stop("quiet RunSpotQC must not inspect images")
    },
    .package = "SeuratObject"
  )
  testthat::local_mocked_bindings(
    spot_qc_has_plottable_coordinates = function(...) {
      coordinate_probes <<- coordinate_probes + 1L
      stop("quiet RunSpotQC must not inspect coordinates")
    },
    .package = "scop"
  )

  run_quiet <- function(verbose) {
    suppressWarnings(RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = verbose
    ))
  }

  expect_s4_class(run_quiet(FALSE), "Seurat")
  old_options <- options(log_message.verbose = FALSE)
  on.exit(options(old_options), add = TRUE)
  expect_s4_class(run_quiet(NULL), "Seurat")
  expect_identical(image_probes, 0L)
  expect_identical(coordinate_probes, 0L)
})

test_that("RunSpatialVariableFeatures reports exactly what the returned object retains", {
  skip_if_not_installed("BiocNeighbors")
  cases <- list(
    both = list(set = TRUE, store = TRUE),
    features_only = list(set = TRUE, store = FALSE),
    tool_only = list(set = FALSE, store = TRUE),
    neither = list(set = FALSE, store = FALSE)
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    messages <- suppressWarnings(testthat::capture_messages(
      out <- RunSpatialVariableFeatures(
        make_receipt_spatial_object(),
        assay = "RNA",
        layer = "counts",
        method = "moran",
        coord.cols = c("col", "row"),
        k = 1,
        nfeatures = 2,
        min_spots = 1,
        set_variable_features = case$set,
        store_results = case$store,
        backend = "r",
        verbose = TRUE
      )
    ))
    plain <- receipt_plain(messages)

    expect_true(grepl(
      "Spatial variable features completed: 3 tested, 2 ranked in the top set",
      plain,
      fixed = TRUE
    ), info = case_name)
    expect_false(grepl("Stored 2 spatial variable features", plain, fixed = TRUE))
    expect_identical(
      !is.null(out@tools[["SpatialVariableFeatures"]]),
      case$store,
      info = case_name
    )
    expect_identical(
      length(SeuratObject::VariableFeatures(out, assay = "RNA")) == 2L,
      case$set,
      info = case_name
    )

    if (isTRUE(case$store)) {
      expect_true(grepl(
        "returned object tool bundle `SpatialVariableFeatures`",
        plain,
        fixed = TRUE
      ), info = case_name)
      expect_true(grepl(
        "Plot returned object `SpatialVariableFeaturePlot(<returned_object>",
        plain,
        fixed = TRUE
      ), info = case_name)
    } else {
      expect_false(grepl("SpatialVariableFeaturePlot", plain, fixed = TRUE))
    }
    if (isTRUE(case$set)) {
      expect_true(grepl("top features as assay `VariableFeatures`", plain, fixed = TRUE))
    }
    if (!isTRUE(case$set) && !isTRUE(case$store)) {
      expect_true(grepl("results not retained", plain, fixed = TRUE))
      expect_false(grepl("Saved", plain, fixed = TRUE))
    }
    expect_false(grepl("srt", plain, fixed = TRUE))
  }
})

test_that("RunSpatialVariableFeatures Plot hint preserves resolved plotting context", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_receipt_svf_multi_image_object()
  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpatialVariableFeatures(
      srt,
      assay = "ALT",
      layer = "counts",
      method = "moran",
      image = "slice1",
      coord.cols = c("ignored_x", "ignored_y"),
      k = 1,
      nfeatures = 2,
      min_spots = 1,
      backend = "r",
      verbose = TRUE
    )
  ))
  plot_call <- paste0(
    "SpatialVariableFeaturePlot(<returned_object>, plot_type = \"combined\", ",
    "assay = \"ALT\", image = \"slice1\", coord.cols = c(\"x\", \"y\"))"
  )
  expect_receipt_plot_hint(messages, plot_call)
  expect_identical(
    out@tools[["SpatialVariableFeatures"]]$parameters$assay,
    "ALT"
  )
  expect_error(
    SpatialVariableFeaturePlot(out, plot_type = "combined"),
    "Multiple spatial images"
  )
  executable_call <- sub("<returned_object>", "out", plot_call, fixed = TRUE)
  expect_no_error(eval(parse(text = executable_call)))
})

test_that("RunSpotQC gives an all-image hint with an explicit result placeholder", {
  srt <- make_receipt_multi_image_object()
  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- RunSpotQC(
      srt,
      assay = "RNA",
      UMI_threshold = 0,
      gene_threshold = 0,
      mito_threshold = 100,
      verbose = TRUE
    )
  ))
  plot_call <- paste0(
    "lapply(SeuratObject::Images(<returned_object>), function(image) ",
    "SpatialSpotPlot(<returned_object>, group.by = \"SpotQC\", image = image))"
  )
  plain <- receipt_plain(messages)
  expect_receipt_plot_hint(messages, plot_call)
  expect_s4_class(out, "Seurat")
  expect_false(grepl("srt", plain, fixed = TRUE))
})

test_that("RunSpatialNetwork emits actual graph state and an explicit result placeholder", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_receipt_spatial_object()
  out <- NULL
  messages <- testthat::capture_messages(
    out <- RunSpatialNetwork(srt, k = 1, verbose = TRUE)
  )
  plain <- receipt_plain(messages)
  expect_true(grepl(
    'Spatial graph completed: "knn", `k` = 1, 4 nodes, 3 edges',
    plain,
    fixed = TRUE
  ))
  expect_true(grepl("Scope raw coordinates; 4 observations", plain, fixed = TRUE))
  expect_true(grepl('Saved graph "knn_k1"', plain, fixed = TRUE))
  expect_receipt_plot_hint(
    messages,
    'SpatialNetworkPlot(<returned_object>, graph.name = "knn_k1")'
  )

  replaced <- NULL
  replaced_messages <- testthat::capture_messages(
    replaced <- RunSpatialNetwork(out, k = 1, overwrite = TRUE, verbose = TRUE)
  )
  replaced_plain <- receipt_plain(replaced_messages)
  expect_true(grepl('Replaced graph "knn_k1"', replaced_plain, fixed = TRUE))
  expect_true("knn_k1" %in% names(replaced@tools[["SpatialNetwork"]]$graphs))
})

test_that("RunSpatialNetwork reports radius and safely quotes a custom graph name", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_receipt_spatial_object()
  graph_name <- 'radius "quoted"\\graph'
  out <- NULL
  messages <- testthat::capture_messages(
    out <- RunSpatialNetwork(
      srt,
      method = "radius",
      radius = 1.5,
      graph.name = graph_name,
      verbose = TRUE
    )
  )
  plain <- receipt_plain(messages)
  expect_true(grepl(
    'Spatial graph completed: "radius", `radius` = 1.5',
    plain,
    fixed = TRUE
  ))
  expect_true(graph_name %in% names(out@tools[["SpatialNetwork"]]$graphs))
  plot_call <- paste0(
    "SpatialNetworkPlot(<returned_object>, graph.name = ",
    deparse1(graph_name, width.cutoff = 500L),
    ")"
  )
  expect_receipt_plot_hint(messages, plot_call)
})

test_that("RunSpatialNetwork respects verbose NULL and emits no receipt on failure", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_receipt_spatial_object()
  old_options <- options(log_message.verbose = FALSE)
  on.exit(options(old_options), add = TRUE)
  quiet <- NULL
  quiet_messages <- testthat::capture_messages(
    quiet <- RunSpatialNetwork(srt, k = 1, verbose = NULL)
  )
  expect_length(quiet_messages, 0L)
  expect_s4_class(quiet, "Seurat")

  failure_messages <- testthat::capture_messages(
    expect_error(
      RunSpatialNetwork(srt, k = ncol(srt), verbose = TRUE),
      "number of nodes minus one"
    )
  )
  expect_false(grepl(
    "Spatial graph completed",
    receipt_plain(failure_messages),
    fixed = TRUE
  ))
})
