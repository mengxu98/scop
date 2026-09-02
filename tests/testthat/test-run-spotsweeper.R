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

with_mock_spotsweeper <- function(
  code,
  artifact_mode = c("success", "error", "missing", "invalid")
) {
  artifact_mode <- match.arg(artifact_mode)
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
  failing_artifact <- function(...) {
    stop("mock artifact failure")
  }
  missing_artifact <- function(spe, ...) {
    spe
  }
  invalid_artifact <- function(spe, ...) {
    cdata <- as.data.frame(SummarizedExperiment::colData(spe))
    cdata$artifact <- rep(NA, nrow(cdata))
    SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(cdata)
    spe
  }
  selected_artifact <- switch(artifact_mode,
    success = fake_artifact,
    error = failing_artifact,
    missing = missing_artifact,
    invalid = invalid_artifact
  )
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      allowed <- c(
        "MicTott/SpotSweeper", "SpatialExperiment",
        "SummarizedExperiment", "S4Vectors"
      )
      expect_true(length(packages) > 0L)
      expect_true(all(packages %in% allowed))
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      if (identical(package, "SpatialExperiment")) {
        return(getExportedValue(package, name))
      }
      if (identical(package, "SpotSweeper")) {
        return(switch(name,
          localOutliers = fake_local,
          findArtifacts = selected_artifact,
          stop("unexpected SpotSweeper function")
        ))
      }
      stop("unexpected namespace function")
    }
  )
  force(code)
}

capture_spotsweeper_logs <- function(code) {
  events <- list()
  testthat::local_mocked_bindings(
    .package = "scop",
    log_message = function(
      ...,
      verbose = NULL,
      message_type = c("info", "success", "warning", "error", "running", "ask")
    ) {
      message_type <- match.arg(message_type)
      text <- paste(vapply(
        list(...),
        function(x) paste(as.character(x), collapse = " "),
        character(1)
      ), collapse = " ")
      if (identical(message_type, "error")) {
        stop(text, call. = FALSE)
      }
      if (!isFALSE(verbose)) {
        events[[length(events) + 1L]] <<- list(
          type = message_type,
          text = text
        )
      }
      invisible(NULL)
    }
  )
  force(code)
  events
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
  expect_identical(out@tools$SpotSweeper$status, "completed")
  expect_true(all(out@tools$SpotSweeper$artifact_status$status == "completed"))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) %in% c("Pass", "Fail")))
  expect_true(all(as.character(out$SpotSweeper_QC) %in% c("Pass", "Fail")))
  expect_identical(out@tools$SpotSweeper$parameters$coordinate_space, "raw")
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
  expect_true(all(is.na(out$SpotSweeper_artifact)))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
  expect_equal(as.character(out$SpotSweeper_QC[1:2]), c("Fail", "Fail"))
  expect_true(all(as.character(out$SpotSweeper_QC[3:6]) == "Partial"))
  expect_identical(out@tools$SpotSweeper$status, "partial")
  expect_true(all(out@tools$SpotSweeper$artifact_status$status == "skipped"))
  expect_true("skip_reason" %in% colnames(out@tools$SpotSweeper$artifacts))
  expect_true(all(!is.na(out@tools$SpotSweeper$artifacts$skip_reason)))
})

test_that("RunSpotSweeper marks artifact backend failures as partial", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  out <- suppressWarnings(with_mock_spotsweeper({
    RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      n_order = 2,
      verbose = FALSE
    )
  }, artifact_mode = "error"))

  expect_identical(out@tools$SpotSweeper$status, "partial")
  expect_true(all(out@tools$SpotSweeper$artifact_status$status == "failed"))
  expect_true(all(is.na(out$SpotSweeper_artifact)))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
  expect_true(all(as.character(out$SpotSweeper_QC[3:6]) == "Partial"))
  expect_true(all(grepl("mock artifact failure", out@tools$SpotSweeper$artifact_status$reason)))
})

test_that("RunSpotSweeper rejects missing or invalid artifact backend output truthfully", {
  skip_if_no_spotsweeper_infra()
  expected_reason <- c(
    missing = "returned no artifact column",
    invalid = "returned invalid artifact values"
  )
  for (artifact_mode in names(expected_reason)) {
    srt <- make_spotsweeper_seurat()
    out <- suppressWarnings(with_mock_spotsweeper({
      RunSpotSweeper(
        srt,
        layer = "counts",
        coord.cols = c("x", "y"),
        sample.by = "sample",
        n_neighbors = 2,
        n_order = 2,
        verbose = FALSE
      )
    }, artifact_mode = artifact_mode))

    status <- out@tools$SpotSweeper$artifact_status
    expect_identical(out@tools$SpotSweeper$status, "partial")
    expect_true(all(status$status == "failed"))
    expect_true(all(grepl(expected_reason[[artifact_mode]], status$reason)))
    expect_true(all(is.na(out$SpotSweeper_artifact)))
    expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
    expect_true(all(as.character(out$SpotSweeper_QC[3:6]) == "Partial"))
  }
})

test_that("RunSpotSweeper preserves mixed per-sample artifact statuses", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  srt$artifact_mito_percent <- c(1, 2, 3, 1, 1, 1)
  srt$artifact_mito_sum <- c(1, 2, 3, 1, 1, 1)
  out <- suppressWarnings(with_mock_spotsweeper({
    RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      mito_percent = "artifact_mito_percent",
      mito_sum = "artifact_mito_sum",
      n_neighbors = 2,
      n_order = 2,
      verbose = FALSE
    )
  }))

  status <- out@tools$SpotSweeper$artifact_status
  expect_equal(status$sample, c("S1", "S2"))
  expect_equal(status$status, c("completed", "skipped"))
  expect_equal(status$n_spots, c(3L, 3L))
  expect_equal(status$n_artifacts, c(0L, NA_integer_))
  expect_identical(out@tools$SpotSweeper$status, "partial")
  expect_true(all(!is.na(out$SpotSweeper_artifact[1:3])))
  expect_true(all(is.na(out$SpotSweeper_artifact[4:6])))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc[1:3]) == "Pass"))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc[4:6]) == "NotEvaluated"))
  expect_equal(
    as.character(out$SpotSweeper_QC),
    c("Fail", "Pass", "Pass", "Partial", "Partial", "Partial")
  )
})

test_that("RunSpotSweeper preserves custom sample order when artifacts are not requested", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  srt$sample <- rep(c("S2", "S1"), each = 3)
  out <- with_mock_spotsweeper({
    RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      run_artifact = FALSE,
      store_results = TRUE,
      verbose = FALSE
    )
  })

  status <- out@tools$SpotSweeper$artifact_status
  expect_identical(out@tools$SpotSweeper$status, "completed")
  expect_equal(status$sample, c("S2", "S1"))
  expect_true(all(status$status == "not_requested"))
  expect_equal(status$n_spots, c(3L, 3L))
  expect_true(all(is.na(status$n_artifacts)))
  expect_true(all(out$SpotSweeper_artifact_status == "not_requested"))
  expect_true(all(is.na(out$SpotSweeper_artifact)))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
})

test_that("RunSpotSweeper refuses filtered output when artifact detection is partial", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat()
  before <- list(
    meta.data = srt@meta.data,
    assays = srt@assays,
    reductions = srt@reductions,
    graphs = srt@graphs,
    tools = srt@tools
  )
  condition <- tryCatch(
    suppressWarnings(with_mock_spotsweeper({
      RunSpotSweeper(
        srt,
        layer = "counts",
        coord.cols = c("x", "y"),
        sample.by = "sample",
        n_neighbors = 2,
        n_order = 2,
        return_filtered = TRUE,
        verbose = FALSE
      )
    }, artifact_mode = "error")),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "fully filtered SpotSweeper object")
  expect_identical(srt@meta.data, before$meta.data)
  expect_identical(srt@assays, before$assays)
  expect_identical(srt@reductions, before$reductions)
  expect_identical(srt@graphs, before$graphs)
  expect_identical(srt@tools, before$tools)
})

test_that("RunSpotSweeper keeps partial metadata when detailed storage is disabled", {
  skip_if_no_spotsweeper_infra()
  srt <- make_spotsweeper_seurat(include_mito = FALSE)
  out <- with_mock_spotsweeper({
    RunSpotSweeper(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      n_order = 2,
      store_results = FALSE,
      verbose = FALSE
    )
  })

  expect_false("SpotSweeper" %in% names(out@tools))
  expect_true(all(out$SpotSweeper_artifact_status == "skipped"))
  expect_true(all(is.na(out$SpotSweeper_artifact)))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
  expect_equal(
    as.character(out$SpotSweeper_QC),
    c("Fail", "Fail", "Partial", "Partial", "Partial", "Partial")
  )
})

test_that("RunSpotSweeper console status is truthful and respects verbose", {
  skip_if_no_spotsweeper_infra()

  partial_events <- capture_spotsweeper_logs(with_mock_spotsweeper({
    RunSpotSweeper(
      make_spotsweeper_seurat(include_mito = FALSE),
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      verbose = TRUE
    )
  }))
  expect_true(any(vapply(
    partial_events,
    function(x) identical(x$type, "running") && grepl("completed partially", x$text),
    logical(1)
  )))
  expect_false(any(vapply(
    partial_events,
    function(x) identical(x$type, "success"),
    logical(1)
  )))

  local_only_events <- capture_spotsweeper_logs(with_mock_spotsweeper({
    RunSpotSweeper(
      make_spotsweeper_seurat(),
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      run_artifact = FALSE,
      verbose = TRUE
    )
  }))
  expect_true(any(vapply(
    local_only_events,
    function(x) {
      identical(x$type, "success") &&
        grepl("local-outlier QC", x$text) &&
        grepl("not requested", x$text)
    },
    logical(1)
  )))

  silent_events <- capture_spotsweeper_logs(with_mock_spotsweeper({
    RunSpotSweeper(
      make_spotsweeper_seurat(),
      layer = "counts",
      coord.cols = c("x", "y"),
      sample.by = "sample",
      n_neighbors = 2,
      verbose = FALSE
    )
  }, artifact_mode = "error"))
  expect_length(silent_events, 0L)
})

test_that("RunSpotSweeper partial completion survives warning escalation", {
  skip_if_no_spotsweeper_infra()
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)

  cases <- list(
    readiness = list(include_mito = FALSE, artifact_mode = "success", status = "skipped"),
    backend_error = list(include_mito = TRUE, artifact_mode = "error", status = "failed"),
    missing_output = list(include_mito = TRUE, artifact_mode = "missing", status = "failed"),
    invalid_output = list(include_mito = TRUE, artifact_mode = "invalid", status = "failed")
  )
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    out <- NULL
    expect_no_error(
      out <- suppressMessages(with_mock_spotsweeper({
        RunSpotSweeper(
          make_spotsweeper_seurat(include_mito = case$include_mito),
          layer = "counts",
          coord.cols = c("x", "y"),
          sample.by = "sample",
          n_neighbors = 2,
          verbose = TRUE
        )
      }, artifact_mode = case$artifact_mode)),
      message = case_name
    )
    expect_s4_class(out, "Seurat")
    expect_identical(out@tools$SpotSweeper$status, "partial", info = case_name)
    expect_true(
      all(out@tools$SpotSweeper$artifact_status$status == case$status),
      info = case_name
    )
  }
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
  expect_true(all(is.na(out$SpotSweeper_artifact)))
  expect_true(all(as.character(out$SpotSweeper_artifact_qc) == "NotEvaluated"))
  expect_true(all(out$SpotSweeper_artifact_status == "not_requested"))
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
