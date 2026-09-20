edge_contract_object <- function() {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(
    matrix(c(0, 0, 0, 1:9), 3, dimnames = list(paste0("g", 1:3), letters[1:4])), sparse = TRUE))
  object$x <- c(0, 1, 0, 1)
  object$y <- c(0, 0, 1, 1)
  object
}

test_that("Cell2location stores NA rows for filtered spots and plots after reload", {
  object <- edge_contract_object()
  signatures <- matrix(1:6, 3, dimnames = list(rownames(object), c("A", "B")))
  result_dir <- tempfile("cell2location_filtered_")
  dir.create(result_dir)
  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) TRUE,
    conda_python = function(...) file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "R"),
    resolve_conda = function(...) "mamba",
    runner_script_path = function(...) "cell2location.py",
    runner_write_json = function(...) invisible(NULL),
    runner_read_json = function(...) list(status = "complete"),
    srt_to_h5ad = function(object, path, ...) { file.create(path); invisible(path) },
    runner_system2 = function(command, args, env, stdout, stderr) {
      files <- cell2location_result_files(result_dir)
      dir.create(dirname(files$abundance), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$signatures), recursive = TRUE, showWarnings = FALSE)
      abundance <- matrix(c(2, 3, 4, 8, 7, 6), 3,
        dimnames = list(c("b", "c", "d"), c("A", "B")))
      utils::write.csv(abundance, files$abundance)
      utils::write.csv(abundance / rowSums(abundance), files$proportions)
      utils::write.csv(signatures, files$signatures)
      file.create(files$manifest, stdout, stderr)
      0L
    }
  )
  result <- RunCell2location(object, reference_signatures = signatures,
    result_dir = result_dir, assay = "RNA", verbose = FALSE)
  expect_identical(colnames(result), colnames(object))
  expect_identical(result@tools$Cell2location$cells, c("b", "c", "d"))
  expect_identical(result@tools$Cell2location$input_summary$dropped_spots, "a")
  for (name in c("abundance", "proportions")) {
    expect_identical(rownames(result@tools$Cell2location[[name]]), colnames(object))
    expect_true(all(is.na(result@tools$Cell2location[[name]]["a", ])))
  }
  file <- tempfile(fileext = ".rds")
  saveRDS(result, file)
  reloaded <- readRDS(file)
  for (mode in c("proportion", "abundance")) {
    plots <- Cell2locationPlot(reloaded, plot_type = mode, combine = FALSE, show_na = TRUE)
    expect_named(plots, c("A", "B"))
    expect_true(is.na(plots$A$data$.value[1]))
    expect_equal(plots$A$data$.value[-1], if (mode == "proportion") c(.2, .3, .4) else c(2, 3, 4))
    expect_no_error(ggplot2::ggplotGrob(plots$A))
  }
  expect_s3_class(SpatialDeconvolutionPlot(reloaded, tool_name = "Cell2location"), "patchwork")
})

test_that("historical Cell2location rows require complete modeled/dropped provenance", {
  object <- edge_contract_object()
  values <- matrix(c(.4, .2, .3, .6, .8, .7), 3,
    dimnames = list(c("d", "b", "c"), c("A", "B")))
  object@tools$Cell2location <- list(proportions = values, abundance = values * 10,
    cells = c("b", "c", "d"), input_summary = list(dropped_spots = "a"))
  before <- object
  p <- Cell2locationPlot(object, combine = FALSE, show_na = TRUE)
  expect_equal(p$A$data$.value, c(NA, .2, .3, .4))
  expect_identical(object, before)
  for (change in c("undeclared", "unknown", "duplicate", "missing_model", "overlap")) {
    bad <- object
    if (change == "undeclared") bad@tools$Cell2location$input_summary <- NULL
    if (change == "unknown") rownames(bad@tools$Cell2location$proportions)[1] <- "other"
    if (change == "duplicate") rownames(bad@tools$Cell2location$proportions)[1] <- "b"
    if (change == "missing_model") bad@tools$Cell2location$proportions <- values[-2, , drop = FALSE]
    if (change == "overlap") bad@tools$Cell2location$input_summary$dropped_spots <- c("a", "b")
    expect_error(Cell2locationPlot(bad), "invalid|unknown|incomplete")
  }
})

test_that("long points and jitter honor explicit scales and named-list returns", {
  object <- edge_contract_object()
  long <- data.frame(spot = colnames(object), score = c(1, 2, 40, 100))
  for (geom in c("point", "jitter")) {
    for (combine in c(TRUE, FALSE)) {
      p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "score",
        geom = geom, combine = combine, upper_cutoff = 10)
      if (!combine) { expect_named(p, "score"); p <- p[[1]] }
      expect_equal(p$scales$get_scales("colour")$limits, c(1, 10))
    }
    p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "score", geom = geom)
    expect_equal(p$scales$get_scales("colour")$limits, c(1, 100))
    p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "score", geom = geom, upper_quantile = .5)
    expect_equal(p$scales$get_scales("colour")$limits, c(1, 21))
  }
  # An ID column named x must not be overwritten before resolving y.
  names(long)[1] <- "x"
  p <- SpatialSpotPlot(object, plot.data = long, spot.by = "x", color.by = "score")
  expect_equal(unname(p$data$y), unname(object$y))
})

test_that("generic dominant maps honor combine and preserve full input", {
  object <- edge_contract_object()
  object@tools$Demo <- list(proportions = matrix(c(.2, .8, .3, .7, .8, .2, .7, .3), 4,
    dimnames = list(colnames(object), c("A", "B"))))
  before <- object
  p <- SpatialDeconvolutionPlot(object, "Demo", plot_type = "dominant", combine = FALSE)
  expect_named(p, "value")
  expect_s3_class(p[[1]], "ggplot")
  expect_identical(object, before)
  expect_s3_class(SpatialDeconvolutionPlot(object, "Demo", plot_type = "dominant"), "ggplot")
})

test_that("PRECAST light storage preserves numbers and roundtrip while reducing payload", {
  object <- edge_contract_object()
  object$sample <- c("S1", "S1", "S2", "S2")
  cells <- colnames(object)
  result <- list(domains = stats::setNames(c("D1", "D1", "D2", "D2"), cells),
    embedding = matrix(1:8, 4, dimnames = list(cells, c("PC1", "PC2"))),
    aligned_coords = data.frame(x = object$x + 1, y = object$y, row.names = cells),
    features = rownames(object), raw_result = list(trace = sin(seq_len(5000))))
  args <- list(srt = object, result = result, method = "PRECAST", sample.by = "sample",
    assay = "RNA", layer = "counts", image = NULL, coord.cols = c("x", "y"),
    coordinate_space = "raw", reduction.name = "PRECAST", cluster_colname = "domain",
    tool_name = "SpatialIntegration", store_results = TRUE)
  full <- do.call(spatial_integration_apply_result, c(args, list(store_object = TRUE)))
  light <- do.call(spatial_integration_apply_result, c(args, list(store_object = FALSE)))
  expect_identical(full@meta.data, light@meta.data)
  expect_identical(full@reductions, light@reductions)
  for (field in c("embedding", "domains", "aligned_coords", "summary")) {
    expect_identical(full@tools$SpatialIntegration$methods$PRECAST[[field]], light@tools$SpatialIntegration$methods$PRECAST[[field]])
  }
  files <- c(tempfile(fileext = ".rds"), tempfile(fileext = ".rds"))
  saveRDS(full, files[1]); saveRDS(light, files[2])
  expect_lt(unname(file.info(files[2])$size), unname(file.info(files[1])$size))
  reloaded <- readRDS(files[2])
  expect_identical(reloaded@reductions, light@reductions)
  expect_identical(reloaded@tools$SpatialIntegration, light@tools$SpatialIntegration)
  expect_no_error(ggplot2::ggplotGrob(SpatialIntegrationPlot(reloaded)))
})
