input_contract_object <- function() {
  counts <- matrix(c(1:6, 6:1, 2, 1, 2, 1, 2, 1), 3, byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), paste0("s", 1:6)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "RNA")
  object[["Spatial"]] <- SeuratObject::CreateAssayObject(counts = counts * 10)
  object$x <- 1:6
  object$y <- c(0, 0, 1, 1, 2, 2)
  object$sample <- rep(c("A", "B"), each = 3)
  object
}

test_that("PRECAST consumes the selected counts rather than the default assay", {
  object <- input_contract_object()
  selected <- GetAssayData5(object, assay = "Spatial", layer = "counts")
  before <- object
  captured <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "PRECAST")
      switch(name,
        CreatePRECASTObject = function(seuList, ...) {
          captured <<- seuList
          stop("captured validated input")
        },
        function(...) stop("training must not run"))
    }
  )
  expect_error(RunSpatialIntegration(object, sample.by = "sample", assay = "Spatial",
    verbose = FALSE), "captured validated input")
  for (sample in names(captured)) {
    cells <- colnames(captured[[sample]])
    expect_identical(SeuratObject::DefaultAssay(captured[[sample]]), "Spatial")
    expect_equal(GetAssayData5(captured[[sample]], layer = "counts"), selected[, cells])
    expect_equal(unname(captured[[sample]]$row), unname(object$y[cells]))
    expect_equal(unname(captured[[sample]]$col), unname(object$x[cells]))
  }
  expect_identical(object, before)
  normalized <- selected / 3.7
  object <- SeuratObject::SetAssayData(object, assay = "Spatial", layer = "data", new.data = normalized)
  expect_error(RunSpatialIntegration(object, sample.by = "sample", assay = "Spatial",
    layer = "data", verbose = FALSE), "integer counts")
})

test_that("PRECAST rejects coordinate-incompatible and managed backend arguments", {
  object <- input_contract_object()
  testthat::local_mocked_bindings(check_r = function(...) stop("backend touched"))
  expect_error(RunSpatialIntegration(object, sample.by = "sample",
    adj_params = list(type = "fixed_distance"), verbose = FALSE), "array indices")
  expect_error(RunSpatialIntegration(object, sample.by = "sample",
    create_params = list(seuList = list()), verbose = FALSE), "managed")
  expect_error(RunSpatialIntegration(object, sample.by = "sample",
    adj_params = list(number = -1), verbose = FALSE), "positive integer")
})

test_that("invalid integration payloads fail before object mutation", {
  object <- input_contract_object()
  before <- object
  cells <- colnames(object)
  bad <- list(
    list(domains = setNames(rep(NA_character_, length(cells)), cells)),
    list(domains = setNames(rep(" ", length(cells)), cells)),
    list(embedding = matrix(NA_real_, length(cells), 2, dimnames = list(cells, NULL))),
    list(embedding = matrix(Inf, length(cells), 2, dimnames = list(cells, NULL))),
    list(aligned_coords = data.frame(x = rep(NA_real_, length(cells)), y = 1, row.names = cells))
  )
  for (payload in bad) {
    testthat::local_mocked_bindings(spatial_integration_run_backend = function(...) payload)
    expect_error(RunSpatialIntegration(object, sample.by = "sample", verbose = FALSE),
      "non-missing|finite")
    expect_identical(object, before)
  }
  labels <- data.frame(domain = rep(c("a", "b"), 3), row.names = rev(cells))
  aligned <- spatial_integration_standardize_named_vector(labels, cells)
  expect_identical(unname(aligned), rev(labels$domain))
})

test_that("SVG surfaces inherit the saved assay, layer and custom coordinates", {
  object <- input_contract_object()
  object$custom_x <- object$x + 100
  object$custom_y <- object$y + 200
  object <- RunSpatialVariableFeatures(object, assay = "Spatial", layer = "counts",
    coord.cols = c("custom_x", "custom_y"), min_spots = 1, nfeatures = 2,
    k = 2, backend = "r", verbose = FALSE)
  SeuratObject::DefaultAssay(object) <- "RNA"
  path <- tempfile(fileext = ".rds")
  saveRDS(object, path)
  object <- readRDS(path)
  p <- SpatialVariableFeaturePlot(object, plot_type = "surface", features = "g1",
    overlay_image = FALSE, combine = FALSE)[[1]]
  expect_equal(p$data$.value, 10 * (1:6))
  expect_equal(unname(p$data$x), unname(object$custom_x))
  p <- SpatialVariableFeaturePlot(object, plot_type = "surface", features = "g1",
    assay = "RNA", coord.cols = c("x", "y"), overlay_image = FALSE, combine = FALSE)[[1]]
  expect_equal(p$data$.value, as.numeric(1:6))
  expect_equal(unname(p$data$x), unname(object$x))
})

test_that("SVG backend identity and statistics are checked before storage", {
  object <- input_contract_object()
  before <- object
  result <- data.frame(feature = rownames(object), p_value = c(.01, .03, .7))
  malformed <- list(
    transform(result, feature = c("g1", "g1", "g3")),
    transform(result, feature = c("g1", "unknown", "g3")),
    transform(result, p_value = c(-.1, .1, .2)),
    transform(result, q_value = c(.1, Inf, .3)),
    transform(result, p_value = c("bad", ".1", ".2"))
  )
  for (payload in malformed) {
    testthat::local_mocked_bindings(spatial_variable_run_sparkx = function(...) payload)
    expect_error(RunSpatialVariableFeatures(object, method = "SPARKX", layer = "counts",
      min_spots = 1, verbose = FALSE), "feature IDs|\\[0, 1\\]|non-numeric")
    expect_identical(object, before)
  }
  expr <- GetAssayData5(object, layer = "counts")
  tested <- setNames(rep(6L, nrow(expr)), rownames(expr))
  result$q_value <- c(.9, NA, NA)
  out <- spatial_variable_finalize_result(result, expr, tested, "SPARKX")
  expect_equal(out$q_value[match("g2", out$feature)], p.adjust(result$p_value, "BH")[2])
  expect_equal(out$q_value[match("g1", out$feature)], .9)
  expect_error(spatial_variable_result_features(data.frame(pval = c(.1, .2)), rownames(expr)),
    "explicit feature IDs")
})

test_that("SVG surfaces resolve the saved selection in a multiple-image object", {
  data(visium_human_pancreas_sub, package = "scop")
  object <- suppressWarnings(visium_human_pancreas_sub[1:8, 1:12])
  object[["slice2"]] <- object[["slice1"]]
  expect_error(RunSpatialVariableFeatures(object, assay = "Spatial", layer = "counts",
    min_spots = 1, backend = "r", verbose = FALSE), "Multiple spatial images")
  out <- RunSpatialVariableFeatures(object, assay = "Spatial", layer = "counts",
    image = "slice2", min_spots = 1, k = 3, backend = "r", verbose = FALSE)
  expect_identical(out@tools$SpatialVariableFeatures$parameters$image, "slice2")
  expect_s3_class(SpatialVariableFeaturePlot(out, plot_type = "surface",
    features = rownames(out)[1], overlay_image = FALSE), "ggplot")
})

test_that("Visium HD never inherits ordinary Visium physical calibration", {
  data(visium_human_pancreas_sub, package = "scop")
  object <- suppressWarnings(visium_human_pancreas_sub[, 1:6])
  object@misc$scop_spatial_input <- list(slice1 = list(technology = "visium_hd"))
  expect_identical(spatialcellchat_detect_technology(object, "auto", "slice1"), "visium_hd")
  expect_identical(spatialcellchat_detect_level(object, "auto", "visium_hd", NULL, "slice1"), "spot")
  object[["slice1"]]@scale.factors$spot <- 4
  args <- list(srt = object, cells = colnames(object), image = "slice1",
    coord.cols = c("x", "y"), technology = "visium_hd", coordinate.unit = "pixel", ratio = NULL, tol = NULL)
  expect_error(do.call(spatialcellchat_metric_coordinates, args), "explicit calibrated")
  args$ratio <- .5
  expect_error(do.call(spatialcellchat_metric_coordinates, args), "explicit.*tol")
  args$tol <- 1
  metric <- do.call(spatialcellchat_metric_coordinates, args)
  expect_equal(metric$data$x, metric$data$x_raw * .5)
  expect_equal(metric$source$scale_to_micron, .5)
  expect_equal(metric$source$tol_um, 1)
  args$technology <- "visium"
  args["ratio"] <- list(NULL)
  args["tol"] <- list(NULL)
  ordinary <- do.call(spatialcellchat_metric_coordinates, args)
  expect_equal(ordinary$source$scale_to_micron, 65 / 4)
  expect_equal(ordinary$source$tol_um, 32.5)
  object@misc$scop_spatial_input <- NULL
  object[["Spatial.008um"]] <- object[["Spatial"]]
  SeuratObject::DefaultAssay(object[["slice1"]]) <- "Spatial.008um"
  expect_identical(spatialcellchat_detect_technology(object, "auto", "slice1"), "visium_hd")
})
