make_spatial_gradient_seurat <- function() {
  counts <- matrix(
    c(
      1, 0, 2, 3,
      0, 4, 1, 0,
      5, 1, 0, 2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(1, 2, 1, 2)
  srt$y <- c(1, 1, 2, 2)
  srt
}

make_spatial_gradient_result <- function() {
  list(
    screening = data.frame(
      variable = rep(c("Gene1", "Gene2"), each = 3),
      distance = rep(c(0, 1, 2), times = 2),
      value = c(1, 2, 3, 3, 2, 1),
      estimate = c(1.1, 1.9, 3.1, 3.0, 2.1, 0.9),
      reference = "trajectory",
      mode = "trajectory",
      stringsAsFactors = FALSE
    ),
    significance = data.frame(
      variable = c("Gene1", "Gene2"),
      tot_var = c(0.8, 0.7),
      p_value = c(0.001, 0.01),
      fdr = c(0.002, 0.02),
      stringsAsFactors = FALSE
    ),
    model_fits = data.frame(
      variable = c("Gene1", "Gene1", "Gene2"),
      model = c("ascending", "descending", "descending"),
      mae = c(0.1, 0.4, 0.2),
      rmse = c(0.2, 0.5, 0.3),
      r2 = c(0.95, 0.6, 0.9),
      stringsAsFactors = FALSE
    ),
    top_variables = data.frame(
      variable = c("Gene1", "Gene2"),
      rank = c(1, 2),
      fdr = c(0.002, 0.02),
      best_model = c("ascending", "descending"),
      mae = c(0.1, 0.2),
      rmse = c(0.2, 0.3),
      stringsAsFactors = FALSE
    ),
    parameters = data.frame(
      key = c("reference", "assay", "seed"),
      value = c("trajectory", "RNA", "123"),
      stringsAsFactors = FALSE
    )
  )
}

test_that("SPATA2 output tables are standardized to stable column names", {
  significance <- sgf_standardize_significance(data.frame(
    variables = "Gene1",
    p.value = 0.01,
    q_value = 0.02,
    tot_var = 0.8
  ))
  expect_named(significance[1:4], c("variable", "tot_var", "p_value", "fdr"))
  expect_equal(significance$variable, "Gene1")
  expect_equal(significance$fdr, 0.02)

  model_fits <- sgf_standardize_model_fits(data.frame(
    genes = "Gene1",
    models = "ascending",
    mae = 0.1,
    rmse = 0.2
  ))
  expect_named(model_fits[1:4], c("variable", "model", "mae", "rmse"))
  expect_equal(model_fits$model, "ascending")

  screening <- sgf_standardize_screening_df(data.frame(
    feature = c("Gene1", "Gene1"),
    dist = c(0, 1),
    expr = c(1, 2),
    fitted = c(1.1, 1.9)
  ), reference = "trajectory")
  expect_named(screening[1:6], c("variable", "distance", "value", "estimate", "reference", "mode"))
  expect_equal(unique(screening$reference), "trajectory")
})

test_that("numeric annotation thresholds are converted for SPATA2", {
  expect_equal(sgf_format_annotation_threshold(0), ">0")
  expect_equal(sgf_format_annotation_threshold("0"), ">0")
  expect_equal(sgf_format_annotation_threshold(" > 0 "), ">0")
  expect_equal(sgf_format_annotation_threshold("<= 1"), "<=1")
  expect_equal(sgf_format_annotation_threshold("kmeans_high"), "kmeans_high")
  expect_error(sgf_format_annotation_threshold(NA_real_), "annotation.threshold")
})

test_that("SPATA2 trajectory preparation does not forward unsupported verbose", {
  testthat::local_mocked_bindings(
    sgf_spata_fun = function(fun, required = TRUE) {
      force(fun)
      force(required)
      function(...) list(...)
    }
  )

  args <- sgf_prepare_trajectory(
    object = list(id = "mock"),
    trajectory_id = "traj",
    start = c(1, 2),
    end = c(3, 4),
    traj_df = NULL,
    width = NULL,
    verbose = TRUE
  )

  expect_equal(args$id, "traj")
  expect_true(isTRUE(args$overwrite))
  expect_false("verbose" %in% names(args))
})

test_that("cpp backend coordinates prefer metadata coord.cols", {
  srt <- make_spatial_gradient_seurat()
  coords <- sgf_cpp_coords(srt, image = NULL, coord.cols = c("x", "y"))

  expect_equal(rownames(coords), colnames(srt))
  expect_equal(unname(coords$x), unname(srt$x))
  expect_equal(unname(coords$y), unname(srt$y))
})

test_that("top variable table merges significance and best model fits deterministically", {
  significance <- data.frame(
    variable = c("Gene3", "Gene1", "Gene2"),
    p_value = c(0.03, 0.001, 0.01),
    fdr = c(0.06, 0.002, 0.02),
    stringsAsFactors = FALSE
  )
  model_fits <- data.frame(
    variable = c("Gene1", "Gene1", "Gene2"),
    model = c("ascending", "descending", "descending"),
    mae = c(0.1, 0.4, 0.2),
    rmse = c(0.2, 0.5, 0.3),
    stringsAsFactors = FALSE
  )

  top <- sgf_top_variables(
    screening_out = list(),
    significance = significance,
    model_fits = model_fits,
    nfeatures = 2,
    sign_var = "fdr",
    sign_threshold = 0.05
  )

  expect_equal(top$variable, c("Gene1", "Gene2"))
  expect_equal(top$rank, c(1L, 2L))
  expect_equal(top$best_model, c("ascending", "descending"))
  expect_equal(top$rmse, c(0.2, 0.3))
})

test_that("spatial gradient storage keeps only plain data frames", {
  srt <- make_spatial_gradient_seurat()
  result <- make_spatial_gradient_result()

  srt <- sgf_store_result(
    srt = srt,
    result_name = "mock",
    result = result,
    assay = "RNA",
    set_variable_features = FALSE
  )

  stored <- srt@tools[["SpatialGradientFeatures"]][["mock"]]
  expect_equal(names(stored), c("screening", "significance", "model_fits", "top_variables", "parameters"))
  expect_true(all(vapply(stored, is.data.frame, logical(1))))
  expect_false(any(vapply(stored, methods::is, logical(1), class2 = "SPATA2")))
  expect_equal(srt@misc[["SpatialGradientFeatures"]], c("Gene1", "Gene2"))
})

test_that("cpp backend stores annotation gradient result tables", {
  srt <- make_spatial_gradient_seurat()
  srt$source_score <- c(1, 0, 0, 0)

  srt <- RunSpatialGradientFeatures(
    srt,
    reference = "annotation",
    backend = "cpp",
    result_name = "cpp_annotation",
    variables = c("Gene1", "Gene2"),
    annotation.variable = "source_score",
    annotation.threshold = 0,
    layer = "counts",
    coord.cols = c("x", "y"),
    n_random = 0,
    n_bins = 3,
    min_spots = 1,
    sign_threshold = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  stored <- srt@tools[["SpatialGradientFeatures"]][["cpp_annotation"]]
  expect_equal(names(stored), c("screening", "significance", "model_fits", "top_variables", "parameters"))
  expect_true(all(vapply(stored, is.data.frame, logical(1))))
  expect_false(any(vapply(stored, methods::is, logical(1), class2 = "SPATA2")))
  expect_true(nrow(stored$screening) > 0)
  expect_true(nrow(stored$top_variables) > 0)
  expect_equal(stored$parameters$value[match("backend", stored$parameters$key)], "cpp")
  expect_true(all(c("norm_var", "rel_var", "linear_r2") %in% colnames(stored$significance)))
  expect_true(all(c("ascending", "descending", "peak", "valley", "linear") %in% stored$model_fits$model))
  expect_true(any(stored$significance$tot_var != stored$significance$linear_r2, na.rm = TRUE))
})

test_that("cpp backend stores trajectory gradient result tables", {
  srt <- make_spatial_gradient_seurat()

  srt <- RunSpatialGradientFeatures(
    srt,
    reference = "trajectory",
    backend = "cpp",
    result_name = "cpp_trajectory",
    variables = c("Gene1", "Gene2"),
    start = c(1, 1),
    end = c(2, 2),
    layer = "counts",
    coord.cols = c("x", "y"),
    n_random = 0,
    n_bins = 3,
    min_spots = 1,
    sign_threshold = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  stored <- srt@tools[["SpatialGradientFeatures"]][["cpp_trajectory"]]
  expect_true(all(vapply(stored, is.data.frame, logical(1))))
  expect_true(nrow(stored$significance) > 0)
  expect_true(nrow(stored$model_fits) > 0)
  expect_equal(unique(stored$screening$mode), "trajectory")
})

test_that("SpatialGradientPlot handles stored results and missing tables clearly", {
  srt <- make_spatial_gradient_seurat()
  result <- make_spatial_gradient_result()
  srt <- sgf_store_result(srt, "mock", result, assay = "RNA", set_variable_features = FALSE)

  expect_s3_class(
    SpatialGradientPlot(srt, result_name = "mock", plot_type = "line", theme_use = NULL),
    "ggplot"
  )
  p_line_lm <- SpatialGradientPlot(
    srt,
    result_name = "mock",
    plot_type = "line",
    line_fit = "lm",
    theme_use = NULL
  )
  expect_s3_class(p_line_lm, "ggplot")
  expect_s3_class(p_line_lm$layers[[2]]$geom, "GeomSmooth")
  expect_s3_class(
    SpatialGradientPlot(srt, result_name = "mock", plot_type = "summary", theme_use = NULL),
    "ggplot"
  )
  expect_s3_class(
    SpatialGradientPlot(srt, result_name = "mock", plot_type = "model", theme_use = NULL),
    "ggplot"
  )
  expect_s3_class(
    SpatialGradientPlot(
      srt,
      result_name = "mock",
      plot_type = "surface",
      features = "Gene1",
      layer = "counts",
      coord.cols = c("x", "y"),
      overlay_image = FALSE,
      theme_use = NULL
    ),
    "ggplot"
  )

  no_screening <- result
  no_screening$screening <- data.frame()
  srt <- sgf_store_result(srt, "no_screening", no_screening, assay = "RNA", set_variable_features = FALSE)
  expect_error(
    SpatialGradientPlot(srt, result_name = "no_screening", plot_type = "line", theme_use = NULL),
    "screening data"
  )

  empty_top <- result
  empty_top$top_variables <- empty_top$top_variables[0, , drop = FALSE]
  srt <- sgf_store_result(srt, "empty_top", empty_top, assay = "RNA", set_variable_features = FALSE)
  expect_error(
    SpatialGradientPlot(srt, result_name = "empty_top", plot_type = "summary", theme_use = NULL),
    "No top spatial gradient variables"
  )
})

test_that("SpatialGradientPlot combined returns a patchwork object when patchwork is available", {
  testthat::skip_if_not_installed("patchwork")
  srt <- make_spatial_gradient_seurat()
  result <- make_spatial_gradient_result()
  srt <- sgf_store_result(srt, "mock", result, assay = "RNA", set_variable_features = FALSE)

  p <- SpatialGradientPlot(
    srt,
    result_name = "mock",
    plot_type = "combined",
    features = "Gene1",
    layer = "counts",
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = NULL
  )
  expect_s3_class(p, "patchwork")
})

test_that("SpatialGradientPlot reuses stored assay layer for surface plots", {
  testthat::skip_if_not_installed("patchwork")
  srt <- make_spatial_gradient_seurat()
  result <- make_spatial_gradient_result()
  result$parameters <- data.frame(
    key = c("assay", "layer"),
    value = c("RNA", "counts"),
    stringsAsFactors = FALSE
  )
  srt <- sgf_store_result(srt, "mock_counts", result, assay = "RNA", set_variable_features = FALSE)

  p <- expect_warning(
    SpatialGradientPlot(
      srt,
      result_name = "mock_counts",
      plot_type = "combined",
      features = "Gene1",
      coord.cols = c("x", "y"),
      overlay_image = FALSE,
      theme_use = NULL
    ),
    regexp = NA
  )
  expect_s3_class(p, "patchwork")
})

test_that("RunSpatialGradientFeatures has a clear optional SPATA2 dependency error", {
  testthat::skip_if(requireNamespace("SPATA2", quietly = TRUE))
  srt <- make_spatial_gradient_seurat()

  expect_error(
    RunSpatialGradientFeatures(
      srt,
      reference = "trajectory",
      backend = "spata2",
      variables = "Gene1",
      start = c(1, 1),
      end = c(2, 2),
      layer = "counts",
      n_random = 1,
      verbose = FALSE
    ),
    "install SPATA2"
  )
})

test_that("SPATA2 integration stores only stable result tables", {
  testthat::skip_if_not_installed("SPATA2")
  testthat::skip_if_not(
    identical(Sys.getenv("SCOP_RUN_SPATA2_INTEGRATION"), "true"),
    "Set SCOP_RUN_SPATA2_INTEGRATION=true to run the slow SPATA2 integration fixture"
  )

  data(visium_human_pancreas_sub, package = "scop")
  srt <- subset(
    visium_human_pancreas_sub,
    cells = colnames(visium_human_pancreas_sub)[1:120]
  )
  coords <- srt@meta.data[, c("x", "y")]
  start <- as.numeric(coords[which.min(coords[["x"]]), c("x", "y")])
  end <- as.numeric(coords[which.max(coords[["x"]]), c("x", "y")])
  srt <- RunSpatialGradientFeatures(
    srt,
    reference = "trajectory",
    result_name = "tiny_trajectory",
    variables = rownames(srt)[1],
    start = start,
    end = end,
    assay = "Spatial",
    layer = "counts",
    platform = "VisiumSmall",
    n_random = 5,
    nfeatures = 1,
    verbose = FALSE
  )

  stored <- srt@tools[["SpatialGradientFeatures"]][["tiny_trajectory"]]
  expect_equal(names(stored), c("screening", "significance", "model_fits", "top_variables", "parameters"))
  expect_true(all(vapply(stored, is.data.frame, logical(1))))
  expect_false(any(vapply(stored, methods::is, logical(1), class2 = "SPATA2")))
})
