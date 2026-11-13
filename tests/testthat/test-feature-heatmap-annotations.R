make_feature_annotation_srt <- function() {
  counts <- matrix(
    rep(seq_len(8L), each = 4L), nrow = 4L,
    dimnames = list(paste0("gene", 1:4), paste0("cell", 1:8))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$group <- rep(c("A", "B"), each = 4L)
  object$time <- rep(c(0, 1, 2, 3), 2)
  object
}

capture_feature_annotations <- function(object, ...) {
  annotations <- list()
  legends <- list()
  original_annotation <- ComplexHeatmap::anno_simple
  original_legend <- ComplexHeatmap::Legend
  testthat::local_mocked_bindings(
    anno_simple = function(...) {
      args <- list(...)
      annotations[[length(annotations) + 1L]] <<- args
      do.call(original_annotation, args)
    },
    Legend = function(...) {
      args <- list(...)
      legends[[as.character(args$title)]] <<- args
      do.call(original_legend, args)
    },
    .package = "ComplexHeatmap"
  )
  result <- FeatureHeatmap(
    object, features = rownames(object), cell_annotation = "time",
    cell_annotation_palette = "viridis", lib_normalize = FALSE,
    exp_method = "raw", use_raster = FALSE, verbose = FALSE, ...
  )
  list(result = result, annotations = annotations, legends = legends)
}

test_that("numeric cell annotations retain their type when also used for grouping", {
  skip_if_not_installed("ComplexHeatmap")
  object <- make_feature_annotation_srt()
  for (args in list(list(group.by = "time"), list(group.by = "group", split.by = "time"))) {
    captured <- suppressWarnings(do.call(capture_feature_annotations, c(list(object), args)))
    annotation <- captured$annotations[[1L]]
    expect_type(captured$result$cell_metadata$time, "double")
    expect_true(is.function(annotation$col))
    expect_equal(unname(annotation$x), unname(object$time[match(names(annotation$x), colnames(object))]))
    expect_false(anyNA(annotation$x))
    expect_true(is.function(captured$legends$time$col_fun))
  }
})

test_that("all-missing numeric cell annotations draw a missing-value legend", {
  skip_if_not_installed("ComplexHeatmap")
  object <- make_feature_annotation_srt()
  for (value in list(rep(NA_real_, 8L), rep(Inf, 8L), rep(NaN, 8L))) {
    object$time <- value
    captured <- capture_feature_annotations(object, group.by = "group")
    annotation <- captured$annotations[[1L]]
    expect_true(all(is.na(annotation$x)))
    expect_identical(annotation$na_col, "grey90")
    expect_identical(captured$legends$time$labels, "Missing")
    expect_null(captured$legends$time$col_fun)
  }
})

test_that("missingness is evaluated among selected cells", {
  skip_if_not_installed("ComplexHeatmap")
  object <- make_feature_annotation_srt()
  object$time[1:4] <- NA_real_
  captured <- capture_feature_annotations(object, group.by = "group", cells = colnames(object)[1:4])
  expect_identical(captured$legends$time$labels, "Missing")
  expect_length(captured$annotations[[1L]]$x, 4L)
})

test_that("partial missingness preserves finite values and continuous colors", {
  skip_if_not_installed("ComplexHeatmap")
  object <- make_feature_annotation_srt()
  object$time <- c(NA, 1, Inf, 3, 4, NaN, 6, 7)
  captured <- capture_feature_annotations(object, group.by = "group", cell_order = rev(colnames(object)))
  annotation <- captured$annotations[[1L]]
  expect_true(is.function(annotation$col))
  expect_equal(sum(is.na(annotation$x)), 3)
  expected <- object$time[match(names(annotation$x), colnames(object))]
  expected[!is.finite(expected)] <- NA_real_
  expect_equal(unname(annotation$x), unname(expected))
  expect_identical(annotation$na_col, "transparent")
})

test_that("factor annotations stay discrete", {
  skip_if_not_installed("ComplexHeatmap")
  object <- make_feature_annotation_srt()
  object$time <- factor(object$time)
  captured <- capture_feature_annotations(object, group.by = "group")
  expect_false(is.function(captured$annotations[[1L]]$col))
  expect_identical(captured$legends$time$labels, levels(object$time))
})
