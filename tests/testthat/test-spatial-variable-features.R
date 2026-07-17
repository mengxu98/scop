make_spatial_variable_seurat <- function() {
  counts <- matrix(
    c(
      5, 4, 0, 0, 0,
      0, 0, 4, 5, 4,
      1, 1, 1, 1, 1,
      3, 0, 3, 0, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:5))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(1, 2, 3, 4, 5)
  srt$y <- c(1, 1, 2, 2, 3)
  srt
}

test_that("native spatial variable feature results keep normalized columns", {
  srt <- RunSpatialVariableFeatures(
    make_spatial_variable_seurat(),
    method = "moran",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 3,
    verbose = FALSE
  )

  result <- srt@tools[["SpatialVariableFeatures"]][["result"]]
  expect_identical(
    srt@tools$SpatialVariableFeatures$source$coordinate_space,
    "raw"
  )
  expect_true(all(c(
    "feature", "rank", "method", "statistic", "score",
    "p_value", "q_value", "mean", "variance", "n_spots"
  ) %in% colnames(result)))
  expect_equal(unique(result$method), "moran")
  expect_equal(result$rank, seq_len(nrow(result)))
  expect_equal(length(srt@misc[["SpatialVariableFeatures"]]), 3)
  expect_named(srt@tools[["SpatialVariableFeatures"]]$summary, c("n_features", "top_features"))
})

test_that("SPARKX backend output is normalized without storing backend objects", {
  srt <- make_spatial_variable_seurat()
  testthat::local_mocked_bindings(
    spatial_variable_get_fun = function(pkg, fun) {
      expect_identical(pkg, "SPARK")
      expect_identical(fun, "sparkx")
      function(count_in, locus_in, ...) {
        expect_equal(rownames(count_in), paste0("Gene", 1:4))
        expect_equal(nrow(locus_in), ncol(count_in))
        list(
          res_mtest = data.frame(
            combinedPval = c(0.2, 0.01, 0.05, 0.3),
            adjustedPval = c(0.2, 0.04, 0.08, 0.3),
            row.names = rownames(count_in)
          )
        )
      }
    }
  )

  out <- RunSpatialVariableFeatures(
    srt,
    method = "SPARKX",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  result <- out@tools[["SpatialVariableFeatures"]][["result"]]
  expect_equal(unique(result$method), "SPARKX")
  expect_equal(result$feature[[1]], "Gene2")
  expect_equal(out@misc[["SpatialVariableFeatures"]], c("Gene2", "Gene3"))
  expect_false(any(vapply(out@tools[["SpatialVariableFeatures"]], methods::is, logical(1), class2 = "SPARK")))
})

test_that("nnSVG backend output is normalized through lightweight helpers", {
  srt <- make_spatial_variable_seurat()
  testthat::local_mocked_bindings(
    spatial_variable_require_package = function(pkg) invisible(TRUE),
    spatial_variable_make_spe = function(expr, coords, assay = NULL) {
      list(expr = expr, coords = coords, assay = assay)
    },
    spatial_variable_get_fun = function(pkg, fun) {
      expect_identical(pkg, "nnSVG")
      expect_identical(fun, "nnSVG")
      function(spe, assay_name = "counts", ...) {
        expect_identical(assay_name, "counts")
        list(
          row_data = data.frame(
            LR_stat = c(2, 9, 1, 4),
            pval = c(0.2, 0.001, 0.5, 0.03),
            padj = c(0.2, 0.004, 0.5, 0.06),
            row.names = rownames(spe$expr)
          )
        )
      }
    },
    spatial_variable_row_data = function(x) x$row_data
  )

  out <- RunSpatialVariableFeatures(
    srt,
    method = "nnSVG",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  result <- out@tools[["SpatialVariableFeatures"]][["result"]]
  expect_equal(unique(result$method), "nnSVG")
  expect_equal(result$feature[[1]], "Gene2")
  expect_equal(out@misc[["SpatialVariableFeatures"]], c("Gene2", "Gene4"))
})

test_that("SpatialVariableFeaturePlot uses stored result and SCOP spatial plotting", {
  srt <- RunSpatialVariableFeatures(
    make_spatial_variable_seurat(),
    method = "moran",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 3,
    verbose = FALSE
  )

  p_summary <- SpatialVariableFeaturePlot(
    srt,
    plot_type = "summary",
    nfeatures = 2,
    theme_use = NULL
  )
  expect_s3_class(p_summary, "ggplot")
  expect_null(p_summary$labels$size)
  expect_false("size" %in% names(p_summary$layers[[2]]$mapping))

  p_surface <- SpatialVariableFeaturePlot(
    srt,
    plot_type = "surface",
    nfeatures = 2,
    layer = "counts",
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = NULL
  )
  expect_true(inherits(p_surface, "ggplot") || inherits(p_surface, "patchwork"))
})

test_that("SpatialVariableFeaturePlot combined returns patchwork when available", {
  testthat::skip_if_not_installed("patchwork")
  srt <- RunSpatialVariableFeatures(
    make_spatial_variable_seurat(),
    method = "moran",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  p <- SpatialVariableFeaturePlot(
    srt,
    plot_type = "combined",
    nfeatures = 2,
    layer = "counts",
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = NULL
  )
  expect_s3_class(p, "patchwork")
})
