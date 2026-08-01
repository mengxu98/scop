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
  skip_if_not_installed("BiocNeighbors")
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
  expect_equal(length(srt@tools[["SpatialVariableFeatures"]]$summary$top_features), 3)
  expect_named(
    srt@tools[["SpatialVariableFeatures"]]$summary,
    c("n_features", "top_features", "top_feature_summary")
  )
})

test_that("C++ spatial scores agree with the R reference", {
  score_matrix <- getFromNamespace("spatial_variable_score_matrix", "scop")
  run_knn <- getFromNamespace("spatial_variable_run_knn", "scop")
  expr <- matrix(
    c(
      1, 0, 2, 0, 4, 1,
      3, NA, 1, 2, 0, 4,
      2, 2, 2, 2, 2, 2,
      0, 1, 0, 3, 0, 5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:6))
  )
  edges <- data.frame(
    from = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 3L),
    to = c(2L, 3L, 4L, 5L, 6L, 1L, 4L, 6L)
  )

  for (method in c("moran", "geary")) {
    reference <- score_matrix(expr, edges, method = method)
    dense <- run_knn(expr, edges, method = method, backend = "cpp")
    sparse <- run_knn(
      methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix"),
      edges,
      method = method,
      backend = "cpp"
    )

    expect_equal(dense$statistic, reference$statistic, tolerance = 1e-12)
    expect_equal(dense$score, reference$score, tolerance = 1e-12)
    expect_equal(sparse$statistic, reference$statistic, tolerance = 1e-12)
    expect_equal(sparse$score, reference$score, tolerance = 1e-12)
    expect_true(all(is.na(dense$p_value)))
  }
})

test_that("native backend selection is explicit and reproducible", {
  skip_if_not_installed("BiocNeighbors")
  args <- list(
    srt = make_spatial_variable_seurat(),
    method = "geary",
    layer = "counts",
    coord.cols = c("x", "y"),
    min_spots = 1,
    nfeatures = 3,
    verbose = FALSE
  )
  cpp <- do.call(RunSpatialVariableFeatures, c(args, backend = "cpp"))
  reference <- do.call(RunSpatialVariableFeatures, c(args, backend = "r"))
  cpp_result <- cpp@tools$SpatialVariableFeatures$result
  reference_result <- reference@tools$SpatialVariableFeatures$result
  cpp_result <- cpp_result[order(cpp_result$feature), ]
  reference_result <- reference_result[order(reference_result$feature), ]

  expect_equal(cpp_result$statistic, reference_result$statistic, tolerance = 1e-12)
  expect_identical(cpp@tools$SpatialVariableFeatures$parameters$backend, "cpp")
  expect_identical(reference@tools$SpatialVariableFeatures$parameters$backend, "r")
  permutation_args <- c(args, nperm = 20, seed = 20260730)
  cpp_permutation <- do.call(
    RunSpatialVariableFeatures,
    c(permutation_args, backend = "cpp")
  )
  r_permutation <- do.call(
    RunSpatialVariableFeatures,
    c(permutation_args, backend = "r")
  )
  cpp_p <- cpp_permutation@tools$SpatialVariableFeatures$result
  r_p <- r_permutation@tools$SpatialVariableFeatures$result
  cpp_p <- cpp_p[order(cpp_p$feature), ]
  r_p <- r_p[order(r_p$feature), ]
  expect_equal(
    cpp_p$p_value,
    r_p$p_value,
    tolerance = 0
  )
})

test_that("large sparse permutation boundaries agree exactly", {
  run_knn <- getFromNamespace("spatial_variable_run_knn", "scop")
  set.seed(20260731)
  expr <- Matrix::rsparsematrix(12, 5000, density = 0.01)
  expr@x <- abs(round(expr@x * 5))
  rownames(expr) <- paste0("Gene", seq_len(nrow(expr)))
  edges <- data.frame(
    from = rep(seq_len(ncol(expr)), each = 6L),
    to = as.integer((rep(seq_len(ncol(expr)), each = 6L) +
      rep(seq_len(6L), times = ncol(expr)) - 1L) %% ncol(expr) + 1L)
  )

  for (method in c("moran", "geary")) {
    set.seed(91)
    reference <- run_knn(
      expr, edges, method = method, nperm = 10L, backend = "r"
    )
    set.seed(91)
    candidate <- run_knn(
      expr, edges, method = method, nperm = 10L, backend = "cpp"
    )
    expect_identical(candidate$statistic, reference$statistic)
    expect_identical(candidate$p_value, reference$p_value)
  }
})

test_that("SPARKX backend output is normalized without storing backend objects", {
  srt <- make_spatial_variable_seurat()
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(pkg, fun) {
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
  expect_equal(
    out@tools[["SpatialVariableFeatures"]]$summary$top_features,
    c("Gene2", "Gene3")
  )
  expect_false(any(vapply(out@tools[["SpatialVariableFeatures"]], methods::is, logical(1), class2 = "SPARK")))
})

test_that("nnSVG backend output is normalized through lightweight helpers", {
  srt <- make_spatial_variable_seurat()
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    spatial_variable_make_spe = function(expr, coords, assay = NULL) {
      list(expr = expr, coords = coords, assay = assay)
    },
    get_namespace_fun = function(pkg, fun) {
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
  expect_equal(
    out@tools[["SpatialVariableFeatures"]]$summary$top_features,
    c("Gene2", "Gene4")
  )
})

test_that("SpatialVariableFeaturePlot uses stored result and SCOP spatial plotting", {
  skip_if_not_installed("BiocNeighbors")
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
  skip_if_not_installed("BiocNeighbors")
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
