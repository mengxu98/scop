make_meringue_seurat <- function() {
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

mock_meringue_get_fun <- function(fun) {
  switch(fun,
    getSpatialNeighbors = function(pos, filterDist = NA, binary = TRUE, verbose = FALSE) {
      out <- matrix(1, nrow = nrow(pos), ncol = nrow(pos))
      diag(out) <- 0
      rownames(out) <- rownames(pos)
      colnames(out) <- rownames(pos)
      out
    },
    moranTest = function(x, weight, alternative = "greater") {
      statistic <- mean(as.numeric(x))
      c(
        observed = statistic,
        expected = 0,
        sd = 1,
        p.value = 1 / (1 + statistic)
      )
    },
    moranPermutationTest = function(z, w, alternative = "greater", N = 10000, seed = 0, ncores = 1, plot = FALSE, ...) {
      statistic <- mean(as.numeric(z))
      c(
        observed = statistic,
        expected = 0,
        sd = 1,
        p.value = 1 / (1 + statistic),
        n = N
      )
    },
    spatialCrossCorMatrix = function(mat, weight) {
      out <- stats::cor(t(as.matrix(mat)))
      out[!is.finite(out)] <- 0
      out
    },
    spatialCrossCorTest = function(x, y, w, n = 1000, ncores = 1, plot = FALSE, ...) {
      0.5
    },
    groupSigSpatialPatterns = function(pos, mat, scc, plot = FALSE, verbose = TRUE, ...) {
      stats::setNames(rep(c("module1", "module2"), length.out = nrow(mat)), rownames(mat))
    },
    stop("Unexpected MERINGUE function: ", fun)
  )
}

test_that("RunMERINGUE stores normalized autocorrelation, cross-correlation, and module tables", {
  testthat::local_mocked_bindings(
    meringue_require_package = function(pkg) invisible(TRUE),
    meringue_get_fun = mock_meringue_get_fun,
    meringue_package_version = function(pkg) "mock"
  )

  out <- RunMERINGUE(
    make_meringue_seurat(),
    layer = "counts",
    coord.cols = c("x", "y"),
    mode = c("autocorrelation", "cross_correlation", "modules"),
    min_spots = 1,
    nfeatures = 3,
    pairwise_features = c("Gene1", "Gene2", "Gene4"),
    verbose = FALSE
  )

  stored <- out@tools[["MERINGUE"]]
  expect_true(all(c(
    "autocorrelation", "cross_correlation", "modules",
    "coords", "weight", "features", "pairwise_features", "parameters"
  ) %in% names(stored)))
  expect_true(all(c(
    "feature", "rank", "statistic", "expected", "sd",
    "p_value", "q_value", "score", "mean", "variance", "n_spots"
  ) %in% colnames(stored$autocorrelation)))
  expect_equal(stored$autocorrelation$rank, seq_len(nrow(stored$autocorrelation)))
  expect_equal(stored$autocorrelation$feature[[1]], "Gene2")
  expect_equal(out@misc[["MERINGUEFeatures"]], c("Gene2", "Gene1", "Gene4"))
  expect_true(all(c("feature1", "feature2", "rank", "correlation", "p_value") %in% colnames(stored$cross_correlation)))
  expect_true(all(c("feature", "module", "module_size", "rank") %in% colnames(stored$modules)))
  expect_true(all(vapply(stored[c("autocorrelation", "cross_correlation", "modules", "parameters")], is.data.frame, logical(1))))
  expect_false(any(vapply(stored, methods::is, logical(1), class2 = "MERINGUE")))
})

test_that("MERINGUE module output with groups is normalized", {
  out <- list(
    hc = list(labels = c("Gene1", "Gene2")),
    groups = stats::setNames(factor(c("0", "1")), c("Gene1", "Gene2")),
    prs = list()
  )

  modules <- meringue_normalize_modules(out, features = c("Gene1", "Gene2"))
  expect_equal(modules$feature, c("Gene1", "Gene2"))
  expect_equal(modules$module, c("0", "1"))
  expect_equal(modules$module_size, c(1L, 1L))
})

test_that("pairwise_features controls cross-correlation scope without implicit O(n^2) expansion", {
  testthat::local_mocked_bindings(
    meringue_require_package = function(pkg) invisible(TRUE),
    meringue_get_fun = mock_meringue_get_fun,
    meringue_package_version = function(pkg) "mock"
  )

  out <- RunMERINGUE(
    make_meringue_seurat(),
    layer = "counts",
    coord.cols = c("x", "y"),
    mode = "cross_correlation",
    min_spots = 1,
    pairwise_features = c("Gene1", "Gene2"),
    verbose = FALSE
  )

  stored <- out@tools[["MERINGUE"]]
  expect_equal(stored$pairwise_features, c("Gene1", "Gene2"))
  expect_equal(nrow(stored$cross_correlation), 1L)
  expect_equal(stored$cross_correlation$feature1, "Gene1")
  expect_equal(stored$cross_correlation$feature2, "Gene2")
  expect_equal(nrow(stored$autocorrelation), 0L)
})

test_that("RunMERINGUE can compute optional cross-correlation p values", {
  testthat::local_mocked_bindings(
    meringue_require_package = function(pkg) invisible(TRUE),
    meringue_get_fun = mock_meringue_get_fun,
    meringue_package_version = function(pkg) "mock"
  )

  out <- RunMERINGUE(
    make_meringue_seurat(),
    layer = "counts",
    coord.cols = c("x", "y"),
    mode = "cross_correlation",
    min_spots = 1,
    pairwise_features = c("Gene1", "Gene2"),
    cross_cor_params = list(test = TRUE, n = 5),
    verbose = FALSE
  )

  expect_equal(out@tools[["MERINGUE"]]$cross_correlation$p_value, 0.5)
})

test_that("RunMERINGUE top features are directly plottable with SpatialSpotPlot", {
  testthat::local_mocked_bindings(
    meringue_require_package = function(pkg) invisible(TRUE),
    meringue_get_fun = mock_meringue_get_fun,
    meringue_package_version = function(pkg) "mock"
  )

  out <- RunMERINGUE(
    make_meringue_seurat(),
    layer = "counts",
    coord.cols = c("x", "y"),
    mode = "autocorrelation",
    min_spots = 1,
    nfeatures = 2,
    verbose = FALSE
  )

  p <- SpatialSpotPlot(
    out,
    features = out@misc[["MERINGUEFeatures"]],
    layer = "counts",
    coord.cols = c("x", "y"),
    overlay_image = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("RunMERINGUE has a clear optional dependency error", {
  testthat::skip_if(requireNamespace("MERINGUE", quietly = TRUE))

  expect_error(
    RunMERINGUE(
      make_meringue_seurat(),
      layer = "counts",
      coord.cols = c("x", "y"),
      min_spots = 1,
      verbose = FALSE
    ),
    "RunMERINGUE"
  )
})
