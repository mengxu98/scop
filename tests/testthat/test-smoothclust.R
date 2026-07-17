make_smoothclust_seurat <- function() {
  counts <- matrix(
    c(
      5, 4, 0, 0, 0,
      0, 0, 4, 5, 4,
      1, 1, 1, 1, 1,
      3, 0, 3, 0, 3,
      0, 2, 1, 4, 2
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("Spot", 1:5))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(0, 1, 0, 1, 2)
  srt$y <- c(0, 0, 1, 1, 2)
  srt
}

with_mock_smoothclust <- function(code) {
  fake_smooth <- function(
    input,
    spatial_coords,
    method,
    bandwidth,
    k,
    truncate,
    n_threads,
    ...
  ) {
    expect_true(is.matrix(input))
    expect_equal(ncol(input), nrow(spatial_coords))
    expect_true(method %in% c("uniform", "kernel", "knn"))
    input
  }
  fake_smoothness <- function(spatial_coords, labels, k = 6, n_threads = 1) {
    expect_equal(length(labels), nrow(spatial_coords))
    list(
      n_discordant = rep(0L, length(labels)),
      mean_discordant = 0
    )
  }
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      expect_identical(packages, "lmweber/smoothclust")
      invisible(TRUE)
    },
    smoothclust_get_fun = function(fun) {
      switch(fun,
        smoothclust = fake_smooth,
        smoothness_metric = fake_smoothness,
        stop("unexpected smoothclust function")
      )
    }
  )
  force(code)
}

test_that("RunSmoothClust requires explicit n_clusters", {
  expect_error(
    RunSmoothClust(
      make_smoothclust_seurat(),
      layer = "counts",
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "n_clusters"
  )
})

test_that("RunSmoothClust stores clusters and normalized schema", {
  srt <- make_smoothclust_seurat()
  with_mock_smoothclust({
    out <- RunSmoothClust(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      features = c("Gene1", "Gene2", "Gene4", "missing"),
      min_spots = 1,
      smooth_method = "knn",
      k = 2,
      n_clusters = 2,
      n_pcs = 2,
      store_smoothed = TRUE,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("SmoothClust_cluster" %in% colnames(out@meta.data))
  expect_false(any(is.na(out$SmoothClust_cluster)))
  expect_true("SmoothClust" %in% names(out@tools))
  store <- out@tools$SmoothClust
  expect_true(all(c(
    "clusters", "coords", "features", "feature_selection",
    "pca", "kmeans", "smoothness", "parameters", "smoothed"
  ) %in% names(store)))
  expect_equal(store$features, c("Gene1", "Gene2", "Gene4"))
  expect_equal(nrow(store$clusters), ncol(out))
  expect_equal(colnames(store$smoothed), colnames(out))
  expect_equal(store$parameters$smooth_method, "knn")
  expect_identical(store$source$coordinate_space, "raw")
  expect_equal(store$parameters$n_clusters, 2)
})

test_that("RunSmoothClust rejects non-finite coordinates without mutating input", {
  srt <- make_smoothclust_seurat()
  srt$x[5] <- NA
  before_metadata <- srt@meta.data
  before_tools <- srt@tools
  before_assays <- srt@assays
  before_reductions <- srt@reductions
  before_graphs <- srt@graphs
  with_mock_smoothclust({
    expect_error(
      RunSmoothClust(
        srt,
        layer = "counts",
        coord.cols = c("x", "y"),
        nfeatures = 3,
        min_spots = 1,
        smooth_method = "knn",
        k = 2,
        n_clusters = 2,
        n_pcs = 2,
        verbose = FALSE
      ),
      "non-finite"
    )
  })

  expect_identical(srt@meta.data, before_metadata)
  expect_identical(srt@tools, before_tools)
  expect_identical(srt@assays, before_assays)
  expect_identical(srt@reductions, before_reductions)
  expect_identical(srt@graphs, before_graphs)
})

test_that("RunSmoothClust validates feature and coordinate filters clearly", {
  srt <- make_smoothclust_seurat()
  expect_error(
    RunSmoothClust(
      srt,
      layer = "counts",
      coord.cols = c("missing_x", "y"),
      n_clusters = 2,
      verbose = FALSE
    ),
    "Spatial coordinates"
  )
  expect_error(
    RunSmoothClust(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      features = "missing",
      n_clusters = 2,
      verbose = FALSE
    ),
    "features"
  )
})

test_that("SmoothClust clusters are directly plottable with SpatialSpotPlot", {
  srt <- make_smoothclust_seurat()
  with_mock_smoothclust({
    out <- RunSmoothClust(
      srt,
      layer = "counts",
      coord.cols = c("x", "y"),
      min_spots = 1,
      smooth_method = "knn",
      k = 2,
      n_clusters = 2,
      n_pcs = 2,
      verbose = FALSE
    )
  })

  p <- SpatialSpotPlot(
    out,
    group.by = "SmoothClust_cluster",
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = "theme_blank"
  )
  expect_s3_class(p, "ggplot")
})
