make_giotto_cluster_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0, 4, 1,
      0, 3, 0, 4, 1, 0,
      2, 1, 5, 0, 0, 3,
      0, 2, 0, 5, 1, 4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(1, 2, 3, 1, 2, 3)
  srt$y <- c(1, 1, 1, 2, 2, 2)
  srt
}

test_that("Giotto cluster extraction aligns metadata by Seurat cell names", {
  giotto_metadata <- data.frame(
    cell_ID = c("Spot3", "Spot1", "Spot2"),
    leiden_clus = c("C", "A", "B"),
    stringsAsFactors = FALSE
  )
  clusters <- giotto_extract_clusters(
    giotto_metadata = giotto_metadata,
    cells = c("Spot1", "Spot2", "Spot3"),
    cluster_name = "leiden_clus",
    method = "leiden"
  )
  expect_equal(names(clusters), c("Spot1", "Spot2", "Spot3"))
  expect_equal(as.character(clusters), c("A", "B", "C"))
})

test_that("Giotto result helper keeps the full Giotto object outside Seurat", {
  result <- giotto_result(
    result_type = "cluster",
    giotto = list(mock = TRUE),
    clusters = data.frame(cluster = c("1", "2"), row.names = c("Spot1", "Spot2"))
  )

  expect_s3_class(result, "giotto2_result")
  expect_s3_class(result, "giotto2_cluster")
  expect_equal(result$giotto, list(mock = TRUE))
  expect_s3_class(result$clusters, "data.frame")
})

test_that("Giotto parameter helpers reject unsafe pass-through arguments", {
  expect_error(
    giotto_call(function(x) x, list(x = 1, bad = 2)),
    "Unsupported Giotto argument"
  )
  expect_error(
    giotto_merge_args(
      defaults = list(return_gobject = TRUE),
      extra = list(return_gobject = FALSE),
      arg_name = "cluster_params",
      reserved = "return_gobject"
    ),
    "cannot override"
  )
})

test_that("Giotto input preparation filters non-finite coordinates", {
  srt <- make_giotto_cluster_seurat()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$x[1] <- NA
  input <- giotto_prepare_input(
    srt = srt,
    assay = "RNA",
    layer = "data",
    features = rownames(srt),
    coord.cols = c("x", "y")
  )

  expect_equal(input$cells, colnames(srt)[-1])
  expect_true(all(is.finite(input$spatial_locs$sdimx)))
  expect_true(all(is.finite(input$spatial_locs$sdimy)))
})

test_that("RunGiottoCluster errors before mutating non-Seurat input", {
  expect_error(
    RunGiottoCluster(matrix(1, nrow = 2, ncol = 2), verbose = FALSE),
    "Seurat"
  )
})

test_that("RunGiottoCluster runs with installed Giotto", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_cluster_seurat()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  old_options <- options(
    giotto.has_conda = TRUE,
    giotto.use_conda = TRUE,
    giotto.update_param = TRUE,
    giotto.no_python_warn = FALSE
  )
  on.exit(options(old_options), add = TRUE)

  before_meta <- srt@meta.data
  before_tools <- srt@tools
  before_misc <- srt@misc

  out <- RunGiottoCluster(
    srt,
    method = "leiden",
    features = rownames(srt),
    k = 20,
    dims = 1:2,
    preprocess_params = list(method = "exact", name = "custom_pca"),
    store_giotto = FALSE,
    verbose = FALSE
  )
  expect_s3_class(out, "giotto2_result")
  expect_s3_class(out, "giotto2_cluster")
  expect_false(inherits(out, "Seurat"))
  expect_false(is.null(out$giotto))
  expect_s3_class(out$clusters, "data.frame")
  expect_equal(rownames(out$clusters), colnames(srt))
  expect_equal(length(out$cluster_vector), ncol(srt))
  expect_equal(out$cells, colnames(srt))
  expect_equal(out$features, rownames(srt))
  expect_equal(out$parameters$k, ncol(srt) - 1L)
  expect_false("Giotto_cluster" %in% colnames(srt@meta.data))
  expect_equal(srt@meta.data, before_meta)
  expect_equal(srt@tools, before_tools)
  expect_equal(srt@misc, before_misc)
  expect_true(getOption("giotto.has_conda"))
  expect_true(getOption("giotto.use_conda"))
  expect_true(getOption("giotto.update_param"))
  expect_false(getOption("giotto.no_python_warn"))
})

test_that("standard_spatial_scop does not dispatch to standalone Giotto results", {
  srt <- make_giotto_cluster_seurat()

  expect_error(
    standard_scop(
      srt,
      workflow = "spatial",
      assay = "RNA",
      do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE,
      do_spatial_cluster = FALSE,
      spatial_cluster_method = "Giotto",
      linear_reduction_dims = 3,
      nonlinear_reduction_dims = 2,
      verbose = FALSE
    ),
    "BayesSpace"
  )
})
