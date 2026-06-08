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

test_that("Giotto cluster metadata is written without storing full Giotto object", {
  srt <- make_giotto_cluster_seurat()
  clusters <- factor(c("1", "1", "2", "2", "1", "2"))
  names(clusters) <- colnames(srt)

  srt <- giotto_add_cluster_metadata(
    srt = srt,
    clusters = clusters,
    cluster_colname = "Giotto_cluster"
  )
  expect_true("Giotto_cluster" %in% colnames(srt@meta.data))
  expect_equal(as.character(srt$Giotto_cluster), as.character(clusters))

  srt@tools[["GiottoCluster"]] <- list(
    clusters = data.frame(cluster = clusters, row.names = names(clusters)),
    giotto_metadata = data.frame(cell_ID = names(clusters), leiden_clus = clusters),
    parameters = list(method = "leiden"),
    features = rownames(srt),
    cells = colnames(srt)
  )
  expect_null(srt@tools[["GiottoCluster"]][["giotto"]])
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

  out <- RunGiottoCluster(
    srt,
    method = "leiden",
    features = rownames(srt),
    k = 2,
    dims = 1:2,
    preprocess_params = list(method = "exact"),
    store_giotto = FALSE,
    verbose = FALSE
  )
  expect_s4_class(out, "Seurat")
  expect_true("Giotto_cluster" %in% colnames(out@meta.data))
  expect_equal(length(out$Giotto_cluster), ncol(out))
  expect_null(out@tools[["GiottoCluster"]][["giotto"]])
  expect_true(getOption("giotto.has_conda"))
  expect_true(getOption("giotto.use_conda"))
  expect_true(getOption("giotto.update_param"))
  expect_false(getOption("giotto.no_python_warn"))
  expect_s3_class(
    SpatialSpotPlot(
      out,
      group.by = "Giotto_cluster",
      coord.cols = c("x", "y"),
      overlay_image = FALSE
    ),
    "ggplot"
  )
})

test_that("standard_spatial_scop dispatches to Giotto clustering", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_cluster_seurat()

  out <- standard_scop(
    srt,
    workflow = "spatial",
    assay = "RNA",
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE,
    spatial_cluster_method = "Giotto",
    giotto_cluster_params = list(
      features = rownames(srt),
      k = 2,
      dims = 1:2,
      preprocess_params = list(method = "exact"),
      store_giotto = FALSE
    ),
    linear_reduction_dims = 3,
    nonlinear_reduction_dims = 2,
    verbose = FALSE
  )
  expect_true("Giotto_cluster" %in% colnames(out@meta.data))
  expect_equal(out@tools[["standard_spatial_scop"]][["cluster_col"]], "Giotto_cluster")
})
