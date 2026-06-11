make_cluster_tree_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0, 4, 1,
      0, 3, 0, 4, 1, 0,
      2, 1, 5, 0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt[["RNA_snn_res.1"]] <- factor(c("0", "1", "1", "2", "3", "4"))
  srt[["RNA_snn_res.0.2"]] <- factor(c("0", "0", "1", "1", "1", "2"))
  srt[["RNA_snn_res.0.8"]] <- factor(c("0", "0", "3", "1", "1", "2"))
  srt[["SCT_snn_res.0.2"]] <- factor(c("A", "A", "B", "B", "B", "C"))
  srt
}

test_that("ClusterTreePlot detects and sorts Seurat resolution columns", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(srt, prefix = "RNA_snn")

  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1"))
  expect_equal(info$resolution, c(0.2, 0.8, 1))
})

test_that("ClusterTreePlot auto-detection keeps one clustering prefix", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(srt)

  expect_equal(unique(info$prefix), "RNA_snn")
  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1"))
})

test_that("ClusterTreePlot auto-detection respects uppercase default SNN prefix", {
  srt <- make_cluster_tree_seurat()
  srt@meta.data[c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1")] <- NULL
  srt[["RNA_SNN_res.0.2"]] <- factor(c("0", "0", "1", "1", "1", "2"))
  srt[["RNA_SNN_res.0.8"]] <- factor(c("0", "0", "3", "1", "1", "2"))
  srt[["SCT_snn_res.0.8"]] <- factor(c("A", "A", "C", "B", "B", "C"))
  srt[["SCT_snn_res.1"]] <- factor(c("A", "D", "C", "B", "B", "C"))

  info <- clustertree_resolve_cluster_cols(srt)

  expect_equal(unique(info$prefix), "RNA_SNN")
  expect_equal(info$column, c("RNA_SNN_res.0.2", "RNA_SNN_res.0.8"))
})

test_that("ClusterTreePlot respects manual cluster column order", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(
    srt,
    cluster_cols = c("RNA_snn_res.1", "RNA_snn_res.0.2")
  )

  expect_equal(info$column, c("RNA_snn_res.1", "RNA_snn_res.0.2"))
})

test_that("ClusterTreePlot filters prefix and resolutions", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(
    srt,
    prefix = "RNA_snn",
    resolutions = c(0.2, 1)
  )

  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.1"))
})

test_that("ClusterTreePlot builds node and edge statistics", {
  srt <- make_cluster_tree_seurat()
  info <- clustertree_resolve_cluster_cols(srt, prefix = "RNA_snn")
  tree_data <- clustertree_build_data(srt@meta.data, info)

  node <- tree_data$nodes[tree_data$nodes$node_id == "RNA_snn_res.0.2::1", ]
  expect_equal(node$size, 3)

  edge <- tree_data$all_edges[
    tree_data$all_edges$from_node == "RNA_snn_res.0.2::1" &
      tree_data$all_edges$to_node == "RNA_snn_res.0.8::1",
  ]
  expect_equal(edge$count, 2)
  expect_equal(edge$in_prop, 1)
  expect_equal(edge$out_prop, 2 / 3)
})

test_that("ClusterTreePlot returns ggplot for cluster tree", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn")

  expect_s3_class(plot, "ggplot")
})

test_that("ClusterTreePlot returns ggplot for a single marker overlay", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn", features = "Gene1")

  expect_s3_class(plot, "ggplot")
})

test_that("ClusterTreePlot accepts numeric metadata feature overlays", {
  srt <- make_cluster_tree_seurat()
  srt$G2M_score <- seq_len(ncol(srt))

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn", features = "G2M_score")
  feature_values <- clustertree_feature_values(
    srt,
    features = "G2M_score",
    cluster_cols = c("RNA_snn_res.0.2", "RNA_snn_res.0.8")
  )

  expect_s3_class(plot, "ggplot")
  expect_named(feature_values, "G2M_score")
})

test_that("ClusterTreePlot combines or lists multiple marker overlays", {
  srt <- make_cluster_tree_seurat()

  plots <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = c("Gene1", "Gene2"),
    combine = FALSE
  )
  expect_type(plots, "list")
  expect_named(plots, c("Gene1", "Gene2"))

  testthat::skip_if_not_installed("patchwork")
  plot <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = c("Gene1", "Gene2"),
    combine = TRUE
  )
  expect_s3_class(plot, "patchwork")
})

test_that("ClusterTreePlot preserves named feature-list groups", {
  srt <- make_cluster_tree_seurat()

  plots <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = list(Beta = c("Gene1", "Gene2"), Alpha = "Gene1"),
    combine = FALSE
  )

  expect_type(plots, "list")
  expect_named(plots, c("Beta", "Alpha"))
})

test_that("ClusterTreePlot aligns feature values to metadata cells", {
  srt <- make_cluster_tree_seurat()
  values <- c(Cell3 = 10, Cell1 = 2, Cell5 = 4)

  feature_nodes <- clustertree_aggregate_feature_values(
    values = values,
    meta_data = srt@meta.data,
    cluster_cols = "RNA_snn_res.0.2"
  )

  expect_equal(
    feature_nodes$feature_value[feature_nodes$node_id == "RNA_snn_res.0.2::0"],
    2
  )
  expect_equal(
    feature_nodes$feature_value[feature_nodes$node_id == "RNA_snn_res.0.2::1"],
    7
  )
})

test_that("ClusterTreePlot errors without multi-resolution columns", {
  counts <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(paste0("Gene", 1:2), paste0("Cell", 1:3))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    ClusterTreePlot(srt),
    "No multi-resolution Seurat clustering columns"
  )
})
