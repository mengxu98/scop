make_pr306_spatial_object <- function() {
  counts <- matrix(
    c(3, 1, 0, 2, 0, 4, 1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  srt
}

make_pr306_multi_image_object <- function() {
  srt <- make_pr306_spatial_object()
  slice1 <- data.frame(x = c(0, 1), y = c(0, 1), row.names = c("spot1", "spot2"))
  slice2 <- data.frame(x = c(2, 3), y = c(2, 3), row.names = c("spot3", "spot4"))
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(slice1, type = "centroids", assay = assay, key = "s1_")
  srt[["slice2"]] <- SeuratObject::CreateFOV(slice2, type = "centroids", assay = assay, key = "s2_")
  srt
}

test_that("PR 306 spatial graph naming and overwrite contract is stable", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_pr306_spatial_object()

  srt <- RunSpatialNetwork(srt, k = 1, verbose = FALSE)
  expect_equal(srt@tools$SpatialNetwork$active_graph, "knn_k1")
  expect_named(srt@tools$SpatialNetwork$graphs, "knn_k1")
  expect_named(
    srt@tools$SpatialNetwork$graphs$knn_k1$edges,
    c("from", "to", "distance", "weight")
  )
  expect_true(all(srt@tools$SpatialNetwork$graphs$knn_k1$edges$weight == 1))

  srt <- RunSpatialNetwork(srt, k = 2, verbose = FALSE)
  expect_setequal(names(srt@tools$SpatialNetwork$graphs), c("knn_k1", "knn_k2"))
  expect_error(
    RunSpatialNetwork(srt, k = 1, verbose = FALSE),
    "already exists"
  )
  srt <- RunSpatialNetwork(srt, k = 1, overwrite = TRUE, verbose = FALSE)
  expect_equal(srt@tools$SpatialNetwork$active_graph, "knn_k1")

  srt <- RunSpatialNetwork(srt, method = "radius", radius = 1.5, verbose = FALSE)
  expect_equal(srt@tools$SpatialNetwork$active_graph, "radius_r1p5")
})

test_that("PR 306 spatial graph plotting supports all result input modes", {
  skip_if_not_installed("BiocNeighbors")
  srt <- RunSpatialNetwork(make_pr306_spatial_object(), k = 1, verbose = FALSE)
  res <- srt@tools$SpatialNetwork

  p <- SpatialNetworkPlot(object = srt)
  expect_s3_class(p, "ggplot")
  expect_null(p$labels$subtitle)
  expect_equal(p$layers[[2]]$aes_params$size, 6)
  expect_s3_class(ggplot2::calc_element("axis.text.x", p$theme), "element_blank")
  expect_s3_class(SpatialNetworkPlot(res = res), "ggplot")
  expect_s3_class(SpatialNetworkPlot(object = srt, res = res), "ggplot")
})

test_that("PR 306 new spatial entry points require explicit multi-image selection", {
  srt <- make_pr306_multi_image_object()

  expect_error(RunSpatialNetwork(srt, k = 1, verbose = FALSE), "Multiple spatial images")
  expect_error(srt_to_giotto(srt), "Multiple spatial images")
  expect_error(SpatialCellPlot(object = srt), "Multiple spatial images")
})

test_that("PR 306 cell plotting rejects spot centers as polygons", {
  centers <- data.frame(
    cell_id = c("spot1", "spot2"),
    x = c(0, 1),
    y = c(0, 1)
  )
  expect_error(
    SpatialCellPlot(boundaries = centers),
    "at least three distinct vertices"
  )
})

test_that("spatial boundary plots use the shared axis-free spatial theme", {
  boundaries <- data.frame(
    cell_id = rep(c("cell1", "cell2"), each = 4),
    polygon_id = rep(c("p1", "p2"), each = 4),
    ring_id = 1,
    vertex_order = rep(1:4, 2),
    x = c(0, 1, 1, 0, 1.2, 2.2, 2.2, 1.2),
    y = c(0, 0, 1, 1, 0, 0, 1, 1),
    cell_type = rep(c("A", "B"), each = 4)
  )
  p <- SpatialCellPlot(boundaries = boundaries, group.by = "cell_type")
  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::calc_element("axis.text.x", p$theme), "element_blank")
})
