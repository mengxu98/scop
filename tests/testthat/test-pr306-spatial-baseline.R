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

  p <- SpatialNetworkPlot(srt = srt)
  expect_s3_class(p, "ggplot")
  expect_null(p$labels$subtitle)
  expect_equal(p$layers[[2]]$aes_params$size, 6)
  expect_s3_class(ggplot2::calc_element("axis.text.x", p$theme), "element_blank")
  expect_s3_class(SpatialNetworkPlot(res = res), "ggplot")
  expect_s3_class(SpatialNetworkPlot(srt = srt, res = res), "ggplot")
})

test_that("PR 306 new spatial entry points require explicit multi-image selection", {
  srt <- make_pr306_multi_image_object()

  expect_error(RunSpatialNetwork(srt, k = 1, verbose = FALSE), "Multiple spatial images")
  expect_error(srt_to_giotto(srt), "Multiple spatial images")
  expect_error(SpatialCellPlot(srt = srt), "Multiple spatial images")
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

test_that("spatial boundary feature values use the requested assay layer", {
  counts <- matrix(
    c(1, 2, 3, 4, 4, 3, 2, 1),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("cell", 1:4))
  )
  srt <- SeuratObject::CreateSeuratObject(counts)
  srt[["RNA"]] <- SeuratObject::CreateAssay5Object(
    counts = counts,
    data = counts * 10
  )
  boundaries <- data.frame(
    cell_id = rep(paste0("cell", 1:4), each = 4),
    polygon_id = rep(paste0("p", 1:4), each = 4),
    ring_id = 1,
    vertex_order = rep(1:4, 4),
    x = rep(c(0, 1, 1, 0), 4) + rep(0:3, each = 4),
    y = rep(c(0, 0, 1, 1), 4)
  )
  p <- SpatialCellPlot(
    srt = srt,
    boundaries = boundaries,
    features = "Gene1",
    assay = "RNA",
    layer = "data"
  )
  plotted <- p$layers[[1L]]$data
  expect_equal(
    sort(unique(plotted$.value)),
    sort(as.numeric(c(10, 20, 30, 40)))
  )
  expect_error(
    SpatialCellPlot(
      srt = srt,
      boundaries = boundaries,
      features = "Gene1",
      assay = "RNA",
      layer = "scale.data"
    ),
    "not available"
  )
  bad_boundaries <- boundaries
  bad_boundaries$cell_id[bad_boundaries$cell_id == "cell1"] <- "missing-cell"
  expect_error(
    SpatialCellPlot(
      srt = srt,
      boundaries = bad_boundaries,
      features = "Gene1",
      assay = "RNA",
      layer = "data"
    ),
    "Boundary cell IDs"
  )
  p2 <- SpatialCellPlot(
    srt = srt,
    boundaries = boundaries,
    features = c("Gene1", "Gene2"),
    assay = "RNA",
    layer = "counts"
  )
  expect_s3_class(p2, "patchwork")
})
