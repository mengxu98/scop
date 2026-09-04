make_spatial_object_alias_fixture <- function() {
  counts <- matrix(
    c(
      5, 4, 0, 1, 5, 4, 0, 1, 5,
      0, 1, 4, 5, 0, 1, 4, 5, 0,
      1, 1, 1, 1, 1, 1, 1, 1, 1,
      3, 0, 3, 0, 3, 0, 3, 0, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("Gene", 1:4),
      paste0("Spot", 1:9)
    )
  )
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  srt$x <- rep(1:3, 3)
  srt$y <- rep(1:3, each = 3)
  srt
}

test_that("native spatial APIs accept object as an input alias", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_spatial_object_alias_fixture()

  qc <- suppressWarnings(RunSpotQC(
    object = srt,
    assay = "RNA",
    qc_metrics = c("umi", "gene"),
    UMI_threshold = 0,
    gene_threshold = 0,
    verbose = FALSE
  ))
  expect_s4_class(qc, "Seurat")
  expect_s3_class(
    SpatialSpotPlot(
      object = qc,
      group.by = "SpotQC",
      overlay_image = FALSE,
      theme_use = NULL
    ),
    "ggplot"
  )

  network <- RunSpatialNetwork(
    object = qc,
    k = 1,
    verbose = FALSE
  )
  expect_s4_class(network, "Seurat")
  expect_s3_class(
    SpatialNetworkPlot(
      object = network,
      graph.name = "knn_k1",
      theme_use = NULL
    ),
    "ggplot"
  )

  svf <- RunSpatialVariableFeatures(
    object = network,
    assay = "RNA",
    layer = "counts",
    method = "moran",
    backend = "r",
    coord.cols = c("x", "y"),
    nfeatures = 2,
    min_spots = 1,
    verbose = FALSE
  )
  expect_s4_class(svf, "Seurat")
  expect_s3_class(
    SpatialVariableFeaturePlot(
      object = svf,
      plot_type = "summary",
      theme_use = NULL
    ),
    "ggplot"
  )
})

test_that("spatial object aliases reject ambiguous or missing input", {
  srt <- make_spatial_object_alias_fixture()
  expect_error(
    RunSpotQC(srt = srt, object = srt, verbose = FALSE),
    "only one of.*srt.*object"
  )
  expect_error(
    SpatialSpotPlot(verbose = FALSE),
    "through.*srt.*object"
  )
  expect_error(
    RunSpatialNetwork(object = list(), verbose = FALSE),
    "object.*Seurat"
  )
})
