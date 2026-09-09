make_spatial_srt_fixture <- function() {
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

test_that("native spatial APIs take srt as the Seurat argument", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_spatial_srt_fixture()

  qc <- suppressWarnings(RunSpotQC(
    srt = srt,
    assay = "RNA",
    qc_metrics = c("umi", "gene"),
    UMI_threshold = 0,
    gene_threshold = 0,
    verbose = FALSE
  ))
  expect_s4_class(qc, "Seurat")
  expect_s3_class(
    SpatialSpotPlot(
      srt = qc,
      group.by = "SpotQC",
      overlay_image = FALSE,
      theme_use = NULL
    ),
    "ggplot"
  )

  network <- RunSpatialNetwork(
    srt = qc,
    k = 1,
    verbose = FALSE
  )
  expect_s4_class(network, "Seurat")
  expect_s3_class(
    SpatialNetworkPlot(
      srt = network,
      graph.name = "knn_k1",
      theme_use = NULL
    ),
    "ggplot"
  )

  svf <- RunSpatialVariableFeatures(
    srt = network,
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
      srt = svf,
      plot_type = "summary",
      theme_use = NULL
    ),
    "ggplot"
  )
})

test_that("spatial APIs reject object= and invalid srt input", {
  srt <- make_spatial_srt_fixture()
  expect_error(
    RunSpotQC(object = srt, verbose = FALSE),
    "unused argument"
  )
  expect_error(
    SpatialSpotPlot(verbose = FALSE),
    "srt"
  )
  expect_error(
    RunSpatialNetwork(srt = list(), verbose = FALSE),
    "srt.*Seurat"
  )
})

test_that("non-spatial wrappers take srt and reject object=", {
  srt <- make_spatial_srt_fixture()
  expect_error(
    RunVECTOR(object = srt, verbose = FALSE),
    "unused argument"
  )
  expect_error(
    VECTORPlot(object = srt),
    "srt"
  )
  expect_error(
    RunDEtest(object = srt, verbose = FALSE),
    "srt"
  )
  expect_error(
    PrepareSCExplorer(object = srt, verbose = FALSE),
    "unused argument"
  )
})
