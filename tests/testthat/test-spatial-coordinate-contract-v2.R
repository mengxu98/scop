make_coordinate_contract_object <- function(image_class = c("VisiumV1", "VisiumV2")) {
  image_class <- match.arg(image_class)
  counts <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  raw <- data.frame(
    x = c(5, 10, 15, 20),
    y = c(4, 8, 12, 16),
    row.names = colnames(srt)
  )
  scale_factors <- structure(
    list(spot = 1, fiducial = 1, hires = 0.5, lowres = 0.25),
    class = "scalefactors"
  )
  raster <- array(1, dim = c(20, 30, 3))
  if (identical(image_class, "VisiumV1")) {
    image <- methods::new(
      "VisiumV1",
      image = raster,
      scale.factors = scale_factors,
      coordinates = data.frame(
        imagerow = raw$y,
        imagecol = raw$x,
        row.names = rownames(raw)
      ),
      spot.radius = 0.1,
      assay = "RNA",
      key = "v1_"
    )
  } else {
    fov <- SeuratObject::CreateFOV(
      data.frame(
        x = raw$y,
        y = raw$x,
        row.names = rownames(raw)
      ),
      type = "centroids",
      assay = "RNA",
      key = "v2_"
    )
    image <- methods::new(
      "VisiumV2",
      image = raster,
      scale.factors = scale_factors,
      molecules = fov@molecules,
      boundaries = fov@boundaries,
      assay = fov@assay,
      key = fov@key
    )
  }
  srt[["slice"]] <- image
  srt$group <- factor(c("A", "A", "B", "B"))
  list(object = srt, raw = raw)
}

test_that("VisiumV1 and VisiumV2 share full-resolution raw coordinates", {
  for (image_class in c("VisiumV1", "VisiumV2")) {
    fixture <- make_coordinate_contract_object(image_class)
    low_raw <- SpatialCoordinates(
      fixture$object,
      image = "slice",
      space = "raw",
      image.scale = "lowres"
    )
    high_raw <- SpatialCoordinates(
      fixture$object,
      image = "slice",
      space = "raw",
      image.scale = "hires"
    )
    expect_equal(low_raw$data$x, fixture$raw$x, info = image_class)
    expect_equal(low_raw$data$y, fixture$raw$y, info = image_class)
    expect_equal(high_raw$data[, c("x", "y")], low_raw$data[, c("x", "y")])
    expect_identical(low_raw$source$coordinate_contract_version, 2L)
    expect_identical(low_raw$source$image_class[[1L]], image_class)
  }
})

test_that("display coordinates scale and flip exactly once", {
  fixture <- make_coordinate_contract_object("VisiumV1")
  low <- SpatialCoordinates(
    fixture$object,
    image = "slice",
    space = "display",
    image.scale = "lowres"
  )
  high <- SpatialCoordinates(
    fixture$object,
    image = "slice",
    space = "display",
    image.scale = "hires"
  )
  expect_equal(low$data$x, fixture$raw$x * 0.25)
  expect_equal(low$data$y, 20 - fixture$raw$y * 0.25)
  expect_equal(high$data$x, fixture$raw$x * 0.5)
  expect_equal(high$data$y, 20 - fixture$raw$y * 0.5)
  expect_identical(low$source$scale_name, "lowres")
  expect_identical(high$source$scale_name, "hires")
})

test_that("VisiumV2 x/y names follow Seurat row-column display semantics", {
  fixture <- make_coordinate_contract_object("VisiumV2")
  image <- fixture$object[["slice"]]
  seurat_coords <- SeuratObject::GetTissueCoordinates(
    image,
    scale = "lowres",
    full = FALSE
  )
  display <- SpatialCoordinates(
    fixture$object,
    image = "slice",
    space = "display",
    image.scale = "lowres"
  )$data
  expect_equal(display$x, seurat_coords$y)
  expect_equal(display$y, dim(image@image)[1L] - seurat_coords$x)
})

test_that("FOV and metadata coordinates use identity transforms", {
  counts <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt[["fov"]] <- SeuratObject::CreateFOV(
    data.frame(x = 1:4, y = 5:8, row.names = colnames(srt)),
    type = "centroids",
    assay = "RNA",
    key = "fov_"
  )
  fov <- SpatialCoordinates(srt, image = "fov", space = "display")
  expect_equal(fov$data$x, 1:4)
  expect_equal(fov$data$y, 5:8)
  expect_identical(fov$transform$scale_name, "identity")

  srt@images <- list()
  srt$col <- 1:4
  srt$row <- 5:8
  metadata <- SpatialCoordinates(srt, coord.cols = c("col", "row"), space = "display")
  expect_equal(metadata$data[, c("x", "y")], fov$data[, c("x", "y")])
  expect_false(metadata$transform$y_flip)
})

test_that("invalid selected scale fails without mutating the object", {
  skip_if_not_installed("BiocNeighbors")
  fixture <- make_coordinate_contract_object("VisiumV1")
  fixture$object[["slice"]]@scale.factors$lowres <- NA_real_
  before <- serialize(fixture$object, NULL)
  raw <- SpatialCoordinates(
    fixture$object,
    image = "slice",
    space = "raw",
    image.scale = "lowres"
  )
  expect_equal(raw$data$x, fixture$raw$x)
  expect_equal(raw$data$y, fixture$raw$y)
  expect_error(
    SpatialCoordinates(
      fixture$object,
      image = "slice",
      space = "display",
      image.scale = "lowres"
    ),
    "lowres.*missing or invalid"
  )
  network <- RunSpatialNetwork(
    fixture$object,
    image = "slice",
    k = 1,
    verbose = FALSE
  )
  expect_identical(
    network@tools$SpatialNetwork$graphs[[1L]]$source$coordinate_space,
    "raw"
  )
  expect_s3_class(
    SpatialNetworkPlot(res = network@tools$SpatialNetwork),
    "ggplot"
  )
  expect_identical(serialize(fixture$object, NULL), before)
})

test_that("point pie and network use the same display coordinates", {
  skip_if_not_installed("BiocNeighbors")
  skip_if_not_installed("scatterpie")
  fixture <- make_coordinate_contract_object("VisiumV1")
  srt <- fixture$object
  expected <- SpatialCoordinates(
    srt,
    image = "slice",
    space = "display",
    image.scale = "lowres"
  )$data
  point <- SpatialSpotPlot(
    srt,
    group.by = "group",
    image = "slice",
    image.scale = "lowres",
    crop = FALSE
  )
  values <- matrix(
    rep(c(0.7, 0.3), each = ncol(srt)),
    nrow = ncol(srt),
    dimnames = list(colnames(srt), c("A", "B"))
  )
  pie <- SpatialSpotPlot(
    srt,
    values = values,
    plot_type = "pie",
    image = "slice",
    image.scale = "lowres",
    crop = FALSE
  )
  network_object <- RunSpatialNetwork(srt, image = "slice", k = 1, verbose = FALSE)
  network <- SpatialNetworkPlot(
    network_object,
    image.scale = "lowres",
    group.by = "group"
  )
  node_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    network$layers
  )
  expect_length(node_layers, 1L)
  node_data <- node_layers[[1L]]$data
  node_data <- node_data[match(expected$cell_id, node_data$cell_id), , drop = FALSE]
  expect_equal(point$data[, c("x", "y")], expected[, c("x", "y")])
  expect_equal(pie$data[, c("x", "y")], expected[, c("x", "y")])
  expect_equal(
    unname(as.matrix(node_data[, c("x", "y")])),
    unname(as.matrix(expected[, c("x", "y")]))
  )
})

test_that("raw distances do not depend on display image scale", {
  fixture <- make_coordinate_contract_object("VisiumV2")
  low <- SpatialCoordinates(fixture$object, image = "slice", image.scale = "lowres")$data
  high <- SpatialCoordinates(fixture$object, image = "slice", image.scale = "hires")$data
  expect_equal(
    unname(as.vector(stats::dist(low[, c("x", "y")]))),
    unname(as.vector(stats::dist(high[, c("x", "y")])))
  )
})

test_that("crop changes limits but not transformed spot coordinates", {
  fixture <- make_coordinate_contract_object("VisiumV2")
  cropped <- SpatialSpotPlot(
    fixture$object,
    group.by = "group",
    image = "slice",
    crop = TRUE
  )
  full <- SpatialSpotPlot(
    fixture$object,
    group.by = "group",
    image = "slice",
    crop = FALSE
  )
  expect_equal(cropped$data[, c("x", "y")], full$data[, c("x", "y")])
  expect_false(is.null(cropped$coordinates$limits$x))
  expect_true(is.null(full$coordinates$limits$x))
})

test_that("subsetting and RenameCells preserve coordinate identity and order", {
  for (image_class in c("VisiumV1", "VisiumV2")) {
    fixture <- make_coordinate_contract_object(image_class)
    subsetted <- suppressWarnings(fixture$object[, c("Spot4", "Spot2")])
    expect_identical(
      SpatialCoordinates(subsetted, image = "slice")$data$cell_id,
      colnames(subsetted)
    )
    renamed <- SeuratObject::RenameCells(
      fixture$object,
      new.names = paste0("Renamed", seq_len(ncol(fixture$object)))
    )
    expect_identical(
      SpatialCoordinates(renamed, image = "slice")$data$cell_id,
      colnames(renamed)
    )
  }
})

test_that("public HE overlay wrappers expose image.scale", {
  wrappers <- c(
    "SpatialCoordinates", "SpatialSpotPlot", "SpatialNetworkPlot",
    "SpatialDeconvolutionPlot", "Cell2locationPlot", "STdeconvolvePlot",
    "SpatialEcoTyperSpatialPlot", "SpatialDMPlot", "SpatialIntegrationPlot",
    "SpatialGradientPlot", "SpatialNeighborhoodPlot",
    "SpatialVariableFeaturePlot"
  )
  expect_true(all(vapply(wrappers, function(name) {
    "image.scale" %in% names(formals(getExportedValue("scop", name)))
  }, logical(1))))
  expect_true(all(vapply(wrappers, function(name) {
    args <- setdiff(names(formals(getExportedValue("scop", name))), "...")
    identical(utils::tail(args, 1L), "image.scale")
  }, logical(1))))
})

test_that("analysis code cannot bypass the shared raw extractor", {
  files <- list.files(
    testthat::test_path("..", "..", "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  skip_if(length(files) == 0L, "package source R files are not available")
  allowed <- c("SpatialCore.R", "SpatialCellPlot.R", "data.R")
  bypass <- vapply(files, function(file) {
    any(grepl("GetTissueCoordinates\\s*\\(", readLines(file, warn = FALSE)))
  }, logical(1))
  expect_setequal(basename(files[bypass]), allowed)
})

test_that("old coordinate-dependent network results are rejected", {
  skip_if_not_installed("BiocNeighbors")
  fixture <- make_coordinate_contract_object("VisiumV1")
  out <- RunSpatialNetwork(fixture$object, image = "slice", k = 1, verbose = FALSE)
  graph_name <- out@tools$SpatialNetwork$active_graph
  out@tools$SpatialNetwork$graphs[[graph_name]]$source$coordinate_contract_version <- NULL
  out@tools$SpatialNetwork$graphs[[graph_name]]$parameters$coordinate_contract_version <- NULL
  expect_error(GetSpatialGraph(out), "rerun\\s+RunSpatialNetwork")
  expect_error(SpatialNetworkPlot(out), "rerun\\s+RunSpatialNetwork")
})
