raw_default_spatial_methods <- c(
  "RunSpaNorm",
  "RunSpotSweeper",
  "RunRCTD",
  "RunCARD",
  "RunSpatialDWLS",
  "RunBANKSY",
  "RunCytoSPACE",
  "RunSmoothClust",
  "RunMERINGUE",
  "RunSpatialVariableFeatures",
  "RunSpatialGradientFeatures",
  "RunSpatialNeighborhood",
  "RunStatialKontextual",
  "RunSpatialIntegration",
  "RunMistyR"
)

test_that("distance-sensitive spatial producers default to raw coordinates", {
  defaults <- vapply(raw_default_spatial_methods, function(method) {
    fun <- getExportedValue("scop", method)
    choices <- eval(formals(fun)$coordinate_space, envir = environment(fun))
    expect_identical(
      choices,
      c("raw", "legacy_display"),
      info = method
    )
    choices[[1L]]
  }, character(1))

  expect_identical(unname(defaults), rep("raw", length(defaults)))
})

test_that("registry advertises the migrated raw coordinate contract", {
  registry <- scop:::spatial_method_registry()
  rows <- registry[match(raw_default_spatial_methods, registry$method), , drop = FALSE]

  expect_false(anyNA(rows$method))
  expect_identical(rows$method, raw_default_spatial_methods)
  expect_identical(
    rows$coordinate_space_current,
    rep("raw", length(raw_default_spatial_methods))
  )
  expect_identical(
    rows$coordinate_space_target,
    rep("raw", length(raw_default_spatial_methods))
  )
  expect_true(all(rows$coordinate_requirement == "distance_sensitive"))
})

test_that("raw is the real public default on an image-backed Visium object", {
  skip_if_not_installed("BiocNeighbors")
  data("visium_human_pancreas_sub", package = "scop")
  srt <- visium_human_pancreas_sub
  image <- SeuratObject::Images(srt)[[1L]]

  raw <- scop:::spatial_analysis_coords(
    srt,
    image = image,
    coordinate_space = "raw"
  )
  legacy <- scop:::spatial_analysis_coords(
    srt,
    image = image,
    coordinate_space = "legacy_display"
  )

  expect_identical(raw$source$coordinate_space, "raw")
  expect_identical(legacy$source$coordinate_space, "legacy_display")
  expect_identical(rownames(raw$data), rownames(legacy$data))
  expect_false(isTRUE(all.equal(
    raw$data[, c("x", "y")],
    legacy$data[, c("x", "y")]
  )))

  # A numeric radius has the units of the selected coordinate space. On this
  # real image, 50 raw acquisition units are smaller than the nearest spacing,
  # while 50 legacy display units connect many of the same first 50 spots.
  cells <- rownames(raw$data)[seq_len(50L)]
  raw_distance <- as.matrix(stats::dist(raw$data[cells, c("x", "y")]))
  legacy_distance <- as.matrix(stats::dist(legacy$data[cells, c("x", "y")]))
  expect_identical(sum(raw_distance > 0 & raw_distance <= 50), 0L)
  expect_gt(sum(legacy_distance > 0 & legacy_distance <= 50), 0L)

  cytospace_coords <- scop:::cytospace_get_spatial_coords(
    srt,
    spot_ids = cells,
    image = image
  )
  expect_identical(
    attr(cytospace_coords, "spatial_source")$coordinate_space,
    "raw"
  )

  subset <- suppressWarnings(srt[, cells])
  subset$coda_label <- factor(rep(c("A", "B"), length.out = ncol(subset)))
  out <- RunSpatialNeighborhood(
    subset,
    group.by = "coda_label",
    image = image,
    k = 2,
    verbose = FALSE
  )
  expect_identical(
    out@tools$SpatialNeighborhood$source$coordinate_space,
    "raw"
  )
})

test_that("legacy display coordinates remain an explicit compatibility path", {
  skip_if_not_installed("BiocNeighbors")
  data("visium_human_pancreas_sub", package = "scop")
  srt <- visium_human_pancreas_sub
  image <- SeuratObject::Images(srt)[[1L]]
  cells <- colnames(srt)[seq_len(50L)]
  subset <- suppressWarnings(srt[, cells])
  subset$coda_label <- factor(rep(c("A", "B"), length.out = ncol(subset)))

  out <- RunSpatialNeighborhood(
    subset,
    group.by = "coda_label",
    image = image,
    coordinate_space = "legacy_display",
    k = 2,
    verbose = FALSE
  )
  expect_identical(
    out@tools$SpatialNeighborhood$source$coordinate_space,
    "legacy_display"
  )
})
