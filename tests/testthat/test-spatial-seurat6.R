seurat6_spatial_object <- function(compact = TRUE) {
  skip_if_not("compact" %in% names(formals(getS3method("CreateSegmentation", "data.frame", envir = asNamespace("SeuratObject")))))
  ids <- c("cell.a_b", "cell_a.b", "cell3")
  counts <- Matrix::Matrix(matrix(seq_len(18), 6,
    dimnames = list(paste0("g", 1:6), ids)
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts, assay = "custom.Polygons")
  outer <- rbind(c(0, 0), c(4, 0), c(4, 6), c(0, 6), c(0, 0))
  hole <- rbind(c(1, 1), c(1, 2), c(2, 2), c(2, 1), c(1, 1))
  xy <- rbind(outer, hole, sweep(outer, 2, c(10, 20), "+"), sweep(outer, 2, c(30, 50), "+"))
  polygons <- data.frame(cell = rep(ids, c(10, 5, 5)), x = xy[, 1L], y = xy[, 2L])
  seg <- SeuratObject::CreateSegmentation(polygons, compact = compact)
  if (!compact) seg[["sf.data"]] <- polygons
  centroids <- SeuratObject::CreateCentroids(data.frame(x = c(2, 12, 32), y = c(3, 23, 53), cell = ids))
  image <- methods::new("VisiumV2",
    boundaries = list(centroids = centroids, segmentations = seg),
    assay = "custom.Polygons", key = "slice_", image = array(1, c(80, 60, 3)),
    scale.factors = Seurat::scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1)
  )
  if ("coords_x_orientation" %in% methods::slotNames(image)) {
    methods::slot(image, "coords_x_orientation") <- "horizontal"
  } else {
    attr(image, "coords_x_orientation") <- "horizontal"
  }
  object[["slice_with_underscores"]] <- image
  object$label <- factor(c("A", "B", NA), levels = c("A", "B", "C"))
  object
}

test_that("HD segmentation resolves cell type and official coordinate orientation", {
  object <- seurat6_spatial_object()
  info <- spatial_input_info(object)
  expect_identical(info$data_type, "cell")
  expect_identical(info$coordinate_units, "pixel")
  expect_identical(spatialcellchat_detect_technology(object, "auto", "slice_with_underscores"), "visium_hd")
  expect_identical(spatialcellchat_detect_level(object, "auto", "visium_hd", NULL, "slice_with_underscores"), "cell")
  expect_equal(SpatialCoordinates(object)$data$x, c(2, 12, 32))
  expect_equal(SpatialCoordinates(object)$data$y, c(3, 23, 53))
  object[["other"]] <- SeuratObject::CreateAssay5Object(counts = GetAssayData5(object))
  SeuratObject::DefaultAssay(object) <- "other"
  expect_identical(spatial_input_info(object)$assay, "custom.Polygons")
  renamed <- SeuratObject::RenameCells(object, add.cell.id = "sample")
  expect_identical(spatial_input_info(renamed)$data_type, "cell")
  legacy_image <- object[["slice_with_underscores"]]
  attr(legacy_image, "coords_x_orientation") <- NULL
  methods::slot(object, "images")[["slice_with_underscores"]] <- legacy_image
  expect_identical(spatial_image_x_orientation(object, "slice_with_underscores"), "horizontal")
})

test_that("compact and legacy segmentations preserve holes and ROI alignment", {
  for (compact in c(TRUE, FALSE)) {
    object <- seurat6_spatial_object(compact)
    before <- object
    p <- SpatialCellPlot(object, group.by = "label", xlim = c(0, 5), ylim = c(0, 7), overlay_image = FALSE)
    dat <- p$layers[[1]]$data
    expect_setequal(dat$cell_id, "cell.a_b")
    expect_length(unique(dat$ring_id), 2L)
    expect_equal(range(dat$x), c(0, 4))
    expect_equal(range(dat$y), c(74, 80))
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_identical(object, before)
    expect_error(SpatialCellPlot(object, boundary = "absent"), "not present")
    expect_error(SpatialCellPlot(object, xlim = c(100, 110)), "No cells")
    p2 <- SpatialSpotPlot(object, group.by = "label", xlim = c(0, 5), ylim = c(0, 7), overlay_image = FALSE)
    expect_equal(p2$data$x, 2)
    expect_equal(p2$data$y, 77)
  }
})

test_that("centroid-only cell plots and named boundary selection work", {
  object <- seurat6_spatial_object()
  object[["slice_with_underscores"]][["segmentations"]] <- NULL
  p <- SpatialCellPlot(object, group.by = "label", overlay_image = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 3L)
  expect_no_warning(ggplot2::ggplotGrob(p))
  expect_error(SpatialCellPlot(object, xlim = c(4, 1)), "increasing")
})

test_that("compound boundary identifiers cannot merge distinct polygons", {
  dat <- data.frame(
    cell_id = rep(c("a.b", "a"), each = 3),
    polygon_id = rep(c("c", "b.c"), each = 3), x = rep(c(0, 1, 0), 2), y = rep(c(0, 0, 1), 2)
  )
  expect_length(unique(spatial_boundary_validate(dat)$.polygon_group), 2L)
})

test_that("FOV morphology origin and loaded molecules share tissue coordinates", {
  object <- seurat6_spatial_object()
  object[["slice_with_underscores"]] <- NULL
  coords <- data.frame(x = c(12, 14, 16), y = c(22, 24, 26), row.names = colnames(object))
  object[["fov"]] <- SeuratObject::CreateFOV(coords,
    type = "centroids", assay = "custom.Polygons",
    molecules = data.frame(x = c(12, 15, 80), y = c(22, 25, 90), gene = c("g1", "g1", "g2"))
  )
  object@misc$atera.morphology <- list(
    image = matrix(1:12, 3, 4),
    pixel.size = 2, origin = c(x = 10, y = 20)
  )
  display <- spatial_dim_coords(object, image = "fov")
  expect_equal(display$image$xmin, 10)
  expect_equal(display$image$ymin, 20)
  expect_equal(display$image$width, 16)
  expect_equal(display$image$height, 28)
  expect_equal(dim(display$image$image), c(4L, 3L))
  p <- SpatialCellPlot(object, group.by = "label", molecules = "g1", xlim = c(10, 14), ylim = c(20, 24))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  molecule_data <- p$layers[[length(p$layers)]]$data
  expect_equal(molecule_data$x, 12)
  expect_equal(molecule_data$y, 22)
  expect_error(SpatialCellPlot(object, molecules = "absent"), "have not been loaded")
  expect_equal(SpatialCoordinates(object)$data$x, coords$x)
})

test_that("Atera uses official loader defaults and preserves the returned object", {
  skip_if_not("LoadAtera" %in% getNamespaceExports("Seurat"))
  object <- seurat6_spatial_object()
  object[["slice_with_underscores"]] <- NULL
  object[["atera_fov"]] <- SeuratObject::CreateFOV(data.frame(
    x = 1:3, y = c(2, 5, 9),
    row.names = colnames(object)
  ), type = "centroids", assay = "custom.Polygons")
  object@misc$atera.molecules <- list(handle = "retained", fov = "atera_fov")
  seen <- NULL
  local_mocked_bindings(LoadAtera = function(data.dir, ...) {
    seen <<- list(...)
    object
  }, .package = "Seurat")
  out <- ReadSpatialData(tempdir(), technology = "atera", bpcells.dir = "cache")
  expect_identical(seen, list(bpcells.dir = "cache"))
  expect_identical(out@misc$atera.molecules, object@misc$atera.molecules)
  expect_identical(spatial_input_info(out)$data_type, "cell")
  expect_identical(spatial_input_info(out)$coordinate_units, "micron")
})

test_that("mixed HD loader requests preserve bins and polygons separately", {
  object <- seurat6_spatial_object()
  bin_counts <- Matrix::Matrix(matrix(1:24, 6,
    dimnames = list(paste0("g", 1:6), paste0("bin", 1:4))
  ), sparse = TRUE)
  bins <- SeuratObject::CreateSeuratObject(counts = bin_counts, assay = "custom.008um")
  object <- merge(object, bins)
  object[["slice.008um"]] <- SeuratObject::CreateFOV(data.frame(
    x = 1:4, y = 4:1,
    row.names = colnames(bin_counts)
  ), type = "centroids", assay = "custom.008um")
  seen <- NULL
  local_mocked_bindings(Load10X_Spatial = function(data.dir, bin.size, segmentation.type = NULL, ...) {
    seen <<- bin.size
    object
  }, .package = "Seurat")
  out <- ReadSpatialData(tempdir(), technology = "visium_hd", bin.size = c(8, "polygons"))
  expect_equal(seen, c(8, "polygons"))
  expect_identical(spatial_input_info(out, image = "slice.008um")$data_type, "bin")
  expect_identical(spatial_input_info(out, image = "slice_with_underscores")$data_type, "cell")
  expect_error(spatial_input_info(out, image = "slice.008um", assay = "custom.Polygons"), "absent from assay")
  expect_error(ReadSpatialData(tempdir(), technology = "visium_hd", bin.size = c(8, "08")), "unique")
})
