make_coordinate_routing_object <- function() {
  counts <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$x <- c(0, 1, 2, 3)
  srt$y <- c(0, 0, 1, 1)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 0), row.names = c("Spot1", "Spot2")),
    type = "centroids", assay = assay, key = "cr1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(1, 1), row.names = c("Spot3", "Spot4")),
    type = "centroids", assay = assay, key = "cr2_"
  )
  srt
}

test_that("legacy metadata paths do not bypass strict multi-image selection", {
  srt <- make_coordinate_routing_object()

  expect_error(
    sgf_cpp_coords(
      srt,
      image = NULL,
      coord.cols = c("x", "y"),
      coordinate_space = "legacy_display"
    ),
    "Multiple spatial images"
  )
  expect_error(
    rctd_get_spatial_coords(
      srt,
      spot_ids = colnames(srt),
      image = NULL,
      coord.cols = c("x", "y"),
      coordinate_space = "legacy_display"
    ),
    "Multiple spatial images"
  )
})

test_that("shared coordinate routing rejects partial image spot sets", {
  srt <- make_coordinate_routing_object()

  expect_error(
    rctd_get_spatial_coords(
      srt,
      spot_ids = colnames(srt),
      image = "slice1",
      coord.cols = c("x", "y"),
      coordinate_space = "raw"
    ),
    "missing for one or more requested spots"
  )
})
