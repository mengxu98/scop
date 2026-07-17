make_spatial_dwls_pair <- function() {
  spatial_counts <- matrix(
    c(
      10, 1, 8,
      1, 9, 2,
      7, 1, 8,
      1, 7, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:3))
  )
  ref_counts <- matrix(
    c(
      9, 8, 1, 1,
      1, 1, 8, 9,
      8, 7, 1, 2,
      1, 2, 7, 8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:4))
  )
  spatial <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(spatial_counts, sparse = TRUE), "dgCMatrix")
  )
  spatial$x <- c(1, 2, 3)
  spatial$y <- c(1, 1, 2)
  reference <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(ref_counts, sparse = TRUE), "dgCMatrix")
  )
  reference$celltype <- c("Alpha", "Alpha", "Beta", "Beta")
  list(spatial = spatial, reference = reference)
}

test_that("RunSpatialDWLS writes standard deconvolution metadata and summary", {
  pair <- make_spatial_dwls_pair()
  out <- RunSpatialDWLS(
    pair$spatial,
    reference = pair$reference,
    reference_label = "celltype",
    normalize = FALSE,
    verbose = FALSE
  )

  expect_true(all(c(
    "SpatialDWLS_prop_Alpha",
    "SpatialDWLS_prop_Beta",
    "SpatialDWLS_dominant_type",
    "SpatialDWLS_max_prop"
  ) %in% colnames(out@meta.data)))
  expect_true("SpatialDWLS" %in% names(out@tools))
  expect_equal(colnames(out@tools$SpatialDWLS$weights), c("Alpha", "Beta"))
  expect_identical(rownames(out@tools$SpatialDWLS$proportions), colnames(out))
  expect_identical(out@tools$SpatialDWLS$cells, colnames(out))
  expect_identical(out@tools$SpatialDWLS$source$coordinate_space, "raw")
  expect_named(out@tools$SpatialDWLS$summary, c("n_spots", "n_types", "dominant_counts", "max_prop"))
  expect_equal(out@tools$SpatialDWLS$summary$n_spots, 3)
  expect_s3_class(SpatialDeconvolutionPlot(
    out,
    tool_name = "SpatialDWLS",
    cell_types = "Alpha",
    overlay_image = FALSE,
    coord.cols = c("x", "y")
  ), "ggplot")
})
