make_spatial_network_receipt_image <- function(
  lowres = NA_real_,
  hires = 0.5
) {
  counts <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  scale_factors <- structure(
    list(spot = 1, fiducial = 1, hires = hires, lowres = lowres),
    class = "scalefactors"
  )
  coordinates <- data.frame(
    imagerow = c(4, 8, 12, 16),
    imagecol = c(5, 10, 15, 20),
    row.names = colnames(srt)
  )
  srt[["slice"]] <- methods::new(
    "VisiumV1",
    image = array(1, dim = c(20, 30, 3)),
    scale.factors = scale_factors,
    coordinates = coordinates,
    spot.radius = 0.1,
    assay = SeuratObject::DefaultAssay(srt),
    key = "network_"
  )
  srt
}

test_that("RunSpatialNetwork receipt selects a usable display scale or the raw graph", {
  skip_if_not_installed("BiocNeighbors")

  hires_messages <- testthat::capture_messages(
    hires_out <- RunSpatialNetwork(
      make_spatial_network_receipt_image(),
      image = "slice",
      k = 1,
      verbose = TRUE
    )
  )
  hires_plain <- cli::ansi_strip(paste(hires_messages, collapse = "\n"))
  hires_call <- paste0(
    "SpatialNetworkPlot(<returned_object>, graph.name = \"slice_knn_k1\", ",
    "image.scale = \"hires\")"
  )
  expect_true(grepl(hires_call, hires_plain, fixed = TRUE))
  hires_plot <- NULL
  expect_no_error(
    hires_plot <- eval(parse(text = sub(
      "<returned_object>",
      "hires_out",
      hires_call,
      fixed = TRUE
    )))
  )
  expect_s3_class(hires_plot, "ggplot")

  raw_messages <- testthat::capture_messages(
    raw_out <- RunSpatialNetwork(
      make_spatial_network_receipt_image(hires = NA_real_),
      image = "slice",
      k = 1,
      verbose = TRUE
    )
  )
  raw_plain <- cli::ansi_strip(paste(raw_messages, collapse = "\n"))
  raw_call <- paste0(
    "SpatialNetworkPlot(res = <returned_object>@tools[[\"SpatialNetwork\"]], ",
    "graph.name = \"slice_knn_k1\")"
  )
  expect_true(grepl(raw_call, raw_plain, fixed = TRUE))
  expect_false(grepl("image.scale", raw_plain, fixed = TRUE))
  raw_plot <- NULL
  expect_no_error(
    raw_plot <- eval(parse(text = sub(
      "<returned_object>",
      "raw_out",
      raw_call,
      fixed = TRUE
    )))
  )
  expect_s3_class(raw_plot, "ggplot")
})
