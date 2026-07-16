make_card_seurat_pair <- function() {
  spatial_counts <- matrix(
    c(
      10, 8, 1,
      0, 2, 9,
      6, 0, 1,
      1, 7, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:3))
  )
  ref_counts <- matrix(
    c(
      9, 7, 1, 0,
      0, 1, 8, 7,
      5, 4, 1, 0,
      1, 0, 6, 5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:4))
  )
  spatial <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(spatial_counts, sparse = TRUE), "dgCMatrix")
  )
  spatial$col <- c(1, 2, 3)
  spatial$row <- c(1, 2, 1)
  reference <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(ref_counts, sparse = TRUE), "dgCMatrix")
  )
  reference$celltype <- c("Alpha", "Alpha", "Beta", "Beta")
  reference$sample <- c("S1", "S1", "S2", "S2")
  list(spatial = spatial, reference = reference)
}

with_mock_card <- function(code) {
  create_fun <- function(
    sc_count,
    sc_meta,
    spatial_count,
    spatial_location,
    ct.varname,
    sample.varname,
    ct.select,
    minCountGene,
    minCountSpot
  ) {
    expect_s4_class(sc_count, "dgCMatrix")
    expect_s4_class(spatial_count, "dgCMatrix")
    expect_equal(colnames(spatial_count), paste0("Spot", 1:3))
    expect_equal(rownames(spatial_location), paste0("Spot", 1:3))
    expect_equal(ct.varname, ".scop_cell_type")
    expect_equal(sample.varname, ".scop_sample")
    expect_equal(ct.select, c("Alpha", "Beta"))
    expect_equal(minCountGene, 100)
    expect_equal(minCountSpot, 5)
    list(sc_meta = sc_meta)
  }
  deconv_fun <- function(CARD_object) {
    CARD_object$Proportion_CARD <- matrix(
      c(
        0.80, 0.35, 0.10,
        0.20, 0.65, 0.90
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("Alpha", "Beta"), paste0("Spot", 1:3))
    )
    CARD_object
  }
  testthat::local_mocked_bindings(
    card_resolve_backend_package = function() "CARD",
    get_namespace_fun = function(package, name) {
      expect_identical(package, "CARD")
      switch(name,
        createCARDObject = create_fun,
        CARD_deconvolution = deconv_fun,
        stop("unexpected function")
      )
    }
  )
  force(code)
}

test_that("RunCARD writes proportions and tool results", {
  pair <- make_card_seurat_pair()
  with_mock_card({
    out <- RunCARD(
      pair$spatial,
      reference = pair$reference,
      reference_label = "celltype",
      sample_varname = "sample",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$CARD_prop_Alpha), c(0.80, 0.35, 0.10))
  expect_equal(unname(out$CARD_prop_Beta), c(0.20, 0.65, 0.90))
  expect_equal(unname(out$CARD_dominant_type), c("Alpha", "Beta", "Beta"))
  expect_equal(unname(out$CARD_max_prop), c(0.80, 0.65, 0.90))
  expect_true("CARD" %in% names(out@tools))
  expect_equal(colnames(out@tools$CARD$weights), c("Alpha", "Beta"))
  expect_identical(rownames(out@tools$CARD$proportions), colnames(out))
  expect_identical(out@tools$CARD$cells, colnames(out))
  expect_identical(out@tools$CARD$source$coordinate_space, "raw")
  expect_equal(out@tools$CARD$backend_package, "CARD")
  expect_named(out@tools$CARD$summary, c("n_spots", "n_types", "dominant_counts", "max_prop"))
})

test_that("RunCARD validates inputs before backend work", {
  pair <- make_card_seurat_pair()
  expect_error(
    RunCARD(matrix(1, nrow = 2), pair$reference, "celltype", verbose = FALSE),
    "Seurat"
  )
  with_mock_card({
    expect_error(
      RunCARD(
        pair$spatial,
        pair$reference,
        reference_label = "missing",
        verbose = FALSE
      ),
      "reference_label"
    )
    expect_error(
      RunCARD(
        pair$spatial,
        pair$reference,
        reference_label = "celltype",
        features = "AbsentGene",
        verbose = FALSE
      ),
      "No shared features"
    )
    expect_error(
      RunCARD(
        pair$spatial,
        pair$reference,
        reference_label = "celltype",
        sample_varname = "missing",
        verbose = FALSE
      ),
      "sample_varname"
    )
  })
})

test_that("CARD stored results use SpatialDeconvolutionPlot", {
  pair <- make_card_seurat_pair()
  pair$spatial$CARD_prop_Alpha <- c(0.8, 0.4, 0.1)
  pair$spatial$CARD_prop_Beta <- c(0.2, 0.6, 0.9)
  pair$spatial$CARD_dominant_type <- c("Alpha", "Beta", "Beta")
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      if (identical(packages, "scatterpie")) {
        testthat::skip_if_not_installed("scatterpie")
      }
      invisible(TRUE)
    }
  )
  pair$spatial@tools$CARD <- scop:::spatial_result_build(
    bundle = list(
      proportions = as.matrix(pair$spatial@meta.data[, c("CARD_prop_Alpha", "CARD_prop_Beta")]),
      cells = colnames(pair$spatial)
    ),
    method = "CARD",
    result_type = "deconvolution",
    provenance = list(producer = "RunCARD", backend_id = "card")
  )
  colnames(pair$spatial@tools$CARD$proportions) <- c("Alpha", "Beta")
  p1 <- SpatialDeconvolutionPlot(
    pair$spatial,
    tool_name = "CARD",
    plot_type = "dominant",
    overlay_image = FALSE
  )
  p2 <- SpatialDeconvolutionPlot(
    pair$spatial,
    tool_name = "CARD",
    plot_type = "pie",
    overlay_image = FALSE
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
