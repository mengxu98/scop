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

test_that("CARD argument routing supports Bioconductor underscore formals", {
  args <- list(
    ct.varname = ".scop_cell_type",
    ct_varname = ".scop_cell_type",
    sample.varname = ".scop_sample",
    sample_varname = ".scop_sample",
    ct.select = c("Alpha", "Beta"),
    ct_select = c("Alpha", "Beta")
  )
  fun <- function(ct_varname, sample_varname, ct_select) NULL

  matched <- scop:::card_match_formals(fun, args)
  expect_named(matched, c("ct_varname", "sample_varname", "ct_select"))
  expect_identical(matched$ct_varname, ".scop_cell_type")
  expect_identical(matched$sample_varname, ".scop_sample")
  expect_identical(matched$ct_select, c("Alpha", "Beta"))
})

test_that("CARD uses an installed CARDspa backend without installing CARD", {
  testthat::local_mocked_bindings(
    check_r = function(...) stop("CARD installation must not be attempted"),
    get_namespace_fun = function(package, name) {
      if (identical(package, "CARD")) stop("CARD is not installed")
      expect_identical(package, "CARDspa")
      expect_true(name %in% c("createCARDObject", "CARD_deconvolution"))
      function(...) NULL
    },
    .package = "scop"
  )

  expect_identical(scop:::card_resolve_backend_package(), "CARDspa")
})

test_that("CARD checks its repository only when no installed backend is usable", {
  installed <- FALSE
  testthat::local_mocked_bindings(
    check_r = function(repository, verbose = FALSE) {
      expect_identical(repository, "YingMa0107/CARD")
      expect_false(verbose)
      installed <<- TRUE
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      if (isTRUE(installed) && identical(package, "CARD")) {
        return(function(...) NULL)
      }
      stop("backend unavailable")
    },
    .package = "scop"
  )

  expect_identical(scop:::card_resolve_backend_package(), "CARD")
  expect_true(installed)
})

test_that("CARDspa one-step API removes backend-unsupported informative-zero spots", {
  st_counts <- methods::as(Matrix::Matrix(
    matrix(
      c(1, 2, 3, 4, 0, 6),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
    ),
    sparse = TRUE
  ), "dgCMatrix")
  ref_counts <- methods::as(Matrix::Matrix(
    matrix(1:8, nrow = 2, dimnames = list(c("G1", "G2"), paste0("C", 1:4))),
    sparse = TRUE
  ), "dgCMatrix")
  ref_meta <- data.frame(
    .scop_cell_type = c("A", "A", "B", "B"),
    .scop_sample = "sample1",
    row.names = colnames(ref_counts)
  )
  coords <- data.frame(x = 1:3, y = c(1, 2, 1), row.names = colnames(st_counts))
  one_step <- function(
    sc_count,
    sc_meta,
    spatial_count,
    spatial_location,
    ct_varname,
    ct_select,
    sample_varname,
    mincountgene,
    mincountspot
  ) {
    expect_identical(ct_varname, ".scop_cell_type")
    expect_identical(sample_varname, ".scop_sample")
    expect_identical(colnames(spatial_count), c("S1", "S2"))
    expect_identical(rownames(spatial_location), c("S1", "S2"))
    list(Proportion_CARD = matrix(
      c(0.8, 0.4, 0.2, 0.6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("A", "B"), c("S1", "S2"))
    ))
  }
  create_one_step <- function(spatial_count, ...) {
    list(
      info_parameters = list(
        ct.select = c("A", "B"),
        ct.varname = ".scop_cell_type",
        sample.varname = ".scop_sample"
      ),
      sc_eset = ref_counts,
      spatial_countMat = spatial_count
    )
  }
  create_ref <- function(...) {
    list(basis = matrix(
      c(1, 0, 0, 1),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("A", "B"))
    ))
  }
  testthat::local_mocked_bindings(
    card_resolve_backend_package = function() "CARDspa",
    get_namespace_fun = function(package, name) {
      switch(name,
        createCARDObject = create_one_step,
        CARD_deconvolution = one_step,
        create_ref = create_ref,
        select_info = function(...) "G1",
        stop("unexpected CARDspa function")
      )
    }
  )

  result <- scop:::card_run_backend(
    st_counts = st_counts,
    ref_counts = ref_counts,
    ref_meta = ref_meta,
    coords = coords,
    ct_select = c("A", "B"),
    minCountGene = 100,
    minCountSpot = 5,
    create_card_params = list(),
    card_deconvolution_params = list()
  )
  expect_identical(result$package, "CARDspa")
  expect_equal(dim(result$weights), c(2L, 2L))
  expect_identical(result$dropped_spots, "S3")
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
