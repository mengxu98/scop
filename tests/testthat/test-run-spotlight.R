make_spotlight_seurat_pair <- function() {
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
  list(spatial = spatial, reference = reference)
}

with_mock_spotlight <- function(fake_fun, code) {
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "SPOTlight")
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "SPOTlight")
      expect_identical(name, "SPOTlight")
      fake_fun
    }
  )
  force(code)
}

test_that("RunSPOTlight writes proportions and tool results", {
  pair <- make_spotlight_seurat_pair()
  fake_fun <- function(x, y, groups, mgs, gene_id, group_id, weight_id, min_prop, scale, ...) {
    expect_s4_class(x, "dgCMatrix")
    expect_s4_class(y, "dgCMatrix")
    expect_identical(groups, as.character(pair$reference$celltype))
    expect_true(all(c(gene_id, group_id, weight_id) %in% colnames(mgs)))
    expect_true(all(unique(groups) %in% unique(mgs[[group_id]])))
    expect_identical(min_prop, 0.01)
    expect_true(scale)
    list(
      mat = matrix(
        c(
          0.80, 0.20,
          0.35, 0.65,
          0.10, 0.90
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("Spot", 1:3), c("Alpha", "Beta"))
      ),
      res_ss = stats::setNames(c(0.01, 0.02, 0.03), paste0("Spot", 1:3)),
      NMF = list()
    )
  }

  with_mock_spotlight(fake_fun, {
    out <- RunSPOTlight(
      pair$spatial,
      reference = pair$reference,
      reference_label = "celltype",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SPOTlight_prop_Alpha), c(0.80, 0.35, 0.10))
  expect_equal(unname(out$SPOTlight_prop_Beta), c(0.20, 0.65, 0.90))
  expect_equal(unname(out$SPOTlight_dominant_type), c("Alpha", "Beta", "Beta"))
  expect_equal(unname(out$SPOTlight_max_prop), c(0.80, 0.65, 0.90))
  expect_true("SPOTlight" %in% names(out@tools))
  expect_equal(colnames(out@tools$SPOTlight$weights), c("Alpha", "Beta"))
  expect_equal(out@tools$SPOTlight$parameters$reference_label, "celltype")
})

test_that("RunSPOTlight validates inputs before backend work", {
  pair <- make_spotlight_seurat_pair()
  expect_error(
    RunSPOTlight(matrix(1, nrow = 2), pair$reference, verbose = FALSE),
    "Seurat"
  )
  with_mock_spotlight(function(...) stop("backend should not run"), {
    expect_error(
      RunSPOTlight(
        pair$spatial,
        pair$reference,
        reference_label = "missing",
        verbose = FALSE
      ),
      "reference_label"
    )
    expect_error(
      RunSPOTlight(
        pair$spatial,
        pair$reference,
        reference_label = "celltype",
        features = "AbsentGene",
        verbose = FALSE
      ),
      "No shared features"
    )
    expect_error(
      RunSPOTlight(
        pair$spatial,
        pair$reference,
        reference_label = "celltype",
        mgs = data.frame(gene = "Gene1", cluster = "Alpha"),
        verbose = FALSE
      ),
      "missing required"
    )
  })
})

test_that("RunSPOTlight accepts user marker table and transposed backend matrices", {
  pair <- make_spotlight_seurat_pair()
  mgs <- data.frame(
    gene = c("Gene1", "Gene2"),
    cluster = c("Alpha", "Beta"),
    weight = c(2, 3)
  )
  fake_fun <- function(x, y, groups, mgs, ...) {
    expect_equal(nrow(mgs), 2)
    matrix(
      c(
        0.7, 0.2, 0.1,
        0.3, 0.8, 0.9
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("Alpha", "Beta"), paste0("Spot", 1:3))
    )
  }

  with_mock_spotlight(fake_fun, {
    out <- RunSPOTlight(
      pair$spatial,
      reference = pair$reference,
      reference_label = "celltype",
      mgs = mgs,
      prefix = "SL",
      tool_name = "SL_tool",
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SL_prop_Alpha), c(0.7, 0.2, 0.1))
  expect_equal(unname(out$SL_prop_Beta), c(0.3, 0.8, 0.9))
  expect_true("SL_tool" %in% names(out@tools))
})

test_that("SPOTlight results reuse SCOP SpatialSpotPlot", {
  pair <- make_spotlight_seurat_pair()
  pair$spatial$SPOTlight_prop_Alpha <- c(0.8, 0.4, 0.1)
  pair$spatial$SPOTlight_prop_Beta <- c(0.2, 0.6, 0.9)
  pair$spatial$SPOTlight_dominant_type <- c("Alpha", "Beta", "Beta")
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      if (identical(packages, "scatterpie")) {
        skip_if_not_installed("scatterpie")
      }
      invisible(TRUE)
    }
  )
  p1 <- SpatialSpotPlot(
    pair$spatial,
    group.by = "SPOTlight_dominant_type",
    overlay_image = FALSE
  )
  p2 <- SpatialSpotPlot(
    pair$spatial,
    group.by = "SPOTlight_dominant_type",
    plot_type = "pie",
    overlay_image = FALSE
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("standard spatial workflow dispatches to RunSPOTlight", {
  pair <- make_spotlight_seurat_pair()
  original_standard_scop <- standard_scop
  testthat::local_mocked_bindings(
    RunSpotQC = function(srt, ...) srt,
    standard_scop = function(srt, ...) srt,
    RunSPOTlight = function(srt, reference, reference_label, assay, reference_assay, ...) {
      expect_identical(reference, pair$reference)
      expect_identical(reference_label, "celltype")
      expect_identical(assay, "RNA")
      expect_identical(reference_assay, "RNA")
      srt$SPOTlight_dominant_type <- "Alpha"
      srt
    }
  )

  out <- original_standard_scop(
    pair$spatial,
    workflow = "spatial",
    assay = "RNA",
    reference = pair$reference,
    reference_label = "celltype",
    reference_assay = "RNA",
    do_spot_qc = TRUE,
    do_spatial_variable_features = FALSE,
    do_spatial_cluster = FALSE,
    do_deconvolution = TRUE,
    deconvolution_method = "SPOTlight",
    verbose = FALSE
  )

  expect_equal(unique(out$SPOTlight_dominant_type), "Alpha")
  expect_equal(out@tools$standard_spatial_scop$parameters$deconvolution_method, "SPOTlight")
})
