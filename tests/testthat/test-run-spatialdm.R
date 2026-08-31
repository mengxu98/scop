make_spatialdm_test_object <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  counts <- Matrix::Matrix(
    matrix(c(
      4, 0, 1, 2,
      0, 3, 2, 0,
      2, 1, 0, 4
    ), nrow = 3, byrow = TRUE),
    sparse = TRUE
  )
  rownames(counts) <- c("L1", "R1", "G1")
  colnames(counts) <- paste0("spot", 1:4)
  object <- Seurat::CreateSeuratObject(counts = counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object[["col"]] <- c(0, 1, 0, 1)
  object[["row"]] <- c(0, 0, 1, 1)
  object
}

mock_spatialdm_execution <- function(input, species, lr.database, parameters, result.name, envname, verbose) {
  cells <- input$cells
  global <- data.frame(
    interaction = "L1_R1", Ligand0 = "L1", Receptor0 = "R1",
    moran_r = 0.8, z_score = 4, pvalue = 0.001, fdr = 0.01,
    selected = TRUE, n_spots = 2, score_type = "spatial_association",
    stringsAsFactors = FALSE
  )
  rownames(global) <- global$interaction
  local_matrix <- matrix(c(0.5, 0.1, 0.2, 0.7),
    nrow = 1,
    dimnames = list("L1_R1", cells)
  )
  weights <- Matrix::Diagonal(length(cells))
  dimnames(weights) <- list(cells, cells)
  list(
    global = global,
    local = list(
      local_i = local_matrix, local_i_r = local_matrix,
      local_z = local_matrix, local_p = local_matrix,
      selected_spots = local_matrix > 0.4
    ),
    coordinates = data.frame(cell_id = cells, x = c(0, 1, 0, 1), y = c(0, 0, 1, 1)),
    weights = list(secreted = weights, contact = weights),
    manifest = list(versions = list(SpatialDM = "0.3.1"))
  )
}

test_that("SpatialDM stores spatial association semantics without fake edges", {
  srt <- make_spatialdm_test_object()
  testthat::local_mocked_bindings(
    spatialdm_execute = mock_spatialdm_execution,
    .package = "scop"
  )
  out <- scop::RunSpatialDM(
    srt,
    lr.database = data.frame(
      ligand = "L1", receptor = "R1", annotation = "Secreted Signaling"
    ),
    image = NULL, coord.cols = c("col", "row"), l = 1.2,
    verbose = FALSE
  )
  expect_identical(out@tools$SpatialDM$method, "SpatialDM")
  expect_identical(out@tools$SpatialDM$result_type, "spatial_ligand_receptor_association")
  expect_identical(out@tools$SpatialDM$global$score_type, "spatial_association")
  expect_identical(out@tools$SpatialDM$parameters$coordinate_space, "raw")
  expect_identical(out@tools$SpatialDM$cells, colnames(srt))
  expect_identical(out@tools$SpatialDM$features, rownames(srt))
  expect_identical(out@tools$SpatialDM$source$spatial$coordinate_space, "raw")
  expect_identical(out@tools$SpatialDM$provenance$commit, getFromNamespace(".spatialdm_backend_commit", "scop"))
  expect_true(is.null(out@tools$CCC))
  expect_equal(scop::GetSpatialDMResult(out, type = "global")$moran_r, 0.8)
  expect_s3_class(scop::SpatialDMPlot(out, plot_type = "global"), "ggplot")
  expect_true(inherits(
    scop::SpatialDMPlot(out, plot_type = "weights", spot = "spot1"),
    c("ggplot", "patchwork")
  ))
  expect_true(inherits(
    scop::SpatialDMPlot(out, plot_type = "local", pair = "L1_R1"),
    c("patchwork", "ggplot")
  ))
})

test_that("SpatialDM method does not require group.by in RunCCC", {
  srt <- make_spatialdm_test_object()
  mock_runner <- function(srt, verbose = TRUE, ...) {
    srt@tools$SpatialDM <- list(
      method = "SpatialDM", results = list(default = list()),
      active_result = "default", result = list(),
      global = data.frame(), local = list(), parameters = list(),
      provenance = list()
    )
    srt
  }
  testthat::local_mocked_bindings(
    RunSpatialDM = mock_runner,
    .package = "scop"
  )
  out <- scop::RunCCC(
    srt,
    methods = "SpatialDM",
    method_params = list(SpatialDM = list(coord.cols = c("col", "row"))),
    verbose = FALSE
  )
  expect_identical(out@tools$RunCCC$completed_methods, "SpatialDM")
  expect_identical(out@tools$CCC$long_table, data.frame())
  expect_identical(out@tools$CCC$association_methods, "SpatialDM")
})

test_that("SpatialDM never silently selects one of several images", {
  srt <- make_spatialdm_test_object()
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 1), row.names = c("spot1", "spot2")),
    type = "centroids", assay = assay, key = "sdm1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(1, 2), row.names = c("spot3", "spot4")),
    type = "centroids", assay = assay, key = "sdm2_"
  )
  expect_error(
    scop::RunSpatialDM(srt, l = 1.2, verbose = FALSE),
    "Multiple spatial images"
  )
})

test_that("SpatialDM backend failure leaves the object unchanged", {
  srt <- make_spatialdm_test_object()
  before <- srt@tools
  testthat::local_mocked_bindings(
    spatialdm_execute = function(...) stop("SpatialDM backend failed"),
    .package = "scop"
  )
  expect_error(
    scop::RunSpatialDM(srt, l = 1.2, verbose = FALSE),
    "backend failed"
  )
  expect_identical(srt@tools, before)
})

test_that("SpatialDM is registered as a non-edge CCC method", {
  spec <- getFromNamespace("ccc_method_spec", "scop")("SpatialDM")
  expect_false(spec$requires_group)
  expect_true(spec$requires_spatial)
  expect_false(spec$supports_unified_edges)
  expect_identical(
    getFromNamespace("ccc_bundle_long_table", "scop")(
      structure(list(), class = "Seurat"), "SpatialDM", bundle = list()
    ),
    data.frame()
  )
})

test_that("SpatialDM uses a standalone pinned Python environment", {
  req <- getFromNamespace("env_requirements", "scop")(modules = "spatialdm", verbose = FALSE)
  expect_identical(req$python, "3.10-1")
  expect_match(unname(req$packages[["SpatialDM"]]), "StatBiomed/SpatialDM.git@9b0f559d")
  expect_match(unname(req$packages[["scipy"]]), "scipy==1.11.4")
  expect_match(unname(req$packages[["anndata"]]), "anndata==0.10.8")
  expect_error(
    getFromNamespace("normalize_env_modules", "scop")(c("spatialdm", "scanpy")),
    "standalone"
  )
})
