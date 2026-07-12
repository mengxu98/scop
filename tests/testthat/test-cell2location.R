make_cell2location_pair <- function() {
  spatial_counts <- matrix(
    c(
      5, 2, 0,
      0, 4, 3,
      2, 1, 6,
      1, 0, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:3))
  )
  reference_counts <- matrix(
    c(
      4, 3, 0, 0,
      0, 1, 5, 4,
      3, 2, 1, 0,
      1, 0, 2, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:4))
  )
  spatial <- Seurat::CreateSeuratObject(
    methods::as(Matrix::Matrix(spatial_counts, sparse = TRUE), "dgCMatrix")
  )
  spatial$x <- c(1, 2, 3)
  spatial$y <- c(1, 2, 1)
  spatial$sample <- "slide1"
  reference <- Seurat::CreateSeuratObject(
    methods::as(Matrix::Matrix(reference_counts, sparse = TRUE), "dgCMatrix")
  )
  reference$celltype <- c("Alpha", "Alpha", "Beta", "Beta")
  reference$sample <- c("r1", "r1", "r2", "r2")
  list(spatial = spatial, reference = reference)
}

test_that("cell2location environment module pins compatible versions", {
  req <- env_requirements(modules = "cell2location")
  expect_identical(unname(req$packages[["cell2location"]]), "cell2location==0.1.5")
  expect_identical(unname(req$packages[["scvi-tools"]]), "scvi-tools==1.3.3")
  expect_identical(unname(req$packages[["numba"]]), "numba==0.66.0")
  expect_identical(unname(req$packages[["llvmlite"]]), "llvmlite==0.48.0")
  modules <- getFromNamespace("normalize_env_modules", "scop")("cell2location")
  expect_setequal(modules, c("cell2location", "scanpy", "scvi"))
})

test_that("cell2location input preparation filters labels and aligns genes", {
  pair <- make_cell2location_pair()
  prepared <- getFromNamespace("cell2location_prepare_inputs", "scop")(
    srt = pair$spatial,
    reference = pair$reference,
    reference_label = "celltype",
    reference_signatures = NULL,
    assay = "RNA",
    reference_assay = "RNA",
    layer = "counts",
    reference_layer = "counts",
    features = c("Gene3", "Gene2", "Gene1"),
    spatial_batch = "sample",
    reference_batch = "sample",
    reference_covariates = NULL,
    min_cells = 2L,
    verbose = FALSE
  )
  expect_identical(prepared$features, c("Gene3", "Gene2", "Gene1"))
  expect_identical(colnames(prepared$spatial), paste0("Spot", 1:3))
  expect_identical(colnames(prepared$reference), paste0("Cell", 1:4))
  expect_equal(prepared$summary$n_reference_cells, 4L)
})

test_that("cell2location rejects normalized input and invalid signatures", {
  pair <- make_cell2location_pair()
  pair$spatial[["RNA"]] <- SeuratObject::SetAssayData(
    pair$spatial[["RNA"]],
    layer = "data",
    new.data = SeuratObject::LayerData(pair$spatial, assay = "RNA", layer = "counts") / 3
  )
  counts_fun <- getFromNamespace("cell2location_counts", "scop")
  expect_error(counts_fun(pair$spatial, "RNA", "data", "Spatial"), "raw integer counts")

  signatures_fun <- getFromNamespace("cell2location_read_signatures", "scop")
  bad <- matrix(-1, nrow = 2, dimnames = list(c("g1", "g2"), "A"))
  expect_error(signatures_fun(bad), "finite and non-negative")
})

test_that("RunCell2location writes abundance, proportions, and reproducible tools", {
  pair <- make_cell2location_pair()
  signatures <- matrix(
    c(4, 1, 2, 1, 1, 4, 1, 2),
    nrow = 4,
    dimnames = list(paste0("Gene", 1:4), c("Alpha", "Beta"))
  )
  result_dir <- tempfile("cell2location result with spaces ")
  dir.create(result_dir)
  testthat::local_mocked_bindings(
    .package = "scop",
    cell2location_check_python = function(...) "python",
    cell2location_runner_path = function() "cell2location_runner.py",
    cell2location_write_json = function(...) invisible(NULL),
    cell2location_read_json = function(...) list(status = "complete"),
    srt_to_h5ad = function(srt, path, ...) {
      file.create(path)
      invisible(path)
    },
    cell2location_run_system2 = function(command, args, env, stdout, stderr) {
      files <- getFromNamespace("cell2location_result_files", "scop")(result_dir)
      dir.create(dirname(files$abundance), recursive = TRUE, showWarnings = FALSE)
      dir.create(dirname(files$signatures), recursive = TRUE, showWarnings = FALSE)
      abundance <- data.frame(
        Alpha = c(8, 2, 1),
        Beta = c(2, 6, 9),
        row.names = paste0("Spot", 1:3)
      )
      proportions <- abundance / rowSums(abundance)
      utils::write.csv(abundance, files$abundance)
      utils::write.csv(proportions, files$proportions)
      utils::write.csv(signatures, files$signatures)
      file.create(files$manifest)
      file.create(stdout)
      file.create(stderr)
      0L
    }
  )

  out <- RunCell2location(
    pair$spatial,
    result_dir = result_dir,
    reference_signatures = signatures,
    assay = "RNA",
    spatial_batch = "sample",
    verbose = FALSE
  )

  expect_equal(unname(out$Cell2location_abundance_Alpha), c(8, 2, 1))
  expect_equal(unname(out$Cell2location_prop_Alpha), c(0.8, 0.25, 0.1))
  expect_equal(unname(out$Cell2location_dominant_type), c("Alpha", "Beta", "Beta"))
  expect_true("Cell2location" %in% names(out@tools))
  expect_equal(out@tools$Cell2location$manifest$status, "complete")
  expect_identical(colnames(out@tools$Cell2location$reference_signatures), c("Alpha", "Beta"))
})

test_that("RunCell2location does not mutate Seurat when Python fails", {
  pair <- make_cell2location_pair()
  signatures <- matrix(
    1,
    nrow = 4,
    ncol = 2,
    dimnames = list(paste0("Gene", 1:4), c("Alpha", "Beta"))
  )
  result_dir <- tempfile("cell2location_failed_")
  testthat::local_mocked_bindings(
    .package = "scop",
    cell2location_check_python = function(...) "python",
    cell2location_runner_path = function() "cell2location_runner.py",
    cell2location_write_json = function(...) invisible(NULL),
    srt_to_h5ad = function(srt, path, ...) {
      file.create(path)
      invisible(path)
    },
    cell2location_run_system2 = function(command, args, env, stdout, stderr) {
      writeLines("synthetic backend failure", stderr)
      file.create(stdout)
      2L
    }
  )
  expect_error(
    RunCell2location(
      pair$spatial,
      result_dir = result_dir,
      reference_signatures = signatures,
      assay = "RNA",
      verbose = FALSE
    ),
    "synthetic backend\\s+failure"
  )
  expect_false(any(startsWith(colnames(pair$spatial@meta.data), "Cell2location_")))
  expect_false("Cell2location" %in% names(pair$spatial@tools))
})

test_that("Cell2locationPlot exposes all result views", {
  pair <- make_cell2location_pair()
  abundance <- matrix(
    c(8, 2, 2, 6, 1, 9),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Spot", 1:3), c("Alpha", "Beta"))
  )
  proportions <- abundance / rowSums(abundance)
  pair$spatial@tools$Cell2location <- list(
    abundance = abundance,
    proportions = proportions
  )
  pair$spatial$Cell2location_dominant_type <- c("Alpha", "Beta", "Beta")

  expect_s3_class(Cell2locationPlot(pair$spatial, "proportion", overlay_image = FALSE), "ggplot")
  expect_s3_class(Cell2locationPlot(pair$spatial, "abundance", overlay_image = FALSE), "ggplot")
  expect_s3_class(Cell2locationPlot(pair$spatial, "dominant", overlay_image = FALSE), "ggplot")
  if (requireNamespace("scatterpie", quietly = TRUE)) {
    expect_s3_class(Cell2locationPlot(pair$spatial, "pie", overlay_image = FALSE), "ggplot")
  }
})

test_that("standard spatial workflow dispatches cell2location signatures", {
  pair <- make_cell2location_pair()
  signatures <- matrix(
    1,
    nrow = 4,
    ncol = 2,
    dimnames = list(paste0("Gene", 1:4), c("Alpha", "Beta"))
  )
  original <- getFromNamespace("standard_spatial_scop", "scop")
  called <- FALSE
  testthat::local_mocked_bindings(
    .package = "scop",
    standard_scop = function(srt, ...) srt,
    RunCell2location = function(srt, result_dir, reference_signatures, ...) {
      called <<- TRUE
      expect_identical(result_dir, "c2l")
      expect_identical(reference_signatures, signatures)
      srt$Cell2location_dominant_type <- "Alpha"
      srt
    }
  )

  out <- original(
    pair$spatial,
    assay = "RNA",
    reference = NULL,
    reference_label = NULL,
    do_spot_qc = FALSE,
    do_spatial_variable_features = FALSE,
    do_deconvolution = TRUE,
    deconvolution_method = "Cell2location",
    deconvolution_params = list(
      result_dir = "c2l",
      reference_signatures = signatures
    ),
    verbose = FALSE
  )
  expect_true(called)
  expect_equal(unique(out$Cell2location_dominant_type), "Alpha")
  expect_identical(
    out@tools$standard_spatial_scop$parameters$deconvolution_method,
    "Cell2location"
  )
})
