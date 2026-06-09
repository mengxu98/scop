make_giotto_spatial_methods_seurat <- function() {
  set.seed(42)
  counts <- matrix(rpois(20 * 16, lambda = 8), nrow = 20)
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Spot", seq_len(ncol(counts)))
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt <- Seurat::FindVariableFeatures(srt, nfeatures = 12, verbose = FALSE)
  srt$x <- rep(seq_len(4), each = 4)
  srt$y <- rep(seq_len(4), times = 4)
  srt$celltype <- rep(c("A", "B", "C", "D"), length.out = ncol(srt))
  srt
}

test_that("Giotto follow-up methods reject non-Seurat input", {
  expect_error(
    RunGiottoCellProximity(matrix(1, nrow = 2, ncol = 2), group.by = "celltype", verbose = FALSE),
    "Seurat"
  )
  expect_error(
    RunGiottoSpatialGenes(matrix(1, nrow = 2, ncol = 2), verbose = FALSE),
    "Seurat"
  )
  expect_error(
    RunGiottoSpatialModules(matrix(1, nrow = 2, ncol = 2), verbose = FALSE),
    "Seurat"
  )
})

test_that("Giotto cell proximity stores pairwise enrichment tables only in tools", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()
  old_options <- options(
    giotto.has_conda = TRUE,
    giotto.use_conda = TRUE,
    giotto.update_param = TRUE,
    giotto.no_python_warn = FALSE
  )
  on.exit(options(old_options), add = TRUE)

  out <- RunGiottoCellProximity(
    srt,
    group.by = "celltype",
    number_of_simulations = 3,
    network_method = "Delaunay",
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_false("GiottoCellProximity" %in% colnames(out@meta.data))
  expect_s3_class(out@tools[["GiottoCellProximity"]][["enrichment"]], "data.frame")
  expect_gt(nrow(out@tools[["GiottoCellProximity"]][["enrichment"]]), 0)
  expect_true(all(c("group_1", "group_2") %in% colnames(out@tools[["GiottoCellProximity"]][["enrichment"]])))
  expect_equal(out@tools[["GiottoCellProximity"]][["cells"]], colnames(srt))
  expect_null(out@tools[["GiottoCellProximity"]][["giotto"]])
  expect_true(getOption("giotto.has_conda"))
  expect_true(getOption("giotto.use_conda"))
  expect_true(getOption("giotto.update_param"))
  expect_false(getOption("giotto.no_python_warn"))
})

test_that("Giotto spatial genes stores result table and optional variable features", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()

  out <- RunGiottoSpatialGenes(
    srt,
    features = rownames(srt),
    top_n = 5,
    network_method = "Delaunay",
    set_variable_features = TRUE,
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_s3_class(out@tools[["GiottoSpatialGenes"]][["results"]], "data.frame")
  expect_gt(nrow(out@tools[["GiottoSpatialGenes"]][["results"]]), 0)
  expect_lte(length(out@tools[["GiottoSpatialGenes"]][["top_features"]]), 5)
  expect_equal(out@misc[["GiottoSpatialGenes"]], out@tools[["GiottoSpatialGenes"]][["top_features"]])
  expect_true(all(SeuratObject::VariableFeatures(out) %in% rownames(out)))
  expect_null(out@tools[["GiottoSpatialGenes"]][["giotto"]])
})

test_that("Giotto spatial modules stores feature-level module objects in tools", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()

  out <- RunGiottoSpatialModules(
    srt,
    features = rownames(srt)[1:12],
    k = 2,
    network_method = "Delaunay",
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_false("GiottoSpatialModules" %in% colnames(out@meta.data))
  expect_true("spatial_cor" %in% names(out@tools[["GiottoSpatialModules"]]))
  expect_true("modules" %in% names(out@tools[["GiottoSpatialModules"]]))
  expect_type(out@tools[["GiottoSpatialModules"]][["module_tables"]], "list")
  expect_equal(out@tools[["GiottoSpatialModules"]][["features"]], rownames(srt)[1:12])
  expect_null(out@tools[["GiottoSpatialModules"]][["giotto"]])
})
