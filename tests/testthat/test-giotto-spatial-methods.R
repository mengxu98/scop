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

expect_giotto_does_not_modify_seurat <- function(srt, before_meta, before_tools, before_misc, before_variable_features) {
  expect_equal(srt@meta.data, before_meta)
  expect_equal(srt@tools, before_tools)
  expect_equal(srt@misc, before_misc)
  expect_equal(SeuratObject::VariableFeatures(srt), before_variable_features)
}

test_that("Giotto cell proximity returns standalone enrichment result", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()
  before_meta <- srt@meta.data
  before_tools <- srt@tools
  before_misc <- srt@misc
  before_variable_features <- SeuratObject::VariableFeatures(srt)
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
    network_name = "public_network",
    network_params = list(name = "actual_network"),
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "giotto2_result")
  expect_s3_class(out, "giotto2_cell_proximity")
  expect_false(inherits(out, "Seurat"))
  expect_false(is.null(out$giotto))
  expect_s3_class(out$enrichment, "data.frame")
  expect_gt(nrow(out$enrichment), 0)
  expect_true(all(c("group_1", "group_2") %in% colnames(out$enrichment)))
  expect_equal(out$cells, colnames(srt))
  expect_equal(out$parameters$network_name, "actual_network")
  expect_giotto_does_not_modify_seurat(srt, before_meta, before_tools, before_misc, before_variable_features)
  expect_true(getOption("giotto.has_conda"))
  expect_true(getOption("giotto.use_conda"))
  expect_true(getOption("giotto.update_param"))
  expect_false(getOption("giotto.no_python_warn"))
})

test_that("Giotto spatial genes returns standalone feature result", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()
  before_meta <- srt@meta.data
  before_tools <- srt@tools
  before_misc <- srt@misc
  before_variable_features <- SeuratObject::VariableFeatures(srt)

  out <- RunGiottoSpatialGenes(
    srt,
    features = rownames(srt),
    top_n = 5,
    network_method = "Delaunay",
    set_variable_features = TRUE,
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "giotto2_result")
  expect_s3_class(out, "giotto2_spatial_genes")
  expect_false(inherits(out, "Seurat"))
  expect_false(is.null(out$giotto))
  expect_s3_class(out$results, "data.frame")
  expect_gt(nrow(out$results), 0)
  expect_lte(length(out$top_features), 5)
  expect_true(all(out$top_features %in% rownames(srt)))
  expect_giotto_does_not_modify_seurat(srt, before_meta, before_tools, before_misc, before_variable_features)
})

test_that("Giotto spatial modules returns standalone module result", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_spatial_methods_seurat()
  before_meta <- srt@meta.data
  before_tools <- srt@tools
  before_misc <- srt@misc
  before_variable_features <- SeuratObject::VariableFeatures(srt)

  out <- RunGiottoSpatialModules(
    srt,
    features = rownames(srt)[1:12],
    k = 2,
    network_method = "Delaunay",
    store_giotto = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "giotto2_result")
  expect_s3_class(out, "giotto2_spatial_modules")
  expect_false(inherits(out, "Seurat"))
  expect_false(is.null(out$giotto))
  expect_true("spatial_cor" %in% names(out))
  expect_true("modules" %in% names(out))
  expect_type(out$module_tables, "list")
  expect_equal(out$features, rownames(srt)[1:12])
  expect_giotto_does_not_modify_seurat(srt, before_meta, before_tools, before_misc, before_variable_features)
})
