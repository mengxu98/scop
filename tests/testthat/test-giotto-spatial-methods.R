make_giotto_method_seurat <- function(n_genes = 8, n_spots = 12) {
  set.seed(11)
  counts <- matrix(
    rpois(n_genes * n_spots, lambda = 3),
    nrow = n_genes,
    dimnames = list(paste0("Gene", seq_len(n_genes)), paste0("Spot", seq_len(n_spots)))
  )
  counts[counts == 0] <- 1
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- rep(seq_len(4), each = 3)
  srt$y <- rep(seq_len(3), times = 4)
  srt$group <- rep(c("A", "B", "C"), length.out = n_spots)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)
  srt
}

test_that("Giotto standalone wrappers reject non-Seurat input", {
  expect_error(RunGiottoCellProximity(matrix(1), group.by = "x", verbose = FALSE), "Seurat")
  expect_error(RunGiottoSpatialGenes(matrix(1), verbose = FALSE), "Seurat")
  expect_error(RunGiottoHMRF(matrix(1), spatial_genes = "Gene1", verbose = FALSE), "Seurat")
  expect_error(RunGiottoSpatialModules(matrix(1), verbose = FALSE), "Seurat")
})

test_that("RunGiottoCellProximity stores proximity table without full Giotto object", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_method_seurat()
  old_options <- options(giotto.has_conda = TRUE, giotto.use_conda = TRUE)
  on.exit(options(old_options), add = TRUE)

  out <- RunGiottoCellProximity(
    srt,
    group.by = "group",
    number_of_simulations = 10,
    store_giotto = FALSE,
    verbose = FALSE,
    seed = 11
  )
  expect_s4_class(out, "Seurat")
  expect_true("GiottoCellProximity" %in% names(out@tools))
  expect_s3_class(out@tools[["GiottoCellProximity"]][["enrichment"]], "data.frame")
  expect_null(out@tools[["GiottoCellProximity"]][["giotto"]])
  expect_true(getOption("giotto.has_conda"))
  expect_true(getOption("giotto.use_conda"))
})

test_that("RunGiottoSpatialGenes stores top features and can update VariableFeatures", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_method_seurat(n_genes = 10, n_spots = 15)

  out <- RunGiottoSpatialGenes(
    srt,
    features = rownames(srt),
    top_n = 5,
    store_giotto = FALSE,
    binSpect_params = list(do_parallel = FALSE, cores = 1, verbose = FALSE),
    verbose = FALSE,
    seed = 11
  )
  expect_s4_class(out, "Seurat")
  expect_true("GiottoSpatialGenes" %in% names(out@tools))
  expect_s3_class(out@tools[["GiottoSpatialGenes"]][["results"]], "data.frame")
  expect_lte(length(out@misc[["GiottoSpatialGenes"]]), 5)
  expect_equal(SeuratObject::VariableFeatures(out), out@misc[["GiottoSpatialGenes"]])
  expect_null(out@tools[["GiottoSpatialGenes"]][["giotto"]])
})

test_that("RunGiottoSpatialModules stores feature-level module results", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_giotto_method_seurat(n_genes = 10, n_spots = 15)

  out <- RunGiottoSpatialModules(
    srt,
    features = rownames(srt),
    k = 2,
    store_giotto = FALSE,
    verbose = FALSE,
    seed = 11
  )
  expect_s4_class(out, "Seurat")
  expect_true("GiottoSpatialModules" %in% names(out@tools))
  expect_true(length(out@tools[["GiottoSpatialModules"]][["module_tables"]]) > 0L)
  expect_null(out@tools[["GiottoSpatialModules"]][["giotto"]])
})

test_that("RunGiottoHMRF skips without explicit HMRF smoke-test opt-in", {
  testthat::skip_if_not_installed("Giotto")
  testthat::skip_if_not(
    identical(Sys.getenv("SCOP_RUN_GIOTTO_HMRF"), "true"),
    "Set SCOP_RUN_GIOTTO_HMRF=true with a Giotto HMRF Python environment"
  )
  srt <- make_giotto_method_seurat(n_genes = 10, n_spots = 15)
  out <- RunGiottoHMRF(
    srt,
    spatial_genes = rownames(srt)[1:6],
    k = 2,
    betas = 0,
    dimensions_to_use = 1:2,
    output_folder = tempdir(),
    store_giotto = FALSE,
    verbose = FALSE,
    seed = 11
  )
  expect_s4_class(out, "Seurat")
  expect_true("GiottoHMRF" %in% names(out@tools))
  expect_null(out@tools[["GiottoHMRF"]][["giotto"]])
})
