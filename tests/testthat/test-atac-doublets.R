test_that("ATAC scDblFinder uses capture groups and preserves partial-assay cells", {
  skip_if_not_installed("scDblFinder")
  skip_if_not_installed("SingleCellExperiment")
  counts <- withr::with_seed(83, Matrix::Matrix(matrix(
    stats::rpois(200 * 120, 2), 200, 120,
    dimnames = list(paste0("p", 1:200), paste0("cell", 1:120))
  ), sparse = TRUE))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = seq(1, by = 200, length.out = 200),
    end = seq(100, by = 200, length.out = 200)
  ))
  atac <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  rna <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1, 2, 121,
    dimnames = list(c("g1", "g2"), c(colnames(atac), "rna_only"))
  ), sparse = TRUE))
  rna[["peaks"]] <- atac[["peaks"]]
  rna$capture <- c(rep(c("cap/A", "cap_B"), each = 60), NA)
  rna$db.scDblFinder_score <- c(rep(NA_real_, 120), 0.25)
  rna$db.scDblFinder_class <- c(rep(NA_character_, 120), "upstream")
  before <- SeuratObject::LayerData(rna, assay = "RNA", layer = "counts")
  set.seed(11)
  rng <- .Random.seed
  expect_warning(result <- RunDoubletCalling(rna,
    assay = "peaks", sample.by = "capture",
    db_rate = 0.1, seed = 31, nfeatures = 20, artificialDoublets = 60,
    nrounds = 5, iter = 1, verbose = FALSE
  ), "very low numbers")
  expect_identical(.Random.seed, rng)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = SeuratObject::LayerData(atac[["peaks"]], layer = "counts")))
  expect_warning(native <- withr::with_seed(31, scDblFinder::scDblFinder(sce,
    samples = rna$capture[1:120], dbr = 0.1, aggregateFeatures = TRUE,
    processing = "normFeatures", multiSampleMode = "split", returnType = "sce",
    nfeatures = 20, artificialDoublets = 60, nrounds = 5, iter = 1,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 31, progressbar = FALSE),
    verbose = FALSE
  )), "very low numbers")
  expect_equal(unname(result$db.scDblFinder_score[1:120]), native$scDblFinder.score)
  expect_identical(unname(result$db.scDblFinder_class[1:120]), as.character(native$scDblFinder.class))
  expect_equal(unname(tail(result$db.scDblFinder_score, 1)), 0.25)
  expect_identical(unname(tail(result$db.scDblFinder_class, 1)), "upstream")
  expect_identical(SeuratObject::LayerData(result, assay = "RNA", layer = "counts"), before)
  expect_identical(colnames(result), colnames(rna))
  receipt <- SeuratObject::Misc(result[["peaks"]], "scop_atac_doublets")
  expect_identical(unname(receipt$captures), unname(rna$capture[1:120]))
  expect_identical(receipt$parameters$multiSampleMode, "split")
  expect_true(receipt$parameters$aggregateFeatures)
  expect_error(RunscDblFinder(rna, assay = "peaks", verbose = FALSE), "sample.by")
  expect_error(RunscDblFinder(rna, assay = "peaks", sample.by = "capture", aggregateFeatures = FALSE), "fixed")
  rna$capture[1] <- NA
  expect_error(RunscDblFinder(rna, assay = "peaks", sample.by = "capture"), "capture/library ID")
})
