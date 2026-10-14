make_chromvar_fixture <- function() {
  skip_if_not_installed("chromVAR")
  skip_if_not_installed("Biostrings")
  set.seed(321)
  n <- 200L
  counts <- Matrix::Matrix(matrix(stats::rpois(n * 8, lambda = 4), n, 8,
    dimnames = list(paste0("p", seq_len(n)), paste0("cell", 1:8))
  ), sparse = TRUE)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = (seq_len(n) - 1L) * 100L + 1L, end = seq_len(n) * 100L
  ))
  names(ranges) <- rownames(counts)
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  incidence <- Matrix::Matrix(matrix(stats::rbinom(n * 4L, 1, 0.2), n, 4L,
    dimnames = list(rownames(object), paste0("M", 1:4))
  ), sparse = TRUE)
  pwm <- stats::setNames(rep(list(matrix(1, 4, 3,
    dimnames = list(c("A", "C", "G", "T"), NULL)
  )), 4), colnames(incidence))
  Signac::Motifs(object[["peaks"]]) <- Signac::CreateMotifObject(data = incidence, pwm = pwm)
  gc <- sample(10:90, n, replace = TRUE)
  sequence <- paste0(vapply(gc, function(k) paste0(strrep("G", k), strrep("A", 100 - k)), ""), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = sequence))
  background <- matrix(sample.int(n, n * 20, replace = TRUE),
    nrow = n,
    dimnames = list(rownames(object), NULL)
  )
  list(object = object, genome = genome, incidence = incidence, background = background)
}

test_that("chromVAR preserves both native scores and deviations with explicit backgrounds", {
  fixture <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  object <- fixture$object
  original <- serialize(object, NULL)
  ranges <- Signac::granges(object[["peaks"]])
  names(ranges) <- rownames(object)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SeuratObject::LayerData(object[["peaks"]], layer = "counts")),
    rowRanges = ranges
  )
  se <- chromVAR::addGCBias(se, genome = fixture$genome)
  expected <- chromVAR::computeDeviations(se,
    annotations = fixture$incidence,
    background_peaks = fixture$background
  )
  actual <- RunChromVAR(object, fixture$genome,
    genome.id = "toy",
    background.peaks = fixture$background
  )
  expect_identical(serialize(object, NULL), original)
  expect_identical(SeuratObject::DefaultAssay(actual), "peaks")
  expect_equal(
    SeuratObject::LayerData(actual[["chromvar"]], layer = "data"),
    chromVAR::deviationScores(expected)
  )
  receipt <- SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")
  expect_equal(receipt$deviations, expected)
  expect_identical(receipt$background_peaks, fixture$background)
  expect_equal(receipt$gc_bias, stats::setNames(SummarizedExperiment::rowRanges(se)$bias, rownames(object)))
  expect_identical(receipt$background_source, "explicit")
  expect_identical(receipt$nonfinite_scores, 0L)
  expect_error(RunChromVAR(actual, fixture$genome, genome.id = "toy"), "unused name")
})

test_that("sampled chromVAR backgrounds match the native seeded backend", {
  fixture <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  object <- fixture$object
  ranges <- Signac::granges(object[["peaks"]])
  names(ranges) <- rownames(object)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SeuratObject::LayerData(object[["peaks"]], layer = "counts")),
    rowRanges = ranges
  )
  se <- chromVAR::addGCBias(se, genome = fixture$genome)
  expected <- withr::with_seed(19, chromVAR::getBackgroundPeaks(se, niterations = 20, w = 0.2, bs = 10))
  rownames(expected) <- rownames(object)
  set.seed(88)
  rng <- .Random.seed
  actual <- RunChromVAR(object, fixture$genome,
    genome.id = "toy",
    niterations = 20, w = 0.2, bs = 10, seed = 19
  )
  expect_identical(.Random.seed, rng)
  receipt <- SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")
  expect_identical(receipt$background_peaks, expected)
  expect_identical(receipt$background_source, "chromVAR::getBackgroundPeaks")
})

test_that("chromVAR rejects stale indices, unknown references and zero-depth input", {
  f <- make_chromvar_fixture()
  expect_error(RunChromVAR(f$object, "BSgenome.fake", genome.id = "toy"), "loaded BSgenome")
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "other"), "reference label")
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", seed = NA), "[Ss]eed")
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", niterations = 1), "integers")
  bad <- f$background[rev(seq_len(nrow(f$background))), , drop = FALSE]
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = bad), "exact[[:space:]]+source-peak")
  bad <- f$background
  bad[1, 1] <- 0
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = bad), "valid 1-based")
  object <- f$object
  counts <- SeuratObject::LayerData(object[["peaks"]], layer = "counts")
  counts[1, ] <- 0
  object[["peaks"]] <- SeuratObject::SetAssayData(object[["peaks"]], layer = "counts", new.data = counts)
  expect_error(RunChromVAR(object, f$genome, genome.id = "toy"), "positive total counts")
})

test_that("chromVAR alignment checks compare cached ID values and order", {
  skip_if_not_installed("fastmatch")
  fixture <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  cached <- paste0(as.vector(rownames(fixture$object)), "")
  fastmatch::fmatch(cached[1], cached)
  expect_false(is.null(attr(cached, ".match.hash")))
  fixture$object[["peaks"]]@motifs@data@Dimnames[[1L]] <- cached
  rownames(fixture$background) <- cached
  expect_false(identical(rownames(fixture$background), rownames(fixture$object)))
  actual <- RunChromVAR(fixture$object, fixture$genome,
    genome.id = "toy",
    background.peaks = fixture$background
  )
  expect_identical(
    as.vector(colnames(actual[["chromvar"]])),
    as.vector(colnames(fixture$object))
  )
  fixture$object[["peaks"]]@motifs@data@Dimnames[[1L]] <- rev(as.vector(cached))
  expect_error(RunChromVAR(fixture$object, fixture$genome,
    genome.id = "toy",
    background.peaks = fixture$background
  ), "aligned exactly")
})

test_that("chromVAR rejects loaded reference labels and lengths that disagree", {
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  f <- make_chromvar_fixture()
  GenomeInfoDb::genome(f$genome) <- "other"
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background), "reference.*labels.*agree")
  GenomeInfoDb::genome(f$genome) <- "toy"
  GenomeInfoDb::seqlengths(f$object[["peaks"]]) <- c(chr1 = 21000L)
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background), "lengths.*agree")
})

test_that("native pattern sparse motif incidence is supported", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  f$object[["peaks"]]@motifs@data <- as(f$incidence, "nMatrix")
  actual <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background)
  ranges <- Signac::granges(f$object[["peaks"]])
  names(ranges) <- rownames(f$object[["peaks"]])
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SeuratObject::LayerData(f$object[["peaks"]], layer = "counts")),
    rowRanges = ranges
  )
  se <- chromVAR::addGCBias(se, genome = f$genome)
  native <- chromVAR::computeDeviations(se, annotations = f$object[["peaks"]]@motifs@data, background_peaks = f$background)
  expect_equal(SeuratObject::LayerData(actual[["chromvar"]], layer = "data"), chromVAR::deviationScores(native))
})

test_that("explicit stored sparse zeros are not interpreted as motif matches", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  dirty <- f$incidence
  dirty@x[seq_len(10L)] <- 0
  f$object[["peaks"]]@motifs@data <- dirty
  before <- serialize(f$object, NULL)
  actual <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background)
  ranges <- Signac::granges(f$object[["peaks"]])
  names(ranges) <- rownames(f$object[["peaks"]])
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SeuratObject::LayerData(f$object[["peaks"]], layer = "counts")),
    rowRanges = ranges
  )
  se <- chromVAR::addGCBias(se, genome = f$genome)
  native <- chromVAR::computeDeviations(se, annotations = Matrix::drop0(dirty), background_peaks = f$background)
  expect_equal(SeuratObject::LayerData(actual[["chromvar"]], layer = "data"), chromVAR::deviationScores(native))
  expect_identical(serialize(f$object, NULL), before)
  f$object[["peaks"]]@motifs@data <- as(dirty, "lMatrix")
  logical_result <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background)
  expect_equal(SeuratObject::LayerData(logical_result[["chromvar"]], layer = "data"), chromVAR::deviationScores(native))
})

test_that("native assay-key collisions preserve caller RNG state", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  SeuratObject::Key(f$object[["peaks"]]) <- "chromvar_"
  set.seed(902L)
  before <- .Random.seed
  expect_warning(actual <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background), "[Kk]ey")
  expect_true(identical(.Random.seed, before))
  expect_true("chromvar" %in% names(actual))
})

test_that("partial and renamed assays preserve source analyses and metadata", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  cells <- c("cell8", "rna_only", paste0("cell", 1:7))
  counts <- Matrix::Matrix(matrix(seq_len(18), 2, dimnames = list(c("g1", "g2"), cells)), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts, assay = "RNA")
  object[["peaks"]] <- f$object[["peaks"]]
  object[["tiny"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(18), 9, dimnames = list(cells, c("tiny_1", "tiny_2"))), key = "tiny_", assay = "RNA"
  )
  object <- subset(object, cells = c("cell8", "rna_only", "cell3", "cell1"))
  expect_warning(object[["peaks"]] <- subset(object[["peaks"]], features = rownames(object[["peaks"]])[seq(2, 200, 2)]), "Different cells and/or features")
  object <- SeuratObject::RenameCells(object, add.cell.id = "renamed")
  ids <- rownames(object[["peaks"]])
  background <- withr::with_seed(13L, matrix(sample.int(length(ids), length(ids) * 20L, replace = TRUE),
    length(ids),
    dimnames = list(ids, NULL)
  ))
  before <- serialize(object, NULL)
  options_before <- options()
  actual <- RunChromVAR(object, f$genome, genome.id = "toy", assay = "peaks", background.peaks = background)
  expect_identical(serialize(object, NULL), before)
  expect_identical(options(), options_before)
  expect_identical(actual[["RNA"]], object[["RNA"]])
  expect_identical(actual[["peaks"]], object[["peaks"]])
  expect_identical(actual[["tiny"]], object[["tiny"]])
  expect_identical(actual[[]], object[[]])
  expect_identical(SeuratObject::DefaultAssay(actual), "RNA")
  expect_identical(colnames(actual), colnames(object))
  expect_identical(as.vector(colnames(actual[["chromvar"]])), as.vector(colnames(object[["peaks"]])))
  expect_false(any(c("nCount_chromvar", "nFeature_chromvar") %in% colnames(actual[[]])))
  receipt <- SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")
  expect_identical(as.vector(receipt$peaks), as.vector(ids))
  expect_identical(receipt$background_peaks, background)
})

test_that("native degenerate deviations are retained without zero imputation", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  n <- nrow(f$object[["peaks"]])
  background <- matrix(rep(seq_len(n), 5L), n, dimnames = list(rownames(f$object[["peaks"]]), NULL))
  actual <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = background)
  receipt <- SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")
  scores <- SeuratObject::LayerData(actual[["chromvar"]], layer = "data")
  expect_equal(scores, chromVAR::deviationScores(receipt$deviations))
  expect_true(any(!is.finite(scores)))
  expect_equal(receipt$nonfinite_scores, sum(!is.finite(scores)))
})

test_that("single-cell and zero-cell-depth inputs fail before creating an assay", {
  f <- make_chromvar_fixture()
  one <- subset(f$object, cells = colnames(f$object)[1])
  expect_error(RunChromVAR(one, f$genome, genome.id = "toy", background.peaks = f$background), "at least two")
  counts <- SeuratObject::LayerData(f$object[["peaks"]], layer = "counts")
  counts[, 1] <- 0
  f$object[["peaks"]] <- SeuratObject::SetAssayData(f$object[["peaks"]], layer = "counts", new.data = counts)
  expect_error(RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background), "positive total counts")
})

test_that("background sampling preserves an absent caller random seed", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  actual <- RunChromVAR(f$object, f$genome, genome.id = "toy", niterations = 5, w = 0.2, bs = 10)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_equal(ncol(SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")$background_peaks), 5)
})

test_that("explicit FASTA bias agrees with the same small DNA reference", {
  f <- make_chromvar_fixture()
  old <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old, default = TRUE), add = TRUE)
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  directory <- tempfile("chromvar_fasta_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- file.path(directory, "reference.fa")
  Biostrings::writeXStringSet(f$genome, path)
  Rsamtools::indexFa(path)
  actual <- RunChromVAR(f$object, Rsamtools::FaFile(path), genome.id = "toy", background.peaks = f$background)
  native <- RunChromVAR(f$object, f$genome, genome.id = "toy", background.peaks = f$background)
  expect_equal(SeuratObject::LayerData(actual[["chromvar"]], layer = "data"), SeuratObject::LayerData(native[["chromvar"]], layer = "data"))
  expect_equal(SeuratObject::Misc(actual[["chromvar"]], "scop_chromvar")$gc_bias, SeuratObject::Misc(native[["chromvar"]], "scop_chromvar")$gc_bias)
})

test_that("degenerate matching retains the native background diagnostic", {
  f <- make_chromvar_fixture()
  f$genome <- Biostrings::DNAStringSet(c(chr1 = strrep("A", 20000L)))
  ranges <- Signac::granges(f$object[["peaks"]])
  names(ranges) <- rownames(f$object[["peaks"]])
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SeuratObject::LayerData(f$object[["peaks"]], layer = "counts")), rowRanges = ranges
  )
  se <- chromVAR::addGCBias(se, genome = f$genome)
  native <- tryCatch(withr::with_seed(1L, chromVAR::getBackgroundPeaks(se, niterations = 5, w = 0.2, bs = 10)), error = identity)
  expect_s3_class(native, "error")
  before <- serialize(f$object, NULL)
  actual <- tryCatch(RunChromVAR(f$object, f$genome, genome.id = "toy", niterations = 5, w = 0.2, bs = 10), error = identity)
  expect_s3_class(actual, "error")
  expect_identical(conditionMessage(actual), conditionMessage(native))
  expect_identical(serialize(f$object, NULL), before)
})
