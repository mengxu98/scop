make_atac_recount_fixture <- function(directory) {
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1", "chr2"), start = c(1, 101, 1), end = c(50, 150, 100)
  ))
  GenomeInfoDb::seqlevels(ranges) <- c("chr1", "chr2", "chr3")
  cells <- c("S1_A", "S1_B", "S2_A", "S2_B")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1), j = 1:4, x = 1:4,
    dims = c(3, 4), dimnames = list(c("p1", "p2", "p3"), cells)
  )
  paths <- file.path(directory, c("s1.tsv.gz", "s2.tsv.gz"))
  records <- list(
    c(
      "chr1\t10\t20\tA\t5", "chr1\t40\t120\tA\t2",
      "chr1\t130\t140\tB\t3", "chr1\t4100\t4120\tA\t4", "chr2\t10\t30\tB\t4"
    ),
    c("chr1\t15\t25\tA\t8", "chr1\t110\t120\tB\t1", "chr1\t4110\t4120\tB\t2")
  )
  for (i in seq_along(paths)) {
    plain <- sub(".gz$", "", paths[[i]])
    writeLines(records[[i]], plain)
    Rsamtools::bgzip(plain, dest = paths[[i]], overwrite = TRUE)
    Rsamtools::indexTabix(paths[[i]], format = "bed")
  }
  object <- ReadATACData(counts,
    ranges = ranges, fragments = paths,
    fragment.cells = list(c(S1_B = "B", S1_A = "A"), c(S2_B = "B", S2_A = "A")),
    genome = "hg38", verbose = FALSE
  )
  object$sample <- c("S1", "S1", "S2", "S2")
  object
}

recount_test_ranges <- function() {
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr2", "chr1", "chr1", "chr3"),
    start = c(1, 101, 1, 1), end = c(100, 150, 50, 50)
  ))
  names(ranges) <- c("second_chrom", "right_peak", "left_peak", "absent_chrom")
  ranges
}

test_that("recount preserves shared ranges, requested cells and native record counts", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("atac_recount_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  ranges <- recount_test_ranges()
  cells <- c("S2_B", "S1_A", "S2_A", "S1_B")
  original <- serialize(object, NULL)
  result <- RecountATAC(object, ranges, cells = cells, process_n = 1, verbose = FALSE)
  expect_s4_class(result$counts, "dgCMatrix")
  expect_identical(colnames(result$counts), cells)
  expect_identical(rownames(result$counts), c("chr2-1-100", "chr1-101-150", "chr1-1-50", "chr3-1-50"))
  expect_identical(names(result$ranges), rownames(result$counts))
  expect_identical(result$receipt$input_peak_ids, names(ranges))
  expect_true(CheckATACData(ReadATACData(result$counts, ranges = result$ranges, verbose = FALSE))$valid)
  expect_identical(result$receipt$cells, cells)
  expect_equal(unname(as.matrix(result$counts)), matrix(c(
    0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 1, 0, 1, 1, 0, 0
  ), nrow = 4))
  native <- Signac::FeatureMatrix(Signac::Fragments(object[["peaks"]]), ranges,
    cells = cells, keep_all_features = TRUE, process_n = 100, verbose = FALSE
  )
  expect_equal(result$counts, native[rownames(result$counts), colnames(result$counts), drop = FALSE])
  expect_identical(original, serialize(object, NULL))
  expect_false(result$receipt$column5_weighted)
  expect_false(result$receipt$binary)
  expect_true("chr3" %in% result$receipt$declared_reference_contigs)
  expect_identical(result$receipt$missing_chromosomes[[1]], "chr3")
  expect_identical(result$receipt$missing_chromosomes[[2]], c("chr2", "chr3"))
  expect_identical(result$receipt$unmatched_cells, character())
  expect_match(result$receipt$fragments$fragment_md5, "^[a-f0-9]{32}$")
  expect_match(result$receipt$fragments$index_md5, "^[a-f0-9]{32}$")
  expect_identical(result$receipt$parameters$process_n, 1L)
  expect_false(result$receipt$parameters$keep_all_features)
})

test_that("native multi-chunk recount is invariant above the default 2000-feature boundary", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("atac_chunks_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = seq(1, by = 2, length.out = 2105), end = seq(1, by = 2, length.out = 2105)
  ))
  small <- RecountATAC(object, ranges, process_n = 1000, verbose = FALSE)
  default <- RecountATAC(object, ranges, verbose = FALSE)
  single <- Signac::FeatureMatrix(Signac::Fragments(object[["peaks"]]), ranges,
    cells = colnames(object), keep_all_features = FALSE, process_n = 3000, verbose = FALSE
  )
  expect_equal(small$counts, default$counts)
  expect_equal(default$counts, single[rownames(default$counts), colnames(default$counts), drop = FALSE])
  expect_gt(sum(default$counts[2001:2105, , drop = FALSE]), 0)
  expect_identical(default$receipt$parameters$process_n, 2000L)
})

test_that("recount distinguishes real zero overlaps from empty or invalid requests", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("atac_zero_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = "chr1", start = 9000, end = 9001))
  zero <- RecountATAC(object, ranges, cells = "S1_A", verbose = FALSE)
  expect_s4_class(zero$counts, "dgCMatrix")
  expect_identical(dim(zero$counts), c(1L, 1L))
  expect_equal(sum(zero$counts), 0)
  expect_identical(zero$receipt$zero_count_cells, "S1_A")
  zero_ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(9000, 9100), end = c(9001, 9101)
  ))
  zero_assay <- RecountATAC(object, zero_ranges,
    cells = c("S1_A", "S1_B"),
    new.assay = "zero", verbose = FALSE
  )
  expect_equal(sum(SeuratObject::GetAssayData(zero_assay, assay = "zero", layer = "counts")), 0)
  expect_identical(dim(zero_assay[["zero"]]), c(2L, 2L))
  expect_error(RecountATAC(object, ranges, cells = character(), verbose = FALSE), "source-assay cell IDs")
  expect_error(RecountATAC(object, ranges[0], verbose = FALSE), "at least one interval")
  expect_error(RecountATAC(object, ranges, cells = "unknown", verbose = FALSE), "source-assay cell IDs")
  expect_error(RecountATAC(object, ranges, new.assay = "onepeak", verbose = FALSE), "two peaks and two cells")
  missing <- recount_test_ranges()[4]
  expect_error(RecountATAC(object, missing, verbose = FALSE), "No requested chromosomes")
  expect_error(RecountATAC(object, ranges, process_n = 0, verbose = FALSE), "positive integer")
  duplicated <- recount_test_ranges()[c(1, 1)]
  expect_error(RecountATAC(object, duplicated, verbose = FALSE), "Duplicate genomic")
  colliding <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr_a", "chr-a"),
    start = c(1, 1), end = c(10, 10)
  ))
  expect_error(RecountATAC(object, colliding, verbose = FALSE), "IDs collide")
  undeclared <- object
  GenomeInfoDb::seqlevels(undeclared[["peaks"]]@ranges) <- c("chr1", "chr2")
  undeclared[["peaks"]]@seqinfo <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), genome = "hg38")
  expect_error(RecountATAC(undeclared, recount_test_ranges(), verbose = FALSE), "not declared.*chr3")
  mixed_style <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "1"),
    start = c(1, 1), end = c(50, 50)
  ))
  expect_error(RecountATAC(object, mixed_style, verbose = FALSE), "not declared.*1")
  GenomeInfoDb::genome(ranges) <- "hg19"
  expect_error(RecountATAC(object, ranges, verbose = FALSE), "genome labels must agree")
  unknown <- object
  GenomeInfoDb::genome(unknown[["peaks"]]@ranges) <- NA_character_
  GenomeInfoDb::genome(unknown[["peaks"]]@seqinfo) <- NA_character_
  expect_error(RecountATAC(unknown, recount_test_ranges(), verbose = FALSE), "verified source-assay genome")
})

test_that("recount validates only selected files but never fabricates missing barcodes", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_{literal}_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  ranges <- recount_test_ranges()[2:3]
  bad <- object
  bad[["peaks"]]@fragments[[1]]@cells["S1_A"] <- "ABSENT"
  expect_error(RecountATAC(bad, ranges, cells = "S1_A", verbose = FALSE), "barcodes are absent")
  expect_s4_class(RecountATAC(bad, ranges, cells = "S1_B", verbose = FALSE)$counts, "dgCMatrix")
  bad <- object
  bad[["peaks"]]@fragments[[2]]@path <- bad[["peaks"]]@fragments[[1]]@path
  expect_error(RecountATAC(bad, ranges, verbose = FALSE), "different files")
  index <- file.path(directory, "s2.tsv.gz.tbi")
  backup <- readBin(index, what = "raw", n = file.info(index)$size)
  unlink(index)
  expect_error(RecountATAC(object, ranges, verbose = FALSE), "index is missing.*recount_\\{literal\\}_")
  writeBin(backup, index)
  unlink(file.path(directory, "s2.tsv.gz"))
  expect_error(RecountATAC(object, ranges, verbose = FALSE), "files are missing")
  expect_s4_class(RecountATAC(object, ranges, cells = c("S1_A", "S1_B"), verbose = FALSE)$counts, "dgCMatrix")
  con <- file(file.path(directory, "s1.tsv.gz"), "ab")
  writeBin(as.raw(0), con)
  close(con)
  expect_error(RecountATAC(object, ranges, cells = "S1_A", verbose = FALSE), "hash differs")
})

test_that("new-assay recount is atomic and preserves source analyses", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("atac_writeback_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  source_counts <- SeuratObject::GetAssayData(object, layer = "counts")
  object <- SeuratObject::SetAssayData(object, layer = "data", new.data = source_counts * 7)
  SeuratObject::VariableFeatures(object) <- rownames(object)[1:2]
  object[["peaks"]] <- SeuratObject::AddMetaData(object[["peaks"]],
    metadata = data.frame(old_peak_annotation = 1:3, row.names = rownames(object))
  )
  object[["oldlsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(8), ncol = 2, dimnames = list(colnames(object), c("OLD_1", "OLD_2"))),
    key = "OLD_", assay = "peaks"
  )
  motifs <- Signac::CreateMotifObject(data = Matrix::sparseMatrix(
    i = c(1, 3), j = c(1, 2), x = 1, dims = c(3, 2),
    dimnames = list(rownames(object), c("M1", "M2"))
  ))
  Signac::Motifs(object[["peaks"]]) <- motifs
  links <- Signac::granges(object[["peaks"]])[1]
  links$score <- 0.5
  Signac::Links(object[["peaks"]]) <- links
  before <- serialize(object, NULL)
  ranges <- recount_test_ranges()[1:3]
  for (name in c("peaks", "oldlsi", "invalid-name")) {
    expect_error(RecountATAC(object, ranges, new.assay = name, verbose = FALSE), "unused name")
    expect_identical(serialize(object, NULL), before)
  }
  collision <- object
  collision$nCount_new <- 9
  expect_error(RecountATAC(collision, ranges, new.assay = "new", verbose = FALSE), "QC metadata collisions")
  requested <- c("S2_B", "S1_A", "S1_B")
  result <- RecountATAC(object, ranges, cells = requested, new.assay = "recount", verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_identical(serialize(object, NULL), before)
  expect_identical(serialize(result[["peaks"]], NULL), serialize(object[["peaks"]], NULL))
  expect_identical(serialize(result[["oldlsi"]], NULL), serialize(object[["oldlsi"]], NULL))
  expect_equal(result[[]][, colnames(object[[]]), drop = FALSE], object[[]])
  expect_identical(SeuratObject::DefaultAssay(result), SeuratObject::DefaultAssay(object))
  expected_cells <- colnames(object)[colnames(object) %in% requested]
  expect_identical(as.character(colnames(result[["recount"]])), expected_cells)
  receipt <- SeuratObject::Misc(result[["recount"]], "scop_atac_recount")
  expect_identical(receipt$requested_cells, requested)
  expect_identical(receipt$cells, expected_cells)
  expect_length(SeuratObject::VariableFeatures(result[["recount"]]), 0)
  expect_null(Signac::Motifs(result[["recount"]]))
  expect_length(Signac::Links(result[["recount"]]), 0)
  expect_length(Signac::Annotation(result[["recount"]]), 0)
  expect_false("old_peak_annotation" %in% colnames(result[["recount"]][[]]))
  expect_equal(
    SeuratObject::GetAssayData(result, assay = "recount", layer = "data"),
    SeuratObject::GetAssayData(result, assay = "recount", layer = "counts")
  )
  expect_true(CheckATACData(result, assay = "recount")$valid)
  expect_true(CheckATACData(result, assay = "recount")$fragments_ready)
  path <- file.path(directory, "recount.rds")
  saveRDS(result, path)
  restored <- readRDS(path)
  expect_true(CheckATACData(restored, assay = "recount")$fragments_ready)
  expect_equal(
    unname(SeuratObject::GetAssayData(restored, assay = "recount", layer = "counts")),
    unname(SeuratObject::GetAssayData(result, assay = "recount", layer = "counts"))
  )
})

test_that("recount ignores inactive captures after native subset and rename", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_subset_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  object <- subset(object, cells = c("S1_A", "S1_B"))
  object <- SeuratObject::RenameCells(object, new.names = c("renamed_A", "renamed_B"))
  unlink(c(file.path(directory, "s2.tsv.gz"), file.path(directory, "s2.tsv.gz.tbi")))
  expect_true(CheckATACData(object)$fragments_ready)
  ranges <- recount_test_ranges()[2:3]
  result <- RecountATAC(object, ranges, verbose = FALSE)
  native <- Signac::FeatureMatrix(Signac::Fragments(object[["peaks"]])[[1]], ranges,
    cells = colnames(object), keep_all_features = FALSE, verbose = FALSE
  )
  expect_equal(result$counts, native[rownames(result$counts), colnames(result$counts), drop = FALSE])
  expect_equal(nrow(result$receipt$fragments), 1L)
})

test_that("recount rejects conflicting known lengths before assigning a reference", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_reference_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- c(chr1 = 10000L, chr2 = 10000L, chr3 = 10000L)
  ranges <- recount_test_ranges()[2:3]
  GenomeInfoDb::seqlengths(ranges) <- c(chr2 = 10000L, chr1 = 9999L, chr3 = 10000L)[GenomeInfoDb::seqlevels(ranges)]
  expect_error(RecountATAC(object, ranges, verbose = FALSE), "lengths conflict")
})

test_that("recount preserves known source lengths on output ranges and assays", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_lengths_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  lengths <- c(chr1 = 10000L, chr2 = 10000L, chr3 = 10000L)
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- lengths
  ranges <- recount_test_ranges()[2:3]
  result <- RecountATAC(object, ranges, verbose = FALSE)
  expect_identical(GenomeInfoDb::seqlengths(result$ranges), lengths[GenomeInfoDb::seqlevels(result$ranges)])
  expect_true(all(GenomeInfoDb::genome(result$ranges) == "hg38"))
  written <- RecountATAC(object, ranges, new.assay = "common", verbose = FALSE)
  expect_identical(GenomeInfoDb::seqlengths(written[["common"]]), lengths[GenomeInfoDb::seqlevels(written[["common"]])])
})

test_that("native FeatureMatrix counts record overlaps rather than read support or FRiP", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_units_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_recount_fixture(directory)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(1L, 60L, 1L), end = c(150L, 90L, 50L)
  ))
  result <- RecountATAC(object, ranges, cells = c("S1_A", "S1_B"), process_n = 1L, verbose = FALSE)
  native <- Signac::FeatureMatrix(Signac::Fragments(object[["peaks"]])[[1]], ranges,
    cells = c("S1_A", "S1_B"), keep_all_features = FALSE, process_n = 1L, verbose = FALSE
  )
  expect_equal(result$counts, native[rownames(result$counts), colnames(result$counts), drop = FALSE])
  expect_equal(unname(as.matrix(result$counts)), matrix(c(2, 1, 2, 1, 0, 0), nrow = 3))
  expect_gt(sum(result$counts[, "S1_A"]), 2)
  expect_false(result$receipt$column5_weighted)
  expect_false(result$receipt$binary)
})

test_that("native BED endpoints and duplicate records remain explicit", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_endpoints_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c("chr1\t10\t20\tA\t5", "chr1\t10\t20\tA\t7", "chr1\t15\t25\tB\t9"), plain)
  path <- paste0(plain, ".gz")
  Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  peaks <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
  counts <- Matrix::Matrix(matrix(c(99, 0, 0, 99), 2,
    dimnames = list(c("p1", "p2"), c("sample_A", "sample_B"))
  ), sparse = TRUE)
  object <- ReadATACData(counts,
    ranges = peaks, genome = "hg38", fragments = path,
    fragment.cells = c(sample_A = "A", sample_B = "B"), verbose = FALSE
  )
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(10L, 11L, 20L, 21L, 1L), end = c(10L, 11L, 20L, 21L, 50L)
  ))
  original <- serialize(object, NULL)
  result <- RecountATAC(object, ranges, process_n = 1L, verbose = FALSE)
  expect_identical(serialize(object, NULL), original)
  native <- Signac::FeatureMatrix(Signac::Fragments(object[["peaks"]]), ranges,
    cells = vapply(colnames(object), identity, character(1), USE.NAMES = FALSE), keep_all_features = FALSE, process_n = 1L, verbose = FALSE
  )
  expect_equal(result$counts, native[rownames(result$counts), colnames(result$counts), drop = FALSE])
  expect_equal(unname(as.matrix(result$counts)), matrix(c(0, 2, 2, 0, 2, 0, 0, 1, 1, 1), nrow = 5))
  expect_equal(unname(as.matrix(SeuratObject::GetAssayData(object, layer = "counts"))), unname(as.matrix(counts)))
})

test_that("partial-assay recount preserves RNA and explicitly selected cell orders", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("recount_partial_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  atac <- make_atac_recount_fixture(directory)
  rna_cells <- c("S2_B", "rna_only", "S1_A", "S2_A", "S1_B")
  rna <- Matrix::Matrix(matrix(seq_len(10), 2, dimnames = list(c("g1", "g2"), rna_cells)), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  object[["peaks"]] <- atac[["peaks"]]
  object[["tiny"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(10), 5, dimnames = list(rna_cells, c("tiny_1", "tiny_2"))), key = "tiny_", assay = "RNA"
  )
  original <- serialize(object, NULL)
  ranges <- recount_test_ranges()[c(3, 2)]
  cells <- c("S1_B", "S1_A", "S2_B")
  matrix_result <- RecountATAC(object, ranges, assay = "peaks", cells = cells, verbose = FALSE)
  written <- RecountATAC(object, ranges, assay = "peaks", cells = cells, new.assay = "common", verbose = FALSE)
  stored <- colnames(object)[colnames(object) %in% cells]
  expect_identical(colnames(matrix_result$counts), cells)
  expect_identical(as.character(colnames(written[["common"]])), as.character(stored))
  # Native matching may add a .match.hash cache to cell-name vectors.
  actual <- SeuratObject::GetAssayData(written, assay = "common", layer = "counts")
  expect_identical(rownames(actual), rownames(matrix_result$counts))
  expect_equal(unname(actual), unname(matrix_result$counts[, stored, drop = FALSE]))
  expect_identical(written[["peaks"]], object[["peaks"]])
  expect_identical(written[["RNA"]], object[["RNA"]])
  expect_identical(written[["tiny"]], object[["tiny"]])
  expect_identical(colnames(written), colnames(object))
  expect_identical(SeuratObject::DefaultAssay(written), "RNA")
  expect_identical(serialize(object, NULL), original)
  matrix_only <- ReadATACData(SeuratObject::GetAssayData(atac, layer = "counts"),
    ranges = Signac::granges(atac[["peaks"]]), genome = "hg38", verbose = FALSE
  )
  expect_error(RecountATAC(matrix_only, ranges, verbose = FALSE), "attached local indexed fragments")
})

test_that("recount uses native CSI-only fragments without substituting an index", {
  skip_if_not_installed("Rsamtools")
  tabix <- Sys.which("tabix")
  skip_if(!nzchar(tabix), "tabix is required to construct a CSI index")
  directory <- tempfile("recount_csi_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c("chr1\t10\t20\tA\t5", "chr1\t110\t120\tB\t3"), plain)
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"))
  expect_identical(system2(tabix, c("-C", "-p", "bed", shQuote(path))), 0L)
  peaks <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
  counts <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2), dimnames = list(c("p1", "p2"), c("A", "B")))
  object <- ReadATACData(counts, ranges = peaks, fragments = path, genome = "hg38", verbose = FALSE)
  result <- RecountATAC(object, peaks, verbose = FALSE)
  expect_equal(unname(as.matrix(result$counts)), diag(2))
  expect_identical(result$receipt$fragments$index, paste0(path, ".csi"))
  expect_false(file.exists(paste0(path, ".tbi")))
})
