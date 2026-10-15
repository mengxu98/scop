make_atac_qc_object <- function() {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2, 3), j = c(1, 1, 1, 2, 3, 4),
    x = c(2, 3, 5, 4, 6, 8), dims = c(3, 4),
    dimnames = list(c("p1", "p2", "p3"), c("S1_A", "S1_B", "S2_A", "S2_B"))
  )
  peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1", "chr2"), start = c(1, 51, 201), end = c(100, 150, 250)
  ))
  object <- suppressWarnings(ReadATACData(counts, ranges = peaks, genome = "hg38", verbose = FALSE))
  object$TSS.enrichment <- 5
  object
}

add_atac_qc_fragments <- function(object, directory) {
  paths <- file.path(directory, c("s1.tsv.gz", "s2.tsv.gz"))
  records <- list(
    c("chr1\t0\t50\tA\t5", "chr1\t60\t100\tA\t2", "chr1\t100\t300\tB\t3"),
    c("chr1\t0\t50\tA\t7", "chr1\t100\t150\tB\t4", "chr1\t200\t400\tB\t2")
  )
  maps <- list(c(S1_B = "B", S1_A = "A"), c(S2_B = "B", S2_A = "A"))
  fragments <- lapply(seq_along(paths), function(i) {
    plain <- sub(".gz$", "", paths[[i]])
    writeLines(records[[i]], plain)
    Rsamtools::bgzip(plain, dest = paths[[i]], overwrite = TRUE)
    Rsamtools::indexTabix(paths[[i]], format = "bed")
    Signac::CreateFragmentObject(paths[[i]], cells = maps[[i]], tolerance = 0, verbose = FALSE)
  })
  Signac::Fragments(object[["peaks"]]) <- fragments
  object
}

test_that("QC counts mapped fragment records separately from read support", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_fragments_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- add_atac_qc_fragments(make_atac_qc_object(), directory)
  result <- RunATACQC(object, verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_identical(colnames(result), colnames(object))
  expect_equal(unname(result$frequency_count), c(2, 1, 1, 2))
  expect_equal(unname(result$reads_count), c(7, 3, 7, 6))
  expect_equal(unname(result$total.fragments), c(2, 1, 1, 2))
  expect_equal(unname(result$nucleosome_signal), c(0, NA, 0, 1))
  expect_false("pct_reads_in_peaks" %in% colnames(result[[]]))
  object$peak_region_fragments <- 1
  object$total.fragments <- 10
  expect_error(RunATACQC(object,
    frip.columns = c(in_peaks = "peak_region_fragments", total = "total.fragments"),
    verbose = FALSE
  ), "overwritten by fragment counting")
  expect_error(RunATACQC(object, max_nucleosome_signal = 4, verbose = FALSE), "missing or invalid.*nucleosome_signal")
  bad <- object
  bad[["peaks"]]@fragments[[2]]@path <- bad[["peaks"]]@fragments[[1]]@path
  expect_error(RunATACQC(bad, verbose = FALSE), "different files")
  subsetted <- subset(object, cells = c("S1_A", "S2_B"))
  expect_equal(unname(RunATACQC(subsetted, verbose = FALSE)$total.fragments), c(2, 2))
})

test_that("FRiP uses explicit matched columns without guessing old peak counts", {
  object <- make_atac_qc_object()
  object$peak_region_fragments <- c(8, 1, 0, 4)
  object$passed_filters <- c(10, 2, 0, 5)
  result <- RunATACQC(object, verbose = FALSE)
  expect_false("pct_reads_in_peaks" %in% colnames(result[[]]))
  expect_error(RunATACQC(object, min_pct_reads_in_peaks = 0.2, verbose = FALSE), "unavailable metric.*pct_reads")
  columns <- c(in_peaks = "peak_region_fragments", total = "passed_filters")
  result <- RunATACQC(object, frip.columns = columns, verbose = FALSE)
  expect_equal(unname(result$pct_reads_in_peaks), c(0.8, 0.5, NA, 0.8))
  expect_error(RunATACQC(object,
    frip.columns = columns, min_pct_reads_in_peaks = 0.2,
    verbose = FALSE
  ), "missing or invalid.*pct_reads")
  object$peak_region_fragments[3] <- 1
  expect_error(RunATACQC(object, frip.columns = columns, verbose = FALSE), "in_peaks <= total")
  expect_error(RunATACQC(object, frip.columns = unname(columns), verbose = FALSE), "in_peaks and total")
  object$pct_reads_in_peaks <- 80
  expect_error(RunATACQC(object, verbose = FALSE), "fractions in 0--1")
})

test_that("blacklist ratios count each overlapping peak once and keep Seurat objects", {
  object <- make_atac_qc_object()
  blacklist <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = rep("chr1", 3), start = c(10, 20, 70), end = c(11, 21, 71)
  ))
  result <- RunATACQC(object, blacklist = blacklist, verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_equal(unname(result$blacklist_ratio), c(0.5, 1, 1, 0))
  filtered <- RunATACQC(object, blacklist = blacklist, max_blacklist_ratio = 0.6, verbose = FALSE)
  expect_identical(colnames(filtered), c("S1_A", "S2_B"))
  expect_true(CheckATACData(filtered)$valid)
  wrong_genome <- blacklist
  GenomeInfoDb::genome(wrong_genome) <- "hg19"
  expect_error(RunATACQC(object, blacklist = wrong_genome, verbose = FALSE), "genome labels must agree")
  empty <- blacklist[0]
  expect_equal(unname(RunATACQC(object, blacklist = empty, verbose = FALSE)$blacklist_ratio), rep(0, 4))
  object[["peaks"]]@ranges <- object[["peaks"]]@ranges[c(3, 1, 2)]
  expect_error(RunATACQC(object, blacklist = blacklist, verbose = FALSE), "not aligned")
})

test_that("QC never silently bypasses an unavailable or invalid requested filter", {
  object <- make_atac_qc_object()
  object$TSS.enrichment <- NULL
  for (threshold in c("min_TSS_enrichment", "max_nucleosome_signal", "max_blacklist_ratio")) {
    args <- list(object = object, verbose = FALSE)
    args[[threshold]] <- 0.5
    expect_error(do.call(RunATACQC, args), "unavailable metric")
  }
  object$TSS.enrichment <- c(5, 4, NA, 6)
  expect_error(RunATACQC(object, min_TSS_enrichment = 2, verbose = FALSE), "missing or invalid.*TSS")
  expect_error(RunATACQC(object,
    tss.positions = Signac::granges(object[["peaks"]]),
    verbose = FALSE
  ), "requires available fragment counts")
  object$pct_reads_in_peaks <- c(0.8, 0.1, 0.6, 0.7)
  expect_error(RunATACQC(object, min_pct_reads_in_peaks = 1, verbose = FALSE), "No cells satisfy")
  expect_error(RunATACQC(object, min_pct_reads_in_peaks = 20, verbose = FALSE), "fraction thresholds use 0--1")
  expect_error(RunATACQC(object, min_TSS_enrichment = c(1, 2), verbose = FALSE), "Invalid threshold")
  filtered <- RunATACQC(object, min_pct_reads_in_peaks = 0.5, verbose = FALSE)
  expect_identical(colnames(filtered), c("S1_A", "S2_A", "S2_B"))
})

test_that("explicit TSS positions use native fast enrichment and replace imported scores", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_tss_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- file.path(directory, "tss.tsv")
  writeLines(c(
    "chr1\t1040\t1060\tA\t3", "chr1\t1040\t1060\tB\t2",
    "chr1\t1970\t2010\tA\t7", "chr1\t1980\t2020\tA\t5", "chr1\t1980\t2020\tB\t2",
    "chr1\t2940\t2960\tA\t4", "chr1\t2940\t2960\tB\t6"
  ), plain)
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  mapping <- c(sample_B = "B", sample_A = "A")
  peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(1, 1501), end = c(1500, 3500)
  ))
  fragment <- Signac::CreateFragmentObject(path, cells = mapping, verbose = FALSE)
  counts <- Signac::FeatureMatrix(fragment, features = peaks, cells = names(mapping), verbose = FALSE)
  object <- ReadATACData(counts,
    ranges = peaks, fragments = path,
    fragment.cells = mapping, genome = "hg38", verbose = FALSE
  )
  object$TSS.enrichment <- 999
  tss <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = 2001, end = 2001, strand = "+"
  ))
  native <- suppressWarnings(Signac::TSSEnrichment(object,
    assay = "peaks",
    tss.positions = tss, fast = TRUE, verbose = FALSE
  ))
  actual <- RunATACQC(object, assay = "peaks", tss.positions = tss, fast = TRUE, verbose = FALSE)
  wrong_tss <- tss
  GenomeInfoDb::genome(wrong_tss) <- "hg19"
  expect_error(RunATACQC(object, tss.positions = wrong_tss, verbose = FALSE), "TSS positions and assay genome labels")
  expect_equal(actual$TSS.enrichment, native$TSS.enrichment)
  expect_equal(actual$TSS.percentile, native$TSS.percentile)
  expect_true(all(is.finite(actual$TSS.enrichment)))
  expect_true(all(actual$TSS.enrichment != 999))
  expect_true(all(RunATACQC(object, verbose = FALSE)$TSS.enrichment == 999))
  rna <- Matrix::Matrix(matrix(1, 2L, 3L,
    dimnames = list(c("g1", "g2"), c(colnames(object), "rna_only"))
  ), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  partial$TSS.enrichment <- 999
  partial$TSS.percentile <- 0.42
  partial_result <- RunATACQC(partial,
    assay = "peaks", tss.positions = tss, fast = TRUE, verbose = FALSE
  )
  expect_identical(colnames(partial_result), colnames(partial))
  native_slow <- suppressWarnings(Signac::TSSEnrichment(object,
    assay = "peaks", tss.positions = tss, fast = FALSE, verbose = FALSE
  ))
  partial_slow <- RunATACQC(partial,
    assay = "peaks", tss.positions = tss, fast = FALSE, verbose = FALSE
  )
  expect_equal(partial_slow$TSS.enrichment[colnames(object)], native_slow$TSS.enrichment)
  expect_equal(partial_slow$TSS.percentile[colnames(object)], native_slow$TSS.percentile)
  expect_equal(partial_slow[["peaks"]]@positionEnrichment$TSS, native_slow[["peaks"]]@positionEnrichment$TSS)
  expect_equal(partial_result$TSS.enrichment[colnames(object)], native$TSS.enrichment)
  expect_equal(partial_result$TSS.percentile[colnames(object)], native$TSS.percentile)
  expect_equal(unname(partial_result$TSS.enrichment["rna_only"]), 999)
  expect_equal(unname(partial_result$TSS.percentile["rna_only"]), 0.42)
  expect_identical(SeuratObject::GetAssayData(partial_result, assay = "RNA", layer = "counts"), rna)
  filtered_partial <- RunATACQC(partial,
    assay = "peaks", tss.positions = tss, fast = TRUE,
    min_TSS_enrichment = 0.15, verbose = FALSE
  )
  expect_identical(colnames(filtered_partial), c("sample_A", "rna_only"))
  expect_error(RunATACQC(partial,
    assay = "peaks", tss.positions = tss, fast = TRUE,
    min_TSS_enrichment = 1000, verbose = FALSE
  ), "No cells satisfy")
  partial$source_numerator <- c(1, 1, NA)
  partial$source_denominator <- c(2, 4, NA)
  partial$frip <- partial$pct_reads_in_peaks <- 0.9
  source_result <- RunATACQC(partial,
    assay = "peaks",
    frip.columns = c(in_peaks = "source_numerator", total = "source_denominator"),
    min_pct_reads_in_peaks = 0.4, verbose = FALSE
  )
  expect_identical(colnames(source_result), c("sample_B", "rna_only"))
  expect_equal(unname(source_result$frip), c(0.5, 0.9))
  expect_identical(source_result$frip, source_result$pct_reads_in_peaks)
  second <- file.path(directory, "second.tsv.gz")
  expect_true(file.copy(path, second))
  expect_true(file.copy(paste0(path, ".tbi"), paste0(second, ".tbi")))
  paired_counts <- cbind(counts, counts)
  colnames(paired_counts) <- c(colnames(counts), "extra_B", "extra_A")
  paired <- ReadATACData(paired_counts,
    ranges = peaks, genome = "hg38",
    fragments = c(path, second),
    fragment.cells = list(mapping, c(extra_B = "B", extra_A = "A")), verbose = FALSE
  )
  capture <- subset(paired, cells = colnames(object))
  unlink(c(second, paste0(second, ".tbi")))
  expect_true(CheckATACData(capture)$fragments_ready)
  resources <- Signac::Fragments(capture[["peaks"]])
  for (fast_mode in c(TRUE, FALSE)) {
    result <- RunATACQC(capture,
      assay = "peaks", tss.positions = tss,
      fast = fast_mode, verbose = FALSE
    )
    reference <- if (fast_mode) native else native_slow
    expect_equal(result$TSS.enrichment, reference$TSS.enrichment)
    expect_equal(result$TSS.percentile, reference$TSS.percentile)
    expect_identical(Signac::Fragments(result[["peaks"]]), resources)
  }
  unlink(paste0(path, ".tbi"))
  expect_error(RunATACQC(capture, assay = "peaks", tss.positions = tss, verbose = FALSE))
})

test_that("the default matrix-only call retains cells without requiring unavailable FRiP", {
  object <- make_atac_qc_object()
  expect_error(RunATACQC(object, frip = NA, verbose = FALSE), "[Ff]rip must be TRUE or FALSE")
  expect_error(RunATACQC(object, frip = 1, verbose = FALSE), "[Ff]rip must be TRUE or FALSE")
  expect_error(RunATACQC(object, frip = c(TRUE, FALSE), verbose = FALSE), "[Ff]rip must be TRUE or FALSE")
  object$TSS.enrichment <- NULL
  result <- RunATACQC(object, verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_identical(colnames(result), colnames(object))
  expect_identical(
    SeuratObject::GetAssayData(result, layer = "counts"),
    SeuratObject::GetAssayData(object, layer = "counts")
  )
  expect_false(any(c("pct_reads_in_peaks", "TSS.enrichment", "nucleosome_signal") %in% colnames(result[[]])))
  object$pct_reads_in_peaks <- c(0.2, 0.5, NA, 1)
  preserved <- RunATACQC(object, verbose = FALSE)
  expect_equal(preserved$pct_reads_in_peaks, object$pct_reads_in_peaks)
  expect_error(RunATACQC(object, min_pct_reads_in_peaks = 0.1, verbose = FALSE), "missing or invalid.*pct_reads")
})

test_that("explicit fragment FRiP counts each record once across the peak union", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_frip_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- add_atac_qc_fragments(make_atac_qc_object(), directory)
  object$pct_reads_in_peaks <- 0.01 # stale imported values must be replaced
  for (chunk_size in c(0, -1, 1.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    expect_error(RunATACQC(object,
      frip = TRUE, frip.chunk.size = chunk_size, verbose = FALSE
    ), "[Ff]rip.chunk.size must be one positive integer")
  }
  result <- RunATACQC(object, frip = TRUE, frip.chunk.size = 1L, verbose = FALSE)
  expect_equal(unname(result$frip), c(1, 1, 1, 0.5))
  expect_identical(result$frip, result$pct_reads_in_peaks)
  receipt <- SeuratObject::Misc(result[["peaks"]], "scop_frip")
  expect_equal(unname(receipt$numerator), c(2, 1, 1, 1))
  expect_equal(unname(receipt$denominator), c(2, 1, 1, 2))
  expect_identical(receipt$peak_ids, rownames(object))
  expect_error(RunATACQC(object, frip = TRUE, frip.columns = c(in_peaks = "x", total = "y")), "Choose either")
  expect_error(RunATACQC(make_atac_qc_object(), frip = TRUE, verbose = FALSE), "requires available fragments")
})

test_that("fragment FRiP preserves duplicate records and avoids multi-peak inflation", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_frip_overlap_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c(
    "chr1\t0\t5\tA\t7", "chr1\t0\t5\tA\t7",
    "chr1\t0\t40\tA\t2", "chr1\t10\t20\tA\t1", "chr1\t30\t40\tB\t4",
    "chr2\t0\t5\tB\t3"
  ), plain)
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"))
  Rsamtools::indexTabix(path, format = "bed")
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = c(1L, 31L), end = c(5L, 40L)
  ))
  counts <- Matrix::Matrix(matrix(c(3L, 1L, 0L, 1L), 2L, 2L,
    dimnames = list(c("p1", "p2"), c("A", "B"))
  ), sparse = TRUE)
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", fragments = path, verbose = FALSE)
  object$TSS.enrichment <- 5
  result <- RunATACQC(object, frip = TRUE, frip.chunk.size = 2L, verbose = FALSE)
  # A: two duplicate rows + one bridging both peaks, over four records.
  # B: one peak hit plus one record on a contig outside the peak universe.
  expect_equal(unname(result$frip), c(3 / 4, 1 / 2))
  expect_equal(unname(result$total.fragments), c(4, 2))
  expect_equal(unname(result$reads_count), c(17, 7))
  expect_equal(unname(SeuratObject::Misc(result[["peaks"]], "scop_frip")$numerator), c(3, 1))
  expect_true(sum(counts[, "A"]) > SeuratObject::Misc(result[["peaks"]], "scop_frip")$numerator[["A"]])
  bedtools <- Sys.which("bedtools")
  if (nzchar(bedtools)) {
    bed <- file.path(directory, "peaks.bed")
    writeLines(c("chr1\t0\t5", "chr1\t30\t40"), bed)
    overlaps <- system2(bedtools, c("intersect", "-u", "-a", shQuote(plain), "-b", shQuote(bed)), stdout = TRUE)
    expect_null(attr(overlaps, "status"))
    barcodes <- vapply(strsplit(overlaps, "\t", fixed = TRUE), `[[`, "", 4L)
    expect_equal(as.integer(table(factor(barcodes, levels = c("A", "B")))), c(3L, 1L))
  }
})

test_that("fragment FRiP reads a native CSI-only indexed resource", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_frip_csi_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  tabix <- Sys.which("tabix")
  skip_if(!nzchar(tabix), "tabix is required to construct a CSI index")
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c(
    "chr1\t0\t5\tA\t1", "chr1\t20\t25\tB\t1",
    "chr1\t40\t45\tA\t1"
  ), plain)
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"))
  status <- system2(tabix, c("-C", "-p", "bed", shQuote(path)))
  expect_identical(status, 0L)
  expect_true(file.exists(paste0(path, ".csi")))
  expect_false(file.exists(paste0(path, ".tbi")))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = c(1L, 21L), end = c(10L, 30L)
  ))
  counts <- Matrix::Matrix(matrix(c(1, 0, 0, 1), 2L, 2L,
    dimnames = list(c("p1", "p2"), c("A", "B"))
  ), sparse = TRUE)
  object <- ReadATACData(counts,
    ranges = ranges, genome = "toy",
    fragments = path, verbose = FALSE
  )
  object$TSS.enrichment <- 5
  expect_true(CheckATACData(object)$fragments_ready)
  result <- RunATACQC(object, frip = TRUE, frip.chunk.size = 1L, verbose = FALSE)
  expect_equal(unname(result$frip), c(1 / 2, 1))
  expect_identical(result$frip, result$pct_reads_in_peaks)
  receipt <- SeuratObject::Misc(result[["peaks"]], "scop_frip")
  expect_equal(unname(receipt$numerator), c(1, 1))
  expect_equal(unname(receipt$denominator), c(2, 1))
  expect_identical(receipt$fragments$index, paste0(normalizePath(path), ".csi"))
  expect_false(file.exists(paste0(path, ".tbi")))
})


test_that("record overlap selection preserves closed boundaries and duplicate records", {
  peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(1L, 31L), end = c(5L, 40L)
  ))
  records <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = c(1L, 1L, 5L, 6L, 30L, 41L),
    end = c(5L, 5L, 31L, 30L, 31L, 50L)
  ))
  hit <- !is.na(GenomicRanges::findOverlaps(
    records, peaks,
    ignore.strand = TRUE, select = "first"
  ))
  expect_identical(hit, c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))
  expect_identical(hit, GenomicRanges::countOverlaps(records, peaks, ignore.strand = TRUE) > 0L)
  expect_identical(!is.na(GenomicRanges::findOverlaps(
    records, peaks[FALSE],
    ignore.strand = TRUE, select = "first"
  )), rep(FALSE, length(records)))
})


test_that("native capture subsetting ignores inactive fragment resources", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_capture_subset_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- add_atac_qc_fragments(make_atac_qc_object(), directory)
  object <- subset(object, cells = c("S1_A", "S1_B"))
  expect_length(Signac::GetFragmentData(Signac::Fragments(object[["peaks"]])[[2L]], "cells"), 0L)
  unlink(file.path(directory, "s2.tsv.gz"))
  counts <- RunATACQC(object, verbose = FALSE)
  expect_equal(unname(counts$total.fragments), c(2, 1))
  direct <- RunATACQC(object, frip = TRUE, verbose = FALSE)
  expect_equal(unname(direct$frip), c(1, 1))
  expect_equal(nrow(SeuratObject::Misc(direct[["peaks"]], "scop_frip")$fragments), 1L)
})

test_that("explicit source recalculation replaces direct FRiP and invalidates its receipt", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("qc_frip_refresh_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- add_atac_qc_fragments(make_atac_qc_object(), directory)
  direct <- RunATACQC(object, frip = TRUE, verbose = FALSE)
  receipt <- SeuratObject::Misc(direct[["peaks"]], "scop_frip")
  preserved <- RunATACQC(direct, verbose = FALSE)
  expect_identical(preserved$frip, direct$frip)
  expect_identical(SeuratObject::Misc(preserved[["peaks"]], "scop_frip"), receipt)
  direct$source_numerator <- 1
  direct$source_denominator <- 4
  refreshed <- RunATACQC(direct,
    frip.columns = c(in_peaks = "source_numerator", total = "source_denominator"), verbose = FALSE
  )
  expect_equal(unname(refreshed$frip), rep(0.25, 4L))
  expect_identical(refreshed$frip, refreshed$pct_reads_in_peaks)
  expect_null(SeuratObject::Misc(refreshed[["peaks"]], "scop_frip"))
})

test_that("RunCellQC validates ATAC thresholds while retaining mark-only semantics", {
  object <- make_atac_qc_object()
  expect_error(RunCellQC(object,
    assay = "peaks", qc_metrics = "atac",
    atac_args = list(min_pct_reads_in_peaks = 0.2, verbose = FALSE), verbose = FALSE
  ), "unavailable metric.*pct_reads")
  object$blacklist_ratio <- c(0.1, NA, 0.1, 0.1)
  expect_error(RunCellQC(object,
    assay = "peaks", qc_metrics = "atac",
    atac_args = list(max_blacklist_ratio = 0.2, verbose = FALSE), verbose = FALSE
  ), "missing or invalid.*blacklist_ratio")
  object$pct_reads_in_peaks <- c(0.1, 0.5, 0.3, 0.4)
  marked <- RunCellQC(object,
    assay = "peaks", qc_metrics = "atac",
    atac_args = list(min_pct_reads_in_peaks = 0.2, verbose = FALSE), verbose = FALSE
  )
  expect_identical(colnames(marked), colnames(object))
  expect_identical(as.character(marked$atac_qc), c("Fail", "Pass", "Pass", "Pass"))
  expect_identical(as.character(marked$CellQC), c("Fail", "Pass", "Pass", "Pass"))
  rna <- Matrix::Matrix(matrix(1, 2L, 5L,
    dimnames = list(c("g1", "g2"), c(colnames(object), "rna_only"))
  ), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  partial$pct_reads_in_peaks <- c(object$pct_reads_in_peaks, NA)
  partial_marked <- RunCellQC(partial,
    assay = "peaks", qc_metrics = "atac",
    atac_args = list(min_pct_reads_in_peaks = 0.2, verbose = FALSE), verbose = FALSE
  )
  expect_identical(colnames(partial_marked), colnames(partial))
  expect_identical(as.character(partial_marked$atac_qc), c("Fail", "Pass", "Pass", "Pass", "Pass"))
})
