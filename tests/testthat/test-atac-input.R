make_atac_input <- function() {
  counts <- Matrix::sparseMatrix(
    i = c(1L, 3L, 2L, 3L, 1L), j = c(1L, 1L, 2L, 2L, 3L),
    x = c(2, 1, 3, 1, 4), dims = c(3L, 3L),
    dimnames = list(c("peak1", "peak2", "peak3"), c("AA-1", "BB-1", "CC-1"))
  )
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1", "chr2"),
    start = c(1L, 101L, 201L), end = c(50L, 150L, 250L)
  ))
  names(ranges) <- rownames(counts)
  list(counts = counts, ranges = ranges)
}

test_that("ATAC imports preserve raw counts and explicitly align peaks, cells and samples", {
  x <- make_atac_input()
  meta <- data.frame(condition = c("stim", "ctrl", "ctrl"), row.names = c("CC-1", "AA-1", "BB-1"))
  samples <- c("CC-1" = "D2", "BB-1" = "D1", "AA-1" = "D1")
  object <- ReadATACData(x$counts,
    ranges = x$ranges[c(3, 1, 2)], metadata = meta,
    cell.names = paste0("capture1_", colnames(x$counts)), sample = samples, genome = "hg38"
  )
  expect_s4_class(object[["peaks"]], "ChromatinAssay")
  actual <- SeuratObject::GetAssayData(object, assay = "peaks", layer = "counts")
  expect_equal(unname(actual), unname(x$counts))
  expect_identical(rownames(actual), c("chr1-1-50", "chr1-101-150", "chr2-201-250"))
  expect_identical(unname(object$input_barcode), colnames(x$counts))
  expect_identical(unname(object$sample), c("D1", "D1", "D2"))
  expect_identical(unname(object$condition), c("ctrl", "ctrl", "stim"))
  expect_identical(object[["peaks"]][[]]$input_peak_id, rownames(x$counts))
  before <- serialize(object, NULL)
  report <- CheckATACData(object, assay = "peaks", genome = "hg38", sample.by = "sample")
  expect_true(report$valid)
  expect_false(report$fragments_ready)
  expect_equal(as.integer(report$samples), c(2L, 1L))
  expect_identical(before, serialize(object, NULL))
  expect_false(CheckATACData(object, genome = "hg19")$valid)
  zero_counts <- x$counts
  zero_counts[1, ] <- 0
  zero_counts[, 1] <- 0
  zeros <- ReadATACData(zero_counts, ranges = x$ranges)
  expect_identical(dim(zeros), dim(x$counts))
  expect_equal(sum(SeuratObject::GetAssayData(zeros, layer = "counts")), sum(zero_counts))
})

test_that("ATAC import rejects ambiguous IDs and non-raw counts", {
  x <- make_atac_input()
  for (value in c(-1, 0.5, Inf, NA_real_)) {
    bad <- x$counts
    bad@x[1] <- value
    expect_error(ReadATACData(bad, ranges = x$ranges), "non-negative integers")
  }
  expect_error(ReadATACData(as.matrix(x$counts), ranges = x$ranges), "sparse matrix")
  expect_error(ReadATACData(list(x$counts, x$counts), ranges = x$ranges), "sparse matrix")
  expect_error(ReadATACData(x$counts[, 1, drop = FALSE], ranges = x$ranges), "at least two cells")
  expect_error(ReadATACData(x$counts[1, , drop = FALSE], ranges = x$ranges[1]), "two peaks")
  bad <- x$counts
  colnames(bad)[2] <- colnames(bad)[1]
  expect_error(ReadATACData(bad, ranges = x$ranges), "Cell IDs")
  bad <- x$counts
  rownames(bad)[2] <- rownames(bad)[1]
  expect_error(ReadATACData(bad, ranges = x$ranges), "Peak IDs")
  bad_ranges <- x$ranges
  names(bad_ranges)[1] <- "other_peak"
  expect_error(ReadATACData(x$counts, ranges = bad_ranges), "Named ranges")
  expect_error(ReadATACData(x$counts, ranges = x$ranges[1:2]), "one interval")
  bad_ranges <- x$ranges[c(1, 1, 3)]
  expect_error(ReadATACData(x$counts, ranges = bad_ranges), "Duplicate genomic")
  expect_error(ReadATACData(x$counts, ranges = x$ranges, sample = c("S1", "S2", "S2")), "named by all")
  expect_error(ReadATACData(x$counts, ranges = x$ranges, metadata = data.frame(a = 1:3)), "matching input cell")
  expect_error(ReadATACData(x$counts, ranges = x$ranges, cell.names = c("A", "A", "B")), "unique non-empty")
  expect_error(ReadATACData(x$counts, ranges = x$ranges, fragment.cells = c(a = "AA-1")), "requires fragments")
  expect_error(ReadATACData(x$counts, ranges = x$ranges, peaks.file = "anything"), "only used with a Matrix Market path")
})

test_that("ATAC import handles reference labels and chromosome boundaries without downloads", {
  x <- make_atac_input()
  GenomeInfoDb::genome(x$ranges) <- "hg38"
  expect_error(ReadATACData(x$counts, ranges = x$ranges, genome = "hg19"), "Genome labels")
  annotation <- x$ranges
  GenomeInfoDb::genome(annotation) <- "hg19"
  expect_error(ReadATACData(x$counts, ranges = x$ranges, annotation = annotation), "Genome labels")
  bad <- x$ranges
  suppressWarnings(GenomeInfoDb::seqlengths(bad) <- c(chr1 = 100L, chr2 = 500L))
  expect_error(ReadATACData(x$counts, ranges = bad), "chromosome lengths")
  bad <- x$ranges
  GenomicRanges::start(bad)[1] <- 0L
  expect_error(ReadATACData(x$counts, ranges = bad), "positive 1-based")
  x$ranges <- unname(x$ranges)
  object <- ReadATACData(x$counts, ranges = x$ranges)
  expect_true(CheckATACData(object, genome = "hg38")$valid)
})

test_that("BED imports convert zero-based starts and retain all cells and peaks", {
  x <- make_atac_input()
  directory <- tempfile("atac_mtx_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  Matrix::writeMM(x$counts, file.path(directory, "matrix.mtx"))
  writeLines(c("chr1\t0\t50", "chr1\t100\t150", "chr2\t200\t250"), file.path(directory, "peaks.bed"))
  writeLines(colnames(x$counts), file.path(directory, "barcodes.tsv"))
  for (file in c("matrix.mtx", "peaks.bed", "barcodes.tsv")) {
    input <- readBin(file.path(directory, file), what = "raw", n = file.info(file.path(directory, file))$size)
    con <- gzfile(file.path(directory, paste0(file, ".gz")), "wb")
    writeBin(input, con)
    close(con)
  }
  object <- ReadATACData(file.path(directory, "matrix.mtx.gz"),
    peaks.file = file.path(directory, "peaks.bed.gz"), barcodes.file = file.path(directory, "barcodes.tsv.gz")
  )
  expect_identical(GenomicRanges::start(Signac::granges(object[["peaks"]])), c(1L, 101L, 201L))
  expect_equal(unname(SeuratObject::GetAssayData(object, layer = "counts")), unname(x$counts))
  expect_identical(object[["peaks"]][[]]$input_peak_id, c("chr1:0-50", "chr1:100-150", "chr2:200-250"))
  expect_identical(dim(object), c(3L, 3L))
  writeLines(c("chr1\t50\t50", "chr1\t100\t150", "chr2\t200\t250"), file.path(directory, "peaks.bed"))
  expect_error(ReadATACData(file.path(directory, "matrix.mtx"),
    peaks.file = file.path(directory, "peaks.bed"), barcodes.file = file.path(directory, "barcodes.tsv")
  ), "end greater")
})

test_that("H5 import uses the real Seurat reader and selects only the Peaks modality", {
  skip_if_not_installed("hdf5r")
  x <- make_atac_input()
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  matrix <- h5$create_group("matrix")
  features <- matrix$create_group("features")
  counts <- rbind(x$counts, Matrix::Matrix(c(5, 0, 8), nrow = 1L, sparse = TRUE))
  matrix[["data"]] <- counts@x
  matrix[["indices"]] <- counts@i
  matrix[["indptr"]] <- counts@p
  matrix[["shape"]] <- dim(counts)
  matrix[["barcodes"]] <- colnames(counts)
  features[["id"]] <- c("chr1:0-50", "chr1:100-150", "chr2:200-250", "ENSG1")
  features[["name"]] <- c("chr1:0-50", "chr1:100-150", "chr2:200-250", "GENE1")
  features[["feature_type"]] <- c("Peaks", "Peaks", "Peaks", "Gene Expression")
  h5$close_all()
  object <- suppressMessages(ReadATACData(path, genome = "hg38"))
  expect_true(CheckATACData(object, genome = "hg38")$valid)
  expect_identical(dim(object), c(3L, 3L))
  expect_equal(unname(SeuratObject::GetAssayData(object, layer = "counts")), unname(x$counts))
  expect_identical(GenomicRanges::start(Signac::granges(object[["peaks"]])), c(1L, 101L, 201L))
})

test_that("fragment mappings distinguish repeated barcodes across captures and survive reload", {
  skip_if_not_installed("Rsamtools")
  x <- make_atac_input()
  counts <- cbind(x$counts[, 1:2], x$counts[, 1:2])
  colnames(counts) <- c("S1_AA-1", "S1_BB-1", "S2_AA-1", "S2_BB-1")
  directory <- tempfile("atac_fragments_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  paths <- file.path(directory, c("s1.tsv.gz", "s2.tsv.gz"))
  for (path in paths) {
    plain <- sub(".gz$", "", path)
    writeLines(c("chr1\t0\t50\tAA-1\t1", "chr1\t100\t150\tBB-1\t2"), plain)
    Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
    Rsamtools::indexTabix(path, format = "bed")
  }
  maps <- list(
    c("S1_AA-1" = "AA-1", "S1_BB-1" = "BB-1"),
    c("S2_AA-1" = "AA-1", "S2_BB-1" = "BB-1")
  )
  expect_error(ReadATACData(counts,
    ranges = x$ranges,
    fragments = c(paths[1], file.path(directory, ".", "s1.tsv.gz")),
    fragment.cells = maps
  ), "different files")
  expect_error(ReadATACData(counts, ranges = x$ranges, fragments = paths), "explicit named")
  expect_error(ReadATACData(counts,
    ranges = x$ranges, fragments = paths,
    fragment.cells = list(maps[[1]], maps[[1]])
  ), "exactly once")
  absent <- maps
  absent[[2]][1] <- "NOT_IN_FILE-1"
  expect_error(ReadATACData(counts,
    ranges = x$ranges, fragments = paths,
    fragment.cells = absent, verbose = FALSE
  ), "[Cc]ells")
  renamed <- ReadATACData(x$counts[, 1:2],
    ranges = x$ranges,
    cell.names = names(maps[[1]]), fragments = paths[1], verbose = FALSE
  )
  renamed_map <- Signac::GetFragmentData(Signac::Fragments(renamed[["peaks"]])[[1]], "cells")
  expect_identical(unname(renamed_map), unname(maps[[1]]))
  expect_identical(as.character(names(renamed_map)), as.character(names(maps[[1]])))
  object <- ReadATACData(counts, ranges = x$ranges, fragments = paths, fragment.cells = maps, verbose = FALSE)
  report <- CheckATACData(object)
  expect_true(report$valid)
  expect_true(report$fragments_ready)
  expect_equal(report$fragments$n_cells, c(2L, 2L))
  expect_identical(Signac::GetFragmentData(Signac::Fragments(object[["peaks"]])[[2]], "cells"), maps[[2]])
  saveRDS(object, file.path(directory, "object.rds"))
  restored <- readRDS(file.path(directory, "object.rds"))
  expect_equal(CheckATACData(restored), report)
  bad <- restored
  bad[["peaks"]]@fragments[[1]]@cells[] <- "AA-1"
  expect_true(CheckATACData(bad)$valid)
  expect_false(CheckATACData(bad)$fragments_ready)
  capture <- subset(restored, cells = names(maps[[1]]))
  capture_report <- CheckATACData(capture)
  expect_true(capture_report$fragments_ready)
  expect_identical(capture_report$fragments$active, c(TRUE, FALSE))
  expect_equal(capture_report$fragments$n_cells, c(2L, 0L))
  incomplete <- capture
  incomplete[["peaks"]]@fragments[[1]]@cells <- maps[[1]][1L]
  expect_false(CheckATACData(incomplete)$fragments_ready)
  unlink(paths[2])
  inactive_missing <- CheckATACData(capture)
  expect_true(inactive_missing$fragments_ready)
  expect_false(inactive_missing$fragments$file_exists[2L])
  expect_identical(inactive_missing$warnings, capture_report$warnings)
  unlink(paths[1])
  expect_false(CheckATACData(capture)$fragments_ready)
  missing <- CheckATACData(restored)
  expect_true(missing$valid)
  expect_false(missing$fragments_ready)
  expect_match(paste(missing$warnings, collapse = " "), "missing")
  expect_equal(unname(SeuratObject::GetAssayData(restored, layer = "counts")), unname(counts))
})

test_that("ATAC checks return useful invalid reports without altering existing objects", {
  expect_false(CheckATACData(matrix(0, 2, 2))$valid)
  x <- make_atac_input()
  object <- ReadATACData(x$counts, ranges = x$ranges)
  expect_false(CheckATACData(object, assay = "missing")$valid)
  expect_false(CheckATACData(object, sample.by = "absent")$valid)
  object$sample <- c("A", NA, "B")
  expect_false(CheckATACData(object, sample.by = "sample")$valid)
  object[["peaks"]]@ranges <- object[["peaks"]]@ranges[c(3, 1, 2)]
  report <- CheckATACData(object)
  expect_false(report$valid)
  expect_match(paste(report$errors, collapse = " "), "not aligned")
})

test_that("existing native assay seqinfo contributes reference and bounds checks", {
  x <- make_atac_input()
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = c("chr1", "chr2"),
    seqlengths = c(1000L, 1000L), genome = "hg38"
  )
  chromatin <- Signac::CreateChromatinAssay(counts = x$counts, ranges = x$ranges, genome = seqinfo)
  object <- SeuratObject::CreateSeuratObject(chromatin, assay = "peaks")
  expect_true(CheckATACData(object, genome = "hg38")$valid)
  GenomeInfoDb::genome(object[["peaks"]]@ranges) <- "hg19"
  expect_false(CheckATACData(object)$valid)
  GenomeInfoDb::genome(object[["peaks"]]@ranges) <- NA_character_
  GenomeInfoDb::seqlengths(object[["peaks"]]@seqinfo) <- c(100L, 1000L)
  expect_false(CheckATACData(object)$valid)
  expect_match(paste(CheckATACData(object)$errors, collapse = " "), "assay seqinfo")
})

test_that("optional ATAC backends are checked without installation", {
  testthat::local_mocked_bindings(check_r = function(pkg, install, verbose) {
    expect_false(install)
    expect_false(verbose)
    FALSE
  }, .package = "scop")
  expect_error(scop:::atac_backend_fun("missing_atac_backend", "fun"), "install it separately")
})

test_that("native ID cache attributes do not masquerade as an alignment change", {
  skip_if_not_installed("fastmatch")
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = c(1L, 101L), end = c(50L, 150L)
  ))
  counts <- Matrix::Matrix(matrix(c(1L, 0L, 0L, 2L), 2L, 2L,
    dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
  ), sparse = TRUE)
  cells <- c("A-1", "B-1")
  fastmatch::fmatch("A-1", cells)
  expect_false(is.null(attr(cells, ".match.hash")))
  object <- ReadATACData(counts,
    ranges = ranges, cell.names = cells,
    genome = "hg38", verbose = FALSE
  )
  expect_identical(as.vector(colnames(object)), c("A-1", "B-1"))
  expect_identical(rownames(object), c("chr1-1-50", "chr1-101-150"))
  peak_ids <- rownames(object)
  fastmatch::fmatch(peak_ids[[1]], peak_ids)
  expect_false(is.null(attr(peak_ids, ".match.hash")))
  object[["peaks"]]@counts@Dimnames[[1L]] <- peak_ids
  expect_true(CheckATACData(object)$valid)
  object@assays[["peaks"]]@counts@Dimnames[[1L]] <- rev(as.vector(peak_ids))
  expect_false(CheckATACData(object)$valid)
})


test_that("imported QC metadata cannot silently contradict the raw matrix", {
  x <- make_atac_input()
  metadata <- data.frame(
    nCount_peaks = Matrix::colSums(x$counts),
    nFeature_peaks = Matrix::colSums(x$counts > 0), row.names = colnames(x$counts)
  )
  object <- ReadATACData(x$counts, ranges = x$ranges, metadata = metadata, verbose = FALSE)
  expect_identical(object[[]][, names(metadata)], metadata)
  metadata$nCount_peaks[1] <- 0
  expect_error(
    ReadATACData(x$counts, ranges = x$ranges, metadata = metadata, verbose = FALSE),
    "conflicts with the imported raw count matrix"
  )
  names(metadata)[1] <- ""
  expect_error(
    ReadATACData(x$counts, ranges = x$ranges, metadata = metadata, verbose = FALSE),
    "[Mm]etadata must be a data frame"
  )
})


test_that("single-modality H5 requires explicit Peaks rather than coordinate-like RNA IDs", {
  skip_if_not_installed("hdf5r")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(2, 3),
    dims = c(2, 2), dimnames = list(c("chr1:0-50", "chr1:100-150"), c("A-1", "B-1"))
  )
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  file <- hdf5r::H5File$new(path, mode = "w")
  matrix <- file$create_group("matrix")
  features <- matrix$create_group("features")
  matrix[["data"]] <- counts@x
  matrix[["indices"]] <- counts@i
  matrix[["indptr"]] <- counts@p
  matrix[["shape"]] <- dim(counts)
  matrix[["barcodes"]] <- colnames(counts)
  features[["id"]] <- rownames(counts)
  features[["name"]] <- rownames(counts)
  features[["feature_type"]] <- rep("Gene Expression", 2)
  file$close_all()
  expect_error(ReadATACData(path, verbose = FALSE), "explicitly declare a Peaks modality")
  file <- hdf5r::H5File$new(path, mode = "r+")
  file$link_delete("matrix/features/feature_type")
  file$close_all()
  expect_error(ReadATACData(path, verbose = FALSE), "feature_type table")
})


test_that("native barcode caches are not imported as metadata identity", {
  skip_if_not_installed("fastmatch")
  x <- make_atac_input()
  fastmatch::fmatch(colnames(x$counts), colnames(x$counts))
  expect_false(is.null(attr(colnames(x$counts), ".match.hash")))
  metadata <- data.frame(group = c("A", "B", "B"), row.names = as.vector(colnames(x$counts)))
  object <- ReadATACData(x$counts, ranges = x$ranges, metadata = metadata, verbose = FALSE)
  expect_identical(unname(object$group), metadata$group)
  expect_identical(as.vector(object$input_barcode), as.vector(colnames(x$counts)))
})

test_that("ATAC backends resolve public re-exported genomic accessors", {
  for (fun in c("seqnames", "seqinfo", "start")) {
    expect_identical(
      scop:::atac_backend_fun("GenomicRanges", fun),
      getExportedValue("GenomicRanges", fun)
    )
  }
  for (fun in c("seqnames", "seqinfo", "genome", "genome<-", "seqlengths")) {
    expect_identical(
      scop:::atac_backend_fun("GenomeInfoDb", fun),
      getExportedValue("GenomeInfoDb", fun)
    )
  }
})

test_that("ATAC backends report a missing public function", {
  expect_error(
    scop:::atac_backend_fun("GenomicRanges", "missing_atac_function"),
    "not an exported object"
  )
})
test_that("fragment lookup caches do not invalidate original-barcode metadata", {
  skip_if_not_installed("fastmatch")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = 1,
    dims = c(2, 2), dimnames = list(c("p1", "p2"), c("cellA", "cellB"))
  )
  ids <- colnames(counts)
  fastmatch::fmatch(ids, ids)
  colnames(counts) <- ids
  ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  expect_identical(as.vector(object$input_barcode), as.vector(ids))
  expect_identical(as.vector(colnames(object)), as.vector(ids))
})
