make_peak_filter_fixture <- function() {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  counts <- Matrix::Matrix(matrix(c(3, 1, 0, 2, 0, 1, 0, 4, 3, 0, 0, 2),
    nrow = 4L, dimnames = list(paste0("p", 1:4), paste0("cell", 1:3))
  ), sparse = TRUE)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1", "chr2", "chr2"),
    start = c(1L, 101L, 201L, 401L), end = c(50L, 150L, 250L, 500L)
  ))
  ReadATACData(counts, ranges = ranges, genome = "hg38", verbose = FALSE)
}

test_that("explicit peak criteria preserve retained counts and leave source results untouched", {
  object <- make_peak_filter_fixture()
  before <- serialize(object, NULL)
  source <- SeuratObject::GetAssayData(object, layer = "counts")
  result <- FilterATACPeaks(object, min.cells = 2L, min.counts = 2L, max.width = 50L, verbose = FALSE)
  expect_identical(serialize(object, NULL), before)
  expect_identical(rownames(result[["filtered_peaks"]]), rownames(object)[1:2])
  expect_equal(SeuratObject::GetAssayData(result, assay = "filtered_peaks", layer = "counts"), source[1:2, ])
  expect_identical(as.vector(colnames(result[["filtered_peaks"]])), as.vector(colnames(object)))
  expect_identical(SeuratObject::DefaultAssay(result), "peaks")
  expect_length(Signac::Links(result[["filtered_peaks"]]), 0L)
  expect_length(SeuratObject::VariableFeatures(result[["filtered_peaks"]]), 0L)
  receipt <- SeuratObject::Misc(result[["filtered_peaks"]], "scop_peak_filter")
  expect_identical(receipt$peaks$keep, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(receipt$peaks$count, c(6, 2, 0, 8))
  expect_error(FilterATACPeaks(result, verbose = FALSE), "unused")
})

test_that("blacklist filtering uses a verified reference and native overlap semantics", {
  object <- make_peak_filter_fixture()
  blacklist <- Signac::granges(object[["peaks"]])[4]
  result <- FilterATACPeaks(object, blacklist = blacklist, seqlevels = c("chr1", "chr2"), verbose = FALSE)
  expect_identical(rownames(result[["filtered_peaks"]]), rownames(object)[1:3])
  receipt <- SeuratObject::Misc(result[["filtered_peaks"]], "scop_peak_filter")
  expect_identical(receipt$peaks$pass_blacklist, c(TRUE, TRUE, TRUE, FALSE))
  GenomeInfoDb::genome(blacklist) <- "hg19"
  expect_error(FilterATACPeaks(object, blacklist = blacklist, verbose = FALSE), "genome labels")
  expect_error(FilterATACPeaks(object, min.cells = 100, verbose = FALSE), "at least two peaks")
  expect_error(FilterATACPeaks(object, seqlevels = "unknown", verbose = FALSE), "unique chromosome names")
  expect_error(FilterATACPeaks(object, min.width = NA_real_, verbose = FALSE), "Width limits")
  expect_error(FilterATACPeaks(object, min.counts = 0.5, verbose = FALSE), "non-negative integers")
})

test_that("native nearest annotation keeps equal-distance ties and unmatched peaks", {
  object <- make_peak_filter_fixture()
  annotation <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = c(60L, 60L), end = c(90L, 90L), strand = c("+", "-"),
    gene_id = c("g1", "g2"), gene_name = c("G1", "G2"),
    gene_biotype = "protein_coding", type = "gene", tx_id = 1:2
  ), keep.extra.columns = TRUE)
  GenomeInfoDb::genome(annotation) <- "hg38"
  before <- SeuratObject::GetAssayData(object, layer = "counts")
  result <- AnnotateATACPeaks(object, annotation, verbose = FALSE)
  nearest <- SeuratObject::Misc(result[["peaks"]], "scop_peak_annotation")
  expect_identical(nearest$unannotated_peaks, rownames(object)[3:4])
  expect_identical(as.character(nearest$nearest$peak), rep(rownames(object)[1:2], each = 2L))
  expect_identical(as.integer(nearest$nearest$annotation_index), rep(1:2, 2L))
  expect_equal(as.integer(nearest$nearest$distance), c(9L, 9L, 10L, 10L))
  expect_identical(as.character(nearest$nearest$annotation_gene_id), rep(c("g1", "g2"), 2L))
  expect_equal(SeuratObject::GetAssayData(result, layer = "counts"), before)
  expect_identical(Signac::Annotation(result[["peaks"]]), annotation)
  expect_error(AnnotateATACPeaks(result, annotation, verbose = FALSE), "results already exist")
  GenomeInfoDb::genome(annotation) <- "hg19"
  expect_error(AnnotateATACPeaks(object, annotation, verbose = FALSE), "genome labels")
})

make_peak_annotation <- function(start = c(60L, 60L), end = c(90L, 90L), chromosome = "chr1", genome = "hg38") {
  annotation <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = chromosome, start = start, end = end,
    gene_id = paste0("g", seq_along(start)), gene_name = paste0("G", seq_along(start)),
    gene_biotype = "protein_coding", type = "gene", tx_id = seq_along(start)
  ), keep.extra.columns = TRUE)
  GenomeInfoDb::genome(annotation) <- genome
  annotation
}

test_that("peak filtering uses inclusive widths and closed overlap boundaries", {
  object <- make_peak_filter_fixture()
  before <- serialize(object, NULL)
  for (start in c(50L, 51L)) {
    blacklist <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = "chr1", start = start, end = start))
    GenomeInfoDb::genome(blacklist) <- "hg38"
    result <- FilterATACPeaks(object, blacklist = blacklist, verbose = FALSE)
    expected <- if (start == 50L) rownames(object)[-1] else rownames(object)
    expect_identical(rownames(result[["filtered_peaks"]]), expected)
  }
  result <- FilterATACPeaks(object, min.width = 50L, max.width = 50L, verbose = FALSE)
  expect_identical(rownames(result[["filtered_peaks"]]), rownames(object)[1:3])
  expect_identical(serialize(object, NULL), before)
  unknown <- object
  GenomeInfoDb::genome(unknown[["peaks"]]) <- NA_character_
  GenomeInfoDb::genome(unknown[["peaks"]]@ranges) <- NA_character_
  expect_true(CheckATACData(FilterATACPeaks(unknown, verbose = FALSE), assay = "filtered_peaks")$valid)
  expect_error(FilterATACPeaks(unknown, blacklist = blacklist, verbose = FALSE), "verified[[:space:]]+reference")
})

test_that("peak operations preserve partial assay cells and source analysis", {
  object <- make_peak_filter_fixture()
  object <- subset(object, cells = c("cell3", "cell1"))
  object <- SeuratObject::RenameCells(object, new.names = paste0("renamed_", colnames(object)))
  cells <- colnames(object)
  rna <- Matrix::Matrix(matrix(seq_len(6L), 2, 3,
    dimnames = list(c("g1", "g2"), c(cells, "rna_only"))
  ), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  SeuratObject::VariableFeatures(partial[["peaks"]]) <- rownames(partial[["peaks"]])[1:2]
  partial[["tiny"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(1:6, 3, dimnames = list(colnames(partial), c("tiny_1", "tiny_2"))),
    key = "tiny_", assay = "RNA"
  )
  raw <- SeuratObject::GetAssayData(partial, assay = "peaks", layer = "counts")
  partial[["peaks"]] <- SeuratObject::SetAssayData(partial[["peaks"]], layer = "data", new.data = log1p(raw))
  partial[["peaks"]] <- SeuratObject::SetAssayData(partial[["peaks"]], layer = "scale.data", new.data = as.matrix(raw))
  before <- serialize(partial, NULL)
  expected_counts <- SeuratObject::GetAssayData(partial, assay = "peaks", layer = "counts")
  filtered <- FilterATACPeaks(partial, assay = "peaks", verbose = FALSE)
  expect_identical(colnames(filtered), colnames(partial))
  expect_identical(colnames(filtered[["filtered_peaks"]]), cells)
  expect_identical(SeuratObject::GetAssayData(filtered, assay = "filtered_peaks", layer = "counts"), expected_counts)
  expect_equal(Signac::granges(filtered[["filtered_peaks"]]), Signac::granges(partial[["peaks"]]))
  expect_identical(filtered[["peaks"]], partial[["peaks"]])
  expect_identical(filtered[["tiny"]], partial[["tiny"]])
  expect_equal(SeuratObject::GetAssayData(filtered, assay = "filtered_peaks", layer = "data"), expected_counts)
  expect_equal(nrow(SeuratObject::GetAssayData(filtered, assay = "filtered_peaks", layer = "scale.data")), 0L)
  expect_length(SeuratObject::VariableFeatures(filtered[["filtered_peaks"]]), 0L)
  expect_identical(SeuratObject::DefaultAssay(filtered), "RNA")
  annotation <- make_peak_annotation()
  annotated <- AnnotateATACPeaks(partial, annotation, assay = "peaks", verbose = FALSE)
  expect_identical(colnames(annotated), colnames(partial))
  expect_identical(SeuratObject::GetAssayData(annotated, assay = "peaks", layer = "counts"), expected_counts)
  expect_equal(SeuratObject::GetAssayData(annotated, assay = "RNA", layer = "counts"), rna)
  expect_identical(SeuratObject::DefaultAssay(annotated), "RNA")
  expect_identical(annotated[["tiny"]], partial[["tiny"]])
  expect_identical(SeuratObject::GetAssayData(annotated, assay = "peaks", layer = "data"), log1p(expected_counts))
  expect_identical(serialize(partial, NULL), before)
  receipt <- SeuratObject::Misc(annotated[["peaks"]], "scop_peak_annotation")
  expect_identical(receipt$peaks, rownames(partial[["peaks"]]))
})

test_that("annotation preserves native zero-distance, ties and unmatched semantics", {
  object <- make_peak_filter_fixture()
  for (annotation in list(
    make_peak_annotation(c(50L, 51L), c(50L, 51L)),
    make_peak_annotation(1L, 100L, chromosome = "chr3")
  )) {
    before <- serialize(object, NULL)
    result <- suppressWarnings(AnnotateATACPeaks(object, annotation, verbose = FALSE))
    receipt <- SeuratObject::Misc(result[["peaks"]], "scop_peak_annotation")
    native <- suppressWarnings(GenomicRanges::distanceToNearest(Signac::granges(object[["peaks"]]), annotation,
      select = "all", ignore.strand = TRUE
    ))
    order <- order(S4Vectors::queryHits(native), S4Vectors::subjectHits(native))
    expect_identical(as.integer(receipt$nearest$annotation_index), S4Vectors::subjectHits(native)[order])
    expect_equal(as.integer(receipt$nearest$distance), S4Vectors::mcols(native)$distance[order])
    expect_identical(receipt$unannotated_peaks, rownames(object)[setdiff(seq_len(nrow(object)), S4Vectors::queryHits(native))])
    expect_identical(serialize(object, NULL), before)
  }
  ann <- make_peak_annotation()
  Signac::Annotation(object[["peaks"]]) <- ann
  expect_s4_class(AnnotateATACPeaks(object, ann, verbose = FALSE), "Seurat")
  changed <- ann
  S4Vectors::mcols(changed)$gene_name <- c("other1", "other2")
  expect_error(AnnotateATACPeaks(object, changed, verbose = FALSE), "different gene annotation")
  GenomeInfoDb::genome(changed) <- NA_character_
  expect_error(AnnotateATACPeaks(object, changed, verbose = FALSE), "verified[[:space:]]+reference")
  bad <- make_peak_annotation(0L, 5L)
  expect_error(AnnotateATACPeaks(object, bad, verbose = FALSE), "positive 1-based")
  expect_error(FilterATACPeaks(object, blacklist = bad, verbose = FALSE), "positive 1-based")
})

test_that("native transcript_id normalization does not masquerade as new annotation", {
  object <- make_peak_filter_fixture()
  annotation <- make_peak_annotation()
  S4Vectors::mcols(annotation)$transcript_id <- S4Vectors::mcols(annotation)$tx_id
  S4Vectors::mcols(annotation)$tx_id <- NULL
  Signac::Annotation(object[["peaks"]]) <- annotation
  native <- Signac::Annotation(object[["peaks"]])
  expect_true("tx_id" %in% colnames(S4Vectors::mcols(native)))
  before <- serialize(object, NULL)
  result <- AnnotateATACPeaks(object, annotation, verbose = FALSE)
  expect_identical(Signac::Annotation(result[["peaks"]]), native)
  expect_identical(serialize(object, NULL), before)
})

test_that("explicit peak references reject known chromosome-length conflicts", {
  object <- make_peak_filter_fixture()
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- c(chr1 = 1000L, chr2 = 1000L)
  annotation <- make_peak_annotation()
  GenomeInfoDb::seqlengths(annotation) <- 999L
  expect_error(AnnotateATACPeaks(object, annotation, verbose = FALSE), "length")
  expect_error(FilterATACPeaks(object, blacklist = annotation, verbose = FALSE), "length")
  outside <- make_peak_annotation(1001L, 1010L)
  expect_error(AnnotateATACPeaks(object, outside, verbose = FALSE), "length")
  expect_error(FilterATACPeaks(object, blacklist = outside, verbose = FALSE), "length")
  own_outside <- make_peak_annotation(900L, 950L)
  suppressWarnings(GenomeInfoDb::seqlengths(own_outside) <- 800L)
  expect_error(AnnotateATACPeaks(make_peak_filter_fixture(), own_outside, verbose = FALSE), "length")
})

test_that("annotation follows native left-right tie selection without adding hits", {
  object <- make_peak_filter_fixture()
  annotation <- make_peak_annotation(c(70L, 161L), c(90L, 180L))
  result <- AnnotateATACPeaks(object, annotation, verbose = FALSE)
  nearest <- SeuratObject::Misc(result[["peaks"]], "scop_peak_annotation")$nearest
  native <- GenomicRanges::distanceToNearest(Signac::granges(object[["peaks"]]), annotation,
    select = "all", ignore.strand = TRUE
  )
  native_order <- order(S4Vectors::queryHits(native), S4Vectors::subjectHits(native))
  expect_identical(as.character(nearest$peak), rownames(object)[S4Vectors::queryHits(native)[native_order]])
  expect_identical(as.integer(nearest$annotation_index), S4Vectors::subjectHits(native)[native_order])
  expect_equal(as.integer(nearest$distance), S4Vectors::mcols(native)$distance[native_order])
})

test_that("valid custom peak order is preserved and misaligned count matrices are rejected", {
  object <- make_peak_filter_fixture()
  source <- SeuratObject::GetAssayData(object, layer = "counts")
  order <- c(4, 2, 1, 3)
  reordered <- ReadATACData(source[order, ],
    ranges = Signac::granges(object[["peaks"]])[order],
    genome = "hg38", verbose = FALSE
  )
  result <- FilterATACPeaks(reordered, verbose = FALSE)
  expect_identical(rownames(result[["filtered_peaks"]]), rownames(reordered))
  expect_equal(SeuratObject::GetAssayData(result, assay = "filtered_peaks", layer = "counts"), source[order, ])
  expect_equal(Signac::granges(result[["filtered_peaks"]]), Signac::granges(reordered[["peaks"]]))
  broken <- object
  broken@assays[["peaks"]]@counts <- source[order, ]
  expect_error(FilterATACPeaks(broken, verbose = FALSE), "not aligned")
  expect_error(AnnotateATACPeaks(broken, make_peak_annotation(), verbose = FALSE), "not aligned")
  singleton <- subset(object, cells = "cell1")
  expect_error(FilterATACPeaks(singleton, verbose = FALSE), "at least two peaks and two cells")
})

test_that("filtering preserves native fragment maps without reading inactive missing files", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("peak_fragments_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_peak_filter_fixture()
  paths <- file.path(directory, c("s1.tsv.gz", "s2.tsv.gz"))
  records <- list(c("chr1\t0\t50\tA\t5", "chr1\t60\t90\tB\t2"), "chr1\t0\t50\tC\t7")
  mappings <- list(c(cell2 = "B", cell1 = "A"), c(cell3 = "C"))
  fragments <- lapply(seq_along(paths), function(i) {
    plain <- sub(".gz$", "", paths[[i]])
    writeLines(records[[i]], plain)
    Rsamtools::bgzip(plain, dest = paths[[i]], overwrite = TRUE)
    Rsamtools::indexTabix(paths[[i]], format = "bed")
    Signac::CreateFragmentObject(paths[[i]], cells = mappings[[i]], tolerance = 0, verbose = FALSE)
  })
  Signac::Fragments(object[["peaks"]]) <- fragments
  object <- subset(object, cells = c("cell1", "cell2"))
  object <- SeuratObject::RenameCells(object, new.names = c("sample_A", "sample_B"))
  unlink(c(paths[2], paste0(paths[2], ".tbi")))
  rna <- Matrix::Matrix(matrix(1:6, 2, dimnames = list(c("g1", "g2"), c(colnames(object), "rna_only"))), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  original <- serialize(partial, NULL)
  result <- FilterATACPeaks(partial, assay = "peaks", verbose = FALSE)
  actual_fragments <- Signac::Fragments(result[["filtered_peaks"]])
  original_fragments <- Signac::Fragments(partial[["peaks"]])
  expect_length(actual_fragments, length(original_fragments))
  for (i in seq_along(original_fragments)) {
    expect_identical(Signac::GetFragmentData(actual_fragments[[i]], "path"), Signac::GetFragmentData(original_fragments[[i]], "path"))
    expect_identical(Signac::GetFragmentData(actual_fragments[[i]], "hash"), Signac::GetFragmentData(original_fragments[[i]], "hash"))
    expected_map <- Signac::GetFragmentData(original_fragments[[i]], "cells")
    expected_map <- expected_map[names(expected_map) %in% colnames(partial[["peaks"]])]
    actual_map <- Signac::GetFragmentData(actual_fragments[[i]], "cells")
    expect_identical(as.vector(actual_map), as.vector(expected_map))
    expect_identical(as.character(names(actual_map)), as.character(names(expected_map)))
  }
  expect_identical(CheckATACData(result, assay = "filtered_peaks")$fragments$active, c(TRUE, FALSE))
  expect_true(CheckATACData(result, assay = "filtered_peaks")$fragments_ready)
  expect_identical(colnames(result), colnames(partial))
  # Native Signac adds only fastmatch's pointer cache to counts/data cell IDs.
  uncached <- partial
  attr(uncached@assays[["peaks"]]@counts@Dimnames[[2]], ".match.hash") <- NULL
  attr(uncached@assays[["peaks"]]@data@Dimnames[[2]], ".match.hash") <- NULL
  expect_identical(uncached, unserialize(original))
  annotation <- make_peak_annotation()
  result <- AnnotateATACPeaks(partial, annotation, assay = "peaks", verbose = FALSE)
  expect_identical(Signac::Fragments(result[["peaks"]]), Signac::Fragments(partial[["peaks"]]))
  # Native Signac adds only fastmatch's pointer cache to counts/data cell IDs.
  uncached <- partial
  attr(uncached@assays[["peaks"]]@counts@Dimnames[[2]], ".match.hash") <- NULL
  attr(uncached@assays[["peaks"]]@data@Dimnames[[2]], ".match.hash") <- NULL
  expect_identical(uncached, unserialize(original))
})
