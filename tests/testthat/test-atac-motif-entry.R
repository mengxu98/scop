make_atac_motif_entry_fixture <- function() {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("TFBSTools")
  skip_if_not_installed("motifmatchr")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2), j = c(1, 1, 2, 2), x = 1,
    dims = c(2L, 2L), dimnames = list(c("p1", "p2"), c("A", "B"))
  )
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = "chr1", start = c(101, 301), end = c(200, 400)))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  genome <- Biostrings::DNAStringSet(c(chr1 = paste(rep("ACGT", 200), collapse = "")))
  pfm <- TFBSTools::PFMatrixList(M1 = TFBSTools::PFMatrix(
    ID = "M1", name = "ACGT",
    profileMatrix = matrix(
      c(
        10L, 0L, 0L, 0L, 0L, 10L, 0L, 0L,
        0L, 0L, 10L, 0L, 0L, 0L, 0L, 10L
      ),
      nrow = 4,
      dimnames = list(c("A", "C", "G", "T"), NULL)
    )
  ))
  list(object = object, genome = genome, pfm = pfm)
}

test_that("final ranges can be scanned before counts exist and reused after recount", {
  f <- make_atac_motif_entry_fixture()
  ranges <- Signac::granges(f$object[["peaks"]])
  original <- ranges
  scanned <- AddATACMotifs(ranges, f$genome, f$pfm, "toy", verbose = FALSE)
  native <- Signac::AddMotifs(ranges, genome = f$genome, pfm = f$pfm, verbose = FALSE)
  expect_s4_class(scanned$motifs, "Motif")
  expect_equal(scanned$motifs, native)
  expect_identical(ranges, original)
  expect_identical(scanned$receipt$peaks, rownames(f$object[["peaks"]]))
  expect_identical(rownames(scanned$statistics), scanned$receipt$peaks)
  assay <- f$object[["peaks"]]
  Signac::Motifs(assay) <- scanned$motifs
  assay <- SeuratObject::AddMetaData(assay, scanned$statistics)
  expected <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE)
  expect_equal(Signac::Motifs(assay), Signac::Motifs(expected[["peaks"]]))
  expect_equal(assay[[]], expected[["peaks"]][[]])
  expect_error(AddATACMotifs(c(ranges, ranges[1]), f$genome, f$pfm, "toy"), "unique|[Dd]uplicate")
  expect_error(AddATACMotifs(ranges, f$genome, f$pfm, "other"), "genome")
  expect_error(AddATACMotifs(ranges[FALSE], f$genome, f$pfm, "toy"), "non-empty")
})

test_that("motif addition preserves native incidence, sequence statistics and counts", {
  fixture <- make_atac_motif_entry_fixture()
  native <- Signac::AddMotifs(fixture$object, genome = fixture$genome, pfm = fixture$pfm, verbose = FALSE)
  actual <- AddATACMotifs(fixture$object, fixture$genome, fixture$pfm, "toy", verbose = FALSE)
  expect_equal(Signac::GetMotifData(actual, slot = "data"), Signac::GetMotifData(native, slot = "data"))
  expect_equal(Signac::GetMotifData(actual, slot = "positions"), Signac::GetMotifData(native, slot = "positions"))
  expect_equal(actual[["peaks"]][[]], native[["peaks"]][[]])
  expect_equal(
    SeuratObject::LayerData(actual[["peaks"]], layer = "counts"),
    SeuratObject::LayerData(fixture$object[["peaks"]], layer = "counts")
  )
  expect_identical(SeuratObject::Misc(actual[["peaks"]], "scop_atac_motifs")$motifs, "M1")
  expect_error(AddATACMotifs(actual, fixture$genome, fixture$pfm, "toy"), "already contains")
  expect_error(AddATACMotifs(fixture$object, fixture$genome, fixture$pfm, "other"), "genome.id")
  short_reference <- Biostrings::DNAStringSet(c(chr1 = paste(rep("A", 350), collapse = "")))
  expect_error(AddATACMotifs(fixture$object, short_reference, fixture$pfm, "toy"), "reference lengths")
})

test_that("replacing motifs refreshes native statistics without duplicate columns", {
  f <- make_atac_motif_entry_fixture()
  first <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE)
  first[["peaks"]] <- SeuratObject::AddMetaData(first[["peaks"]], data.frame(
    GC.percent = c(-1, -1), custom = c("a", "b"), row.names = rownames(first)
  ))
  expected <- Signac::RegionStats(Signac::granges(first[["peaks"]]), genome = f$genome, verbose = FALSE)
  expect_warning(actual <- AddATACMotifs(first, f$genome, f$pfm, "toy", replace = TRUE, verbose = FALSE), "Overwriting miscellanous data")
  metadata <- actual[["peaks"]][[]]
  expect_false(anyDuplicated(colnames(metadata)) > 0L)
  expect_equal(metadata$GC.percent, expected$GC.percent)
  expect_equal(metadata$custom, c("a", "b"))
  expect_equal(first[["peaks"]][[]]$GC.percent, c(-1, -1))
})

test_that("motif addition rejects conflicting known source reference lengths", {
  f <- make_atac_motif_entry_fixture()
  GenomeInfoDb::seqlengths(f$object[["peaks"]]) <- c(chr1 = 1000L)
  expect_error(AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE), "lengths.*agree|lengths.*conflict")
})

test_that("sequence reference chromosome identifiers must be unambiguous", {
  f <- make_atac_motif_entry_fixture()
  duplicate <- c(f$genome, f$genome)
  expect_error(AddATACMotifs(f$object, duplicate, f$pfm, "toy", verbose = FALSE), "duplicated sequence names")
})

make_atac_sequence_fixture <- function(duplicate.display = FALSE) {
  f <- make_atac_motif_entry_fixture()
  sequence <- withr::with_seed(7L, paste(sample(c("A", "C", "G", "T"), 4800L, replace = TRUE), collapse = ""))
  starts <- seq(1L, by = 120L, length.out = 40L)
  for (i in seq_len(20L)) substr(sequence, starts[i] + 30L, starts[i] + 37L) <- if (i <= 10L) "ACGTACGT" else "TTGGCCAA"
  reference <- Biostrings::DNAStringSet(c(chr1 = sequence))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = "chr1", start = starts, end = starts + 99L))
  counts <- Matrix::Matrix(matrix(1, 40, 3, dimnames = list(paste0("p", seq_len(40)), c("A", "B", "C"))), sparse = TRUE)
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  make_pfm <- function(sequence, id) {
    bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
    profile <- matrix(0L, 4, length(bases), dimnames = list(c("A", "C", "G", "T"), NULL))
    profile[cbind(match(bases, rownames(profile)), seq_along(bases))] <- 100L
    TFBSTools::PFMatrix(ID = id, name = if (duplicate.display) "shared_display_name" else paste0("display_", id), profileMatrix = profile)
  }
  pfm <- TFBSTools::PFMatrixList(M2 = make_pfm("TTGGCCAA", "M2"), M1 = make_pfm("ACGTACGT", "M1"))
  list(object = object, genome = reference, pfm = pfm)
}

test_that("real sequence matching feeds native explicit and matched enrichment", {
  f <- make_atac_sequence_fixture()
  before <- serialize(f$object, NULL)
  matched <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE)
  expect_identical(serialize(f$object, NULL), before)
  native <- Signac::AddMotifs(f$object, genome = f$genome, pfm = f$pfm, verbose = FALSE)
  expect_equal(Signac::GetMotifData(matched, slot = "data"), Signac::GetMotifData(native, slot = "data"))
  expect_identical(as.character(colnames(Signac::GetMotifData(matched, slot = "data"))), c("M2", "M1"))
  ids <- rownames(matched)
  result <- RunMotifEnrichment(matched, ids[1:10], background = ids, verbose = FALSE)
  expect_equal(result$results, Signac::FindMotifs(native, ids[1:10], background = ids, verbose = FALSE))
  expect_true(result$results$observed[result$results$motif == "M1"] > 0)
  metadata <- native[["peaks"]][[]]
  expected <- withr::with_seed(13L, Signac::MatchRegionStats(metadata[ids[-(1:10)], , drop = FALSE],
    query.feature = metadata[ids[1:10], , drop = FALSE], features.match = c("GC.percent", "sequence.length"), n = 15, verbose = FALSE
  ))
  actual <- RunMotifEnrichment(matched, ids[1:10], n.background = 15, seed = 13L, verbose = FALSE)
  expect_identical(actual$receipt$reference, expected)
  expect_equal(actual$results, Signac::FindMotifs(native, ids[1:10], background = c(ids[1:10], expected), verbose = FALSE))
})

test_that("indexed FASTA and DNAStringSet produce identical native motif scans", {
  f <- make_atac_sequence_fixture()
  directory <- tempfile("motif_fasta_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- file.path(directory, "reference.fa")
  Biostrings::writeXStringSet(f$genome, path)
  Rsamtools::indexFa(path)
  fasta <- Rsamtools::FaFile(path)
  actual <- AddATACMotifs(f$object, fasta, f$pfm, "toy", verbose = FALSE)
  expected <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE)
  expect_equal(Signac::GetMotifData(actual, slot = "data"), Signac::GetMotifData(expected, slot = "data"))
  expect_equal(actual[["peaks"]][[]], expected[["peaks"]][[]])
})

test_that("partial-assay motif updates preserve RNA and native subset identity", {
  f <- make_atac_sequence_fixture()
  rna_cells <- c("C", "rna_only", "A", "B")
  rna <- Matrix::Matrix(matrix(seq_len(8), 2, dimnames = list(c("g1", "g2"), rna_cells)), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  object[["peaks"]] <- f$object[["peaks"]]
  before <- serialize(object, NULL)
  actual <- AddATACMotifs(object, f$genome, f$pfm, "toy", assay = "peaks", verbose = FALSE)
  expect_identical(serialize(object, NULL), before)
  expect_identical(actual[["RNA"]], object[["RNA"]])
  expect_identical(colnames(actual), rna_cells)
  expect_identical(SeuratObject::DefaultAssay(actual), "RNA")
  expect_equal(SeuratObject::GetAssayData(actual, assay = "peaks", layer = "counts"), SeuratObject::GetAssayData(object, assay = "peaks", layer = "counts"))
  trimmed <- subset(actual, cells = c("C", "A"))
  expect_warning(trimmed[["peaks"]] <- subset(trimmed[["peaks"]], features = rownames(trimmed[["peaks"]])[seq(2, 40, 2)]), "Different cells and/or features")
  trimmed <- SeuratObject::RenameCells(trimmed, add.cell.id = "renamed")
  ids <- rownames(trimmed[["peaks"]])
  result <- RunMotifEnrichment(trimmed, ids[1:10], assay = "peaks", background = ids, verbose = FALSE)
  expect_equal(result$results, Signac::FindMotifs(trimmed, ids[1:10], assay = "peaks", background = ids, verbose = FALSE))
  expect_identical(result$receipt$universe, ids)
})

test_that("explicit labels cannot override conflicting loaded-reference metadata", {
  f <- make_atac_motif_entry_fixture()
  GenomeInfoDb::genome(f$genome) <- "other"
  expect_identical(unname(GenomeInfoDb::genome(f$genome)), "other")
  expect_error(AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE), "reference genome labels must agree")
  GenomeInfoDb::genome(f$genome) <- "toy"
  matched <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE)
  expect_identical(SeuratObject::Misc(matched[["peaks"]], "scop_atac_motifs")$genome, "toy")
  GenomeInfoDb::genome(f$object[["peaks"]]) <- NA_character_
  expect_error(AddATACMotifs(f$object, f$genome, f$pfm, NA_character_, verbose = FALSE), "genome.id")
})

test_that("motif identifiers and explicit sequence inputs are validated before scanning", {
  f <- make_atac_motif_entry_fixture()
  expect_error(AddATACMotifs(f$object, "toy", f$pfm, "toy", verbose = FALSE), "loaded sequence reference")
  expect_error(AddATACMotifs(f$object, f$genome, list(M1 = matrix(1, 4, 4)), "toy", verbose = FALSE), "PFMatrixList")
  repeated <- c(f$pfm, f$pfm)
  expect_error(AddATACMotifs(f$object, f$genome, repeated, "toy", verbose = FALSE), "unique non-empty motif")
  expect_error(AddATACMotifs(f$object, f$genome, f$pfm, "toy", replace = NA, verbose = FALSE), "[Rr]eplace must")
  wrong_chromosome <- f$genome
  names(wrong_chromosome) <- "1"
  expect_error(AddATACMotifs(f$object, wrong_chromosome, f$pfm, "toy", verbose = FALSE), "every source-peak chromosome")
})

test_that("duplicate display names follow native naming while motif IDs remain exact", {
  f <- make_atac_sequence_fixture(duplicate.display = TRUE)
  expect_warning(actual <- AddATACMotifs(f$object, f$genome, f$pfm, "toy", verbose = FALSE), "Non-unique motif names")
  expect_warning(native <- Signac::AddMotifs(f$object, genome = f$genome, pfm = f$pfm, verbose = FALSE), "Non-unique motif names")
  expect_identical(as.character(colnames(Signac::GetMotifData(actual, slot = "data"))), c("M2", "M1"))
  expect_equal(Signac::GetMotifData(actual, slot = "motif.names"), Signac::GetMotifData(native, slot = "motif.names"))
})

test_that("explicit PWMatrixList uses the same native motif ordering", {
  f <- make_atac_sequence_fixture()
  pwm <- TFBSTools::toPWM(f$pfm)
  actual <- AddATACMotifs(f$object, f$genome, pwm, "toy", verbose = FALSE)
  native <- Signac::AddMotifs(f$object, genome = f$genome, pfm = pwm, verbose = FALSE)
  expect_equal(Signac::GetMotifData(actual, slot = "data"), Signac::GetMotifData(native, slot = "data"))
  expect_identical(as.character(colnames(Signac::GetMotifData(actual, slot = "data"))), names(pwm))
})
