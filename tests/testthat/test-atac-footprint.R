test_that("footprints refuse missing fragments and plots refuse stale receipts", {
  x <- make_atac_test_object()
  expect_error(RunATACFootprint(x, NULL, "toy", "chr1-4-9999", "M1"), "indexed fragments")
  expect_error(ATACFootprintPlot(x, "M1"), "missing or stale")
  SeuratObject::Misc(x[["peaks"]], "scop_footprints") <- list(M1 = list(cells = rev(colnames(x)), peaks = rownames(x)))
  expect_error(ATACFootprintPlot(x, "M1"), "missing or stale")
})

test_that("bias-corrected footprint matches native Signac on indexed synthetic fragments", {
  skip_if_not_installed("Signac")
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("Biostrings")
  withr::local_preserve_seed()
  set.seed(51)
  directory <- tempfile("footprint-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fasta <- file.path(directory, "toy.fa")
  dna <- Biostrings::DNAStringSet(c(chr1 = paste(sample(c("A", "C", "G", "T"), 100000, TRUE), collapse = "")))
  Biostrings::writeXStringSet(dna, fasta)
  Rsamtools::indexFa(fasta)
  genome <- Rsamtools::FaFile(fasta)
  start <- 10:99900
  plain <- file.path(directory, "fragments.tsv")
  utils::write.table(data.frame("chr1", start, start + 50L, rep(c("raw1", "raw2"), length.out = length(start)), 1),
    plain,
    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
  )
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  counts <- Matrix::Matrix(matrix(c(2, 3, 4, 5), 2, dimnames = list(c("p1", "p2"), c("c1", "c2"))), sparse = TRUE)
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1000, 2000), width = 100))
  fragments <- Signac::CreateFragmentObject(path, cells = c(c1 = "raw1", c2 = "raw2"), verbose = FALSE)
  x <- ReadATACData(counts, ranges = ranges, genome = "toy", fragments = path, fragment.cells = c(c1 = "raw1", c2 = "raw2"), verbose = FALSE)
  pwm <- list(M1 = matrix(1, 4, 6, dimnames = list(c("A", "C", "G", "T"), NULL)))
  incidence <- Matrix::Matrix(matrix(1, 2, 1, dimnames = list(rownames(x), "M1")), sparse = TRUE)
  positions <- GenomicRanges::GRangesList(M1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1040, 2040), width = 6), strand = "+"))
  Signac::Motifs(x[["peaks"]]) <- Signac::CreateMotifObject(data = incidence, pwm = pwm, positions = positions)
  original <- serialize(x, NULL)
  bad <- x
  GenomeInfoDb::genome(bad[["peaks"]]@motifs@positions[["M1"]]) <- "other"
  expect_error(RunATACFootprint(bad, genome, "toy", "chr1-4-99997", "M1", verbose = FALSE), "reference labels")
  bad <- x
  bad[["peaks"]]@motifs@positions[["M1"]] <- GenomicRanges::GRanges("chr1", IRanges::IRanges(400, width = 6), strand = "-")
  expect_error(RunATACFootprint(bad, genome, "toy", "chr1-4-99997", "M1", upstream = 50, downstream = 500, verbose = FALSE), "windows must fit")
  bad[["peaks"]]@motifs@positions[["M1"]] <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50000, width = 6), strand = "+")
  expect_error(RunATACFootprint(bad, genome, "toy", "chr1-4-99997", "M1", in.peaks = TRUE, verbose = FALSE), "No motif positions")
  signac_version <- utils::packageVersion("Signac")
  actual <- RunATACFootprint(x, genome, "toy", "chr1-4-99997", "M1", upstream = 50, downstream = 50, verbose = FALSE)
  expect_identical(serialize(x, NULL), original)
  native_input <- unserialize(serialize(x, NULL))
  expected <- Signac::InsertionBias(native_input, assay = "peaks", genome = genome, region = "chr1-4-99997", verbose = FALSE)
  expected <- Signac::Footprint(expected, assay = "peaks", genome = genome, motif.name = "M1", upstream = 50, downstream = 50, compute.expected = TRUE, verbose = FALSE)
  expect_equal(
    SeuratObject::GetAssayData(actual[["peaks"]], layer = "positionEnrichment"),
    SeuratObject::GetAssayData(expected[["peaks"]], layer = "positionEnrichment")
  )
  expect_identical(serialize(x, NULL), original)
  expect_error(RunATACFootprint(actual, genome, "toy", "chr1-4-99997", "M1", verbose = FALSE), "already exist")
  legacy_line_width <- signac_version <= "1.17.1" && utils::packageVersion("ggplot2") >= "3.4.0"
  plotted <- withCallingHandlers(
    withr::with_options(
      list(lifecycle_verbosity = "warning"),
      ATACFootprintPlot(actual, "M1")
    ),
    warning = function(w) {
      if (legacy_line_width && inherits(w, "lifecycle_warning_deprecated")) {
        expect_match(conditionMessage(w), "size.*aesthetic for lines")
        if (grepl("size.*aesthetic for lines", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    }
  )
  expect_s3_class(plotted, "ggplot")
  stale <- actual
  GenomeInfoDb::genome(stale[["peaks"]]@seqinfo) <- "other"
  expect_error(ATACFootprintPlot(stale, "M1"), "missing or stale|conflicting")
})
