make_atac_tiles_fixture <- function(directory) {
  ranges <- GenomicRanges::GRanges(c("chr1:80-120", "chr1:140-180"))
  GenomeInfoDb::seqlengths(ranges) <- c(chr1 = 250)
  counts <- Matrix::Matrix(matrix(c(2, 1, 3, 1, 1, 2), nrow = 2,
    dimnames = list(c("p1", "p2"), c("S1_A", "S1_B", "S2_A"))), sparse = TRUE)
  paths <- file.path(directory, c("one.tsv.gz", "two.tsv.gz"))
  lines <- list(c("chr1\t90\t100\tA\t7", "chr1\t90\t100\tB\t2",
    "chr1\t95\t105\tB\t1", "chr1\t100\t110\tA\t1"),
    c("chr1\t90\t100\tA\t3", "chr1\t150\t160\tA\t1"))
  for (i in seq_along(paths)) {
    path <- sub(".gz$", "", paths[i])
    writeLines(lines[[i]], path)
    Rsamtools::bgzip(path, dest = paths[i], overwrite = TRUE)
    Rsamtools::indexTabix(paths[i], format = "bed")
  }
  object <- ReadATACData(counts, ranges = ranges, fragments = paths,
    fragment.cells = list(c(S1_A = "A", S1_B = "B"), c(S2_A = "A")),
    genome = "toy", verbose = FALSE)
  object$sample <- c("first", "first", "second")
  object
}

test_that("genome tiles use verified lengths and retain sparse recount semantics", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("tiles_"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_tiles_fixture(directory)
  original <- serialize(object, NULL)
  result <- RunGenomeTiles(object, bin.size = 100, process_n = 1, verbose = FALSE)
  expect_s4_class(result$counts, "dgCMatrix")
  expect_equal(GenomicRanges::start(result$ranges), c(1, 101, 201))
  expect_equal(GenomicRanges::end(result$ranges), c(100, 200, 250))
  expect_equal(unname(as.matrix(result$counts)), matrix(c(1, 1, 0, 2, 1, 0, 1, 1, 0), nrow = 3))
  expect_equal(result$receipt$tiles$n_tiles, 3)
  expect_identical(serialize(object, NULL), original)
  expect_error(RunGenomeTiles(object, seqlengths = c(chr1 = 249)), "conflict")
  expect_error(RunGenomeTiles(object, bin.size = 100, max.tiles = 1), "max.tiles")
  expect_error(RunGenomeTiles(object, seqlengths = c(chr1 = NA_real_)), "known positive")
  expect_error(RunGenomeTiles(object, chromosomes = "chr2"), "known contigs")
  written <- RunGenomeTiles(object, bin.size = 100, new.assay = "tiles", verbose = FALSE)
  expect_identical(SeuratObject::DefaultAssay(written), SeuratObject::DefaultAssay(object))
  expect_equal(SeuratObject::Misc(written[["tiles"]], "scop_atac_recount")$tiles$n_tiles, 3)
})
