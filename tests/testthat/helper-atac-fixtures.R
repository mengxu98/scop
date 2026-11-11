make_atac_test_object <- function() {
  skip_if_not_installed("Signac")
  counts <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 2, 3, 4, 5, 3, 2, 5, 4), 3, 4,
    dimnames = list(paste0("p", 1:3), paste0("c", 1:4))
  ), sparse = TRUE)
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1000, 2000, 3000), width = 100))
  x <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  coords <- matrix(c(1, 2, 3, 4, 2, 5, 3, 6), 4, 2,
    dimnames = list(colnames(x), c("LSI_1", "LSI_2"))
  )
  x[["lsi"]] <- SeuratObject::CreateDimReducObject(coords, key = "LSI_", assay = "peaks")
  x
}

make_atac_plot_test_object <- function() {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  counts <- Matrix::Matrix(matrix(c(1, 2, 2, 3, 3, 4, 4, 5), 2L,
    dimnames = list(c("chr1-1-100", "chr1-101-200"), c("A", "B", "C", "D"))
  ), sparse = TRUE)
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 101), width = 100))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  object$group <- c("one", "one", "two", "two")
  object$verified.depth <- c(10, 20, 30, 40)
  embedding <- cbind(LSI_1 = 1:4, LSI_2 = 4:1, LSI_3 = rep(1, 4))
  rownames(embedding) <- colnames(object)
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(embedding, assay = "peaks", key = "LSI_")
  object
}

make_atac_fragment_test_object <- function(directory) {
  ranges <- GenomicRanges::GRanges(c("chr1:80-120", "chr1:140-180"))
  GenomeInfoDb::seqlengths(ranges) <- c(chr1 = 250)
  counts <- Matrix::Matrix(matrix(c(2, 1, 3, 1, 1, 2),
    nrow = 2,
    dimnames = list(c("p1", "p2"), c("S1_A", "S1_B", "S2_A"))
  ), sparse = TRUE)
  paths <- file.path(directory, c("one.tsv.gz", "two.tsv.gz"))
  lines <- list(
    c(
      "chr1\t90\t100\tA\t7", "chr1\t90\t100\tB\t2",
      "chr1\t95\t105\tB\t1", "chr1\t100\t110\tA\t1"
    ),
    c("chr1\t90\t100\tA\t3", "chr1\t150\t160\tA\t1")
  )
  for (i in seq_along(paths)) {
    path <- sub(".gz$", "", paths[i])
    writeLines(lines[[i]], path)
    Rsamtools::bgzip(path, dest = paths[i], overwrite = TRUE)
    Rsamtools::indexTabix(paths[i], format = "bed")
  }
  object <- ReadATACData(counts,
    ranges = ranges, fragments = paths,
    fragment.cells = list(c(S1_A = "A", S1_B = "B"), c(S2_A = "A")),
    genome = "toy", verbose = FALSE
  )
  object$sample <- c("first", "first", "second")
  object
}
