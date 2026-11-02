make_atac_motif_fixture <- function() {
  skip_if_not_installed("Signac")
  counts <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 2, 3, 4, 5, 3, 2, 5, 4), 3, 4,
    dimnames = list(paste0("p", 1:3), paste0("c", 1:4))), sparse = TRUE)
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1000, 2000, 3000), width = 100))
  x <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  coords <- matrix(c(1, 2, 3, 4, 2, 5, 3, 6), 4, 2,
    dimnames = list(colnames(x), c("LSI_1", "LSI_2")))
  x[["lsi"]] <- SeuratObject::CreateDimReducObject(coords, key = "LSI_", assay = "peaks")
  x
}


test_that("motif logos use explicit IDs from the selected assay", {
  x <- make_atac_motif_fixture()
  incidence <- Matrix::Matrix(matrix(1, nrow(x), 1, dimnames = list(rownames(x), "M1")), sparse = TRUE)
  pwm <- list(M1 = matrix(seq_len(20), 4, 5, dimnames = list(c("A", "C", "G", "T"), NULL)))
  Signac::Motifs(x[["peaks"]]) <- Signac::CreateMotifObject(data = incidence, pwm = pwm)
  expect_error(ATACMotifPlot(x, "unknown"), "unique motif IDs")
  expect_error(ATACMotifPlot(x, c("M1", "M1")), "unique motif IDs")
  skip_if_not_installed("ggseqlogo")
  actual <- ATACMotifPlot(x, "M1", use.names = FALSE)
  expected <- Signac::MotifPlot(x, "M1", use.names = FALSE)
  expect_equal(actual$data, expected$data)
})
