atac_consensus_fixture <- function() {
  gr <- function(start, end) {
    x <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = "chr1", start, end))
    GenomeInfoDb::genome(x) <- "hg38"
    x
  }
  list(
    ctrl1 = gr(c(1, 200), c(20, 220)), ctrl2 = gr(10, 30),
    case1 = gr(c(101, 401), c(120, 410)), case2 = gr(110, 130)
  )
}

test_that("peak consensus matches native interval union without inventing counts", {
  peaks <- atac_consensus_fixture()
  out <- CombineATACPeaks(object = peaks)
  expect_identical(names(formals(CombineATACPeaks))[[1L]], "object")
  expect_false("peaks" %in% names(formals(CombineATACPeaks)))
  expected <- GenomicRanges::reduce(do.call(c, unname(peaks)), ignore.strand = TRUE, min.gapwidth = 1L)
  expect_identical(atac_range_ids(out$ranges), atac_range_ids(expected))
  expect_identical(names(out$ranges), atac_range_ids(expected))
  expect_identical(
    as.list(out$ranges$input_sets),
    list(c("ctrl1", "ctrl2"), c("case1", "case2"), "ctrl1", "case1")
  )
  expect_identical(out$receipt$input_sets$n_peaks, c(2L, 1L, 2L, 1L))
  expect_false("counts" %in% names(out))
  expect_identical(out$receipt$genome, "hg38")
})

test_that("support is filtered within groups before union", {
  peaks <- atac_consensus_fixture()
  groups <- c(case2 = "case", ctrl2 = "control", case1 = "case", ctrl1 = "control")
  out <- CombineATACPeaks(peaks, groups = groups, min.support = 2)
  expect_identical(names(out$ranges), c("chr1-1-30", "chr1-101-130"))
  expect_identical(as.list(out$ranges$support_groups), list("control", "case"))
  expect_identical(out$receipt$input_sets$group, c("control", "control", "case", "case"))
  expect_identical(out$receipt$support_inventory$control$kept, c(TRUE, FALSE))
  expect_identical(out$receipt$support_inventory$case$kept, c(TRUE, FALSE))
  expect_error(CombineATACPeaks(peaks, groups, min.support = 3), "No peaks remain")
})

test_that("width limits apply to the final connected union", {
  peaks <- atac_consensus_fixture()
  out <- CombineATACPeaks(peaks, min.width = 21, max.width = 30)
  expect_identical(names(out$ranges), c("chr1-1-30", "chr1-101-130", "chr1-200-220"))
  expect_identical(out$receipt$discarded_by_width, "chr1-401-410")
  expect_error(CombineATACPeaks(peaks, min.width = 31), "width filtering")
  expect_error(CombineATACPeaks(peaks, min.width = 0), "Width limits")
  expect_error(CombineATACPeaks(peaks, max.width = NA_real_), "Width limits")
})

test_that("reference, source identities and empty inputs are explicit", {
  peaks <- atac_consensus_fixture()
  expect_error(CombineATACPeaks(unname(peaks)), "uniquely named")
  expect_error(CombineATACPeaks(peaks, groups = c(ctrl1 = "x")), "every peak-set")
  expect_error(CombineATACPeaks(peaks, min.support = 1.5), "positive integer")
  GenomeInfoDb::genome(peaks$ctrl1) <- "hg19"
  expect_error(CombineATACPeaks(peaks), "conflicting reference")
  GenomeInfoDb::genome(peaks$ctrl1) <- NA_character_
  expect_error(CombineATACPeaks(peaks), "explicit verified genome")
  expect_identical(CombineATACPeaks(peaks, genome = "hg38")$receipt$genome, "hg38")
  peaks$empty <- peaks$ctrl1[FALSE]
  expect_length(CombineATACPeaks(peaks, genome = "hg38")$ranges, 4L)
  expect_error(CombineATACPeaks(list(empty = peaks$empty)), "No input peaks")
  duplicate <- atac_consensus_fixture()
  duplicate$ctrl1 <- c(duplicate$ctrl1, duplicate$ctrl1[1])
  expect_error(CombineATACPeaks(duplicate), "Duplicate genomic")
})

test_that("cross-group connected unions preserve ordered distinct support metadata", {
  gr <- function(start, end) {
    x <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start, end))
    GenomeInfoDb::genome(x) <- "hg38"
    x
  }
  peaks <- list(
    ctrl1 = gr(c(1, 31), c(10, 40)),
    ctrl2 = gr(c(6, 36), c(15, 45)), case1 = gr(10, 36), case2 = gr(20, 30)
  )
  groups <- c(ctrl1 = "control", ctrl2 = "control", case1 = "case", case2 = "case")
  out <- CombineATACPeaks(peaks, groups, min.support = 2L)
  expect_identical(names(out$ranges), "chr1-1-45")
  expect_identical(
    as.list(out$ranges$input_sets),
    list(c("ctrl1", "ctrl2", "case1", "case2"))
  )
  expect_identical(as.list(out$ranges$support_groups), list(c("control", "case")))
  expect_identical(out$receipt$support_inventory$control$support, c(2L, 2L))
  expect_identical(out$receipt$support_inventory$case$support, 2L)
})


test_that("consensus errors preserve literal braces in caller-supplied set names", {
  peak <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 20L))
  peaks <- stats::setNames(list(c(peak, peak)), "sample{literal}")
  expect_error(CombineATACPeaks(peaks), "sample\\{literal\\}: Duplicate genomic")
})

test_that("consensus rejects intervals beyond lengths inherited from another input", {
  a <- GenomicRanges::GRanges("chr1:1-10")
  b <- GenomicRanges::GRanges("chr1:150-160")
  GenomeInfoDb::genome(a) <- GenomeInfoDb::genome(b) <- "hg38"
  GenomeInfoDb::seqlengths(a) <- 100L
  expect_error(suppressWarnings(CombineATACPeaks(list(a = a, b = b))), "chromosome lengths")
  expect_error(suppressWarnings(CombineATACPeaks(list(a = a, b = b), max.width = 10L)), "chromosome lengths")
})

test_that("adjacent intervals join without inventing per-base support", {
  first <- GenomicRanges::GRanges(c("chr1:101-110", "chr1:1-10"))
  second <- GenomicRanges::GRanges("chr1:11-20")
  peaks <- list(later = first, adjacent = second)
  before <- serialize(peaks, NULL)
  out <- CombineATACPeaks(peaks, min.support = 2L, genome = "hg38")
  expect_identical(names(out$ranges), "chr1-1-20")
  expect_identical(as.list(out$ranges$input_sets), list(c("later", "adjacent")))
  expect_equal(out$receipt$support_inventory$all$support, c(2L, 1L))
  expect_identical(serialize(peaks, NULL), before)
  GenomeInfoDb::genome(first) <- GenomeInfoDb::genome(second) <- "hg38"
  GenomeInfoDb::seqlengths(first) <- 1000L
  GenomeInfoDb::seqlengths(second) <- 999L
  expect_error(CombineATACPeaks(list(first = first, second = second)), "seqlengths")
})
