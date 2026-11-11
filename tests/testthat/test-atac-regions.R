atac_regions_test_centers <- function() {
  ranges <- GenomicRanges::GRanges(c("chr1:100-100", "chr1:99-101"), strand = c("+", "-"))
  names(ranges) <- c("forward", "reverse")
  ranges
}

test_that("region endpoint counting preserves duplicates, strand and barcode mapping", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("regions_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_fragment_test_object(directory)
  original <- serialize(object, NULL)
  result <- RunATACRegionMatrix(object, atac_regions_test_centers(),
    group.by = "sample",
    upstream = 10, downstream = 10, bin.size = 1, verbose = FALSE
  )
  expect_s3_class(result, "scop_atac_regions")
  expect_identical(names(result$matrices), c("first", "second"))
  plus <- numeric(21)
  plus[c(2, 7, 11, 12, 16, 21)] <- c(2, 1, 2, 1, 1, 1)
  expect_equal(unname(result$matrices$first["forward", ]), plus)
  expect_equal(unname(result$matrices$first["reverse", ]), rev(plus))
  expect_equal(sum(result$matrices$second["forward", ]), 2)
  expect_equal(as.numeric(result$cells.per.group), c(2, 1))
  expect_equal(result$positions, -10:10)
  expect_false(result$receipt$column5_weighted)
  expect_identical(result$receipt$endpoint_shift, 0L)
  expect_identical(serialize(object, NULL), original)
})

test_that("region bins include boundaries and partial last bins without sampling", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("regions_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_atac_fragment_test_object(directory)
  result <- RunATACRegionMatrix(object, atac_regions_test_centers(),
    group.by = "sample",
    cells = c("S1_B", "S1_A"), upstream = 10, downstream = 10,
    bin.size = 10, verbose = FALSE
  )
  expect_equal(result$bin.widths, c(10, 10, 1))
  expect_equal(result$positions, c(-5.5, 4.5, 10))
  expect_equal(unname(result$matrices$first[1, ]), c(3, 4, 1))
  center <- RunATACRegionMatrix(object, atac_regions_test_centers()[1],
    group.by = "sample",
    upstream = 0, downstream = 0, bin.size = 1, verbose = FALSE
  )
  expect_equal(unname(center$matrices$first[1, 1]), 2)
  expect_error(RunATACRegionMatrix(object, atac_regions_test_centers(), max.entries = 1), "max.entries")
  expect_error(RunATACRegionMatrix(object, GenomicRanges::GRanges("chr1:1-1")), "boundaries")
  expect_error(RunATACRegionMatrix(object, atac_regions_test_centers(), cells = "unknown"), "source-assay")
  object$sample[1] <- NA
  expect_error(RunATACRegionMatrix(object, atac_regions_test_centers(), group.by = "sample"), "non-missing")
})

test_that("region plots expose raw and normalized values and actual bin widths", {
  x <- list(
    matrices = list(a = matrix(c(2, 4, 6, 8), 2,
      dimnames = list(c("r1", "r2"), c("0", "1"))
    )), positions = c(0, 1),
    bin.widths = c(2, 1), cells.per.group = c(a = 2)
  )
  class(x) <- c("scop_atac_regions", "list")
  heat <- ATACRegionHeatmap(x)
  expect_s3_class(heat, "ggplot")
  expect_equal(heat$data$value, c(0.5, 1, 3, 4))
  expect_equal(heat$data$bin.width, c(2, 2, 1, 1))
  curve <- ATACRegionPlot(x)
  expect_equal(curve$data$value, c(0.75, 3.5))
  raw <- ATACRegionHeatmap(x, normalize = FALSE, per.base = FALSE, log1p = TRUE)
  expect_equal(raw$data$value, log1p(c(2, 4, 6, 8)))
  expect_error(ATACRegionPlot(x, groups = "missing"), "groups")
  expect_error(ATACRegionPlot(x, normalize = NA), "TRUE or FALSE")
  bad <- x
  bad$positions[1] <- NA_real_
  expect_error(ATACRegionPlot(bad), "positions")
  bad <- x
  bad$matrices$b <- x$matrices$a[2:1, , drop = FALSE]
  bad$cells.per.group <- c(a = 2, b = 2)
  expect_error(ATACRegionHeatmap(bad), "Malformed")
})
