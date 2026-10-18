test_that("narrowPeak import distinguishes an empty successful file from missing or malformed output", {
  path <- tempfile(fileext = ".narrowPeak")
  on.exit(unlink(path))
  expect_error(atac_read_narrowpeaks(path, "toy"), "did not produce")
  file.create(path)
  expect_length(atac_read_narrowpeaks(path, "toy"), 0)
  writeLines("chr1\t0\t100\tpeak1\t100\t.\t8\t10\t6\t50", path)
  result <- atac_read_narrowpeaks(path, "toy")
  expect_identical(GenomicRanges::start(result), 1L)
  expect_identical(GenomicRanges::end(result), 100L)
  expect_equal(result$fold_change, 8)
  expect_identical(unname(GenomeInfoDb::genome(result)), "toy")
  writeLines("chr1\t100\t99\tpeak1\t100\t.\t8\t10\t6\t50", path)
  expect_error(atac_read_narrowpeaks(path, "toy"), "Malformed")
  writeLines("chr1\t0\t100", path)
  expect_error(atac_read_narrowpeaks(path, "toy"), "Malformed")
})

test_that("peak-call grouping separates metadata field names from output IDs", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("atac_caller_mapping_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c("chr1\t10\t20\tA\t1", "chr1\t110\t120\tB\t1"), plain)
  fragments <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"))
  Rsamtools::indexTabix(fragments, format = "bed")
  counts <- Matrix::sparseMatrix(
    i = 1:2, j = 1:2, x = 1,
    dimnames = list(c("p1", "p2"), c("A", "B"))
  )
  ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", fragments = fragments, verbose = FALSE)
  object$call_id <- c("capture/A", "capture_A")
  object$sample <- c("group/1", "group_1")
  for (cells in list(character(), NA_character_, c("A", "A"), "unknown")) {
    expect_error(
      CallATACPeaks(object, sample.by = "sample", cells = cells),
      "unique selected-assay cell IDs"
    )
  }
  program <- file.path(directory, "macs2")
  file.create(program)
  call <- CallATACPeaks
  execution <- new.env(parent = environment(call))
  # Unix returns NA for a missing path; reproduce that behavior on every OS.
  execution$Sys.readlink <- function(path) NA_character_
  # Mock only the process boundary; fragment splitting and parsing are native.
  # Real MACS2 acceptance is performed separately on external PBMC data.
  execution$system2 <- function(command, args, stdout, stderr) {
    if (identical(args, "--version")) {
      return("macs2 2.2.9.1")
    }
    args <- gsub("^['\"]|['\"]$", "", args)
    id <- args[match("-n", args) + 1L]
    outdir <- args[match("--outdir", args) + 1L]
    writeLines(
      "chr1\t10\t20\tpeak1\t10\t.\t2\t3\t4\t5",
      file.path(outdir, paste0(id, "_peaks.narrowPeak"))
    )
    writeLines("Mock caller output for the grouping contract", stdout)
    0L
  }
  environment(call) <- execution
  result <- call(object,
    sample.by = "call_id", group.by = "sample",
    macs2.path = program, outdir = file.path(directory, "calls"),
    effective.genome.size = 1000, verbose = FALSE
  )
  expect_identical(names(result$groups), c("call_id", "sample", "group"))
  expect_identical(result$groups$call_id, names(result$peaks))
  expect_identical(result$groups$sample, unname(object$call_id))
  expect_identical(result$groups$group, unname(object$sample))
  expect_identical(result$receipt$sample.by, "call_id")
  expect_identical(result$receipt$group.by, "sample")
  expect_error(call(object,
    sample.by = "call_id", group.by = "sample",
    macs2.path = program, outdir = file.path(directory, "calls"),
    effective.genome.size = 1000, verbose = FALSE
  ), "[Oo]utdir")
  execution$Sys.readlink <- function(path) "missing-symlink-target"
  expect_error(call(object,
    sample.by = "call_id", group.by = "sample",
    macs2.path = program, outdir = file.path(directory, "dangling"),
    effective.genome.size = 1000, verbose = FALSE
  ), "[Oo]utdir")
  expect_false(dir.exists(file.path(directory, "dangling")))
})
