make_atac_export_fixture <- function() {
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr2"),
    start = c(1L, 201L), end = c(100L, 300L)
  ))
  counts <- Matrix::Matrix(matrix(1:6, 2L, 3L, dimnames = list(c("p1", "p2"), c("A", "B", "C"))), sparse = TRUE)
  ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
}

test_that("ATAC export roundtrips all sparse counts, cells and explicit tables", {
  object <- make_atac_export_fixture()
  base <- tempfile("export_")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  path <- file.path(base, "output")
  table <- data.frame(peak = rownames(object), p = c(0.02, 0.9))
  original <- serialize(object, NULL)
  receipt <- ExportATACResults(object, path, tables = list(all_tested = table), counts = TRUE)
  expect_identical(serialize(object, NULL), original)
  restored <- readRDS(file.path(path, "object.rds"))
  expect_equal(restored, object)
  expect_true(CheckATACData(restored)$valid)
  expect_false(CheckATACData(restored)$fragments_ready)
  bed <- read.table(file.path(path, "peaks.bed"), sep = "\t")
  expect_equal(bed[[2]], c(0L, 200L))
  expect_equal(bed[[3]], c(100L, 300L))
  expect_equal(
    as.matrix(Matrix::readMM(file.path(path, "counts.mtx"))),
    unname(as.matrix(SeuratObject::LayerData(object[["peaks"]], layer = "counts")))
  )
  expect_identical(readLines(file.path(path, "barcodes.tsv")), colnames(object))
  expect_equal(read.csv(file.path(path, "tables/all_tested.csv"), row.names = 1), table)
  expect_identical(unname(tools::md5sum(file.path(path, receipt$manifest$file))), receipt$manifest$md5)
  expect_error(ExportATACResults(object, path), "already exists")
  expect_error(ExportATACResults(object, file.path(base, "bad"), tables = list("../escape" = table)), "safely named")
  expect_false(dir.exists(file.path(base, "bad")))
})

test_that("failed exports leave neither a destination nor staging remnants", {
  object <- make_atac_export_fixture()
  base <- tempfile("export_atomic_")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  expect_error(ExportATACResults(object, file.path(base, "failed"), compress = "invalid"))
  expect_identical(list.files(base, all.files = TRUE, no.. = TRUE), character())
})

test_that("export refuses dangling symlink destinations without changing them", {
  skip_on_os("windows")
  object <- make_atac_export_fixture()
  base <- tempfile("export_symlink_")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  path <- file.path(base, "output")
  target <- file.path(base, "missing-target")
  expect_true(file.symlink(target, path))
  expect_error(ExportATACResults(object, path), "already exists")
  expect_identical(Sys.readlink(path), target)
  expect_false(file.exists(target))
  expect_identical(list.files(base, all.files = TRUE, no.. = TRUE), "output")
})

test_that("exported Matrix Market companions are readable by ReadATACData", {
  object <- make_atac_export_fixture()
  object <- subset(object, features = rev(rownames(object)))
  path <- tempfile("export_mex_")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  ExportATACResults(object, path, counts = TRUE)
  restored <- ReadATACData(file.path(path, "counts.mtx"),
    peaks.file = file.path(path, "peaks.bed"),
    barcodes.file = file.path(path, "barcodes.tsv"), genome = "toy", verbose = FALSE
  )
  expect_identical(
    SeuratObject::GetAssayData(restored, layer = "counts"),
    SeuratObject::GetAssayData(object, layer = "counts")
  )
  expect_equal(unname(Signac::granges(restored)), unname(Signac::granges(object)))
  bed <- read.table(file.path(path, "peaks.bed"), sep = "\t")
  expect_identical(as.character(bed[[1]]), as.character(GenomicRanges::seqnames(Signac::granges(object))))
  expect_equal(bed[[2]] + 1L, GenomicRanges::start(Signac::granges(object)))
  expect_equal(bed[[3]], GenomicRanges::end(Signac::granges(object)))
})

test_that("export rejects table filenames which collide case-insensitively", {
  object <- make_atac_export_fixture()
  path <- tempfile("export_collision_")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  expect_error(ExportATACResults(object, path,
    tables = list(DAR = data.frame(p = 1), dar = data.frame(p = 2))
  ), "safely named")
  expect_false(file.exists(path))
})

make_atac_export_fragments <- function(object, directory) {
  paths <- file.path(directory, c("s1.tsv.gz", "s2.tsv.gz"))
  records <- list(
    c("chr1\t0\t50\trawA\t5", "chr1\t60\t90\trawB\t2"),
    c("chr1\t0\t50\trawC\t7")
  )
  maps <- list(c(B = "rawB", A = "rawA"), c(C = "rawC"))
  fragments <- lapply(seq_along(paths), function(i) {
    plain <- sub(".gz$", "", paths[[i]])
    writeLines(records[[i]], plain)
    Rsamtools::bgzip(plain, dest = paths[[i]], overwrite = TRUE)
    Rsamtools::indexTabix(paths[[i]], format = "bed")
    Signac::CreateFragmentObject(paths[[i]], cells = maps[[i]], tolerance = 0, verbose = FALSE)
  })
  Signac::Fragments(object[["peaks"]]) <- fragments
  object
}

test_that("subset and renamed fragment mappings retain inactive resource snapshots", {
  skip_if_not_installed("Rsamtools")
  base <- tempfile("export_fragments_")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  object <- make_atac_export_fragments(make_atac_export_fixture(), base)
  object <- subset(object, cells = c("B", "A"))
  object <- SeuratObject::RenameCells(object, new.names = paste0("renamed_", colnames(object)))
  fragments <- Signac::Fragments(object[["peaks"]])
  inactive <- Signac::GetFragmentData(fragments[[2]], "path")
  unlink(c(inactive, paste0(inactive, ".tbi")))
  expect_true(CheckATACData(object)$fragments_ready)
  original <- serialize(object, NULL)
  path <- file.path(base, "output")
  receipt <- ExportATACResults(object, path, counts = TRUE)
  resources <- read.csv(file.path(path, "fragments.csv"))
  expect_identical(resources$active, c(TRUE, FALSE))
  expect_identical(resources$file_exists, c(TRUE, FALSE))
  expect_identical(resources$n_cells, vapply(fragments, function(x) length(Signac::GetFragmentData(x, "cells")), integer(1)))
  expect_identical(resources$fragment, c(1L, 2L))
  expect_equal(nrow(resources), 2L)
  expect_identical(
    resources$stored_signac_fragment_hash,
    vapply(fragments, function(x) Signac::GetFragmentData(x, "hash")[[1]], character(1))
  )
  mapping <- read.csv(file.path(path, "fragment_cells.csv"))
  native <- Signac::GetFragmentData(fragments[[1]], "cells")
  expect_identical(mapping$cell, names(native))
  expect_identical(mapping$barcode, unname(native))
  expect_true(all(mapping$fragment == 1L))
  expect_identical(serialize(object, NULL), original)
  expect_equal(readRDS(file.path(path, "object.rds")), object)
  expect_true(readRDS(file.path(path, "check.rds"))$fragments_ready)
  expect_identical(unname(tools::md5sum(file.path(path, receipt$manifest$file))), receipt$manifest$md5)
  active <- Signac::GetFragmentData(fragments[[1]], "path")
  unlink(paste0(active, ".tbi"))
  unavailable <- ExportATACResults(object, file.path(base, "matrix_only"))
  expect_true(unavailable$check$valid)
  expect_false(unavailable$check$fragments_ready)
  object[["peaks"]]@fragments[[1]]@cells <- character()
  empty <- ExportATACResults(object, file.path(base, "empty_mapping"))
  expect_identical(
    names(read.csv(file.path(empty$directory, "fragment_cells.csv"))),
    c("fragment", "cell", "barcode")
  )
  expect_false(empty$check$fragments_ready)
})

test_that("partial assay export preserves the full object and aligns companion cells", {
  atac <- make_atac_export_fixture()
  rna <- Matrix::Matrix(matrix(1:8, 2, 4,
    dimnames = list(c("g1", "g2"), c("C", "rna_only", "A", "B"))
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  object[["peaks"]] <- atac[["peaks"]]
  object$label <- paste0("meta_", colnames(object))
  original <- serialize(object, NULL)
  path <- tempfile("export_partial_")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  receipt <- ExportATACResults(object, path, assay = "peaks", counts = TRUE)
  cells <- colnames(object[["peaks"]])
  expect_identical(readLines(file.path(path, "barcodes.tsv")), cells)
  metadata <- object[[]][cells, , drop = FALSE]
  metadata[] <- lapply(metadata, function(x) if (is.factor(x)) as.character(x) else x)
  expect_equal(read.csv(file.path(path, "cells.csv"), row.names = 1), metadata)
  expect_equal(
    as.matrix(Matrix::readMM(file.path(path, "counts.mtx"))),
    unname(as.matrix(SeuratObject::GetAssayData(object, assay = "peaks", layer = "counts")))
  )
  restored <- readRDS(file.path(path, "object.rds"))
  expect_equal(restored, object)
  expect_identical(colnames(restored), colnames(object))
  expect_equal(SeuratObject::GetAssayData(restored, assay = "RNA", layer = "counts"), rna)
  expect_identical(serialize(object, NULL), original)
  expect_equal(receipt$check$n_cells, length(cells))
})

test_that("export validates destinations and inputs without overwriting existing files", {
  object <- make_atac_export_fixture()
  base <- tempfile("export_paths_")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  existing <- file.path(base, "existing.txt")
  writeLines("preserve me", existing)
  expect_error(ExportATACResults(object, existing), "already exists")
  expect_identical(readLines(existing), "preserve me")
  for (outdir in list(NA_character_, character(), c("a", "b"), "", 5)) {
    expect_error(ExportATACResults(object, outdir), "[Oo]utdir")
  }
  for (counts in list(NA, 1, c(TRUE, FALSE))) {
    expect_error(ExportATACResults(object, file.path(base, "invalid"), counts = counts), "[Cc]ounts")
  }
  table <- data.frame(p = 1)
  for (tables in list(
    list(table), list("a/b" = table), list(a = 1),
    setNames(list(table, table), c("same", "same"))
  )) {
    expect_error(ExportATACResults(object, file.path(base, "invalid"), tables = tables), "safely named")
  }
  bad <- table
  bad$unsupported <- list(list(value = 1))
  expect_error(ExportATACResults(object, file.path(base, "failed"), tables = list(bad = bad)))
  expect_identical(list.files(base, all.files = TRUE, no.. = TRUE), "existing.txt")
})

test_that("single-cell peak subsets retain their complete RDS representation", {
  object <- make_atac_export_fixture()
  object <- subset(object, cells = "B", features = rownames(object)[2])
  path <- tempfile("export_small_")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  receipt <- ExportATACResults(object, path, counts = TRUE)
  expect_true(receipt$check$valid)
  expect_equal(readRDS(file.path(path, "object.rds")), object)
  expect_equal(
    as.matrix(Matrix::readMM(file.path(path, "counts.mtx"))),
    unname(as.matrix(SeuratObject::GetAssayData(object, layer = "counts")))
  )
})
