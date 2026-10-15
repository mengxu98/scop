make_gene_activity_fixture <- function(directory) {
  skip_if_not_installed("Rsamtools")
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = c(101, 301, 601), end = c(150, 350, 650)
  ))
  GenomeInfoDb::genome(ranges) <- "hg38"
  cells <- c("S_A", "S_B", "S_C")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = 1:3, x = 1,
    dims = c(3, 3), dimnames = list(c("p1", "p2", "p3"), cells)
  )
  plain <- file.path(directory, "records.tsv")
  writeLines(c(
    "chr1\t100\t120\tA\t5", "chr1\t130\t140\tB\t1",
    "chr1\t300\t320\tB\t2", "chr1\t610\t630\tC\t3"
  ), plain)
  path <- paste0(plain, ".gz")
  Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  object <- ReadATACData(counts,
    ranges = ranges, genome = "hg38", fragments = path,
    fragment.cells = c(S_A = "A", S_B = "B", S_C = "C"), verbose = FALSE
  )
  annotation <- ranges
  annotation$gene_id <- c("ENSG1", "ENSG2", "ENSG3")
  annotation$gene_name <- c("GENE1", "GENE2", "GENE3")
  annotation$gene_biotype <- "protein_coding"
  annotation$type <- "gene"
  annotation$tx_id <- annotation$gene_id
  Signac::Annotation(object[["peaks"]]) <- annotation
  object
}

test_that("gene activity preserves native counts, normalization and source assays", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("gene_activity_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  original <- serialize(object, NULL)
  expected <- Signac::GeneActivity(object,
    assay = "peaks", extend.upstream = 0,
    gene.id = TRUE, process_n = 1, verbose = FALSE
  )
  result <- RunGeneActivity(object,
    extend.upstream = 0, gene.id = TRUE,
    process_n = 1, verbose = FALSE
  )
  actual <- SeuratObject::GetAssayData(result, assay = "ACTIVITY", layer = "counts")
  expect_equal(actual, methods::as(expected, "dgCMatrix"))
  expect_equal(unname(as.matrix(actual)), matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow = 3))
  expect_identical(serialize(object, NULL), original)
  expect_identical(result[["peaks"]], object[["peaks"]])
  expect_identical(SeuratObject::DefaultAssay(result), "peaks")
  # Native fragment validation may cache a fastmatch hash on the source ID
  # vector; that cache is not part of the ordered cell identifiers.
  expect_identical(as.vector(colnames(result[["ACTIVITY"]])), as.vector(colnames(object[["peaks"]])))
  normalized <- Seurat::NormalizeData(SeuratObject::CreateAssayObject(counts = expected),
    normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE
  )
  expect_equal(
    SeuratObject::GetAssayData(result, assay = "ACTIVITY", layer = "data"),
    SeuratObject::GetAssayData(normalized, layer = "data")
  )
  receipt <- SeuratObject::Misc(result[["ACTIVITY"]], "scop_gene_activity")
  expect_identical(receipt$genes, rownames(actual))
  expect_identical(receipt$genome, "hg38")
  expect_identical(receipt$parameters$extend.upstream, 0)
  expect_match(receipt$fragments$fragment_md5, "^[a-f0-9]{32}$")
  expect_error(RunGeneActivity(result, verbose = FALSE), "unused name")
})

test_that("activity requires compatible annotation and actual fragment resources", {
  skip_if_not_installed("Rsamtools")
  directory <- tempfile("gene_activity_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  expect_error(RunGeneActivity(object, normalize = NA, verbose = FALSE), "[Nn]ormalize")
  expect_error(RunGeneActivity(object, scale.factor = 0, verbose = FALSE), "[Ss]cale.factor")
  expect_error(RunGeneActivity(object, new.assay = "peaks", verbose = FALSE), "unused name")
  no_annotation <- object
  Signac::Annotation(no_annotation[["peaks"]]) <- Signac::Annotation(object[["peaks"]])[FALSE]
  expect_error(RunGeneActivity(no_annotation, verbose = FALSE), "transcript annotation")
  path <- Signac::GetFragmentData(Signac::Fragments(object[["peaks"]])[[1]], slot = "path")
  unlink(paste0(path, ".tbi"))
  expect_error(RunGeneActivity(object, verbose = FALSE), "index[[:space:]]+is missing")
})


test_that("stable gene IDs distinguish duplicate display symbols natively", {
  directory <- tempfile("activity_symbols_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  annotation <- Signac::Annotation(object[["peaks"]])
  annotation$gene_name[1:2] <- "SAME"
  Signac::Annotation(object[["peaks"]]) <- annotation
  for (by_id in c(TRUE, FALSE)) {
    native <- Signac::GeneActivity(object, gene.id = by_id, extend.upstream = 0, verbose = FALSE)
    actual <- RunGeneActivity(object, gene.id = by_id, extend.upstream = 0, normalize = FALSE, verbose = FALSE)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "ACTIVITY", layer = "counts"), methods::as(native, "dgCMatrix"))
    expect_identical(rownames(actual[["ACTIVITY"]]), if (by_id) c("ENSG1", "ENSG2", "ENSG3") else c("SAME", "SAME.1", "GENE3"))
  }
})

test_that("shared counting windows never silently merge distinct genes", {
  directory <- tempfile("activity_shared_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  annotation <- Signac::Annotation(object[["peaks"]])
  extra <- annotation[1]
  extra$gene_id <- "ENSG4"
  extra$gene_name <- "GENE4"
  extra$tx_id <- "ENSG4"
  Signac::Annotation(object[["peaks"]]) <- c(annotation, extra)
  before <- serialize(object, NULL)
  native <- Signac::GeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
  if (anyDuplicated(rownames(native))) {
    expect_error(RunGeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE), "duplicate gene IDs")
  } else {
    actual <- RunGeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "ACTIVITY", layer = "counts"), methods::as(native, "dgCMatrix"))
  }
  expect_identical(serialize(object, NULL), before)
})

test_that("annotation known lengths and reference identities must agree", {
  directory <- tempfile("activity_reference_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  annotation <- Signac::Annotation(object[["peaks"]])
  GenomeInfoDb::seqlengths(annotation) <- c(chr1 = 1000L)
  Signac::Annotation(object[["peaks"]]) <- annotation
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- c(chr1 = 900L)
  expect_error(RunGeneActivity(object, verbose = FALSE), "lengths must agree")
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- c(chr1 = 500L)
  expect_error(RunGeneActivity(object, verbose = FALSE), "known chromosome lengths|exceed chromosome lengths")
  GenomeInfoDb::seqlengths(object[["peaks"]]) <- c(chr1 = 1000L)
  GenomeInfoDb::genome(annotation) <- "other"
  expect_error(Signac::Annotation(object[["peaks"]]) <- annotation, "Annotation genome does not match")
  object[["peaks"]]@annotation <- annotation
  expect_error(RunGeneActivity(object, verbose = FALSE), "conflicting genome labels")
})

make_activity_capture_fixture <- function(directory) {
  template <- make_gene_activity_fixture(directory)
  paths <- character(2)
  records <- list(
    c("chr1\t100\t120\tA\t5", "chr1\t130\t140\tB\t1"),
    c("chr1\t300\t320\tA\t2", "chr1\t610\t630\tB\t3")
  )
  for (i in 1:2) {
    plain <- file.path(directory, paste0("capture", i, ".tsv"))
    writeLines(records[[i]], plain)
    paths[i] <- paste0(plain, ".gz")
    Rsamtools::bgzip(plain, dest = paths[i], overwrite = TRUE)
    Rsamtools::indexTabix(paths[i], format = "bed")
  }
  counts <- Matrix::Matrix(matrix(seq_len(12), 3, dimnames = list(c("p1", "p2", "p3"), c("S1_A", "S1_B", "S2_A", "S2_B"))), sparse = TRUE)
  ReadATACData(counts,
    ranges = Signac::granges(template[["peaks"]]), genome = "hg38",
    annotation = Signac::Annotation(template[["peaks"]]), fragments = paths,
    fragment.cells = list(c(S1_A = "A", S1_B = "B"), c(S2_A = "A", S2_B = "B")), verbose = FALSE
  )
}

test_that("capture mappings and inactive missing resources remain exact after rename", {
  directory <- tempfile("activity_captures_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_activity_capture_fixture(directory)
  native <- Signac::GeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
  result <- RunGeneActivity(object, gene.id = TRUE, extend.upstream = 0, normalize = FALSE, process_n = 1, verbose = FALSE)
  expect_equal(SeuratObject::GetAssayData(result, assay = "ACTIVITY", layer = "counts"), methods::as(native, "dgCMatrix"))
  expect_equal(unname(as.matrix(SeuratObject::GetAssayData(result, assay = "ACTIVITY", layer = "counts"))), matrix(c(1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1), 3))
  subset_object <- subset(object, cells = c("S1_A", "S1_B"))
  subset_object <- SeuratObject::RenameCells(subset_object, add.cell.id = "renamed")
  inactive <- Signac::GetFragmentData(Signac::Fragments(subset_object[["peaks"]])[[2]], slot = "path")
  unlink(c(inactive, paste0(inactive, ".tbi")))
  before <- serialize(subset_object, NULL)
  result <- RunGeneActivity(subset_object, gene.id = TRUE, extend.upstream = 0, normalize = FALSE, verbose = FALSE)
  expect_identical(serialize(subset_object, NULL), before)
  expect_identical(Signac::Fragments(result[["peaks"]]), Signac::Fragments(subset_object[["peaks"]]))
  expect_identical(as.vector(colnames(result[["ACTIVITY"]])), paste0("renamed_", c("S1_A", "S1_B")))
  expect_equal(nrow(SeuratObject::Misc(result[["ACTIVITY"]], "scop_gene_activity")$fragments), 1)
  expect_equal(unname(as.matrix(SeuratObject::GetAssayData(result, assay = "ACTIVITY", layer = "counts"))), matrix(c(1, 0, 0, 1, 0, 0), 3))
})

test_that("missing gene identity cannot collapse unrelated annotation rows", {
  directory <- tempfile("activity_missing_ids_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  for (missing in c(NA_character_, "")) {
    annotation <- Signac::Annotation(object[["peaks"]])
    annotation$gene_id[1:2] <- missing
    Signac::Annotation(object[["peaks"]]) <- annotation
    expect_error(RunGeneActivity(object, gene.id = FALSE, extend.upstream = 0, verbose = FALSE), "gene_id.*non-empty|non-empty.*gene_id")
  }
})

test_that("partial activity preserves RNA, reductions, metadata and key-collision RNG", {
  directory <- tempfile("activity_partial_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  atac <- make_activity_capture_fixture(directory)
  cells <- c("S2_B", "rna_only", "S1_A", "S2_A", "S1_B")
  rna <- Matrix::Matrix(matrix(seq_len(10), 2, dimnames = list(c("g1", "g2"), cells)), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  object[["peaks"]] <- atac[["peaks"]]
  object[["tiny"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(10), 5, dimnames = list(cells, c("tiny_1", "tiny_2"))), key = "tiny_", assay = "RNA"
  )
  SeuratObject::Key(object[["peaks"]]) <- "activity_"
  before <- serialize(object, NULL)
  set.seed(83L)
  rng <- .Random.seed
  expect_warning(actual <- RunGeneActivity(object, assay = "peaks", gene.id = TRUE, extend.upstream = 0, normalize = FALSE, verbose = FALSE), "[Kk]ey")
  expect_true(identical(.Random.seed, rng))
  expect_identical(serialize(object, NULL), before)
  expect_identical(actual[["RNA"]], object[["RNA"]])
  expect_identical(actual[["peaks"]], object[["peaks"]])
  expect_identical(actual[["tiny"]], object[["tiny"]])
  expect_identical(actual[[]][, colnames(object[[]]), drop = FALSE], object[[]])
  expect_identical(colnames(actual), colnames(object))
  expect_identical(SeuratObject::DefaultAssay(actual), "RNA")
  expect_identical(as.vector(colnames(actual[["ACTIVITY"]])), as.vector(colnames(object[["peaks"]])))
  receipt <- SeuratObject::Misc(actual[["ACTIVITY"]], "scop_gene_activity")
  expect_identical(as.vector(receipt$cells), as.vector(colnames(actual[["ACTIVITY"]])))
  expect_null(receipt$normalization)
  expect_true(is.na(actual[[]]["rna_only", "nCount_ACTIVITY"]))
})

test_that("native transcript spans and strand-aware promoter boundaries are explicit", {
  directory <- tempfile("activity_promoters_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  template <- make_gene_activity_fixture(directory)
  plain <- file.path(directory, "promoters.tsv")
  writeLines(c(
    "chr1\t79\t80\tA\t4", "chr1\t80\t81\tA\t5", "chr1\t100\t101\tA\t7",
    "chr1\t149\t150\tB\t8", "chr1\t150\t151\tB\t9", "chr1\t185\t186\tA\t2",
    "chr1\t289\t290\tC\t1", "chr1\t290\t291\tC\t2", "chr1\t349\t350\tA\t3",
    "chr1\t369\t370\tB\t5", "chr1\t370\t371\tC\t6", "chr1\t610\t630\tC\t3"
  ), plain)
  path <- paste0(plain, ".gz")
  Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  annotation <- Signac::Annotation(template[["peaks"]])
  GenomicRanges::strand(annotation) <- c("+", "-", "+")
  extra <- annotation[1]
  GenomicRanges::end(extra) <- 250L
  GenomicRanges::start(extra) <- 201L
  extra$tx_id <- "TX1b"
  extra$type <- "transcript"
  annotation <- c(annotation, extra)
  object <- ReadATACData(SeuratObject::GetAssayData(template, layer = "counts"),
    ranges = Signac::granges(template[["peaks"]]), genome = "hg38", annotation = annotation,
    fragments = path, fragment.cells = c(S_A = "A", S_B = "B", S_C = "C"), verbose = FALSE
  )
  native <- Signac::GeneActivity(object,
    features = c("ENSG1", "ENSG2"), gene.id = TRUE,
    extend.upstream = 20, extend.downstream = 10, process_n = 2, verbose = FALSE
  )
  actual <- RunGeneActivity(object,
    features = c("ENSG1", "ENSG2"), gene.id = TRUE,
    extend.upstream = 20, extend.downstream = 10, process_n = 1, normalize = FALSE, verbose = FALSE
  )
  counts <- SeuratObject::GetAssayData(actual, assay = "ACTIVITY", layer = "counts")
  expect_equal(counts, methods::as(native, "dgCMatrix"))
  expect_equal(unname(as.matrix(counts)), matrix(c(3, 1, 2, 1, 0, 1), 2))
  # A record within the inter-transcript gap contributes to the collapsed gene span.
  expect_equal(counts["ENSG1", "S_A"], 3)
  filtered <- RunGeneActivity(object,
    gene.id = TRUE, max.width = 150,
    extend.upstream = 20, normalize = FALSE, verbose = FALSE
  )
  expect_identical(rownames(filtered[["ACTIVITY"]]), c("ENSG2", "ENSG3"))
  expect_error(RunGeneActivity(object, gene.id = TRUE, max.width = 1, verbose = FALSE), "No genes remaining")
})

test_that("one gene identity cannot describe incompatible contigs or strands", {
  directory <- tempfile("activity_identity_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  annotation <- Signac::Annotation(object[["peaks"]])
  GenomeInfoDb::seqlevels(annotation) <- c("chr1", "chr2")
  GenomeInfoDb::genome(annotation) <- "hg38"
  extra <- annotation[1]
  extra$tx_id <- "other_transcript"
  GenomicRanges::seqnames(extra) <- factor("chr2", levels = GenomeInfoDb::seqlevels(extra))
  Signac::Annotation(object[["peaks"]]) <- c(annotation, extra)
  expect_error(RunGeneActivity(object, gene.id = FALSE, extend.upstream = 0, verbose = FALSE), "gene_id.*contig|gene identity")
  extra <- annotation[1]
  GenomicRanges::strand(extra) <- "-"
  Signac::Annotation(object[["peaks"]]) <- c(annotation, extra)
  expect_error(RunGeneActivity(object, gene.id = FALSE, extend.upstream = 20, verbose = FALSE), "gene_id.*contig|gene identity")
})

test_that("native single-gene dimension dropping is reported without reconstruction", {
  directory <- tempfile("activity_dimensions_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  native <- Signac::GeneActivity(object, features = "ENSG1", gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
  if (!inherits(native, "sparseMatrix")) {
    expect_error(RunGeneActivity(object, features = "ENSG1", gene.id = TRUE, extend.upstream = 0, verbose = FALSE), "sparse gene-by-cell matrix")
  } else {
    actual <- RunGeneActivity(object, features = "ENSG1", gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "ACTIVITY", layer = "counts"), methods::as(native, "dgCMatrix"))
  }
})

test_that("active fragment hashes and raw barcodes are actually validated", {
  directory <- tempfile("activity_resources_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- make_gene_activity_fixture(directory)
  altered <- object
  altered[["peaks"]]@fragments[[1]]@hash[1] <- strrep("0", 32)
  expect_error(RunGeneActivity(altered, verbose = FALSE), "hash differs")
  altered <- object
  altered[["peaks"]]@fragments[[1]]@cells[1] <- "absent_barcode"
  expect_error(RunGeneActivity(altered, verbose = FALSE), "barcodes are absent")
  empty <- ReadATACData(SeuratObject::GetAssayData(object, layer = "counts"),
    ranges = Signac::granges(object[["peaks"]]), genome = "hg38",
    annotation = Signac::Annotation(object[["peaks"]]), verbose = FALSE
  )
  expect_error(RunGeneActivity(empty, verbose = FALSE), "local indexed fragments")
})

test_that("native single-cell output remains a matrix requirement", {
  directory <- tempfile("activity_single_cell_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  object <- subset(make_gene_activity_fixture(directory), cells = "S_A")
  native <- Signac::GeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
  if (!inherits(native, "sparseMatrix")) {
    expect_error(RunGeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE), "sparse gene-by-cell matrix")
  } else {
    actual <- RunGeneActivity(object, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "ACTIVITY", layer = "counts"), methods::as(native, "dgCMatrix"))
  }
})
