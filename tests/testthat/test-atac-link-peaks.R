make_atac_link_fixture <- function() {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("withr")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("BSgenome")
  counts <- withr::with_seed(83L, Matrix::Matrix(matrix(stats::rpois(60L * 40L, 4) + 1,
    nrow = 60L, dimnames = list(paste0("p", 1:60), paste0("cell", 1:40))
  ), sparse = TRUE))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = rep(c("chr1", "chr2"), each = 30L),
    start = rep(seq(1L, by = 200L, length.out = 30L), 2L),
    end = rep(seq(100L, by = 200L, length.out = 30L), 2L)
  ))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  rna <- counts[c(1, 31), , drop = FALSE]
  rownames(rna) <- c("g1", "g2")
  object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = rna)
  gene_coords <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr2"),
    start = c(201L, 201L), end = c(400L, 400L), strand = "+", gene_id = c("g1", "g2"),
    gene_name = c("g1", "g2"), gene_biotype = "protein_coding", type = "gene"
  ), keep.extra.columns = TRUE)
  GenomeInfoDb::genome(gene_coords) <- "toy"
  dna <- withr::with_seed(89L, Biostrings::DNAStringSet(stats::setNames(
    vapply(seq_len(2), function(chr) {
      paste0(vapply(seq_len(30), function(block) {
        gc <- seq(0.25, 0.75, length.out = 30)[block]
        paste(sample(c("A", "C", "G", "T"), 200,
          replace = TRUE,
          prob = c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
        ), collapse = "")
      }, character(1)), collapse = "")
    }, character(1)), c("chr1", "chr2")
  )))
  object <- Signac::RegionStats(object, assay = "peaks", genome = dna, verbose = FALSE)
  list(object = object, gene_coords = gene_coords, dna = dna)
}

test_that("RunLinkPeaks agrees with native Signac statistics on constructed paired input", {
  x <- make_atac_link_fixture()
  native_input <- x$object
  SeuratObject::DefaultAssay(native_input) <- "RNA"
  native <- withr::with_seed(23L, Signac::LinkPeaks(native_input,
    peak.assay = "peaks", expression.assay = "RNA", expression.slot = "counts",
    gene.coords = x$gene_coords, min.cells = 0L, n_sample = 10L,
    distance = 10000L, score_cutoff = 0.05, pvalue_cutoff = 1, verbose = FALSE
  ))
  set.seed(71L)
  state <- .Random.seed
  before <- serialize(x$object, NULL)
  actual <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0L, n.sample = 10L, distance = 10000L, pvalue.cutoff = 1,
    score.cutoff = 0.05, seed = 23L, verbose = FALSE
  )
  expect_identical(.Random.seed, state)
  expect_identical(serialize(x$object, NULL), before)
  expect_gt(length(Signac::Links(actual[["peaks"]])), 0L)
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]))
  expect_identical(SeuratObject::DefaultAssay(actual), "peaks")
  expect_equal(
    SeuratObject::GetAssayData(actual, assay = "peaks", layer = "counts"),
    SeuratObject::GetAssayData(x$object, assay = "peaks", layer = "counts")
  )
  receipt <- SeuratObject::Misc(actual[["peaks"]], "scop_link_peaks")
  expect_identical(receipt$n_links, length(Signac::Links(native[["peaks"]])))
  expect_identical(receipt$parameters$n_sample, 10L)
  expect_error(RunLinkPeaks(actual, gene.coords = x$gene_coords, verbose = FALSE), "already exist")
})

test_that("empty native LinkPeaks output is recorded without fabricating links", {
  x <- make_atac_link_fixture()
  object <- RunLinkPeaks(x$object,
    gene.coords = x$gene_coords, min.cells = 0,
    score.cutoff = 1, verbose = FALSE
  )
  expect_length(Signac::Links(object[["peaks"]]), 0L)
  expect_identical(SeuratObject::Misc(object[["peaks"]], "scop_link_peaks")$n_links, 0L)
})

test_that("peak-gene linking rejects stale or non-paired biological inputs", {
  x <- make_atac_link_fixture()
  expect_error(RunLinkPeaks(x$object, expression.assay = "peaks", verbose = FALSE), "measured RNA")
  expect_error(RunLinkPeaks(x$object, expression.slot = "missing", verbose = FALSE), "layers must exist")
  expect_error(RunLinkPeaks(x$object, method = "other", verbose = FALSE), "pearson or spearman")
  expect_error(RunLinkPeaks(x$object, n.sample = 1, verbose = FALSE), "n.sample >= 2")
  expect_error(RunLinkPeaks(x$object, genes.use = "missing", verbose = FALSE), "RNA feature IDs")
  bad <- x$object
  SeuratObject::Misc(bad[["RNA"]], "scop_gene_activity") <- list(source_assay = "peaks")
  expect_error(RunLinkPeaks(bad, verbose = FALSE), "cannot substitute")
  bad <- x$object
  bad[["peaks"]]@meta.features$GC.percent[1] <- NA_real_
  expect_error(RunLinkPeaks(bad, gene.coords = x$gene_coords, verbose = FALSE), "RegionStats")
  GenomeInfoDb::genome(x$gene_coords) <- "hg19"
  expect_error(RunLinkPeaks(x$object, gene.coords = x$gene_coords, verbose = FALSE), "genome[[:space:]]+labels")
})

test_that("single-gene measured RNA retains native LinkPeaks results and cell validation", {
  x <- make_atac_link_fixture()
  values <- SeuratObject::GetAssayData(x$object[["RNA"]], layer = "counts")[1L, , drop = FALSE]
  suppressWarnings(x$object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = values))
  native_input <- x$object
  SeuratObject::DefaultAssay(native_input) <- "RNA"
  native <- withr::with_seed(23L, Signac::LinkPeaks(native_input,
    peak.assay = "peaks", expression.assay = "RNA", expression.slot = "counts",
    gene.coords = x$gene_coords, min.cells = 0L, n_sample = 10L,
    distance = 10000L, score_cutoff = 0.05, pvalue_cutoff = 1, verbose = FALSE
  ))
  before <- serialize(x$object, NULL)
  set.seed(71L)
  state <- .Random.seed
  actual <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0L, n.sample = 10L, distance = 10000L, pvalue.cutoff = 1,
    score.cutoff = 0.05, seed = 23L, verbose = FALSE
  )
  expect_gt(length(Signac::Links(native[["peaks"]])), 0L)
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]), tolerance = 0)
  expect_identical(.Random.seed, state)
  expect_identical(serialize(x$object, NULL), before)
  expect_identical(SeuratObject::GetAssayData(actual[["RNA"]], layer = "counts"), values)
  expect_identical(as.vector(colnames(actual)), as.vector(colnames(x$object)))

  paired <- make_atac_link_fixture()$object
  full_values <- SeuratObject::GetAssayData(paired[["RNA"]], layer = "counts")
  incomplete <- SeuratObject::CreateAssay5Object(counts = full_values)
  SeuratObject::LayerData(incomplete, layer = "data") <- log1p(full_values[, -1L, drop = FALSE])
  suppressWarnings(x$object[["RNA"]] <- incomplete)
  expect_error(RunLinkPeaks(x$object, gene.coords = x$gene_coords, verbose = FALSE), "every paired cell")
})

test_that("paired ATAC cells are explicitly selected from a larger RNA assay", {
  x <- make_atac_link_fixture()
  values <- SeuratObject::GetAssayData(x$object, assay = "RNA", layer = "counts")
  larger <- cbind(values, RNAonly = c(1, 1))
  object <- SeuratObject::CreateSeuratObject(larger, assay = "RNA")
  object[["peaks"]] <- x$object[["peaks"]]
  data <- larger
  data[, "RNAonly"] <- NA_real_
  SeuratObject::LayerData(object[["RNA"]], layer = "data") <- data
  before <- serialize(object, NULL)
  actual <- RunLinkPeaks(object,
    assay = "peaks", expression.slot = "data",
    gene.coords = x$gene_coords, min.cells = 0, n.sample = 10, distance = 10000,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expected <- RunLinkPeaks(x$object,
    assay = "peaks", expression.slot = "counts",
    gene.coords = x$gene_coords, min.cells = 0, n.sample = 10, distance = 10000,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(expected[["peaks"]]))
  expect_identical(colnames(actual), colnames(object))
  expect_identical(SeuratObject::DefaultAssay(actual), "RNA")
  expect_identical(serialize(object, NULL), before)
  expect_equal(SeuratObject::GetAssayData(actual, assay = "RNA", layer = "data"), data)
})

test_that("stale accessibility covariates do not change current-cell backgrounds", {
  x <- make_atac_link_fixture()
  fresh <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n.sample = 10, distance = 10000, pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  stale <- x$object
  stale[["peaks"]] <- SeuratObject::AddMetaData(
    stale[["peaks"]],
    data.frame(count = seq(100000, 1, length.out = 60), row.names = rownames(stale[["peaks"]]))
  )
  actual <- RunLinkPeaks(stale,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n.sample = 10, distance = 10000, pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(fresh[["peaks"]]))
  expect_identical(actual[["peaks"]][[]]$count, stale[["peaks"]][[]]$count)
})

test_that("explicit coordinate identities and known lengths are checked", {
  x <- make_atac_link_fixture()
  duplicate <- x$gene_coords
  duplicate$gene_name <- "g1"
  expect_error(RunLinkPeaks(x$object, gene.coords = duplicate, verbose = FALSE), "one interval")
  bad <- x$gene_coords
  bad$gene_id[1] <- NA_character_
  expect_error(RunLinkPeaks(x$object, gene.coords = bad, gene.id = TRUE, verbose = FALSE), "non-empty")
  GenomeInfoDb::seqlengths(x$object[["peaks"]]) <- c(chr1 = 6000, chr2 = 6000)
  GenomeInfoDb::seqlengths(x$gene_coords) <- c(chr1 = 7000, chr2 = 6000)
  expect_error(RunLinkPeaks(x$object, gene.coords = x$gene_coords, verbose = FALSE), "chromosome lengths")
  bad <- make_atac_link_fixture()
  bad$object[["peaks"]]@meta.features$sequence.length[1] <- 99
  expect_error(RunLinkPeaks(bad$object, gene.coords = bad$gene_coords, verbose = FALSE), "RegionStats")
})

test_that("native annotation collapse rejects ambiguous gene identity", {
  x <- make_atac_link_fixture()
  annotation <- x$gene_coords
  annotation$tx_id <- annotation$gene_id
  Signac::Annotation(x$object[["peaks"]]) <- annotation
  actual <- RunLinkPeaks(x$object,
    expression.slot = "counts", min.cells = 0,
    n.sample = 10, distance = 10000, pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  native <- x$object
  SeuratObject::DefaultAssay(native) <- "RNA"
  native <- withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks",
    expression.assay = "RNA", expression.slot = "counts", min.cells = 0,
    n_sample = 10, distance = 10000, pvalue_cutoff = 1, verbose = FALSE
  ))
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]))
  annotation$gene_id <- "same"
  Signac::Annotation(x$object[["peaks"]]) <- annotation
  expect_error(RunLinkPeaks(x$object, verbose = FALSE), "one contig")
  annotation$gene_id[1] <- ""
  Signac::Annotation(x$object[["peaks"]]) <- annotation
  expect_error(RunLinkPeaks(x$object, verbose = FALSE), "non-empty")
})

test_that("distance boundaries and strict cutoffs follow native Signac", {
  x <- make_atac_link_fixture()
  run <- function(distance, min.distance = NULL) {
    RunLinkPeaks(x$object,
      expression.slot = "counts", gene.coords = x$gene_coords, min.cells = 0,
      n.sample = 10, distance = distance, min.distance = min.distance,
      score.cutoff = 0, pvalue.cutoff = 1, seed = 23, verbose = FALSE
    )
  }
  narrow <- run(199)
  boundary <- run(200)
  expect_false("chr1-401-500" %in% Signac::Links(narrow[["peaks"]])$peak)
  expect_true("chr1-401-500" %in% Signac::Links(boundary[["peaks"]])$peak)
  expect_length(Signac::Links(run(200, 199)[["peaks"]]), 0L)
  expect_equal(Signac::Links(run(200, 0)[["peaks"]]), Signac::Links(boundary[["peaks"]]))
  native <- x$object
  SeuratObject::DefaultAssay(native) <- "RNA"
  native <- withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks",
    expression.assay = "RNA", expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n_sample = 10, distance = 200, score_cutoff = 0,
    pvalue_cutoff = 1, verbose = FALSE
  ))
  expect_equal(Signac::Links(boundary[["peaks"]]), Signac::Links(native[["peaks"]]))
  links <- Signac::Links(boundary[["peaks"]])
  expect_equal(links$pvalue, stats::pnorm(-abs(links$zscore)), tolerance = 0)
})

test_that("subset, renamed cells and spearman retain native associations", {
  x <- make_atac_link_fixture()
  object <- subset(x$object, cells = colnames(x$object)[1:30])
  object <- SeuratObject::RenameCells(object, add.cell.id = "paired")
  object[["existing"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(60),
      ncol = 2,
      dimnames = list(colnames(object), c("OLD_1", "OLD_2"))
    ), key = "OLD_", assay = "RNA"
  )
  native <- object
  SeuratObject::DefaultAssay(native) <- "RNA"
  native <- withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks",
    expression.assay = "RNA", expression.slot = "counts", gene.coords = x$gene_coords,
    method = "spearman", min.cells = 0, n_sample = 10, distance = 10000,
    pvalue_cutoff = 1, verbose = FALSE
  ))
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  actual <- RunLinkPeaks(object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    method = "spearman", min.cells = 0, n.sample = 10, distance = 10000,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]))
  expect_identical(actual[["existing"]], object[["existing"]])
  expect_equal(actual[[]], object[[]])
  expect_identical(SeuratObject::Misc(actual[["peaks"]], "scop_link_peaks")$cells, colnames(object))
})

test_that("constant genes and native detection filtering are not reinterpreted", {
  x <- make_atac_link_fixture()
  values <- SeuratObject::GetAssayData(x$object, assay = "RNA", layer = "counts")
  values[1, ] <- 1
  x$object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = values)
  native <- x$object
  SeuratObject::DefaultAssay(native) <- "RNA"
  expected <- tryCatch(withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks", expression.assay = "RNA", expression.slot = "counts",
    gene.coords = x$gene_coords, min.cells = 0, n_sample = 10, distance = 10000,
    pvalue_cutoff = 1, verbose = FALSE
  )), error = identity)
  actual <- tryCatch(RunLinkPeaks(x$object,
    expression.slot = "counts",
    gene.coords = x$gene_coords, min.cells = 0, n.sample = 10, distance = 10000,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  ), error = identity)
  if (inherits(expected, "error")) {
    expect_s3_class(actual, "error")
    expect_identical(conditionMessage(actual), conditionMessage(expected))
  } else {
    expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(expected[["peaks"]]))
  }
  expect_error(RunLinkPeaks(x$object, min.cells = 40, verbose = FALSE), "min.cells < cell count")
  missing_pair <- x$object
  expect_warning(missing_pair[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = values[, -1, drop = FALSE]), "Assay|Different")
  expect_error(RunLinkPeaks(missing_pair, verbose = FALSE), "cover every selected ATAC cell")
})

test_that("strict detection support and gene-ID mode follow native data", {
  x <- make_atac_link_fixture()
  peaks <- SeuratObject::GetAssayData(x$object, assay = "peaks", layer = "counts")
  peaks[1, 11:40] <- 0
  SeuratObject::LayerData(x$object[["peaks"]], layer = "counts") <- peaks
  rna <- peaks[c(1, 31), , drop = FALSE]
  rownames(rna) <- c("g1", "g2")
  x$object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = rna)
  x$gene_coords$gene_name <- "shared-display-name"
  outputs <- lapply(c(9, 10), function(minimum) {
    native <- x$object
    SeuratObject::DefaultAssay(native) <- "RNA"
    native <- withr::with_seed(23, Signac::LinkPeaks(native,
      peak.assay = "peaks",
      expression.assay = "RNA", expression.slot = "counts", gene.coords = x$gene_coords,
      gene.id = TRUE, min.cells = minimum, n_sample = 10, distance = 10000,
      pvalue_cutoff = 1, verbose = FALSE
    ))
    actual <- RunLinkPeaks(x$object,
      expression.slot = "counts", gene.coords = x$gene_coords,
      gene.id = TRUE, min.cells = minimum, n.sample = 10, distance = 10000,
      pvalue.cutoff = 1, seed = 23, verbose = FALSE
    )
    expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]))
    Signac::Links(actual[["peaks"]])
  })
  expect_true("g1" %in% outputs[[1]]$gene)
  expect_false("g1" %in% outputs[[2]]$gene)
})

test_that("native assay subsets preserve pairing for reversed cell requests", {
  x <- make_atac_link_fixture()
  reordered <- x$object
  reordered@assays$RNA <- subset(reordered[["RNA"]], cells = rev(colnames(reordered)))
  expect_true(methods::validObject(reordered, test = TRUE))
  expect_setequal(colnames(reordered[["RNA"]]), colnames(reordered[["peaks"]]))
  # Native Assay subset retains its stored cell order despite a reversed request.
  expect_identical(colnames(reordered[["RNA"]]), colnames(reordered[["peaks"]]))
  expected <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n.sample = 10, distance = 10000, pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  actual <- RunLinkPeaks(reordered,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n.sample = 10, distance = 10000, pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(expected[["peaks"]]))
})

test_that("negative-strand TSS distance is delegated to native Signac", {
  x <- make_atac_link_fixture()
  GenomicRanges::strand(x$gene_coords)[2] <- "-"
  actual <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n.sample = 10, distance = 200, score.cutoff = 0,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  native <- x$object
  SeuratObject::DefaultAssay(native) <- "RNA"
  native <- withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks",
    expression.assay = "RNA", expression.slot = "counts", gene.coords = x$gene_coords,
    min.cells = 0, n_sample = 10, distance = 200, score_cutoff = 0,
    pvalue_cutoff = 1, verbose = FALSE
  ))
  links <- Signac::Links(actual[["peaks"]])
  expect_equal(links, Signac::Links(native[["peaks"]]))
  expect_false("chr2-1-100" %in% links$peak[links$gene == "g2"])
  expect_true("chr2-401-500" %in% links$peak[links$gene == "g2"])
})

test_that("explicit stable gene IDs do not require unused display symbols", {
  x <- make_atac_link_fixture()
  x$gene_coords$gene_name <- NULL
  native <- x$object
  SeuratObject::DefaultAssay(native) <- "RNA"
  native <- withr::with_seed(23, Signac::LinkPeaks(native,
    peak.assay = "peaks",
    expression.assay = "RNA", expression.slot = "counts", gene.coords = x$gene_coords,
    gene.id = TRUE, min.cells = 0, n_sample = 10, distance = 10000,
    pvalue_cutoff = 1, verbose = FALSE
  ))
  actual <- RunLinkPeaks(x$object,
    expression.slot = "counts", gene.coords = x$gene_coords,
    gene.id = TRUE, min.cells = 0, n.sample = 10, distance = 10000,
    pvalue.cutoff = 1, seed = 23, verbose = FALSE
  )
  expect_gt(length(Signac::Links(actual[["peaks"]])), 0L)
  expect_equal(Signac::Links(actual[["peaks"]]), Signac::Links(native[["peaks"]]))
})
