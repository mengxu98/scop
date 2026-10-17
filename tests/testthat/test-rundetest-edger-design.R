make_edger_design_fixture <- function() {
  set.seed(20261010)
  samples <- paste0("D", rep(seq_len(4), each = 2), "_", rep(c("A", "B"), 4))
  metadata <- data.frame(
    condition = factor(rep(c("A", "B"), 4), levels = c("A", "B")),
    donor = factor(rep(paste0("D", seq_len(4)), each = 2)),
    quality = c(8, 6, 5, 9, 7, 5, 6, 8), row.names = samples
  )
  mu <- matrix(60, 301, 8)
  mu[1:30, metadata$condition == "B"] <- 240
  mu[31:60, metadata$condition == "A"] <- 240
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 12), 301, 8,
    dimnames = list(paste0("peak", seq_len(301)), samples)
  )
  counts[301, ] <- 0
  list(counts = counts, metadata = metadata)
}

edger_design_contrast <- function(design, metadata) {
  cols <- colnames(stats::model.matrix(design, metadata))
  out <- stats::setNames(rep(0, length(cols)), cols)
  out["conditionB"] <- 1
  out
}

run_edger_design_fixture <- function(
  fixture, design = ~ donor + condition,
  contrast = edger_design_contrast(design, fixture$metadata), ...
) {
  scop:::RunDEtest_edgeR(fixture$counts, fixture$metadata$condition,
    condition1 = "A", condition2 = "B", sample_data = fixture$metadata,
    design = design, contrast = contrast, ...
  )
}

native_edger_design_reference <- function(fixture, design, contrast) {
  mm <- stats::model.matrix(design, fixture$metadata)
  y <- edgeR::DGEList(fixture$counts)
  keep <- edgeR::filterByExpr(y, design = mm)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design = mm)
  fit <- edgeR::glmQLFit(y, design = mm, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, contrast = unname(contrast[colnames(mm)]))
  edgeR::topTags(test, n = Inf, sort.by = "none")$table
}

for (model in c("independent", "paired", "covariate")) {
  test_that(paste("explicit edgeR", model, "models match the installed native backend"), {
    skip_if_not_installed("edgeR")
    skip_if_not_installed("statmod")
    fixture <- make_edger_design_fixture()
    design <- switch(model,
      independent = ~condition,
      paired = ~ donor + condition,
      covariate = ~ donor + quality + condition
    )
    contrast <- edger_design_contrast(design, fixture$metadata)
    actual <- run_edger_design_fixture(fixture, design, contrast,
      only.pos = FALSE, p.adjust.method = "BH"
    )
    expected <- native_edger_design_reference(fixture, design, contrast)
    receipt <- attr(actual, "edgeR_model")
    expect_identical(rownames(actual), rownames(expected))
    expect_equal(actual$avg_log2FC, expected$logFC, tolerance = 1e-12)
    expect_equal(actual$p_val, expected$PValue, tolerance = 1e-12)
    expect_equal(actual$F, expected$F, tolerance = 1e-12)
    expect_equal(actual$p_val_BH, expected$FDR, tolerance = 1e-12)
    expect_equal(actual$p_val_adj, expected$FDR, tolerance = 1e-12)
    expect_identical(receipt$versions[["edgeR"]], as.character(utils::packageVersion("edgeR")))
    expect_equal(receipt$design, stats::model.matrix(design, fixture$metadata))
    expect_identical(receipt$tested_features, rownames(expected))
    expect_identical(receipt$filtered_features, "peak301")
    expect_identical(receipt$counts, fixture$counts)
    expect_equal(receipt$residual_df, 8L - ncol(receipt$design))
  })
}

test_that("contrast names align by name and reversed contrasts reverse direction", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  fixture <- make_edger_design_fixture()
  contrast <- edger_design_contrast(~ donor + condition, fixture$metadata)
  a <- run_edger_design_fixture(fixture, contrast = rev(contrast), only.pos = FALSE)
  b <- run_edger_design_fixture(fixture, contrast = -contrast, only.pos = FALSE)
  expect_equal(a$avg_log2FC, -b$avg_log2FC, tolerance = 1e-12)
  expect_equal(a$p_val, b$p_val, tolerance = 1e-12)
  expect_identical(attr(a, "edgeR_model")$positive_direction, "1 * conditionB")
  expect_identical(attr(b, "edgeR_model")$positive_direction, "-1 * conditionB")
})

test_that("all tested peaks retain BH correction before display filters", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  fixture <- make_edger_design_fixture()
  displayed <- run_edger_design_fixture(fixture,
    only.pos = TRUE, logfc.threshold = 1,
    p.adjust.method = "bonferroni"
  )
  full <- attr(displayed, "edgeR_model")$full_results
  expect_gt(nrow(displayed), 0)
  expect_lt(nrow(displayed), nrow(full))
  expect_equal(full$p_val_BH, stats::p.adjust(full$p_val, "BH"))
  expect_equal(full$p_val_adj, stats::p.adjust(full$p_val, "bonferroni"))
  expect_equal(displayed$p_val_BH, full[rownames(displayed), "p_val_BH"])
  empty_display <- run_edger_design_fixture(fixture, only.pos = TRUE, logfc.threshold = Inf)
  expect_equal(nrow(empty_display), 0L)
  expect_equal(nrow(attr(empty_display, "edgeR_model")$full_results), nrow(full))
})

test_that("explicit models reject confounding, zero df, and missing true replicates", {
  fixture <- make_edger_design_fixture()
  fixture$metadata$batch <- fixture$metadata$condition
  expect_error(run_edger_design_fixture(fixture, design = ~ batch + condition), "not full rank")
  fixture$metadata$sample <- factor(rownames(fixture$metadata))
  saturated <- stats::setNames(rep(1, 8), colnames(stats::model.matrix(~ 0 + sample, fixture$metadata)))
  expect_error(run_edger_design_fixture(fixture, design = ~ 0 + sample, contrast = saturated), "no residual degrees")
  fixture$counts <- fixture$counts[, 1:2, drop = FALSE]
  fixture$metadata <- fixture$metadata[1:2, , drop = FALSE]
  expect_error(run_edger_design_fixture(fixture, design = ~condition), "at least 2 biological samples")
})

test_that("invalid design inputs fail before reaching an optional backend", {
  fixture <- make_edger_design_fixture()
  expect_error(run_edger_design_fixture(fixture, contrast = c(conditionB = 1)), "match every design column")
  expect_error(run_edger_design_fixture(fixture, contrast = 1), "named numeric")
  expect_error(run_edger_design_fixture(fixture, contrast = c(conditionB = 0)), "nonzero")
  expect_error(run_edger_design_fixture(fixture, design = ~missing, contrast = c(missing = 1)), "Missing design metadata")
  expect_error(scop:::RunDEtestCheckDesign(~condition, c(conditionB = 1), "wilcox"), "only for sample-level edgeR")
  expect_error(scop:::RunDEtestCheckDesign(NULL, c(conditionB = 1)), "one-sided formula")
  expect_error(scop:::RunDEtestCheckDesign(~., c(conditionB = 1)), "explicit sample metadata")
  fixture$metadata$quality[1] <- NA_real_
  expect_error(run_edger_design_fixture(fixture, design = ~ quality + condition), "missing values")
  fixture$counts[1, 1] <- 0.5
  expect_error(run_edger_design_fixture(fixture), "integer raw counts")
})

test_that("aggregate metadata is constant within each true sample", {
  metadata <- data.frame(
    sample = c("D1_A", "D1_A", "D1_B", "D1_B"),
    donor = "D1", condition = c("A", "A", "B", "B"), quality = c(1, 2, 3, 3)
  )
  expect_error(
    scop:::RunDEtestModelSampleData(metadata, "sample", "condition", ~ quality + condition),
    "quality.*constant within each aggregate sample"
  )
  expect_error(
    scop:::RunDEtestModelSampleData(metadata, "donor", "condition", ~ donor + condition),
    "distinct donor-by-cell-type sample IDs"
  )
  metadata$quality[2] <- 1
  actual <- scop:::RunDEtestModelSampleData(metadata, "sample", "condition", ~ quality + condition)
  expect_identical(rownames(actual), c("D1_A", "D1_B"))
})

make_edger_paired_atac <- function(uneven = FALSE) {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  fixture <- make_edger_design_fixture()
  n_cells <- if (uneven) c(12, 10, 14, 9, 11, 0, 8, 3) else rep(10, 8)
  sample_index <- rep(seq_len(8), n_cells)
  metadata <- fixture$metadata[sample_index, , drop = FALSE]
  metadata$sample <- rownames(fixture$metadata)[sample_index]
  rownames(metadata) <- paste0("cell", seq_len(nrow(metadata)))
  set.seed(17)
  counts <- matrix(stats::rpois(300 * nrow(metadata), lambda = 10), 300,
    dimnames = list(paste0("input", seq_len(300)), rownames(metadata))
  )
  counts[1:30, metadata$condition == "B"] <- counts[1:30, metadata$condition == "B"] * 4
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = seq(1, by = 200, length.out = 300), end = seq(100, by = 200, length.out = 300)
  ))
  names(ranges) <- rownames(counts)
  ReadATACData(Matrix::Matrix(counts, sparse = TRUE),
    ranges = ranges, metadata = metadata,
    genome = "hg19", verbose = FALSE
  )
}

run_edger_paired_atac <- function(object, ...) {
  contrast <- c("(Intercept)" = 0, donorD2 = 0, donorD3 = 0, donorD4 = 0, conditionB = 1)
  RunDEtest(object,
    test.use = "edgeR", assay = "peaks", layer = "counts", feature_type = "peak",
    sample_col = "sample", condition_col = "condition", group1 = "A", group2 = "B",
    design = ~ donor + condition, contrast = contrast, min.cells.sample = 5,
    only.pos = FALSE, fc.threshold = 1, p.adjust.method = "BH", verbose = FALSE, ...
  )
}

test_that("paired ATAC models retain only observed samples and report genomic coordinates", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  object <- make_edger_paired_atac(uneven = TRUE)
  out <- suppressWarnings(run_edger_paired_atac(object))
  result <- out@tools$DEtest_pseudobulk
  model <- result$edgeR_models$All
  expect_identical(colnames(model$counts), c("D1_A", "D1_B", "D2_A", "D2_B", "D3_A", "D4_A"))
  expect_false("D3_B" %in% result$sample_inventory$sample)
  expect_false(result$sample_inventory$kept[result$sample_inventory$sample == "D4_B"])
  expect_equal(as.integer(model$sample_counts), c(4L, 2L))
  expect_equal(model$residual_df, 1L)
  expect_equal(nrow(model$cell_samples), 64L)
  expect_true(all(model$full_results$coordinate_system == "1-based closed"))
  expect_true(all(model$full_results$genome == "hg19"))
  expect_equal(model$full_results$start, seq(1, by = 200, length.out = 300))
  expect_identical(model$input$layer, "counts")
  expect_identical(model$sample_unit$biological_independence, "not verified from metadata")
  expect_identical(result$AllMarkers_edgeR$group2, rep("1 * conditionB", nrow(result$AllMarkers_edgeR)))
})

test_that("paired ATAC rejects donor-only aggregation and variable sample covariates", {
  object <- make_edger_paired_atac()
  object$sample <- as.character(object$donor)
  expect_error(run_edger_paired_atac(object), "distinct donor-by-cell-type sample IDs")
  object <- make_edger_paired_atac()
  object$donor[1] <- "D2"
  expect_error(run_edger_paired_atac(object), "donor.*constant within each aggregate sample")
})

test_that("bulk explicit models persist full tests and propagate errors", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  fixture <- make_edger_design_fixture()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fixture$counts),
    colData = fixture$metadata
  )
  contrast <- edger_design_contrast(~ donor + condition, fixture$metadata)
  out <- RunDEtest(se,
    condition_col = "condition", group1 = "A", group2 = "B",
    design = ~ donor + condition, contrast = contrast, fc.threshold = Inf,
    verbose = FALSE
  )
  store <- S4Vectors::metadata(out)$DEtest
  expect_identical(store$status$status, "success")
  expect_equal(nrow(store$results), 0L)
  model <- store$methods$de_edgeR_qlf$details$groups$all$comparisons$B_vs_A$details$edgeR_model
  expect_equal(nrow(model$full_results), 300L)
  expect_equal(model$full_results$p_val_BH, stats::p.adjust(model$full_results$p_val, "BH"))
  expect_error(RunDEtest(se,
    condition_col = "condition", design = ~ donor + condition,
    contrast = c(conditionB = 1), verbose = FALSE
  ), "match every design column")
})

test_that("legacy edgeR defaults remain compatible with the independent native model", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  fixture <- make_edger_design_fixture()
  legacy <- scop:::RunDEtest_edgeR(fixture$counts, fixture$metadata$condition,
    condition1 = "A", condition2 = "B", only.pos = FALSE
  )
  model <- run_edger_design_fixture(fixture, design = ~condition, only.pos = FALSE)
  expect_null(attr(legacy, "edgeR_model"))
  expect_equal(legacy$avg_log2FC, model$avg_log2FC, tolerance = 1e-12)
  expect_equal(legacy$p_val, model$p_val, tolerance = 1e-12)
  expect_equal(legacy$p_val_adj, stats::p.adjust(legacy$p_val, "bonferroni"))
})

test_that("per-group sample inventories never invent absent combinations", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  object <- make_edger_paired_atac()
  object$cluster <- rep(c("first", "second"), length.out = ncol(object))
  object$cluster[object$sample == "D4_B"] <- "first"
  out <- RunDEtest(object,
    test.use = "edgeR", assay = "peaks", layer = "counts",
    feature_type = "peak", group.by = "cluster", sample_col = "sample",
    condition_col = "condition", group1 = "A", group2 = "B", min.cells.sample = 5,
    design = ~condition, contrast = c("(Intercept)" = 0, conditionB = 1),
    only.pos = FALSE, fc.threshold = 1, verbose = FALSE
  )
  models <- out@tools$DEtest_cluster$edgeR_models
  expect_equal(ncol(models$first$counts), 8L)
  expect_equal(ncol(models$second$counts), 7L)
  expect_false("D4_B" %in% colnames(models$second$counts))
  expect_false("D4_B" %in% models$second$sample_inventory$sample)
})

test_that("bulk peak coordinates follow the requested feature order", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  skip_if_not_installed("GenomicRanges")
  fixture <- make_edger_design_fixture()
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = seq(1, by = 200, length.out = 301), end = seq(100, by = 200, length.out = 301)
  ))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fixture$counts),
    colData = fixture$metadata, rowRanges = ranges
  )
  features <- rev(rownames(fixture$counts))
  out <- RunDEtest(se,
    condition_col = "condition", group1 = "A", group2 = "B",
    feature_type = "peak", features = features, design = ~condition,
    contrast = c("(Intercept)" = 0, conditionB = 1), only.pos = FALSE,
    fc.threshold = 1, verbose = FALSE
  )
  store <- S4Vectors::metadata(out)$DEtest
  expected_features <- features[features != "peak301"]
  expect_identical(store$results$gene, expected_features)
  expect_equal(store$results$start, rev(seq(1, by = 200, length.out = 300)))
  expect_true(all(store$results$coordinate_system == "1-based closed"))
})

test_that("optional edgeR checks never request installation", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  testthat::local_mocked_bindings(check_r = function(package, install, verbose) {
    expect_identical(package, "edgeR")
    expect_false(install)
    FALSE
  }, .package = "scop")
  expect_error(scop:::RunDEtestEdgeRBackend(), "[Ee]dgeR is required")
})

test_that("explicit models reject sample metadata contradictions across groups", {
  object <- make_edger_paired_atac()
  object$cluster <- rep(c("first", "second"), length.out = ncol(object))
  args <- list(
    test.use = "edgeR", assay = "peaks", layer = "counts", feature_type = "peak",
    group.by = "cluster", sample_col = "sample", condition_col = "condition",
    group1 = "A", group2 = "B", min.cells.sample = 5, only.pos = FALSE,
    fc.threshold = 1, verbose = FALSE
  )
  conflicting <- object$sample == "D1_A" & object$cluster == "second"
  object$condition[conflicting] <- "B"
  expect_error(
    do.call(RunDEtest, c(
      list(object = object), args,
      list(design = ~condition, contrast = c("(Intercept)" = 0, conditionB = 1))
    )),
    "constant[[:space:]]+within[[:space:]]+each[[:space:]]+aggregate[[:space:]]+sample"
  )

  object$condition[conflicting] <- "A"
  object$quality[conflicting] <- object$quality[conflicting] + 1
  expect_error(
    do.call(RunDEtest, c(
      list(object = object), args,
      list(
        design = ~ quality + condition,
        contrast = c("(Intercept)" = 0, quality = 0, conditionB = 1)
      )
    )),
    "quality.*constant within each aggregate sample"
  )
})

test_that("nonparticipating conditions and groups do not need design metadata", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  object <- make_edger_paired_atac()
  object$cluster <- rep(c("first", "second"), length.out = ncol(object))
  excluded <- object$sample == "D4_B"
  object$condition <- as.character(object$condition)
  object$condition[excluded] <- "unused"
  object$cluster[excluded] <- "excluded_group"
  object$quality[excluded] <- NA_real_
  object$sample[excluded] <- ""
  out <- RunDEtest(object,
    test.use = "edgeR", assay = "peaks", layer = "counts",
    feature_type = "peak", group.by = "cluster", sample_col = "sample",
    condition_col = "condition", group1 = "A", group2 = "B", min.cells.sample = 5,
    design = ~ quality + condition, contrast = c("(Intercept)" = 0, quality = 0, conditionB = 1),
    only.pos = FALSE, fc.threshold = 1, verbose = FALSE
  )
  models <- out@tools$DEtest_cluster$edgeR_models
  expect_setequal(names(models), c("first", "second"))
  expect_true(all(vapply(models, function(x) ncol(x$counts) == 7L, logical(1))))
  expect_false(any(vapply(models, function(x) anyNA(x$sample_data), logical(1))))
})

test_that("peak genome labels map by seqlevel instead of recycling across peaks", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")
  ids <- c("p1", "p2", "p3", "p4", "p5")
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr2", "chr1", "chr2", "chr1", "chr2"),
    start = seq(1, by = 200, length.out = 5),
    end = seq(100, by = 200, length.out = 5)
  ))
  GenomeInfoDb::genome(ranges) <- "hg19"
  counts <- matrix(1L, 5, 2, dimnames = list(ids, c("s1", "s2")))
  object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), rowRanges = ranges
  )
  features <- rev(ids)
  coordinates <- scop:::RunDEtestFeatureCoordinates(object, "counts", features, "peak")
  expect_identical(rownames(coordinates), features)
  expect_identical(coordinates$genome, rep("hg19", 5))
  expect_identical(coordinates$chromosome, rev(as.character(GenomicRanges::seqnames(ranges))))
  expect_equal(coordinates$start, rev(GenomicRanges::start(ranges)))

  # Peak filtering can leave only one chromosome while retaining unused
  # seqlevels in Seqinfo; genome labels must still have exactly one per row.
  features <- c("p5", "p3", "p1")
  coordinates <- scop:::RunDEtestFeatureCoordinates(object, "counts", features, "peak")
  expect_identical(coordinates$genome, rep("hg19", 3))
  expect_identical(coordinates$chromosome, rep("chr2", 3))
})

test_that("explicit edgeR formulas reject offsets rather than silently discarding them", {
  contrast <- c("(Intercept)" = 0, conditionB = 1)
  expect_error(
    scop:::RunDEtestCheckDesign(~ condition + offset(log(exposure)), contrast),
    "[Oo]ffset\\(\\) terms are not supported"
  )
  expect_error(
    scop:::RunDEtestCheckDesign(~ condition + stats::offset(log(exposure)), contrast),
    "[Oo]ffset\\(\\) terms are not supported"
  )
  # A metadata column named offset is still an ordinary covariate, not offset().
  expect_true(scop:::RunDEtestCheckDesign(
    ~ condition + offset,
    c("(Intercept)" = 0, conditionB = 1, offset = 0)
  ))
})

test_that("bulk explicit models fail when all grouping metadata is missing", {
  skip_if_not_installed("SummarizedExperiment")
  fixture <- make_edger_design_fixture()
  fixture$metadata$cluster <- NA_character_
  object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fixture$counts), colData = fixture$metadata
  )
  expect_error(
    RunDEtest(object,
      condition_col = "condition", group.by = "cluster",
      group1 = "A", group2 = "B", design = ~condition,
      contrast = c("(Intercept)" = 0, conditionB = 1), verbose = FALSE
    ),
    "No bulk sample groups participate in the explicit edgeR comparison"
  )
})


test_that("explicit design errors preserve literal braces in metadata identifiers", {
  metadata <- data.frame(sample = c("s1", "s2"), condition = c("A", "B"))
  expect_error(
    scop:::RunDEtestModelSampleData(metadata, "sample", "condition", ~`bad{literal}`),
    "Missing design metadata columns: bad\\{literal\\}"
  )
})

test_that("explicit Seurat models require a selected pair among three conditions", {
  object <- make_edger_paired_atac()
  object$condition <- as.character(object$condition)
  object$condition[object$sample == "D4_B"] <- "C"
  expect_error(RunDEtest(object,
    assay = "peaks", test.use = "edgeR", layer = "counts",
    sample_col = "sample", condition_col = "condition", min.cells.sample = 1,
    design = ~condition, contrast = c("(Intercept)" = 0, conditionB = 1),
    verbose = FALSE
  ), "one selected condition pair")
})

test_that("new design parameters do not capture existing positional dots", {
  for (method in c("RunDEtest.Seurat", "RunDEtest.SummarizedExperiment")) {
    args <- names(formals(get(method, asNamespace("scop"))))
    expect_true(match("design", args) > match("...", args))
    expect_true(match("contrast", args) > match("...", args))
  }
})

test_that("partial assay and renamed cells retain true sample aggregation", {
  object <- make_edger_paired_atac()
  peaks <- object[["peaks"]]
  rna_counts <- Matrix::Matrix(matrix(3, 3, ncol(object) + 1L,
    dimnames = list(paste0("gene", 1:3), c(colnames(object), "RNAonly"))
  ), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna_counts, assay = "RNA")
  partial[["peaks"]] <- peaks
  meta <- object[[]]
  partial <- SeuratObject::AddMetaData(partial, meta)
  partial$condition <- as.character(partial$condition)
  partial$condition[is.na(partial$condition)] <- "RNA_condition"
  partial <- SeuratObject::RenameCells(partial, add.cell.id = "renamed")
  SeuratObject::DefaultAssay(partial) <- "RNA"
  partial[["existing"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(ncol(partial) * 2),
      ncol = 2,
      dimnames = list(colnames(partial), c("OLD_1", "OLD_2"))
    ), key = "OLD_", assay = "RNA"
  )
  before_rna <- SeuratObject::GetAssayData(partial, assay = "RNA", layer = "counts")
  before_peaks <- SeuratObject::GetAssayData(partial, assay = "peaks", layer = "counts")
  out <- run_edger_paired_atac(partial)
  model <- out@tools$DEtest_pseudobulk$edgeR_models$All
  expected <- sapply(unique(meta$sample), function(sample) {
    Matrix::rowSums(before_peaks[, paste0("renamed_", rownames(meta)[meta$sample == sample]), drop = FALSE])
  })
  expect_equal(model$counts, expected[, colnames(model$counts), drop = FALSE])
  expect_identical(model$cell_samples$cell, colnames(before_peaks))
  expect_identical(colnames(out), colnames(partial))
  expect_identical(SeuratObject::DefaultAssay(out), "RNA")
  expect_equal(SeuratObject::GetAssayData(out, assay = "RNA", layer = "counts"), before_rna)
  expect_equal(SeuratObject::GetAssayData(out, assay = "peaks", layer = "counts"), before_peaks)
  expect_equal(out[[]], partial[[]])
  expect_identical(out[["existing"]], partial[["existing"]])
  # RNA-only condition must not make the selected ATAC assay a three-condition design.
  auto <- RunDEtest(partial,
    assay = "peaks", test.use = "edgeR", layer = "counts",
    sample_col = "sample", condition_col = "condition", min.cells.sample = 5,
    design = ~condition, contrast = c("(Intercept)" = 0, conditionB = 1), verbose = FALSE
  )
  expect_equal(as.integer(auto@tools$DEtest_pseudobulk$edgeR_models$All$sample_counts), c(4L, 4L))
})

test_that("explicit RNA pseudobulk uses the same native sample model", {
  atac <- make_edger_paired_atac()
  counts <- SeuratObject::GetAssayData(atac, assay = "peaks", layer = "counts")
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  object <- SeuratObject::CreateSeuratObject(counts, assay = "RNA", meta.data = atac[[]])
  out <- RunDEtest(object,
    assay = "RNA", test.use = "edgeR", layer = "counts",
    sample_col = "sample", condition_col = "condition", min.cells.sample = 5,
    design = ~ donor + condition,
    contrast = c("(Intercept)" = 0, donorD2 = 0, donorD3 = 0, donorD4 = 0, conditionB = 1),
    only.pos = FALSE, fc.threshold = 1, verbose = FALSE
  )
  model <- out@tools$DEtest_pseudobulk$edgeR_models$All
  expected <- native_edger_design_reference(
    list(counts = model$counts, metadata = model$sample_data),
    ~ donor + condition, model$contrast
  )
  expect_equal(model$full_results$p_val, expected$PValue, tolerance = 1e-12)
  expect_equal(model$full_results$avg_log2FC, expected$logFC, tolerance = 1e-12)
  expect_true(all(model$full_results$feature_type == "gene"))
  expect_false("start" %in% names(model$full_results))
  expect_equal(SeuratObject::GetAssayData(out, layer = "counts"), counts)
})

test_that("native filtering preserves empty tests and rejects post-filter zero libraries", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("statmod")
  fixture <- make_edger_design_fixture()
  fixture$counts[, ] <- 1
  result <- run_edger_design_fixture(fixture)
  expect_equal(nrow(result), 0L)
  expect_length(attr(result, "edgeR_model")$tested_features, 0L)
  expect_identical(attr(result, "edgeR_model")$filtered_features, rownames(fixture$counts))
  fixture$counts[, ] <- 50
  fixture$counts[, 1] <- 0
  fixture$counts[301, ] <- 0
  fixture$counts[301, 1] <- 100
  expect_error(run_edger_design_fixture(fixture), "zero library size after filterByExpr")
})

test_that("legacy SCT markers retain native results without explicit design", {
  skip_if_not_installed("sctransform")
  set.seed(7)
  counts <- Matrix::Matrix(matrix(rpois(80 * 40, 5), 80,
    dimnames = list(paste0("gene", 1:80), paste0("cell", 1:40))
  ), sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts)
  object <- Seurat::SCTransform(object,
    vst.flavor = "v1", variable.features.n = 30,
    return.only.var.genes = FALSE, verbose = FALSE
  )
  SeuratObject::Idents(object) <- rep(c("A", "B"), each = 20)
  expected <- Seurat::FindMarkers(object,
    assay = "SCT", ident.1 = "A", ident.2 = "B",
    features = rownames(object), min.pct = 0, logfc.threshold = 0, only.pos = FALSE,
    verbose = FALSE
  )
  out <- RunDEtest(object,
    assay = "SCT", cells1 = colnames(object)[1:20],
    cells2 = colnames(object)[21:40], features = rownames(object), min.pct = 0,
    fc.threshold = 1, only.pos = FALSE, verbose = FALSE
  )
  actual <- out@tools$DEtest_custom$AllMarkers_wilcox
  expected <- expected[actual$gene, , drop = FALSE]
  expect_equal(actual$p_val, expected$p_val, tolerance = 1e-12)
  expect_equal(actual$avg_log2FC, expected$avg_log2FC, tolerance = 1e-12)
  expect_equal(SeuratObject::GetAssayData(out, assay = "RNA", layer = "counts"), counts)
  expect_identical(out[["SCT"]]@SCTModel.list, object[["SCT"]]@SCTModel.list)
})
