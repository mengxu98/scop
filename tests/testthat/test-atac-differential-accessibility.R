make_test_atac_da_object <- function(normalize = TRUE, add_samples = FALSE) {
  skip_if_not_installed("Signac")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(20260809)
  n_features <- 60L
  n_cells <- 32L
  group <- factor(rep(c("A", "B"), each = n_cells / 2L))
  counts <- matrix(
    stats::rbinom(n_features * n_cells, size = 1, prob = 0.2),
    nrow = n_features,
    ncol = n_cells
  )
  counts[1:10, group == "A"] <- stats::rbinom(10 * 16, 1, 0.8)
  counts[1:10, group == "B"] <- stats::rbinom(10 * 16, 1, 0.1)
  counts[11:20, group == "A"] <- stats::rbinom(10 * 16, 1, 0.1)
  counts[11:20, group == "B"] <- stats::rbinom(10 * 16, 1, 0.8)
  counts[21:25, group == "A"] <- 0
  counts[21:25, group == "B"] <- stats::rbinom(5 * 16, 1, 0.7)
  counts[26:30, group == "A"] <- stats::rbinom(5 * 16, 1, 0.7)
  counts[26:30, group == "B"] <- 0
  counts <- counts * matrix(
    sample.int(3L, n_features * n_cells, replace = TRUE),
    nrow = n_features,
    ncol = n_cells
  )
  if (isTRUE(add_samples)) {
    counts <- counts * 10L
  }

  starts <- seq.int(1000L, by = 200L, length.out = n_features)
  rownames(counts) <- paste0("chr1-", starts, "-", starts + 99L)
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  chromatin <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    sep = c("-", "-")
  )
  srt <- Seurat::CreateSeuratObject(counts = chromatin, assay = "peaks")
  srt$group <- group
  SeuratObject::Idents(srt) <- group

  if (isTRUE(add_samples)) {
    srt$sample <- rep(paste0("sample", seq_len(8L)), each = 4L)
    srt$condition <- rep(c("ctrl", "case"), each = n_cells / 2L)
  }
  if (isTRUE(normalize)) {
    srt <- Signac::RunTFIDF(srt, assay = "peaks", verbose = FALSE)
  }
  srt
}

seurat_find_markers_reference <- function(object, ...) {
  getFromNamespace("FindMarkers.Seurat", "Seurat")(object = object, ...)
}

test_that("ChromatinAssay Wilcoxon markers bypass RNA sparse summaries", {
  srt <- make_test_atac_da_object()
  features <- rownames(srt)[1:30]
  args <- list(
    assay = "peaks",
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    test.use = "wilcox",
    layer = "data",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  reference_args <- args
  reference_args$slot <- reference_args$layer
  reference_args$layer <- NULL

  expected <- do.call(
    seurat_find_markers_reference,
    c(list(object = srt), reference_args)
  )
  actual <- do.call(scop::FindMarkers, c(list(object = srt), args))

  expect_setequal(rownames(actual), rownames(expected))
  actual <- actual[rownames(expected), colnames(expected), drop = FALSE]
  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("ChromatinAssay FindAllMarkers bypasses RNA all-in-one summaries", {
  srt <- make_test_atac_da_object()
  testthat::local_mocked_bindings(
    parallel_all_in_one_dgc = function(...) {
      stop("RNA all-in-one summary reached")
    },
    .package = "scop"
  )

  expect_no_error(markers <- scop::FindAllMarkers(
    srt,
    assay = "peaks",
    features = rownames(srt)[1:20],
    test.use = "wilcox",
    slot = "data",
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = Inf,
    only.pos = FALSE,
    verbose = FALSE
  ))
  expect_s3_class(markers, "data.frame")
  expect_gt(nrow(markers), 0)
})

test_that("RunDEtest LR supports cell vectors and depth covariates", {
  srt <- make_test_atac_da_object()
  features <- rownames(srt)[1:20]
  expected <- suppressWarnings(seurat_find_markers_reference(
    object = srt,
    assay = "peaks",
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    test.use = "LR",
    latent.vars = "nCount_peaks",
    slot = "data",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  ))

  out <- expect_no_error(suppressWarnings(RunDEtest(
    srt,
    assay = "peaks",
    layer = "data",
    group.by = "group",
    group1 = "A",
    group2 = "B",
    features = features,
    test.use = "LR",
    latent.vars = "nCount_peaks",
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )))
  actual <- out@tools$DEtest_custom$AllMarkers_LR
  expect_s3_class(actual, "data.frame")
  expect_setequal(actual$gene, rownames(expected))
  actual <- actual[match(rownames(expected), actual$gene), , drop = FALSE]
  expect_equal(actual$avg_log2FC, expected$avg_log2FC, tolerance = 1e-8)
  expect_equal(actual$p_val, expected$p_val, tolerance = 1e-8)

  all_groups <- expect_no_error(suppressWarnings(RunDEtest(
    srt,
    assay = "peaks",
    layer = "data",
    group.by = "group",
    markers_type = "all",
    features = features,
    test.use = "LR",
    latent.vars = "nCount_peaks",
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )))
  expect_gt(nrow(all_groups@tools$DEtest_group$AllMarkers_LR), 0)

  expect_error(
    suppressWarnings(RunDEtest(
      srt,
      assay = "peaks",
      layer = "data",
      group.by = "group",
      markers_type = "all",
      features = features,
      test.use = "LR",
      latent.vars = "missing_depth_covariate",
      fc.threshold = 1,
      min.pct = 0,
      only.pos = FALSE,
      verbose = FALSE
    )),
    regexp = "missing_depth_covariate|requested variables|not found"
  )
})

test_that("RunDEtest normalizes raw ChromatinAssay data with TF-IDF", {
  srt <- make_test_atac_da_object(normalize = FALSE)
  expected <- Signac::RunTFIDF(srt, assay = "peaks", verbose = FALSE)

  out <- suppressWarnings(RunDEtest(
    srt,
    assay = "peaks",
    layer = "data",
    group.by = "group",
    markers_type = "all",
    features = rownames(srt)[1:20],
    test.use = "wilcox",
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  ))

  expect_equal(
    SeuratObject::LayerData(out[["peaks"]], layer = "data"),
    SeuratObject::LayerData(expected[["peaks"]], layer = "data"),
    tolerance = 1e-12
  )
})

sample_method_packages <- c(
  limma = "limma",
  edgeR = "edgeR",
  DESeq2 = "DESeq2",
  dream = "variancePartition"
)

for (method in names(sample_method_packages)) {
  test_that(paste("RunDEtest reaches the", method, "pseudobulk method"), {
    skip_if_not_installed(sample_method_packages[[method]])
    srt <- make_test_atac_da_object(normalize = FALSE, add_samples = TRUE)
    out <- expect_no_error(suppressWarnings(RunDEtest(
      srt,
      assay = "peaks",
      layer = "counts",
      group1 = "ctrl",
      group2 = "case",
      features = rownames(srt),
      feature_type = "peak",
      test.use = method,
      sample_col = "sample",
      condition_col = "condition",
      min.cells.sample = 1,
      fc.threshold = 1,
      only.pos = FALSE,
      verbose = FALSE
    )))
    markers <- out@tools$DEtest_pseudobulk[[paste0("AllMarkers_", method)]]
    expect_s3_class(markers, "data.frame")
    expect_gt(nrow(markers), 0)
  })
}
