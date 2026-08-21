# End-to-end consistency between the scop spatial wrappers and the original
# backend methods, run with real datasets (visium_human_pancreas_sub,
# panc8_sub). No simulated data and no mocked backends: the wrappers and the
# original pipelines receive the same real input and must return the same
# scientific results.

real_visium_subset <- function(n = 200, seed = 42) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  set.seed(seed)
  cells <- sample(colnames(srt), n)
  srt <- srt[, cells]
  srt
}

visium_assay <- function(srt) {
  SeuratObject::DefaultAssay(srt)
}

real_human_pancreas_visium_subset <- function(n = 150, seed = 7) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  set.seed(seed)
  cells <- sample(colnames(srt), n)
  srt <- srt[, cells]
  srt
}

real_panc8_reference <- function(n = 200, seed = 7) {
  data(panc8_sub)
  set.seed(seed)
  panc8_sub[, sample(colnames(panc8_sub), n)]
}

fast_row_vars <- function(mat) {
  apply(mat, 1, stats::var)
}

test_that("RunBayesSpace clusters match the original BayesSpace pipeline", {
  skip_on_cran()
  skip_if_not_installed("BayesSpace")
  skip_if_not_installed("SingleCellExperiment")

  srt <- real_visium_subset()
  set.seed(2024)
  wrapped <- RunBayesSpace(
    srt,
    q = 3,
    platform = "Visium",
    image = NULL,
    coord.cols = c("x", "y"),
    preprocess = TRUE,
    n.PCs = 10,
    n.HVGs = 500,
    store_sce = FALSE,
    verbose = FALSE
  )

  sce <- Seurat::as.SingleCellExperiment(srt, assay = visium_assay(srt))
  coordinate_input <- getFromNamespace("bayesspace_add_spatial_coords", "scop")(
    srt = srt,
    sce = sce,
    image = NULL,
    platform = "Visium",
    coord.cols = c("x", "y")
  )
  sce <- coordinate_input$sce
  set.seed(2024)
  sce <- BayesSpace::spatialPreprocess(
    sce,
    platform = "Visium",
    n.PCs = 10,
    n.HVGs = 500,
    skip.PCA = FALSE
  )
  set.seed(2024)
  sce <- BayesSpace::spatialCluster(
    sce,
    q = 3,
    use.dimred = "PCA",
    d = 10,
    platform = "Visium"
  )
  original_clusters <- as.character(
    SummarizedExperiment::colData(sce)$spatial.cluster
  )

  expect_identical(unname(wrapped$BayesSpace_cluster), unname(original_clusters))
  expect_identical(
    rownames(wrapped@tools$BayesSpace$colData),
    colnames(srt)
  )
})

test_that("RunMERINGUE autocorrelation matches original moranTest per feature", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- real_visium_subset()
  expr <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "data")
  set.seed(42)
  features <- head(
    rownames(expr)[order(fast_row_vars(expr), decreasing = TRUE)],
    100
  )

  wrapped <- RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    filterDist = NA_real_,
    binary = TRUE,
    alternative = "greater",
    nperm = 0,
    ncores = 1,
    set_variable_features = FALSE,
    verbose = FALSE
  )
  auto <- wrapped@tools$MERINGUE$autocorrelation
  coords <- wrapped@tools$MERINGUE$coords

  weight <- MERINGUE::getSpatialNeighbors(
    pos = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    filterDist = NA_real_,
    binary = TRUE,
    verbose = FALSE
  )
  original <- lapply(features, function(feature) {
    out <- MERINGUE::moranTest(
      expr[feature, ],
      weight = weight,
      alternative = "greater"
    )
    data.frame(
      feature = feature,
      statistic = out[["observed"]],
      expected = out[["expected"]],
      sd = out[["sd"]],
      p_value = out[["p.value"]]
    )
  })
  original <- do.call(rbind, original)

  merged <- merge(
    auto[, c("feature", "statistic", "expected", "sd", "p_value")],
    original,
    by = "feature"
  )
  expect_equal(nrow(merged), length(features))
  expect_equal(merged$statistic.x, merged$statistic.y, tolerance = 1e-12)
  expect_equal(merged$expected.x, merged$expected.y, tolerance = 1e-12)
  expect_equal(merged$sd.x, merged$sd.y, tolerance = 1e-12)
  expect_true(all(is.na(auto$p_value)))
  expect_true(all(is.na(auto$q_value)))
  expect_identical(order(-auto$statistic), seq_len(nrow(auto)))
})

test_that("RunMERINGUE permutation test matches original moranPermutationTest", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- real_visium_subset(n = 120)
  expr <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "data")
  set.seed(42)
  features <- head(
    rownames(expr)[order(fast_row_vars(expr), decreasing = TRUE)],
    20
  )

  wrapped <- RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    filterDist = NA_real_,
    binary = TRUE,
    alternative = "greater",
    nperm = 50,
    ncores = 1,
    seed = 11,
    set_variable_features = FALSE,
    verbose = FALSE
  )
  auto <- wrapped@tools$MERINGUE$autocorrelation
  coords <- wrapped@tools$MERINGUE$coords

  weight <- MERINGUE::getSpatialNeighbors(
    pos = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    filterDist = NA_real_,
    binary = TRUE,
    verbose = FALSE
  )
  original <- lapply(features, function(feature) {
    out <- MERINGUE::moranPermutationTest(
      z = expr[feature, ],
      w = weight,
      alternative = "greater",
      N = 50,
      seed = 11,
      ncores = 1,
      plot = FALSE
    )
    data.frame(
      feature = feature,
      statistic = out[["observed"]],
      expected = out[["expected"]],
      sd = out[["sd"]],
      p_value = out[["p.value"]]
    )
  })
  original <- do.call(rbind, original)

  merged <- merge(
    auto[, c("feature", "statistic", "expected", "sd", "p_value")],
    original,
    by = "feature"
  )
  expect_equal(nrow(merged), length(features))
  expect_equal(merged$statistic.x, merged$statistic.y, tolerance = 1e-12)
  expect_equal(merged$p_value.x, merged$p_value.y, tolerance = 1e-12)
})

test_that("RunSPOTlight deconvolution matches original SPOTlight with same markers", {
  skip_on_cran()
  skip_if_not_installed("SPOTlight")

  srt <- real_human_pancreas_visium_subset()
  reference <- real_panc8_reference()
  genes_use <- intersect(rownames(srt), rownames(reference))
  srt <- srt[genes_use, ]

  wrapped <- RunSPOTlight(
    srt,
    reference = reference,
    reference_label = "celltype",
    marker_top_n = 30,
    min_prop = 0.01,
    scale = TRUE,
    verbose = FALSE
  )
  wrapped_weights <- wrapped@tools$SPOTlight$weights
  mgs <- wrapped@tools$SPOTlight$marker_genes

  ref_counts <- Seurat::GetAssayData(reference, assay = visium_assay(reference), layer = "counts")
  st_counts <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "counts")
  labels <- reference$celltype
  original <- SPOTlight::SPOTlight(
    x = ref_counts,
    y = st_counts,
    groups = as.character(labels),
    mgs = mgs,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    min_prop = 0.01,
    scale = TRUE,
    verbose = FALSE
  )$mat

  original <- original[rownames(wrapped_weights), colnames(wrapped_weights), drop = FALSE]
  expect_identical(rownames(original), rownames(wrapped_weights))
  expect_identical(colnames(original), colnames(wrapped_weights))
  expect_equal(original, wrapped_weights, tolerance = 1e-6)
  expect_equal(
    unname(rowSums(wrapped_weights)),
    rep(1, nrow(wrapped_weights)),
    tolerance = 1e-6
  )
})

test_that("RunSpaNorm normalized counts match the original SpaNorm pipeline", {
  skip_on_cran()
  skip_if_not_installed("SpaNorm")
  skip_if_not_installed("SpatialExperiment")

  srt <- real_visium_subset(n = 500)
  expr <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "data")
  set.seed(3)
  genes <- head(rownames(expr)[order(apply(expr, 1, stats::var), decreasing = TRUE)], 2000)
  srt <- srt[genes, ]
  set.seed(7)
  wrapped <- RunSpaNorm(
    srt,
    layer = "counts",
    coord.cols = c("x", "y"),
    store_results = FALSE,
    store_spe = FALSE,
    verbose = FALSE
  )

  counts <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "counts")
  coords <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = srt,
    image = NULL,
    coord.cols = c("x", "y"),
    coordinate_space = "raw"
  )$data
  cells <- intersect(colnames(srt), rownames(coords))
  coords <- coords[cells, , drop = FALSE]
  counts <- counts[, cells, drop = FALSE]
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE])
  )
  set.seed(7)
  original <- SpaNorm::SpaNorm(spe, verbose = FALSE)
  original_counts <- SummarizedExperiment::assay(original, "logcounts")

  wrapped_counts <- Seurat::GetAssayData(
    wrapped,
    assay = "SpaNorm",
    layer = "data"
  )
  wrapped_counts <- wrapped_counts[rownames(original_counts), cells, drop = FALSE]
  expect_identical(dim(wrapped_counts), dim(original_counts))
  expect_equal(
    as.matrix(wrapped_counts),
    as.matrix(original_counts),
    tolerance = 1e-6
  )
})

test_that("RunSpatialDWLS weights match independent NNLS fitting", {
  skip_on_cran()

  srt <- real_human_pancreas_visium_subset()
  reference <- real_panc8_reference()
  genes_use <- intersect(rownames(srt), rownames(reference))
  srt <- srt[genes_use, ]

  wrapped <- RunSpatialDWLS(
    srt,
    reference = reference,
    reference_label = "celltype",
    min_cells = 2,
    normalize = TRUE,
    verbose = FALSE
  )
  weights <- wrapped@tools$SpatialDWLS$weights
  signatures <- wrapped@tools$SpatialDWLS$signatures

  shared <- intersect(rownames(srt), rownames(reference))
  st_expr <- Seurat::GetAssayData(srt, assay = visium_assay(srt), layer = "counts")
  st_expr <- st_expr[shared, rownames(weights), drop = FALSE]
  st_expr <- getFromNamespace("spatial_dwls_normalize_matrix", "scop")(st_expr)
  st_expr <- st_expr[rownames(signatures), , drop = FALSE]

  qr_sig <- qr(signatures)
  expected <- matrix(
    0,
    nrow = nrow(weights),
    ncol = ncol(weights),
    dimnames = dimnames(weights)
  )
  for (spot in rownames(weights)) {
    coef <- qr.coef(qr_sig, st_expr[, spot])
    coef[!is.finite(coef) | coef < 0] <- 0
    expected[spot, ] <- coef
  }
  expected <- getFromNamespace("spatial_normalize_weights", "scop")(expected)

  expect_equal(weights, expected, tolerance = 1e-12)
})
