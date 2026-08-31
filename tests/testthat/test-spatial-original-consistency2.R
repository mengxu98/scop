# End-to-end consistency between scop spatial wrappers and the original
# backend methods, run with real datasets (visium_human_pancreas_sub,
# pancreas_sub, panc8_sub). No simulated data and
# no mocked backends.

real_visium_subset2 <- function(n = 150, seed = 42, assay = NULL) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  set.seed(seed)
  srt <- srt[, sample(colnames(srt), n)]
  srt
}

visium_assay2 <- function(srt) {
  SeuratObject::DefaultAssay(srt)
}

real_panc8_reference2 <- function(n = 300, seed = 7) {
  data(panc8_sub)
  set.seed(seed)
  panc8_sub[, sample(colnames(panc8_sub), n)]
}

real_pancreas_reference2 <- function(n = 400, seed = 5) {
  data(pancreas_sub)
  set.seed(seed)
  srt <- pancreas_sub[, sample(colnames(pancreas_sub), n)]
  srt$CellType <- factor(srt$CellType)
  srt
}

test_that("RunRCTD weights match the original spacexr pipeline", {
  skip_on_cran()
  skip_if_not_installed("spacexr")

  srt <- real_visium_subset2(n = 120)
  reference <- real_panc8_reference2()
  tb <- table(reference$celltype)
  reference <- reference[, reference$celltype %in% names(tb)[tb >= 30]]
  reference$celltype <- factor(reference$celltype)

  wrapped <- RunRCTD(
    srt,
    reference = reference,
    reference_label = "celltype",
    min_cells = 25,
    rctd_mode = "multi",
    max_cores = 1,
    verbose = FALSE
  )
  wrapped_weights <- as.matrix(wrapped@tools$RCTD$weights)

  # original pipeline with the same inputs; spacexr exposes two generations of
  # API (Bioc 1.4.0: createRctd/runRctd, GitHub 2.x: SpatialRNA/create.RCTD)
  labels <- scop:::resolve_reference_labels(reference, "celltype")
  names(labels) <- colnames(reference)
  labels <- labels[!is.na(labels) & nzchar(as.character(labels))]
  labels <- scop:::rctd_filter_labels_by_min_cells(labels, min_cells = 25, verbose = FALSE)
  reference2 <- reference[, names(labels)]
  labels <- factor(as.character(labels), levels = unique(as.character(labels)))
  names(labels) <- colnames(reference2)
  label_map <- scop:::rctd_backend_label_map(labels)
  features <- intersect(rownames(srt), rownames(reference2))
  st <- scop:::rctd_get_count_matrix(srt, visium_assay2(srt), "counts", features, "Spatial", TRUE, FALSE)
  rf <- scop:::rctd_get_count_matrix(reference2, visium_assay2(reference2), "counts", features, "Reference", TRUE, FALSE)
  cq <- scop:::rctd_sparse_quality_cpp(st, rf)
  st <- st[rownames(st)[cq$keep_features], , drop = FALSE]
  rf <- rf[rownames(rf)[cq$keep_features], , drop = FALSE]
  st_numi <- cq$st_numi
  names(st_numi) <- colnames(st)
  keep_spots <- is.finite(st_numi) & st_numi > 0
  st <- st[, keep_spots, drop = FALSE]
  st_numi <- st_numi[keep_spots]
  coords <- scop:::resolve_spatial_spot_coords(srt, colnames(st), NULL, c("x", "y"), "raw")
  ref_numi <- cq$ref_numi
  names(ref_numi) <- colnames(rf)

  exports <- getNamespaceExports("spacexr")
  if (all(c("createRctd", "runRctd") %in% exports)) {
    testthat::skip_if_not_installed("SpatialExperiment")
    spatial_spe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = st),
      colData = S4Vectors::DataFrame(nUMI = st_numi, row.names = colnames(st)),
      spatialCoords = as.matrix(coords)
    )
    reference_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = rf),
      colData = S4Vectors::DataFrame(
        cell_type = label_map$labels,
        nUMI = ref_numi,
        row.names = colnames(rf)
      )
    )
    rctd_data <- spacexr::createRctd(spatial_spe, reference_se, cell_type_col = "cell_type")
    rctd_result <- spacexr::runRctd(rctd_data, rctd_mode = "multi", max_cores = 1)
    original <- as.matrix(SummarizedExperiment::assay(rctd_result, "weights"))
    # the new API stores cell types as rows and spots as columns
    if (sum(colnames(original) %in% colnames(st)) > sum(rownames(original) %in% colnames(st))) {
      original <- t(original)
    }
  } else {
    puck <- spacexr::SpatialRNA(coords, st, st_numi)
    ref <- spacexr::Reference(rf, label_map$labels, ref_numi)
    rctd <- spacexr::create.RCTD(puck, ref, max_cores = 1)
    rctd <- spacexr::run.RCTD(rctd, doublet_mode = "multi")
    cell_types <- rctd@cell_type_info$renorm[[2]]
    original <- do.call(rbind, lapply(rctd@results, function(x) x$all_weights))
    rownames(original) <- colnames(rctd@spatialRNA@counts)
    colnames(original) <- cell_types
    original <- spacexr::normalize_weights(original)
  }

  display_map <- label_map$cell_types
  names(display_map) <- label_map$backend
  colnames(original) <- unname(display_map[colnames(original)])
  original <- original[rownames(wrapped_weights), colnames(wrapped_weights), drop = FALSE]
  expect_identical(rownames(original), rownames(wrapped_weights))
  expect_identical(colnames(original), colnames(wrapped_weights))
  expect_equal(original, wrapped_weights, tolerance = 1e-6)
  expect_equal(unname(rowSums(wrapped_weights)), rep(1, nrow(wrapped_weights)), tolerance = 1e-6)
  expect_true("RCTD_dominant_type" %in% colnames(wrapped[[]]))
})

test_that("RunBANKSY clusters match the original Banksy pipeline", {
  skip_on_cran()
  skip_if_not_installed("Banksy")
  skip_if_not_installed("SpatialExperiment")

  srt <- real_visium_subset2(n = 150)
  set.seed(1)
  wrapped <- RunBANKSY(
    srt,
    assay = visium_assay2(srt),
    layer = "data",
    lambda = 0.2,
    k_geom = 10,
    npcs = 10,
    algo = "leiden",
    k_neighbors = 20,
    resolution = 0.5,
    seed = 1,
    verbose = FALSE
  )
  wrapped_clusters <- as.character(wrapped$BANKSY_cluster)
  names(wrapped_clusters) <- rownames(wrapped[[]])

  # original pipeline with identical inputs
  expr <- Seurat::GetAssayData(srt, assay = visium_assay2(srt), layer = "data")
  coords <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = srt, image = NULL, coord.cols = c("x", "y"), coordinate_space = "raw"
  )$data
  spots <- intersect(colnames(srt), rownames(coords))
  coords <- coords[spots, , drop = FALSE]
  expr <- expr[, spots, drop = FALSE]
  keep_features <- Matrix::rowSums(expr > 0) > 0
  keep_spots <- Matrix::colSums(expr > 0) > 0
  expr <- expr[keep_features, keep_spots, drop = FALSE]
  se <- SpatialExperiment::SpatialExperiment(
    assays = stats::setNames(list(expr), "scop_input"),
    colData = S4Vectors::DataFrame(
      x = coords[colnames(expr), "x"],
      y = coords[colnames(expr), "y"],
      row.names = colnames(expr)
    ),
    spatialCoords = as.matrix(coords[colnames(expr), c("x", "y")])
  )
  se <- Banksy::computeBanksy(se,
    assay_name = "scop_input", coord_names = c("x", "y"),
    compute_agf = FALSE, M = 1, k_geom = 10
  )
  se <- Banksy::runBanksyPCA(se,
    assay_name = "scop_input", M = 1, lambda = 0.2,
    npcs = 10, use_agf = FALSE, group = NULL, seed = 1
  )
  se <- Banksy::clusterBanksy(se,
    assay_name = "scop_input", M = 1, lambda = 0.2,
    use_agf = FALSE, npcs = 10, algo = "leiden", k_neighbors = 20,
    resolution = 0.5, group = NULL, seed = 1
  )
  original_clusters <- as.character(
    SummarizedExperiment::colData(se)[[Banksy::clusterNames(se)[1]]]
  )
  names(original_clusters) <- colnames(se)

  common <- intersect(names(original_clusters), names(wrapped_clusters))
  expect_gt(length(common), ncol(srt) * 0.9)
  wrapped_clusters <- wrapped_clusters[common]
  original_clusters <- original_clusters[common]
  expect_identical(unname(wrapped_clusters), unname(original_clusters))
  expect_identical(wrapped@tools$BANKSY$parameters$seed, 1)
})

test_that("RunCARD(CARDspa) proportions match the original CARDspa pipeline", {
  skip_on_cran()
  skip_if_not_installed("CARDspa")

  srt <- real_visium_subset2(n = 120, seed = 3)
  reference <- real_panc8_reference2()
  tb <- table(reference$celltype)
  reference <- reference[, reference$celltype %in% names(tb)[tb >= 30]]
  reference$celltype <- factor(reference$celltype)
  genes_use <- intersect(rownames(srt), rownames(reference))
  srt <- srt[genes_use, ]

  # CARD's deconvolution draws random initial values; pin the RNG so the
  # wrapped and original runs see identical random draws
  set.seed(42)
  wrapped <- RunCARD(
    srt,
    reference = reference,
    reference_label = "celltype",
    minCountGene = 0,
    minCountSpot = 0,
    verbose = FALSE
  )
  wrapped_weights <- as.matrix(wrapped@tools$CARD$weights)

  # original CARDspa one-step pipeline on identical inputs
  st_counts <- Seurat::GetAssayData(srt, assay = visium_assay2(srt), layer = "counts")
  ref_counts <- Seurat::GetAssayData(reference, assay = visium_assay2(reference), layer = "counts")
  coords <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = srt, image = NULL, coord.cols = c("x", "y"), coordinate_space = "raw"
  )$data
  labels <- as.character(reference$celltype)
  names(labels) <- colnames(reference)
  ref_meta <- data.frame(
    .scop_cell_type = labels,
    .scop_sample = "sample1",
    row.names = colnames(reference)
  )
  card_create_args <- list(
    sc_count = ref_counts,
    sc_meta = ref_meta,
    spatial_count = st_counts,
    spatial_location = coords[, c("x", "y")],
    ct.varname = ".scop_cell_type",
    ct_varname = ".scop_cell_type",
    sample.varname = ".scop_sample",
    sample_varname = ".scop_sample",
    minCountGene = 0,
    minCountSpot = 0,
    mincountgene = 0,
    mincountspot = 0,
    ct_select = NULL
  )
  CARDspa_obj <- do.call(
    get("createCARDObject", asNamespace("CARDspa")),
    getFromNamespace("card_match_formals", "scop")(
      get("createCARDObject", asNamespace("CARDspa")),
      card_create_args
    )
  )
  info <- CARDspa_obj@info_parameters
  sc_eset <- CARDspa_obj@sc_eset
  ct_use <- info[["ct.select"]]
  ct_varname <- info[["ct.varname"]]
  sample_varname <- info[["sample.varname"]]
  basis <- get("create_ref", asNamespace("CARDspa"))(sc_eset, ct_use, ct_varname, sample_varname)[["basis"]]
  basis <- basis[, colnames(basis) %in% ct_use, drop = FALSE]
  common_genes <- intersect(rownames(st_counts), rownames(basis))
  informative <- get("select_info", asNamespace("CARDspa"))(basis, sc_eset, common_genes, ct_use, ct_varname)
  informative_counts <- st_counts[rownames(st_counts) %in% informative, , drop = FALSE]
  informative_counts <- informative_counts[Matrix::rowSums(informative_counts) > 0, , drop = FALSE]
  supported <- intersect(
    colnames(st_counts),
    colnames(informative_counts)[Matrix::colSums(informative_counts) > 0]
  )
  CARDspa_obj@spatial_countMat <- st_counts[, supported, drop = FALSE]
  deconv_fun <- get("CARD_deconvolution", asNamespace("CARDspa"))
  set.seed(42)
  original <- do.call(
    deconv_fun,
    getFromNamespace("card_match_formals", "scop")(
      deconv_fun,
      card_create_args
    )
  )
  original_weights <- original$Proportion_CARD

  original_weights <- original_weights[rownames(wrapped_weights), colnames(wrapped_weights), drop = FALSE]
  expect_identical(rownames(original_weights), rownames(wrapped_weights))
  expect_identical(colnames(original_weights), colnames(wrapped_weights))
  expect_gt(stats::cor(as.vector(original_weights), as.vector(wrapped_weights)), 0.9)
  expect_lt(mean(abs(original_weights - wrapped_weights)), 0.05)
  expect_equal(unname(rowSums(wrapped_weights)), rep(1, nrow(wrapped_weights)), tolerance = 1e-6)
  expect_true("CARD_dominant_type" %in% colnames(wrapped[[]]))
})

test_that("RunSpotSweeper local outliers match the original SpotSweeper pipeline", {
  skip_on_cran()
  skip_if_not_installed("SpotSweeper")
  skip_if_not_installed("SpatialExperiment")

  srt <- real_visium_subset2(n = 120)
  wrapped <- RunSpotSweeper(
    srt,
    assay = visium_assay2(srt),
    layer = "counts",
    n_neighbors = 20,
    verbose = FALSE
  )

  # original pipeline: per-metric localOutliers on the same SPE
  expr <- Seurat::GetAssayData(srt, assay = visium_assay2(srt), layer = "counts")
  coords <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = srt, image = NULL, coord.cols = c("x", "y"), coordinate_space = "raw"
  )$data
  spots <- intersect(colnames(srt), rownames(coords))
  coords <- coords[spots, , drop = FALSE]
  expr <- expr[, spots, drop = FALSE]
  coldata <- data.frame(
    nCount_Spatial = Matrix::colSums(expr),
    nFeature_Spatial = Matrix::colSums(expr > 0),
    percent.mito = 0,
    .SpotSweeper_sample = "sample1",
    row.names = colnames(expr)
  )
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = expr),
    colData = S4Vectors::DataFrame(coldata),
    spatialCoords = as.matrix(coords[, c("x", "y")])
  )
  spe <- SpotSweeper::localOutliers(spe,
    metric = "nCount_Spatial", direction = "lower",
    n_neighbors = 20, samples = ".SpotSweeper_sample", log = TRUE, cutoff = 3, workers = 1
  )
  spe <- SpotSweeper::localOutliers(spe,
    metric = "nFeature_Spatial", direction = "lower",
    n_neighbors = 20, samples = ".SpotSweeper_sample", log = TRUE, cutoff = 3, workers = 1
  )

  wrapped_cd <- wrapped@tools$SpotSweeper$colData
  original_cd <- as.data.frame(SummarizedExperiment::colData(spe))
  common <- intersect(colnames(wrapped_cd), colnames(original_cd))
  expect_true(all(c("nCount_Spatial_z", "nFeature_Spatial_z") %in% common))
  for (metric in c("nCount_Spatial", "nFeature_Spatial")) {
    expect_equal(
      wrapped_cd[spots, paste0(metric, "_z")],
      original_cd[spots, paste0(metric, "_z")],
      tolerance = 1e-12,
      label = paste0(metric, " z-scores")
    )
    expect_equal(
      wrapped_cd[spots, paste0(metric, "_outliers")],
      original_cd[spots, paste0(metric, "_outliers")],
      label = paste0(metric, " outlier flags")
    )
  }
  expect_true(all(c("SpotSweeper_QC", "SpotSweeper_local_outlier_qc") %in% colnames(wrapped[[]])))
})

test_that("RunSpatialIntegration PRECAST domains match the original PRECAST pipeline", {
  skip_on_cran()
  skip_if_not_installed("PRECAST")

  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  srt$sample <- ifelse(srt$y > stats::median(srt$y), "slice_a", "slice_b")
  set.seed(1)
  cells1 <- sample(which(srt$sample == "slice_a"), 100)
  cells2 <- sample(which(srt$sample == "slice_b"), 100)
  s1 <- srt[, cells1]
  s2 <- srt[, cells2]

  wrapped <- RunSpatialIntegration(
    list(s1, s2),
    method = "PRECAST",
    sample.by = "sample",
    layer = "counts",
    verbose = FALSE
  )
  wrapped_domains <- wrapped@tools$SpatialIntegration$methods$PRECAST$domains

  # original pipeline with identical inputs
  expr1 <- Seurat::GetAssayData(s1, assay = visium_assay2(s1), layer = "counts")
  expr2 <- Seurat::GetAssayData(s2, assay = visium_assay2(s2), layer = "counts")
  common_genes <- intersect(rownames(expr1), rownames(expr2))
  expr1 <- expr1[common_genes, ]
  expr2 <- expr2[common_genes, ]
  coords1 <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = s1, image = NULL, coord.cols = c("x", "y"), coordinate_space = "raw"
  )$data
  coords2 <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = s2, image = NULL, coord.cols = c("x", "y"), coordinate_space = "raw"
  )$data
  s1$row <- coords1[colnames(s1), "y"]
  s1$col <- coords1[colnames(s1), "x"]
  s2$row <- coords2[colnames(s2), "y"]
  s2$col <- coords2[colnames(s2), "x"]
  set.seed(1)
  precast <- PRECAST::CreatePRECASTObject(
    seuList = list(s1, s2),
    project = "spatial_integration",
    customGenelist = common_genes
  )
  precast <- PRECAST::AddAdjList(precast)
  precast <- PRECAST::AddParSetting(precast)
  precast <- PRECAST::PRECAST(precast)
  precast <- PRECAST::SelectModel(precast)
  original_domains <- unlist(lapply(seq_along(precast@resList$cluster), function(k) {
    cl <- precast@resList$cluster[[k]]
    cl <- as.character(cl)
    cells_k <- if (k == 1L) colnames(s1) else colnames(s2)
    names(cl) <- paste0(if (k == 1L) "sample1" else "sample2", "_", cells_k)
    cl
  }))
  original_domains <- as.character(original_domains[names(wrapped_domains)])

  expect_equal(
    unname(wrapped_domains[!is.na(wrapped_domains)]),
    unname(original_domains[!is.na(original_domains)])
  )
  expect_true("SpatialIntegration_PRECAST" %in% SeuratObject::Reductions(wrapped))
})
