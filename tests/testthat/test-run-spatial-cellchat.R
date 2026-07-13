make_spatialcellchat_test_object <- function(samples = "s1") {
  cells <- paste0("spot", seq_len(8))
  counts <- matrix(
    rep(seq_len(24), length.out = 6 * length(cells)),
    nrow = 6,
    dimnames = list(paste0("gene", seq_len(6)), cells)
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$celltype <- rep(c("A", "B"), each = 4)
  srt$sample <- rep(samples, length.out = ncol(srt))
  srt$col <- rep(0:3, 2)
  srt$row <- rep(0:1, each = 4)
  srt
}

mock_spatialcellchat_table <- function(sample, analysis.level) {
  data.frame(
    sender = c("A", "B"),
    receiver = c("B", "A"),
    ligand = c("L1", "L2"),
    receptor = c("R1", "R2"),
    interaction_name = c("L1_R1", "L2_R2"),
    pathway_name = c("PATH1", "PATH2"),
    classification = c("PATH1", "PATH2"),
    pair_lr = c("L1-R1", "L2-R2"),
    interaction_label = c("L1 - R1", "L2 - R2"),
    interaction_display = c("L1 - R1", "L2 - R2"),
    ligand_display = c("L1", "L2"),
    receptor_display = c("R1", "R2"),
    score = c(0.8, 0.4),
    pvalue = c(0.01, 0.04),
    significant = TRUE,
    neglog10_pvalue = -log10(c(0.01, 0.04)),
    method = "SpatialCellChat",
    sample = sample,
    modality = "spatial",
    analysis_level = analysis.level,
    spatially_constrained = TRUE,
    significance_basis = "spatial_permutation",
    passes_filter = TRUE,
    stringsAsFactors = FALSE
  )
}

local_mock_spatialcellchat_backend <- function(fail = FALSE) {
  testthat::local_mocked_bindings(
    spatialcellchat_check_r = function(...) invisible(TRUE),
    spatialcellchat_database = function(...) list(mock = TRUE),
    spatialcellchat_run_one = if (isTRUE(fail)) {
      function(...) stop("backend failed")
    } else {
      function(...) structure(list(mock = TRUE), class = "mock_spatialcellchat")
    },
    spatialcellchat_extract_table = function(chat, sample, analysis.level, do.permutation) {
      mock_spatialcellchat_table(sample, analysis.level)
    },
    spatialcellchat_remote_info = function() {
      list(package_version = "0.1.0", remote_sha = "abc123", remote_repo = "jinworks/SpatialCellChat")
    },
    .package = "scop",
    .env = parent.frame()
  )
}

test_that("coordinate conversion is explicit and micron based", {
  srt <- make_spatialcellchat_test_object()
  metric <- scop:::spatialcellchat_metric_coordinates(
    srt,
    cells = colnames(srt),
    image = NULL,
    coord.cols = c("col", "row"),
    technology = "generic",
    coordinate.unit = "pixel",
    ratio = 0.5,
    tol = 5
  )
  expect_equal(unname(metric$data$x), unname(srt$col * 0.5))
  expect_equal(unname(metric$data$x_raw), unname(srt$col))
  expect_equal(unname(metric$data$x_display), unname(srt$col))
  expect_equal(unname(metric$data$y_display), unname(srt$row))
  expect_identical(metric$source$coordinate_unit, "pixel")
  expect_equal(metric$source$scale_to_micron, 0.5)
  expect_equal(metric$spatial.factors, list(ratio = 1, tol = 5))
  expect_error(
    scop:::spatialcellchat_metric_coordinates(
      srt, colnames(srt), NULL, c("col", "row"), "generic", "pixel", NULL, 5
    ),
    "explicit positive"
  )
})

test_that("strict auto detection refuses ambiguous generic data", {
  srt <- make_spatialcellchat_test_object()
  expect_error(
    scop:::spatialcellchat_detect_technology(srt, "auto"),
    "Supply.*technology"
  )
  expect_error(
    scop:::spatialcellchat_detect_level(srt, "auto", "generic", NULL),
    "Cannot determine"
  )
})

test_that("Visium scale evidence and multi-image selection are strict", {
  srt <- make_spatialcellchat_test_object()
  srt@misc$scalefactors_json <- list(
    slice1 = list(spot_diameter_fullres = 130)
  )
  expect_equal(scop:::spatialcellchat_visium_spot_diameter(srt, "slice1"), 130)
  testthat::local_mocked_bindings(
    Images = function(object) c("slice1", "slice2"),
    .package = "SeuratObject"
  )
  expect_error(
    scop:::spatialcellchat_image_map(srt, "ALL", image = NULL),
    "Multiple spatial images"
  )
})

test_that("composition validation aligns and normalizes spots", {
  cells <- paste0("spot", 1:3)
  composition <- matrix(
    c(2, 1, 1, 3, 4, 1),
    nrow = 3,
    dimnames = list(cells, c("A", "B"))
  )
  out <- scop:::spatialcellchat_validate_composition(
    composition[, c("B", "A")],
    cells,
    groups = c("A", "B"),
    normalize = TRUE
  )
  expect_true(out$normalized)
  expect_equal(unname(rowSums(out$data)), rep(1, 3))
  expect_identical(colnames(out$data), c("A", "B"))
  expect_error(
    scop:::spatialcellchat_validate_composition(composition[-1, ], cells),
    "missing selected spots"
  )
})

test_that("expression validation rejects scientifically invalid backend input", {
  srt <- make_spatialcellchat_test_object()
  expression <- GetAssayData5(srt, assay = "RNA", layer = "data")
  expect_true(scop:::spatialcellchat_validate_expression(expression, colnames(srt)))

  negative <- as.matrix(expression)
  negative[[1L]] <- -1
  expect_error(
    scop:::spatialcellchat_validate_expression(negative, colnames(srt)),
    "non-negative"
  )

  non_finite <- as.matrix(expression)
  non_finite[[1L]] <- Inf
  expect_error(
    scop:::spatialcellchat_validate_expression(non_finite, colnames(srt)),
    "non-finite"
  )

  expect_error(
    scop:::spatialcellchat_validate_expression(expression, rev(colnames(srt))),
    "do not align"
  )
})

test_that("RunSpatialCellChat stores truthful schema and unified CCC rows", {
  local_mock_spatialcellchat_backend()
  srt <- make_spatialcellchat_test_object()
  out <- RunSpatialCellChat(
    srt,
    group.by = "celltype",
    technology = "generic",
    analysis.level = "cell",
    coordinate.unit = "micron",
    tol = 5,
    min.cells = 2,
    min.links = 1,
    nboot = 2,
    database = "custom",
    custom.db = list(mock = TRUE),
    backend = "r",
    verbose = FALSE
  )
  expect_identical(out@tools$SpatialCellChat$schema_version, 1L)
  expect_identical(out@tools$SpatialCellChat$provenance$backend_id, "spatialcellchat")
  expect_identical(out@tools$SpatialCellChat$provenance$remote_sha, "abc123")
  expect_identical(out@tools$SpatialCellChat$source$analysis_unit, "micron")
  expect_equal(out@tools$SpatialCellChat$source$samples$ALL$interaction_range_um, 250)
  expect_equal(out@tools$SpatialCellChat$source$samples$ALL$scale_distance, 0.2)
  expect_null(out@tools$SpatialCellChat$results$default$ALL$native_object)
  stored_coords <- out@tools$SpatialCellChat$results$default$ALL$coordinates
  expect_true(all(c(
    "x", "y", "x_raw", "y_raw", "x_display", "y_display"
  ) %in% colnames(stored_coords)))
  expect_identical(out@tools$SpatialCellChat$results$default$ALL$diagnostics$interpretation, "cell-level communication")
  spatial_rows <- out@tools$CCC$long_table$method == "SpatialCellChat"
  expect_true(any(spatial_rows))
  expect_true(all(out@tools$CCC$long_table$spatially_constrained[spatial_rows]))
  expect_identical(GetSpatialResult(out, method = "RunSpatialCellChat")$method, "SpatialCellChat")
  expect_identical(SpatialResultInfo(out)$plot_function, "SpatialCellChatPlot")
  expect_error(GetCCCObject(out, method = "SpatialCellChat"), "not stored")
  expect_s3_class(SpatialCellChatPlot(out, plot_type = "incoming"), "ggplot")
  expect_s3_class(SpatialCellChatPlot(out, plot_type = "spatial_network"), "ggplot")
  expect_error(SpatialCellChatPlot(out, plot_type = "pathway"), "signaling.*required")
  expect_error(SpatialCellChatPlot(out, plot_type = "lr_spatial"), "pairLR.use.*required")
  expect_s3_class(
    SpatialCellChatPlot(out, plot_type = "pathway", signaling = "PATH1"),
    "ggplot"
  )
  expect_s3_class(
    SpatialCellChatPlot(out, plot_type = "lr_spatial", pairLR.use = "L1_R1"),
    "ggplot"
  )
  expect_s3_class(
    CCCNetworkPlot(out, method = "SpatialCellChat", plot_type = "circle"),
    "recordedplot"
  )
  expect_s3_class(
    CCCHeatmap(out, method = "SpatialCellChat", plot_type = "bubble"),
    "ggplot"
  )
  expect_s3_class(
    CCCStatPlot(out, method = "SpatialCellChat", plot_type = "bar"),
    "ggplot"
  )
})

test_that("cell, spot, and composition routes stay semantically distinct", {
  local_mock_spatialcellchat_backend()
  srt <- make_spatialcellchat_test_object()
  common <- list(
    srt = srt,
    group.by = "celltype",
    technology = "generic",
    coordinate.unit = "micron",
    tol = 5,
    min.cells = 2,
    min.links = 1,
    nboot = 2,
    database = "custom",
    custom.db = list(mock = TRUE),
    backend = "r",
    verbose = FALSE
  )
  cell <- do.call(RunSpatialCellChat, c(common, list(analysis.level = "cell")))
  spot <- do.call(RunSpatialCellChat, c(common, list(analysis.level = "spot")))
  composition <- matrix(
    rep(c(0.7, 0.3), ncol(srt)),
    nrow = ncol(srt), byrow = TRUE,
    dimnames = list(colnames(srt), c("A", "B"))
  )
  composed <- do.call(RunSpatialCellChat, c(common, list(
    analysis.level = "composition",
    composition = composition
  )))
  expect_identical(cell@tools$SpatialCellChat$summary$analysis_level, "cell")
  expect_identical(spot@tools$SpatialCellChat$summary$analysis_level, "spot")
  expect_identical(composed@tools$SpatialCellChat$summary$analysis_level, "composition")
  expect_match(spot@tools$SpatialCellChat$results$default$ALL$diagnostics$interpretation, "spot/domain")
  expect_error(
    do.call(RunSpatialCellChat, c(common, list(
      analysis.level = "spot", contact.dependent = TRUE, contact.range = 10
    ))),
    "not cells"
  )
})

test_that("CellChat and SpatialCellChat rows coexist without method theft", {
  local_mock_spatialcellchat_backend()
  srt <- make_spatialcellchat_test_object()
  cellchat <- mock_spatialcellchat_table("legacy", "cell")
  cellchat$method <- "CellChat"
  cellchat$modality <- NA_character_
  cellchat$analysis_level <- NA_character_
  cellchat$spatially_constrained <- NA
  cellchat$significance_basis <- NA_character_
  cellchat$passes_filter <- NA
  srt@tools$CCC <- list(
    method = "CCC",
    methods = "CellChat",
    long_table = cellchat,
    metadata = list(schema = "scop_ccc_unified_v1")
  )
  out <- RunSpatialCellChat(
    srt,
    group.by = "celltype",
    technology = "generic",
    analysis.level = "cell",
    coordinate.unit = "micron",
    tol = 5,
    min.cells = 2,
    min.links = 1,
    nboot = 2,
    database = "custom",
    custom.db = list(mock = TRUE),
    backend = "r",
    verbose = FALSE
  )
  expect_setequal(out@tools$CCC$methods, c("CellChat", "SpatialCellChat"))
  expect_true(all(c("CellChat", "SpatialCellChat") %in% out@tools$CCC$long_table$method))
  legacy <- out@tools$CCC$long_table$method == "CellChat"
  expect_true(all(is.na(out@tools$CCC$long_table$spatially_constrained[legacy])))
})

test_that("multi-sample runs remain isolated and require sample-aware access", {
  local_mock_spatialcellchat_backend()
  srt <- make_spatialcellchat_test_object(samples = c("s1", "s2"))
  out <- RunSpatialCellChat(
    srt,
    group.by = "celltype",
    sample.by = "sample",
    technology = "generic",
    analysis.level = "cell",
    coordinate.unit = "micron",
    tol = 5,
    min.cells = 1,
    min.links = 1,
    nboot = 2,
    database = "custom",
    custom.db = list(mock = TRUE),
    backend = "r",
    verbose = FALSE
  )
  expect_identical(names(out@tools$SpatialCellChat$results$default), c("s1", "s2"))
  expect_true(all(vapply(
    out@tools$SpatialCellChat$results$default,
    function(x) length(unique(x$interactions$sample)) == 1L,
    logical(1)
  )))
  expect_error(SpatialCellChatPlot(out), "select.*sample")
  expect_s3_class(SpatialCellChatPlot(out, sample = "s1", plot_type = "outgoing"), "ggplot")
})

test_that("backend failure leaves the object unchanged", {
  local_mock_spatialcellchat_backend(fail = TRUE)
  srt <- make_spatialcellchat_test_object()
  before <- srt@tools
  expect_error(
    RunSpatialCellChat(
      srt,
      group.by = "celltype",
      technology = "generic",
      analysis.level = "cell",
      coordinate.unit = "micron",
      tol = 5,
      database = "custom",
      custom.db = list(mock = TRUE),
      verbose = FALSE
    ),
    "backend failed"
  )
  expect_identical(srt@tools, before)
})

test_that("full storage exposes native objects and protects result names", {
  local_mock_spatialcellchat_backend()
  srt <- make_spatialcellchat_test_object()
  out <- RunSpatialCellChat(
    srt,
    group.by = "celltype",
    technology = "generic",
    analysis.level = "cell",
    coordinate.unit = "micron",
    tol = 5,
    min.cells = 2,
    min.links = 1,
    nboot = 2,
    database = "custom",
    custom.db = list(mock = TRUE),
    store.object = "full",
    backend = "r",
    verbose = FALSE
  )
  expect_s3_class(GetCCCObject(out, method = "SpatialCellChat"), "mock_spatialcellchat")
  expect_error(
    RunSpatialCellChat(
      out,
      group.by = "celltype",
      technology = "generic",
      analysis.level = "cell",
      coordinate.unit = "micron",
      tol = 5,
      database = "custom",
      custom.db = list(mock = TRUE),
      verbose = FALSE
    ),
    "already exists"
  )
})

test_that("SpatialCellChat is discoverable through the spatial registries", {
  methods <- ListSpatialMethods(pattern = "SpatialCellChat")
  expect_true(all(c("RunSpatialCellChat", "SpatialCellChatPlot") %in% methods$method))
  backend <- SpatialBackendStatus(backend = "spatialcellchat", api_check = FALSE, refresh = TRUE)
  expect_identical(backend$repository, "jinworks/SpatialCellChat")
})
