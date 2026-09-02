make_fast_path_object <- function(seed = 7) {
  set.seed(seed)
  mat <- Matrix::rsparsematrix(120, 80, density = 0.08)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 60, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  RunPCA(obj, features = SeuratObject::VariableFeatures(obj), npcs = 10, verbose = FALSE)
}

test_that("internal workflows use scop Seurat-compatible entry points", {
  has_namespaced_call <- function(expr, entry_point) {
    if (
      is.call(expr) && length(expr) >= 3L &&
        identical(expr[[1L]], quote(`::`)) &&
        identical(expr[[2L]], quote(Seurat)) &&
        identical(as.character(expr[[3L]]), entry_point)
    ) {
      return(TRUE)
    }
    is.recursive(expr) && any(vapply(
      as.list(expr),
      has_namespaced_call,
      logical(1),
      entry_point = entry_point
    ))
  }
  routes <- list(
    NormalizeData = c(
      "scpagwas_prepare_seurat", "EnrichmentHeatmap", "scissor_heatmap_plot"
    ),
    FindVariableFeatures = c(
      "srt_reorder", "RunKNNMap", "RunKNNPredict", "RunMonocle2",
      "RunDynamicFeatures", "CheckDataList"
    ),
    ScaleData = c(
      "Uncorrected_integrate", "Seurat_integrate", "MNN_integrate",
      "Harmony_integrate", "BBKNN_integrate", "CSS_integrate",
      "Conos_integrate", "ComBat_integrate"
    ),
    FindNeighbors = c(
      "find_neighbors_and_clusters", "RunKNNMap", "RunKNNPredict",
      "RunFR.Seurat", "RunRareQ"
    ),
    FindClusters = c("find_neighbors_and_clusters", "RunStandardWorkflow"),
    FindMultiModalNeighbors = "WNN_integrate"
  )
  ns <- asNamespace("scop")
  for (entry_point in names(routes)) {
    for (workflow in routes[[entry_point]]) {
      code <- body(get(workflow, envir = ns, inherits = FALSE))
      expect_true(
        entry_point %in% all.names(code, functions = TRUE),
        info = paste(workflow, entry_point)
      )
      expect_false(has_namespaced_call(code, entry_point), info = paste(workflow, entry_point))
    }
  }
})

test_that("NormalizeData uses the fast path for Assay5 objects", {
  set.seed(11)
  mat <- Matrix::rsparsematrix(20, 10, density = 0.2)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)

  expect_true(inherits(obj[["RNA"]], "StdAssay"))
  expect_false("data" %in% SeuratObject::Layers(obj[["RNA"]], search = NA))
  out <- NormalizeData(obj, verbose = FALSE)

  expect_true(inherits(out[["RNA"]], "StdAssay"))
  expect_true("data" %in% SeuratObject::Layers(out[["RNA"]], search = NA))
  expect_equal(
    as.matrix(GetAssayData5(out, assay = "RNA", layer = "data")),
    as.matrix(seurat_reference_method(
      "NormalizeData", "default", mat, verbose = FALSE
    )),
    tolerance = 1e-8
  )
})

test_that("GetAssayData5 bypasses JoinLayers for an exact single layer", {
  obj <- make_fast_path_object()
  expected <- SeuratObject::GetAssayData(obj[["RNA"]], layer = "data")
  testthat::local_mocked_bindings(
    JoinLayers = function(...) stop("JoinLayers should not run for a single exact layer"),
    .package = "SeuratObject"
  )

  actual <- GetAssayData5(obj, assay = "RNA", layer = "data")
  expect_identical(actual, expected)
})

test_that("GetAssayData5 still joins split Assay5 layers", {
  obj <- make_fast_path_object()
  obj[["RNA"]] <- split(obj[["RNA"]], f = rep(c("A", "B"), length.out = ncol(obj)))
  expected <- SeuratObject::GetAssayData(
    SeuratObject::JoinLayers(obj[["RNA"]]),
    layer = "data"
  )

  actual <- GetAssayData5(obj, assay = "RNA", layer = "data")
  expect_identical(actual, expected)
})

test_that("FindNeighbors fast path stores Seurat-compatible graphs", {
  obj <- make_fast_path_object()
  out <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:10,
    k.param = 10,
    graph.name = c("RNA_nn", "RNA_snn"),
    verbose = FALSE
  )

  expect_true(all(c("RNA_nn", "RNA_snn") %in% names(out@graphs)))
  expect_s4_class(out@graphs$RNA_nn, "Graph")
  expect_s4_class(out@graphs$RNA_snn, "Graph")
  expect_equal(dim(out@graphs$RNA_snn), rep(ncol(out), 2))
  expect_equal(out@graphs$RNA_snn@assay.used, "RNA")
})

test_that("RunPCA fast path stores total variance and cell embeddings", {
  obj <- make_fast_path_object()
  expect_true("total.variance" %in% names(obj@reductions$pca@misc))
  expect_true(is.finite(obj@reductions$pca@misc$total.variance))
  expect_equal(nrow(SeuratObject::Embeddings(obj[["pca"]])), ncol(obj))
})

test_that("NormalizeData supports legacy Assay objects", {
  mat <- Matrix::rsparsematrix(20, 10, density = 0.2)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  suppressWarnings(
    obj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = mat)
  )

  expect_s4_class(obj[["RNA"]], "Assay")
  expect_no_error(out <- NormalizeData(obj, verbose = FALSE))
  expect_s4_class(out[["RNA"]], "Assay")
  expect_gt(length(out[["RNA"]]@data@x), 0)
  expect_equal(
    as.matrix(out[["RNA"]]@data),
    as.matrix(seurat_reference_method(
      "NormalizeData", "default", mat, verbose = FALSE
    )),
    tolerance = 1e-8
  )
})

test_that("Seurat fast paths support legacy Assay preprocessing", {
  set.seed(12)
  mat <- Matrix::rsparsematrix(80, 40, density = 0.15)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  suppressWarnings(
    obj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = mat)
  )
  SeuratObject::Idents(obj) <- rep(c("A", "B"), length.out = ncol(obj))

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(
    obj,
    features = SeuratObject::VariableFeatures(obj),
    npcs = 10,
    verbose = FALSE
  )
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:10, k.param = 10, verbose = FALSE)
  markers <- FindMarkers(
    obj,
    ident.1 = "A",
    ident.2 = "B",
    verbose = FALSE
  )
  all_markers <- suppressWarnings(FindAllMarkers(obj, verbose = FALSE))

  expect_s4_class(obj[["RNA"]], "Assay")
  expect_equal(nrow(obj[["RNA"]]@scale.data), length(SeuratObject::VariableFeatures(obj)))
  expect_true("pca" %in% SeuratObject::Reductions(obj))
  expect_true(all(c("RNA_nn", "RNA_snn") %in% names(obj@graphs)))
  expect_s3_class(markers, "data.frame")
  expect_s3_class(all_markers, "data.frame")
})

test_that("FindAllMarkers supports feature subsets with Seurat parity", {
  obj <- make_fast_path_object(seed = 13)
  SeuratObject::Idents(obj) <- rep(c("A", "B", "C", "D"), each = 20)
  features <- rownames(obj)[1:40]

  expected <- suppressWarnings(seurat_reference_find_all_markers(
    obj,
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    return.thresh = Inf,
    verbose = FALSE
  ))
  actual <- suppressWarnings(FindAllMarkers(
    obj,
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    return.thresh = Inf,
    verbose = FALSE
  ))

  # Seurat and the scop implementation may order tied rows differently
  # across compilers and platforms; align by (cluster, gene) so the
  # parity comparison is row-order independent.
  expected <- expected[order(expected$cluster, expected$gene), ]
  actual <- actual[order(actual$cluster, actual$gene), ]
  rownames(expected) <- NULL
  rownames(actual) <- NULL
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("FindAllMarkers native path covers metadata groups and FC controls", {
  skip_if(
    is.null(presto_get_fun(install = FALSE, error_on_missing = FALSE)),
    "The runtime-optional Presto backend is unavailable"
  )
  set.seed(20260830)
  counts <- matrix(stats::rpois(150L * 180L, lambda = 1), nrow = 150L)
  groups <- factor(
    rep(c("zeta", "alpha", "mu"), each = 60L),
    levels = c("zeta", "alpha", "mu")
  )
  counts[1:12, groups == "zeta"] <- counts[1:12, groups == "zeta"] + 4L
  counts[13:24, groups == "alpha"] <- counts[13:24, groups == "alpha"] + 4L
  counts[25:36, groups == "mu"] <- counts[25:36, groups == "mu"] + 4L
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object$marker_group <- groups
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  args <- list(
    object = object,
    group.by = "marker_group",
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = 0.05,
    base = exp(1),
    fc.name = "effect",
    return.thresh = Inf,
    only.pos = FALSE,
    verbose = FALSE
  )
  expected <- suppressWarnings(do.call(seurat_reference_find_all_markers, args))
  native_calls <- 0L
  native <- get("marker_all_from_context", asNamespace("scop"))
  testthat::local_mocked_bindings(
    marker_all_from_context = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .package = "scop"
  )
  actual <- suppressWarnings(do.call(FindAllMarkers, args))
  expected <- expected[order(expected$cluster, expected$gene), ]
  actual <- actual[order(actual$cluster, actual$gene), ]
  rownames(expected) <- rownames(actual) <- NULL
  expect_identical(colnames(actual), colnames(expected))
  expect_equal(actual$effect, expected$effect, tolerance = 1e-12)
  expect_equal(actual$pct.1, expected$pct.1, tolerance = 0)
  expect_equal(actual$pct.2, expected$pct.2, tolerance = 0)
  expect_equal(actual$p_val_adj, expected$p_val_adj, tolerance = 5e-10)
  expect_equal(actual$p_val, expected$p_val, tolerance = 5e-7)
  expect_identical(native_calls, 1L)
})

test_that("FindAllMarkers accelerates Seurat-compatible sampling", {
  skip_if(
    is.null(presto_get_fun(install = FALSE, error_on_missing = FALSE)),
    "The runtime-optional Presto backend is unavailable"
  )
  object <- make_fast_path_object(seed = 14)
  SeuratObject::Idents(object) <- rep(c("A", "B", "C", "D"), each = 20)
  args <- list(
    object = object,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    max.cells.per.ident = 12L,
    random.seed = 42L,
    return.thresh = Inf,
    verbose = FALSE
  )
  expected <- suppressWarnings(do.call(seurat_reference_find_all_markers, args))
  native_calls <- 0L
  native <- get("marker_all_from_context", asNamespace("scop"))
  testthat::local_mocked_bindings(
    marker_all_from_context = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .package = "scop"
  )
  actual <- suppressWarnings(do.call(FindAllMarkers, args))
  expected <- expected[order(expected$cluster, expected$gene), ]
  actual <- actual[order(actual$cluster, actual$gene), ]
  rownames(expected) <- rownames(actual) <- NULL
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_identical(native_calls, 1L)
})

test_that("sparse Wilcoxon returns finite p-values for zero-variance genes", {
  counts <- Matrix::Matrix(1, nrow = 5, ncol = 12, sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  SeuratObject::Idents(obj) <- rep(c("A", "B"), each = 6)

  markers <- FindMarkers(
    obj,
    ident.1 = "A",
    ident.2 = "B",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )

  expect_false(any(is.nan(markers$p_val)))
  expect_true(all(markers$p_val == 1))
})

test_that("FindAllMarkers handles a single identity like Seurat", {
  object <- make_fast_path_object(seed = 140)
  SeuratObject::Idents(object) <- rep("A", ncol(object))
  args <- list(
    object = object,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = Inf,
    verbose = FALSE
  )
  expected <- suppressWarnings(do.call(seurat_reference_find_all_markers, args))
  actual <- suppressWarnings(do.call(FindAllMarkers, args))
  expect_identical(actual, expected)
})

test_that("AddModuleScore sparse path preserves Seurat sampling and scores", {
  set.seed(20260831)
  counts <- Matrix::rsparsematrix(600, 240, density = 0.08)
  counts@x <- stats::rpois(length(counts@x), lambda = 3) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  features <- list(
    alpha = paste0("g", 1:25),
    beta = paste0("g", 101:130),
    gamma = paste0("g", 301:320)
  )
  args <- list(
    object = object,
    features = features,
    pool = paste0("g", 1:550),
    nbin = 12,
    ctrl = 10,
    name = "program_",
    seed = 42
  )
  expected <- do.call(seurat_reference_add_module_score, args)
  native_calls <- 0L
  native <- get("module_score_native", asNamespace("scop"))
  testthat::local_mocked_bindings(
    module_score_native = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .package = "scop"
  )
  actual <- do.call(AddModuleScore, args)
  columns <- paste0("program_", 1:3)
  expect_equal(actual[[columns]], expected[[columns]], tolerance = 1e-14)
  expect_identical(native_calls, 1L)
})

test_that("AddModuleScore preserves NULL-seed RNG behavior and fallback branches", {
  object <- make_fast_path_object(seed = 15)
  features <- list(rownames(object)[1:8], rownames(object)[21:28])
  set.seed(99)
  expected <- seurat_reference_add_module_score(
    object,
    features = features,
    nbin = 4,
    ctrl = 3,
    seed = NULL,
    name = "nullseed"
  )
  set.seed(99)
  actual <- AddModuleScore(
    object,
    features = features,
    nbin = 4,
    ctrl = 3,
    seed = NULL,
    name = "nullseed"
  )
  expect_equal(
    actual[[c("nullseed1", "nullseed2")]],
    expected[[c("nullseed1", "nullseed2")]],
    tolerance = 1e-14
  )

  fallback_cases <- list(
    list(features = list(c(features[[1]], "missing"))),
    list(features = list(c(features[[1]], "missing")), search = TRUE)
  )
  for (case in fallback_cases) {
    expected_fallback <- suppressWarnings(do.call(
      seurat_reference_add_module_score,
      c(list(object = object, nbin = 4, ctrl = 3), case)
    ))
    actual_fallback <- suppressWarnings(do.call(
      AddModuleScore,
      c(list(object = object, nbin = 4, ctrl = 3), case)
    ))
    score_columns <- grep("Cluster", colnames(expected_fallback[[]]), value = TRUE)
    expect_equal(
      actual_fallback[[score_columns]],
      expected_fallback[[score_columns]],
      tolerance = 1e-14
    )
  }
})

test_that("AddModuleScore deduplicates features and delegates fractional nbin", {
  object <- make_fast_path_object(seed = 151)
  features <- list(c(rownames(object)[1:8], rownames(object)[1:2]))
  native_calls <- 0L
  native <- get("module_score_native", asNamespace("scop"))
  testthat::local_mocked_bindings(
    module_score_native = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .package = "scop"
  )

  for (nbin in c(4, 4.9)) {
    args <- list(
      object = object,
      features = features,
      nbin = nbin,
      ctrl = 3,
      seed = 1,
      name = "edge"
    )
    expected <- do.call(seurat_reference_add_module_score, args)
    actual <- do.call(AddModuleScore, args)
    expect_equal(actual$edge1, expected$edge1, tolerance = 1e-14)
  }
  expect_identical(native_calls, 1L)
})

test_that("AddModuleScore accelerates inert search/dots and split layers", {
  object <- make_fast_path_object(seed = 17)
  object$batch <- rep(c("a", "b"), each = ncol(object) / 2)
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
  features <- list(rownames(object)[1:8], rownames(object)[21:28])
  args <- list(
    object = object,
    features = features,
    nbin = 4,
    ctrl = 3,
    search = TRUE,
    custom.argument = TRUE,
    seed = 42,
    name = "split"
  )
  expected <- seurat_reference_add_module_score(
    object = object,
    features = features,
    nbin = 4,
    ctrl = 3,
    search = TRUE,
    custom.argument = TRUE,
    seed = 42,
    name = "split"
  )
  native_calls <- 0L
  native <- get("module_score_native", asNamespace("scop"))
  testthat::local_mocked_bindings(
    module_score_native = function(...) {
      native_calls <<- native_calls + 1L
      native(...)
    },
    .package = "scop"
  )
  actual <- AddModuleScore(
    object = object,
    features = features,
    nbin = 4,
    ctrl = 3,
    search = TRUE,
    custom.argument = TRUE,
    seed = 42,
    name = "split"
  )
  expect_equal(actual[[c("split1", "split2")]], expected[[c("split1", "split2")]])
  expect_identical(native_calls, 1L)
})

test_that("CellCycleScoring preserves Seurat scores, phases, and identities", {
  set.seed(20260901)
  counts <- Matrix::rsparsematrix(500, 180, density = 0.1)
  counts@x <- stats::rpois(length(counts@x), lambda = 4) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  s_features <- paste0("g", 1:30)
  g2m_features <- paste0("g", 101:135)
  expected <- seurat_reference_cell_cycle_scoring(
    object,
    s.features = s_features,
    g2m.features = g2m_features,
    ctrl = 10,
    set.ident = TRUE,
    seed = 7,
    nbin = 10
  )
  actual <- CellCycleScoring(
    object,
    s.features = s_features,
    g2m.features = g2m_features,
    ctrl = 10,
    set.ident = TRUE,
    seed = 7,
    nbin = 10
  )
  expect_equal(
    actual[[c("S.Score", "G2M.Score")]],
    expected[[c("S.Score", "G2M.Score")]],
    tolerance = 1e-14
  )
  expect_identical(actual$Phase, expected$Phase)
  expect_identical(as.character(SeuratObject::Idents(actual)), as.character(SeuratObject::Idents(expected)))
  expect_identical(as.character(actual$old.ident), as.character(expected$old.ident))
})

test_that("SelectIntegrationFeatures vectorizes Seurat's exact rank contract", {
  skip_if_not_installed("matrixStats")
  set.seed(20260903)
  genes <- paste0("g", seq_len(300))
  objects <- lapply(seq_len(8), function(i) {
    counts <- Matrix::rsparsematrix(300, 30, density = 0.04)
    counts@x <- stats::rpois(length(counts@x), lambda = 2) + 1
    dimnames(counts) <- list(genes, paste0("c", i, "_", seq_len(30)))
    object <- Seurat::CreateSeuratObject(counts)
    SeuratObject::VariableFeatures(object) <- sample(genes, 80)
    object
  })

  for (nfeatures in c(1, 25, 80, 200, 1000)) {
    expected <- Seurat::SelectIntegrationFeatures(
      objects,
      nfeatures = nfeatures,
      verbose = FALSE
    )
    actual <- SelectIntegrationFeatures(
      objects,
      nfeatures = nfeatures,
      verbose = FALSE
    )
    expect_identical(actual, expected)
  }
  expect_identical(
    SelectIntegrationFeatures(
      objects,
      nfeatures = 80,
      assay = rep("RNA", length(objects)),
      verbose = FALSE
    ),
    Seurat::SelectIntegrationFeatures(
      objects,
      nfeatures = 80,
      assay = rep("RNA", length(objects)),
      verbose = FALSE
    )
  )
})

test_that("module and cell-cycle scoring cover remaining parameter branches", {
  object <- make_fast_path_object(seed = 18)
  features <- list(rownames(object)[1:8], rownames(object)[21:28])
  args <- list(
    object = object,
    features = features,
    assay = "RNA",
    slot = "counts",
    nbin = 4,
    ctrl = 3,
    name = "countscore",
    seed = 8,
    search = FALSE
  )
  expected <- do.call(get("AddModuleScore.Seurat", asNamespace("Seurat")), args)
  actual <- do.call(AddModuleScore, args)
  expect_equal(actual[[c("countscore1", "countscore2")]], expected[[c("countscore1", "countscore2")]], tolerance = 1e-14)

  expect_error(AddModuleScore(object, features = features, k = TRUE))
  expect_error(get("AddModuleScore.Seurat", asNamespace("Seurat"))(
    object,
    features = features,
    k = TRUE
  ))

  cc_expected <- seurat_reference_cell_cycle_scoring(
    object,
    s.features = features[[1]],
    g2m.features = features[[2]],
    ctrl = NULL,
    set.ident = FALSE,
    nbin = 4,
    seed = 9
  )
  cc_actual <- CellCycleScoring(
    object,
    s.features = features[[1]],
    g2m.features = features[[2]],
    ctrl = NULL,
    set.ident = FALSE,
    nbin = 4,
    seed = 9
  )
  expect_equal(cc_actual[[c("S.Score", "G2M.Score")]], cc_expected[[c("S.Score", "G2M.Score")]], tolerance = 1e-14)
  expect_identical(cc_actual$Phase, cc_expected$Phase)
})

test_that("clustering and multimodal wrappers delegate all arguments", {
  set.seed(20260902)
  cells <- paste0("c", seq_len(60L))
  adjacency <- Matrix::rsparsematrix(60L, 60L, density = 0.08)
  adjacency@x[] <- 1
  adjacency <- as(adjacency + Matrix::t(adjacency) > 0, "dgCMatrix")
  diag(adjacency) <- 0
  rownames(adjacency) <- colnames(adjacency) <- cells
  graph <- SeuratObject::as.Graph(adjacency)
  cluster_args <- list(
    object = graph,
    resolution = c(0.2, 0.8),
    algorithm = 1,
    n.start = 1,
    n.iter = 2,
    random.seed = 7,
    verbose = FALSE
  )
  expect_identical(
    do.call(FindClusters, cluster_args),
    do.call(Seurat::FindClusters, cluster_args)
  )

  counts <- Matrix::rsparsematrix(100L, 60L, density = 0.1)
  counts@x <- abs(counts@x) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- cells
  object <- Seurat::CreateSeuratObject(counts)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(stats::rnorm(60L * 12L), nrow = 60L, dimnames = list(cells, paste0("PC_", 1:12))),
    key = "PC_",
    assay = "RNA"
  )
  object[["apca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(stats::rnorm(60L * 10L), nrow = 60L, dimnames = list(cells, paste0("APC_", 1:10))),
    key = "APC_",
    assay = "RNA"
  )
  multimodal_args <- list(
    object = object,
    reduction.list = c("pca", "apca"),
    dims.list = list(1:10, 1:8),
    k.nn = 5,
    knn.range = 20,
    modality.weight.name = c("pca.weight", "apca.weight"),
    verbose = FALSE
  )
  expected <- do.call(Seurat::FindMultiModalNeighbors, multimodal_args)
  actual <- do.call(FindMultiModalNeighbors, multimodal_args)
  expect_equal(as(actual[["wknn"]], "dgCMatrix"), as(expected[["wknn"]], "dgCMatrix"), tolerance = 0)
  expect_equal(as(actual[["wsnn"]], "dgCMatrix"), as(expected[["wsnn"]], "dgCMatrix"), tolerance = 0)
  expect_equal(actual$pca.weight, expected$pca.weight, tolerance = 0)
  expect_equal(actual$apca.weight, expected$apca.weight, tolerance = 0)
})

test_that("clustering delegation covers modularity and Leiden branches", {
  set.seed(20260907)
  cells <- paste0("c", seq_len(40L))
  adjacency <- Matrix::sparseMatrix(
    i = rep(seq_len(40L), 2L),
    j = c(c(2:40, 1), c(3:40, 1, 2)),
    x = 1,
    dims = c(40L, 40L)
  )
  adjacency <- pmax(adjacency, Matrix::t(adjacency))
  rownames(adjacency) <- colnames(adjacency) <- cells
  graph <- SeuratObject::as.Graph(adjacency)
  cases <- list(
    list(algorithm = 1, modularity.fxn = 2, resolution = c(0.3, 0.7), group.singletons = FALSE),
    list(algorithm = 2, resolution = 0.4),
    list(
      algorithm = "leiden",
      leiden_method = "igraph",
      leiden_objective_function = "CPM",
      resolution = 0.4
    )
  )
  for (case in cases) {
    args <- c(
      list(object = graph, n.start = 1, n.iter = 2, random.seed = 9, verbose = FALSE),
      case
    )
    expect_identical(do.call(FindClusters, args), do.call(Seurat::FindClusters, args))
  }
})

test_that("clustering forwards Leidenbase-only parameters unchanged", {
  captured <- NULL
  testthat::local_mocked_bindings(
    FindClusters = function(object, ...) {
      captured <<- c(list(object = object), list(...))
      "delegated"
    },
    .package = "Seurat"
  )
  initial <- c(0L, 1L)
  sizes <- c(2, 3)
  expect_identical(
    FindClusters(
      object = "graph",
      algorithm = 4,
      leiden_method = "leidenbase",
      leiden_objective_function = "modularity",
      initial.membership = initial,
      node.sizes = sizes,
      temp.file.location = "tmp",
      edge.file.name = "edges",
      verbose = FALSE
    ),
    "delegated"
  )
  expect_identical(captured$initial.membership, initial)
  expect_identical(captured$node.sizes, sizes)
  expect_identical(captured$temp.file.location, "tmp")
  expect_identical(captured$edge.file.name, "edges")
})

test_that("multimodal delegation covers optional weight branches", {
  set.seed(20260908)
  cells <- paste0("c", seq_len(50L))
  counts <- Matrix::sparseMatrix(
    i = sample.int(60L, 250L, replace = TRUE),
    j = sample.int(50L, 250L, replace = TRUE),
    x = 1,
    dims = c(60L, 50L),
    dimnames = list(paste0("g", 1:60), cells)
  )
  object <- Seurat::CreateSeuratObject(counts)
  for (spec in list(c("pca", "PC_", 8L), c("apca", "APC_", 7L))) {
    object[[spec[[1L]]]] <- SeuratObject::CreateDimReducObject(
      embeddings = matrix(
        stats::rnorm(50L * as.integer(spec[[3L]])),
        nrow = 50L,
        dimnames = list(cells, paste0(spec[[2L]], seq_len(as.integer(spec[[3L]]))))
      ),
      key = spec[[2L]],
      assay = "RNA"
    )
  }
  args <- list(
    object = object,
    reduction.list = c("pca", "apca"),
    dims.list = list(1:8, 1:7),
    k.nn = 5,
    l2.norm = FALSE,
    knn.graph.name = "custom_knn",
    snn.graph.name = "custom_snn",
    weighted.nn.name = "custom_nn",
    modality.weight.name = c("w1", "w2"),
    knn.range = 15,
    prune.SNN = 0.1,
    sd.scale = 0.8,
    cross.contant.list = list(1e-3, 2e-3),
    smooth = FALSE,
    return.intermediate = TRUE,
    verbose = FALSE
  )
  expected <- do.call(Seurat::FindMultiModalNeighbors, args)
  actual <- do.call(FindMultiModalNeighbors, args)
  expect_equal(actual[["custom_knn"]], expected[["custom_knn"]], tolerance = 0)
  expect_equal(actual[["custom_snn"]], expected[["custom_snn"]], tolerance = 0)
  weights <- SeuratObject::Misc(expected, slot = "modality.weight")
  args$object <- object
  args$modality.weight <- weights
  args$return.intermediate <- FALSE
  expect_equal(
    do.call(FindMultiModalNeighbors, args)[["custom_snn"]],
    do.call(Seurat::FindMultiModalNeighbors, args)[["custom_snn"]],
    tolerance = 0
  )
})

test_that("multimodal delegation forwards smooth and supplied weights", {
  captured <- NULL
  weights <- structure(list(), class = "ModalityWeights")
  testthat::local_mocked_bindings(
    FindMultiModalNeighbors = function(object, ...) {
      captured <<- c(list(object = object), list(...))
      "delegated"
    },
    .package = "Seurat"
  )
  expect_identical(
    FindMultiModalNeighbors(
      object = "object",
      reduction.list = c("pca", "apca"),
      dims.list = list(1:2, 1:2),
      smooth = TRUE,
      modality.weight = weights,
      return.intermediate = TRUE,
      verbose = FALSE
    ),
    "delegated"
  )
  expect_true(captured$smooth)
  expect_identical(captured$modality.weight, weights)
  expect_true(captured$return.intermediate)
})
