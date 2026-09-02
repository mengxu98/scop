test_that("RunDEtest uses scop FindMarkers-compatible sparse Wilcoxon path", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::Matrix(rpois(80 * 30, lambda = 1), nrow = 80, sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  SeuratObject::Idents(srt) <- rep(c("A", "B"), each = 15)

  cells1 <- colnames(srt)[1:15]
  cells2 <- colnames(srt)[16:30]
  features <- rownames(srt)[1:50]

  seurat_fm <- get("FindMarkers.Seurat", asNamespace("Seurat"))(
    srt,
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  scop_fm <- scop::FindMarkers(
    srt,
    ident.1 = "A",
    ident.2 = "B",
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  cell_fm <- scop::FindMarkers(
    srt,
    cells.1 = cells1,
    cells.2 = cells2,
    features = features,
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE
  )
  out <- RunDEtest(
    srt,
    cells1 = cells1,
    cells2 = cells2,
    features = features,
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )
  markers <- out@tools$DEtest_custom$AllMarkers_wilcox

  expect_equal(scop_fm, seurat_fm)
  expect_equal(cell_fm, scop_fm)
  expect_identical(rownames(cell_fm), markers$gene)
  expect_equal(cell_fm$p_val, markers$p_val, tolerance = 1e-12)
  expect_equal(cell_fm$avg_log2FC, markers$avg_log2FC, tolerance = 1e-12)
  expect_equal(cell_fm$pct.1, markers$pct.1, tolerance = 0)
  expect_equal(cell_fm$pct.2, markers$pct.2, tolerance = 0)
})

test_that("FindMarkers native branches exactly match the Seurat method", {
  set.seed(20260909)
  counts <- Matrix::rsparsematrix(
    300,
    150,
    density = 0.1,
    rand.x = function(n) stats::rpois(n, 2) + 1
  )
  dimnames(counts) <- list(paste0("g", 1:300), paste0("c", 1:150))
  object <- Seurat::CreateSeuratObject(counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  SeuratObject::Idents(object) <- factor(rep(c("A", "B", "C"), each = 50))
  seurat_method <- get("FindMarkers.Seurat", asNamespace("Seurat"))
  cases <- list(
    list(ident.1 = "A"),
    list(ident.1 = "A", ident.2 = c("B", "C")),
    list(ident.1 = "A", min.diff.pct = 0.05),
    list(ident.1 = "A", max.cells.per.ident = 25, random.seed = 8),
    list(ident.1 = "A", min.cells.group = 5, min.cells.feature = 5),
    list(ident.1 = "A", only.pos = TRUE),
    list(ident.1 = "A", features = paste0("g", 1:80)),
    list(ident.1 = "A", slot = "counts")
  )
  for (case in cases) {
    args <- c(
      list(object = object, logfc.threshold = 0, min.pct = 0, verbose = FALSE),
      case
    )
    expect_equal(do.call(scop::FindMarkers, args), do.call(seurat_method, args), tolerance = 1e-12)
  }
})

test_that("RunDEtest exports internal marker workers to parallel processes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(2)
  counts <- Matrix::Matrix(rpois(50 * 30, lambda = 2), nrow = 50, sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- rep(c("A", "B", "C"), each = 10)

  out <- RunDEtest(
    srt,
    group.by = "group",
    markers_type = "paired",
    features = rownames(srt)[1:30],
    cores = 2,
    fc.threshold = 1,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  )

  markers <- out@tools$DEtest_group$PairedMarkers_wilcox
  expect_s3_class(markers, "data.frame")
  expect_true(all(c("gene", "group1", "group2") %in% colnames(markers)))
})

test_that("supported RunDEtest Wilcoxon branches stay on the scop backend", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if(
    is.null(presto_get_fun(install = FALSE, error_on_missing = FALSE)),
    "The runtime-optional Presto backend is unavailable"
  )

  set.seed(3)
  counts <- Matrix::rsparsematrix(
    100,
    48,
    density = 0.15,
    rand.x = function(n) stats::rpois(n, 2) + 1
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 24))
  srt$batch <- factor(rep(rep(c("x", "y"), each = 12), 2))
  SeuratObject::Idents(srt) <- srt$group
  cells1 <- colnames(srt)[srt$group == "A"]
  cells2 <- colnames(srt)[srt$group == "B"]
  features <- rownames(srt)[1:60]

  testthat::local_mocked_bindings(
    FindMarkers = function(...) stop("Seurat marker backend was reached"),
    .package = "Seurat"
  )

  expect_no_error(pair <- scop::FindMarkers(
    object = srt,
    assay = "RNA",
    layer = "data",
    cells.1 = cells1,
    cells.2 = cells2,
    features = features,
    test.use = "wilcox",
    logfc.threshold = 0,
    base = 2,
    min.pct = 0,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    min.cells.feature = 3,
    min.cells.group = 3,
    latent.vars = NULL,
    only.pos = FALSE,
    norm.method = "LogNormalize",
    pseudocount.use = 1,
    mean.fxn = NULL,
    verbose = FALSE
  ))
  expect_s3_class(pair, "data.frame")

  expect_no_error(conserved <- get("FindConservedMarkers2", asNamespace("scop"))(
    object = srt,
    grouping.var = "batch",
    cells.1 = cells1,
    cells.2 = cells2,
    features = features,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = FALSE
  ))
  expect_s3_class(conserved, "data.frame")

  for (markers_type in c("all", "paired", "conserved", "disturbed")) {
    args <- list(
      object = srt,
      group.by = "group",
      markers_type = markers_type,
      features = features,
      fc.threshold = 1,
      min.pct = 0,
      only.pos = FALSE,
      cores = 1,
      verbose = FALSE
    )
    if (markers_type %in% c("conserved", "disturbed")) {
      args[["grouping.var"]] <- "batch"
    }
    expect_no_error(do.call(
      get("RunDEtest.Seurat", asNamespace("scop")),
      args
    ))
  }
})

test_that("unsupported FindMarkers arguments still use the Seurat fallback", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(4)
  counts <- Matrix::rsparsematrix(20, 12, density = 0.2)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  reached_seurat <- FALSE
  testthat::local_mocked_bindings(
    FindMarkers.Seurat = function(...) {
      reached_seurat <<- TRUE
      data.frame(row.names = rownames(srt)[1])
    },
    .package = "Seurat"
  )

  scop::FindMarkers(
    object = srt,
    assay = "RNA",
    layer = "counts",
    cells.1 = colnames(srt)[1:6],
    cells.2 = colnames(srt)[7:12],
    features = rownames(srt),
    test.use = "roc",
    logfc.threshold = 0,
    base = 2,
    min.pct = 0,
    min.diff.pct = -Inf,
    max.cells.per.ident = Inf,
    min.cells.feature = 3,
    min.cells.group = 3,
    latent.vars = NULL,
    only.pos = FALSE,
    norm.method = "LogNormalize",
    pseudocount.use = 1,
    mean.fxn = NULL,
    verbose = FALSE
  )

  expect_true(reached_seurat)
})

test_that("RunDEtest all-in-one fallback never installs missing Presto", {
  install_requested <- NULL
  native_all_markers <- get("RunDEtestFindAllMarkers", asNamespace("scop"))

  testthat::local_mocked_bindings(
    marker_assay_is_chromatin = function(...) FALSE,
    marker_context = function(...) list(labels = c("A", "B")),
    presto_get_fun = function(
      fun = "wilcoxauc",
      install = FALSE,
      error_on_missing = TRUE
    ) {
      install_requested <<- install
      NULL
    },
    .package = "scop"
  )

  expect_null(native_all_markers(
    srt = NULL,
    assay = "RNA",
    layer = "data",
    cell_group = factor(c("A", "B")),
    features = NULL,
    test.use = "wilcox",
    logfc.threshold = 0,
    base = 2,
    min.pct = 0,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    min.cells.group = 3,
    latent.vars = NULL,
    only.pos = FALSE,
    norm.method = "LogNormalize",
    pseudocount.use = 1,
    mean.fxn = NULL,
    p.adjust.method = "bonferroni"
  ))
  expect_false(install_requested)
})

test_that("RunDEtest all-in-one markers match the pairwise scop backend", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if(
    is.null(presto_get_fun(install = FALSE, error_on_missing = FALSE)),
    "The runtime-optional Presto backend is unavailable"
  )

  set.seed(5)
  counts <- Matrix::rsparsematrix(
    120,
    60,
    density = 0.15,
    rand.x = function(n) stats::rpois(n, 2) + 1
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  cell_group <- factor(rep(c("A", "B", "C"), each = 20))
  names(cell_group) <- colnames(srt)
  features <- rownames(srt)[1:80]

  native_all_markers <- get("RunDEtestFindAllMarkers", asNamespace("scop"))
  native_args <- list(
    srt = srt,
    assay = "RNA",
    layer = "data",
    cell_group = cell_group,
    features = features,
    test.use = "wilcox",
    logfc.threshold = 0,
    base = 2,
    min.pct = 0,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    min.cells.group = 3,
    latent.vars = NULL,
    only.pos = FALSE,
    norm.method = "LogNormalize",
    pseudocount.use = 1,
    mean.fxn = NULL,
    p.adjust.method = "bonferroni"
  )
  all_in_one <- do.call(native_all_markers, native_args)
  pairwise <- lapply(levels(cell_group), function(group) {
    markers <- scop::FindMarkers(
      object = srt,
      assay = "RNA",
      layer = "data",
      cells.1 = names(cell_group)[cell_group == group],
      cells.2 = names(cell_group)[cell_group != group],
      features = features,
      test.use = "wilcox",
      logfc.threshold = 0,
      base = 2,
      min.pct = 0,
      min.diff.pct = -Inf,
      max.cells.per.ident = Inf,
      min.cells.feature = 3,
      min.cells.group = 3,
      latent.vars = NULL,
      only.pos = FALSE,
      norm.method = "LogNormalize",
      pseudocount.use = 1,
      mean.fxn = NULL,
      verbose = FALSE
    )
    markers$gene <- rownames(markers)
    markers$group1 <- factor(group, levels = levels(cell_group))
    markers$group2 <- "others"
    markers
  })
  pairwise <- do.call(rbind, pairwise)
  pairwise <- pairwise[, colnames(all_in_one), drop = FALSE]
  # The all-in-one and pairwise paths may order tied rows differently
  # across compilers; align by (group1, gene) so the comparison is
  # row-order independent.
  all_in_one <- all_in_one[order(all_in_one$group1, all_in_one$gene), ]
  pairwise <- pairwise[order(pairwise$group1, pairwise$gene), ]
  rownames(pairwise) <- NULL
  rownames(all_in_one) <- NULL

  expect_equal(all_in_one, pairwise, tolerance = 1e-12)
})
