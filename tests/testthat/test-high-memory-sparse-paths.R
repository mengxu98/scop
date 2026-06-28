test_that("scFEA comparison statistics match dense row SDs", {
  mat <- Matrix::Matrix(
    c(
      0, 1, 0, 3,
      2, 0, 0, 4,
      0, 5, 6, 0
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(mat) <- paste0("m", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))

  sparse_sd <- MatrixGenerics::rowSds(mat)
  dense_sd <- apply(as_matrix(mat), 1, stats::sd)

  expect_equal(sparse_sd, dense_sd)
})

test_that("SCENICPlus sparse coercion preserves values and dimnames", {
  mat <- Matrix::rsparsematrix(8, 5, density = 0.25)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))

  out <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")

  expect_s4_class(out, "dgCMatrix")
  expect_equal(dimnames(out), dimnames(mat))
  expect_equal(as_matrix(out), as_matrix(mat))
})

test_that("spatial sparse matrix coercion preserves matrix values", {
  mat <- matrix(c(1, 0, NA, 4, Inf, 0), nrow = 2)
  rownames(mat) <- c("g1", "g2")
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))

  out <- scop:::spatial_integration_sparse_matrix(mat)

  expect_s4_class(out, "dgCMatrix")
  expect_equal(dimnames(out), dimnames(mat))
  expect_true(all(is.finite(out@x)))
  expect_equal(as_matrix(out), matrix(
    c(1, 0, 0, 4, 0, 0),
    nrow = 2,
    dimnames = dimnames(mat)
  ))
})

test_that("spatial integration prepares merged objects without SplitObject", {
  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 4,
      5, 0, 0, 6
    ),
    nrow = 3,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$sample <- rep(c("s1", "s2"), each = 2)
  srt$col <- c(1, 2, 1, 2)
  srt$row <- c(1, 1, 2, 2)

  testthat::local_mocked_bindings(
    SplitObject = function(...) stop("SplitObject should not be called"),
    .package = "Seurat"
  )

  input <- scop:::spatial_integration_prepare_input(
    object = srt,
    sample.by = "sample",
    assay = "RNA",
    layer = "counts",
    features = NULL,
    image = NULL,
    coord.cols = c("col", "row")
  )

  expect_s4_class(input$expr, "dgCMatrix")
  expect_equal(input$samples, c("s1", "s2"))
  expect_equal(names(input$expr_list), c("s1", "s2"))
})

test_that("RunUMAP2 crops requested dims to available embeddings", {
  counts <- Matrix::Matrix(
    c(
      1, 0, 2, 0, 3, 0,
      0, 2, 0, 3, 0, 4,
      5, 0, 4, 0, 3, 0,
      0, 1, 0, 2, 0, 3,
      4, 0, 3, 0, 2, 0
    ),
    nrow = 5,
    byrow = TRUE,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  emb <- matrix(
    stats::rnorm(ncol(srt) * 3),
    nrow = ncol(srt),
    dimnames = list(colnames(srt), paste0("PC_", 1:3))
  )
  srt[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = emb,
    key = "PC_",
    assay = "RNA"
  )

  out <- RunUMAP2(
    srt,
    reduction = "pca",
    dims = 1:10,
    n.components = 2,
    n.neighbors = 3,
    n.epochs = 10,
    reduction.name = "crop_umap",
    reduction.key = "cropUMAP_",
    verbose = FALSE
  )

  expect_true("crop_umap" %in% SeuratObject::Reductions(out))
  expect_equal(ncol(Seurat::Embeddings(out, "crop_umap")), 2)
})

test_that("CIBERSORT matrix validation handles sparse matrices once", {
  mat <- Matrix::Matrix(
    c(1, 0, 3, 4, 5, 6),
    nrow = 2,
    sparse = TRUE
  )
  rownames(mat) <- c("GeneA", "GeneB")
  colnames(mat) <- c("Sample1", "Sample2", "Sample3")

  out <- scop:::cibersort_check_matrix(mat, "count_matrix")

  expect_type(out, "double")
  expect_equal(dimnames(out), dimnames(mat))
  expect_equal(out, as_matrix(mat))
})

test_that("scTenifold network plot subsets sparse networks before dense conversion", {
  net <- Matrix::Matrix(
    c(
      0, 1, 0, 0,
      2, 0, 3, 0,
      0, 4, 0, 5,
      0, 0, 6, 0
    ),
    nrow = 4,
    byrow = TRUE,
    sparse = TRUE
  )
  dimnames(net) <- list(paste0("g", 1:4), paste0("g", 1:4))
  bundle <- list(
    result = list(tensorNetworks = list(WT = net)),
    parameters = list(gKO = "g1")
  )
  dr <- data.frame(
    gene = paste0("g", 1:4),
    FC = c(4, 3, 2, 1),
    p.adj = c(0.01, 0.02, 0.5, 0.8),
    status = factor(c("Knockout", "Significant", "NS", "NS"),
      levels = c("Knockout", "Significant", "NS")
    )
  )

  p <- scop:::sctenifold_plot_network(
    bundle = bundle,
    dr = dr,
    top_n = 3,
    features = NULL,
    label = FALSE,
    pt.size = 2,
    label.size = 3,
    edge_top_n = 5,
    edge_threshold = NULL,
    cols.sig = "#D7301F",
    cols.ns = "grey70",
    cols.ko = "#4575B4",
    title = NULL,
    theme_use = ggplot2::theme_void(),
    theme_args = list()
  )

  expect_s3_class(p, "ggplot")
})

test_that("heatmaps skip temporary row-order heatmap when labels are disabled", {
  skip_if_not_installed("ComplexHeatmap")
  counts <- Matrix::Matrix(
    c(
      1, 0, 3, 0,
      2, 1, 4, 0,
      0, 5, 1, 0
    ),
    nrow = 4,
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(c("A", "A", "B"))

  testthat::local_mocked_bindings(
    row_order = function(...) stop("row_order should not be called"),
    .package = "ComplexHeatmap"
  )

  fh <- FeatureHeatmap(
    srt,
    features = rownames(srt),
    group.by = "group",
    nlabel = 0,
    verbose = FALSE
  )
  gh <- GroupHeatmap(
    srt,
    features = rownames(srt),
    group.by = "group",
    nlabel = 0,
    verbose = FALSE
  )

  expect_s3_class(fh$plot, "ggplot")
  expect_s3_class(gh$plot, "ggplot")
})

test_that("Palantir boundary distance shortcuts match full distance matrices", {
  ms_data <- matrix(
    c(
      0, 0,
      1, 0,
      0, 2,
      2, 2,
      3, 1
    ),
    ncol = 2,
    byrow = TRUE
  )
  early_cell <- 1L
  dm_boundaries <- c(2L, 3L, 5L)

  full_ec <- as.matrix(stats::dist(rbind(
    ms_data[early_cell, , drop = FALSE],
    ms_data[dm_boundaries, , drop = FALSE]
  )))
  shortcut_ec <- sqrt(rowSums(
    sweep(ms_data[dm_boundaries, , drop = FALSE], 2, ms_data[early_cell, ], "-")^2
  ))

  expect_equal(unname(shortcut_ec), unname(full_ec[1, -1]))
  expect_identical(dm_boundaries[which.min(shortcut_ec)], dm_boundaries[which.min(full_ec[1, -1])])

  high_rank_ms <- ms_data[c(1L, 4L), , drop = FALSE]
  boundary_ms <- ms_data[c(2L, 3L, 5L), , drop = FALSE]
  full_boundary <- as.matrix(stats::dist(rbind(high_rank_ms, boundary_ms)))
  dist_part_full <- full_boundary[
    seq_len(nrow(high_rank_ms)),
    nrow(high_rank_ms) + seq_len(nrow(boundary_ms)),
    drop = FALSE
  ]
  dist_part_shortcut <- vapply(
    seq_len(nrow(boundary_ms)),
    function(j) sqrt(rowSums(sweep(high_rank_ms, 2, boundary_ms[j, ], "-")^2)),
    numeric(nrow(high_rank_ms))
  )

  expect_equal(unname(dist_part_shortcut), unname(dist_part_full))
  expect_equal(
    unname(apply(dist_part_shortcut, 1, which.min)),
    unname(apply(dist_part_full, 1, which.min))
  )
})
