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
