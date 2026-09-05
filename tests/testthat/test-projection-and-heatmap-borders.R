projection_theme_with_annotation <- function(...) {
  list(
    list(ggplot2::annotation_custom(grid::nullGrob())),
    list(ggplot2::theme_void())
  )
}

make_projection_test_srt <- function(prefix) {
  counts <- matrix(
    1L,
    nrow = 2L,
    ncol = 4L,
    dimnames = list(c("gene1", "gene2"), paste0(prefix, seq_len(4L)))
  )
  srt <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 2L))
  embedding <- cbind(
    UMAP_1 = c(0, 0.2, 2, 2.2),
    UMAP_2 = c(0, 0.2, 2, 2.2)
  )
  rownames(embedding) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  srt
}

make_heatmap_test_srt <- function() {
  counts <- matrix(
    rep(seq_len(8L), each = 4L),
    nrow = 4L,
    dimnames = list(paste0("gene", seq_len(4L)), paste0("cell", seq_len(8L)))
  )
  srt <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$group <- factor(rep(c("A", "B"), each = 4L))
  srt$batch <- factor(rep(c("X", "Y"), 4L))
  srt
}

capture_heatmap_legends <- function(code) {
  calls <- list()
  original <- ComplexHeatmap::Legend
  testthat::local_mocked_bindings(
    Legend = function(...) {
      args <- list(...)
      calls[[length(calls) + 1L]] <<- args
      do.call(original, args)
    },
    .package = "ComplexHeatmap"
  )
  force(code)
  stats::setNames(
    calls,
    vapply(calls, function(x) as.character(x$title), character(1L))
  )
}

test_that("discrete heatmap legends retain custom border color and width", {
  gp <- heatmap_discrete_legend_gp("red", TRUE, "navy", 2.5)
  expect_identical(gp$fill, "red")
  expect_identical(gp$col, "navy")
  expect_identical(gp$lwd, 2.5)
})

test_that("ProjectionPlot finds query points after annotation layers", {
  query <- make_projection_test_srt("query")
  reference <- make_projection_test_srt("reference")
  assign(
    ".scop_projection_test_theme",
    projection_theme_with_annotation,
    envir = .GlobalEnv
  )
  on.exit(rm(".scop_projection_test_theme", envir = .GlobalEnv), add = TRUE)

  expect_no_error(ProjectionPlot(
    query,
    reference,
    query_group = "group",
    ref_group = "group",
    query_reduction = "umap",
    ref_reduction = "umap",
    query_param = list(
      palette = "Set1",
      cells.highlight = TRUE,
      raster = FALSE,
      theme_use = ".scop_projection_test_theme"
    ),
    ref_param = list(
      palette = "Set1",
      raster = FALSE,
      theme_use = ".scop_projection_test_theme"
    ),
    verbose = FALSE
  ))
})

test_that("GroupHeatmap border settings also control discrete legends", {
  srt <- make_heatmap_test_srt()
  legends <- capture_heatmap_legends(GroupHeatmap(
    srt,
    features = rownames(srt),
    group.by = "group",
    cell_annotation = "batch",
    feature_split = factor(c("M1", "M1", "M2", "M2")),
    exp_method = "raw",
    lib_normalize = FALSE,
    border = FALSE,
    heatmap_border = FALSE,
    cell_annotation_border = FALSE,
    feature_annotation_border = FALSE,
    cell_annotation_border_size = 2.5,
    feature_annotation_border_size = 3.5,
    use_raster = FALSE,
    verbose = FALSE
  ))

  for (title in c("group", "batch")) {
    expect_identical(legends[[title]]$border, FALSE)
    expect_true(is.na(legends[[title]]$legend_gp$col))
    expect_identical(legends[[title]]$legend_gp$lwd, 2.5)
  }
  expect_identical(legends[["Cluster"]]$border, FALSE)
  expect_true(is.na(legends[["Cluster"]]$legend_gp$col))
  expect_identical(legends[["Cluster"]]$legend_gp$lwd, 3.5)
})

test_that("CellCorHeatmap border settings also control group legends", {
  query <- make_projection_test_srt("query")
  reference <- make_projection_test_srt("reference")
  query$group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  reference$group <- factor(c("C", "C", "D", "D"), levels = c("C", "D"))
  other_params <- list(
    query_group = "group",
    query_reduction = "umap",
    query_assay = "RNA",
    query_dims = 1:2,
    query_collapsing = TRUE,
    ref_group = "group",
    ref_reduction = "umap",
    ref_assay = "RNA",
    ref_dims = 1:2,
    ref_collapsing = TRUE
  )
  query@tools[["KNNPredict_classification"]] <- list(
    distance_matrix = matrix(
      c(0.1, 0.8, 0.7, 0.2),
      nrow = 2L,
      dimnames = list(c("C", "D"), c("A", "B"))
    ),
    distance_metric = "cosine",
    other_params = other_params
  )

  legends <- capture_heatmap_legends(CellCorHeatmap(
    query,
    reference,
    query_group = "group",
    ref_group = "group",
    query_reduction = "umap",
    ref_reduction = "umap",
    query_dims = 1:2,
    ref_dims = 1:2,
    border = FALSE,
    heatmap_border = FALSE,
    cell_annotation_border = FALSE,
    feature_annotation_border = FALSE,
    cell_annotation_border_size = 2.5,
    use_raster = FALSE,
    verbose = FALSE
  ))

  for (title in c("Query:group", "Ref:group")) {
    expect_identical(legends[[title]]$border, FALSE)
    expect_true(is.na(legends[[title]]$legend_gp$col))
    expect_identical(legends[[title]]$legend_gp$lwd, 2.5)
  }
})
