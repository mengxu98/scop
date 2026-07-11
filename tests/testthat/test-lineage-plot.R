make_lineage_plot_srt <- function(ncells = 3000L) {
  cells <- paste0("cell", seq_len(ncells))
  counts <- Matrix::sparseMatrix(
    i = rep.int(1L, ncells),
    j = seq_len(ncells),
    x = 1,
    dims = c(1L, ncells),
    dimnames = list("gene1", cells)
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  embedding <- cbind(
    UMAP_1 = seq(-1, 1, length.out = ncells),
    UMAP_2 = sin(seq(-pi, pi, length.out = ncells))
  )
  rownames(embedding) <- cells
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  srt$Lineage1 <- seq(0, 1, length.out = ncells)
  srt
}

test_that("LineagePlot adapts large LOESS fit resolution while preserving trajectory extent", {
  srt <- make_lineage_plot_srt()

  layers <- LineagePlot(
    srt,
    lineages = "Lineage1",
    reduction = "umap",
    return_layer = TRUE
  )

  path_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPath"),
    layers$curve_layer
  )
  fitted <- path_layers[[1]]$data

  expect_lt(nrow(fitted), 2940L)
  expect_gt(nrow(fitted), 500L)
  expect_equal(range(fitted$index), c(31L, 2970L))
})

test_that("LineagePlot keeps small trajectories intact", {
  srt <- make_lineage_plot_srt(ncells = 100L)

  layers <- LineagePlot(
    srt,
    lineages = "Lineage1",
    reduction = "umap",
    return_layer = TRUE
  )
  path_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPath"),
    layers$curve_layer
  )

  expect_equal(nrow(path_layers[[1]]$data), 98L)
  expect_equal(range(path_layers[[1]]$data$index), c(2L, 99L))
})

test_that("LineagePlot handles NA-heavy lineages with whiskers", {
  srt <- make_lineage_plot_srt()
  srt$Lineage2 <- c(seq(0, 1, length.out = 1000L), rep(NA_real_, 2000L))

  layers <- LineagePlot(
    srt,
    lineages = c("Lineage1", "Lineage2"),
    reduction = "umap",
    whiskers = TRUE,
    return_layer = TRUE
  )
  whisker_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomSegment"),
    layers$curve_layer
  )

  expect_length(whisker_layers, 2L)
  expect_true(all(vapply(
    whisker_layers,
    function(layer) all(stats::complete.cases(layer$data)),
    logical(1)
  )))
})

test_that("LineagePlot increases fitting resolution with trajectory size", {
  small <- lineage_plot_fit_index(
    pseudotime = seq(0, 1, length.out = 3000L),
    trim = c(0.01, 0.99),
    span = 0.75
  )
  large <- lineage_plot_fit_index(
    pseudotime = seq(0, 1, length.out = 50000L),
    trim = c(0.01, 0.99),
    span = 0.75
  )

  expect_gt(length(large), length(small))
  expect_lt(length(large), 49000L)
})
