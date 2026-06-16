make_giotto_plot_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0, 4, 1,
      0, 3, 0, 4, 1, 0,
      2, 1, 5, 0, 0, 3,
      0, 2, 0, 5, 1, 4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$x <- c(1, 2, 3, 1, 2, 3)
  srt$y <- c(1, 1, 1, 2, 2, 2)
  srt
}

test_that("GiottoPlot cluster uses scop spatial plot without mutating Seurat", {
  srt <- make_giotto_plot_seurat()
  before_meta <- srt@meta.data
  result <- giotto_result(
    result_type = "cluster",
    giotto = list(mock = TRUE),
    clusters = data.frame(
      cluster = c("A", "A", "B", "B", "C", "C"),
      row.names = colnames(srt),
      stringsAsFactors = FALSE
    ),
    parameters = list(
      cluster_colname = "Giotto_cluster",
      coord.cols = c("x", "y"),
      image = NULL,
      k = 2,
      resolution = 0.5
    )
  )

  p <- GiottoPlot(result, srt = srt, overlay_image = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(srt@meta.data, before_meta)
  expect_false("Giotto_cluster" %in% colnames(srt@meta.data))
  expect_s3_class(plot(result, srt = srt, overlay_image = FALSE), "ggplot")
})

test_that("GiottoPlot proximity returns a scop-style heatmap", {
  result <- giotto_result(
    result_type = "cell_proximity",
    giotto = list(mock = TRUE),
    enrichment = data.frame(
      unified_int = c("A--A", "A--B", "B--B"),
      type_int = c("homo", "hetero", "homo"),
      enrichm = c(2, -1, 1.5),
      group_1 = c("A", "A", "B"),
      group_2 = c("A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    parameters = list(network_method = "Delaunay", number_of_simulations = 10)
  )

  p <- GiottoPlot(result)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot(result), "ggplot")
})

test_that("GiottoPlot spatial genes supports ranking and feature plots", {
  srt <- make_giotto_plot_seurat()
  result <- giotto_result(
    result_type = "spatial_genes",
    giotto = list(mock = TRUE),
    results = data.frame(
      feats = paste0("Gene", 1:4),
      score = c(9, 5, 3, 1),
      stringsAsFactors = FALSE
    ),
    top_features = paste0("Gene", 1:2),
    parameters = list(assay = "RNA", layer = "data", coord.cols = c("x", "y"))
  )

  expect_s3_class(GiottoPlot(result, top_n = 3), "ggplot")
  expect_s3_class(
    GiottoPlot(result, srt = srt, plot_type = "feature", overlay_image = FALSE),
    "ggplot"
  )
})

test_that("GiottoPlot spatial modules uses extracted spatial correlations", {
  result <- giotto_result(
    result_type = "spatial_modules",
    giotto = list(mock = TRUE),
    module_tables = list(
      result.cor_DT = expand.grid(
        feat_ID = paste0("Gene", 1:3),
        variable = paste0("Gene", 1:3),
        stringsAsFactors = FALSE
      )
    ),
    features = paste0("Gene", 1:3),
    parameters = list()
  )
  result$module_tables$result.cor_DT$spat_cor <- c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1)

  expect_s3_class(GiottoPlot(result, top_n = 3), "ggplot")
  expect_s3_class(plot(result, top_n = 3), "ggplot")
})
