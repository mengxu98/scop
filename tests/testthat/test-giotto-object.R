make_scop_giotto_seurat <- function(with_sct = FALSE) {
  counts <- matrix(
    c(
      5, 0, 1, 0, 4, 1,
      0, 3, 0, 4, 1, 0,
      2, 1, 5, 0, 0, 3,
      0, 2, 0, 5, 1, 4,
      3, 0, 2, 1, 0, 4,
      0, 4, 1, 3, 2, 0
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:6), paste0("Spot", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$x <- c(1, 2, 3, 1, 2, 3)
  srt$y <- c(1, 1, 1, 2, 2, 2)
  srt$celltype <- rep(c("A", "B"), each = 3)
  if (isTRUE(with_sct)) {
    srt[["SCT"]] <- Seurat::CreateAssayObject(counts = counts)
    srt <- Seurat::NormalizeData(srt, assay = "SCT", verbose = FALSE)
  }
  srt
}

make_mock_scop_giotto <- function() {
  cells <- paste0("Spot", 1:4)
  scop_giotto(
    giotto = list(mock = TRUE),
    source = list(
      cells = cells,
      features = paste0("Gene", 1:4),
      coordinates = data.frame(
        cell_ID = cells,
        sdimx = c(1, 2, 1, 2),
        sdimy = c(1, 1, 2, 2),
        stringsAsFactors = FALSE
      )
    ),
    results = list(
      cluster = list(
        table = data.frame(
          cell = cells,
          cluster = c("A", "A", "B", "B"),
          row.names = cells,
          stringsAsFactors = FALSE
        )
      ),
      spatial_genes = list(
        table = data.frame(
          feats = paste0("Gene", 1:4),
          score = c(4, 3, 2, 1),
          stringsAsFactors = FALSE
        )
      ),
      cell_proximity = list(
        table = data.frame(
          group_1 = c("A", "A", "B"),
          group_2 = c("A", "B", "B"),
          enrichm = c(1.5, -0.4, 2),
          type_int = c("homo", "hetero", "homo"),
          stringsAsFactors = FALSE
        )
      )
    ),
    active = "cluster",
    parameters = list(spat_unit = "cell", feat_type = "rna")
  )
}

test_that("scop_giotto constructor and print keep object state", {
  g <- make_mock_scop_giotto()

  expect_s3_class(g, "scop_giotto")
  expect_equal(g$active, "cluster")
  expect_output(print(g), "scop_giotto")
})

test_that("GiottoPlot supports scop_giotto without original Seurat", {
  g <- make_mock_scop_giotto()

  expect_s3_class(GiottoPlot(g, plot_type = "spatial"), "ggplot")
  expect_s3_class(GiottoPlot(g, plot_type = "cluster"), "ggplot")
  expect_s3_class(GiottoPlot(g, plot_type = "spatial_genes", top_n = 2), "ggplot")
  expect_s3_class(GiottoPlot(g, plot_type = "cell_proximity"), "ggplot")
  expect_s3_class(plot(g, plot_type = "cluster"), "ggplot")
})

test_that("SeuratToScopGiotto returns standalone object without modifying Seurat", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_scop_giotto_seurat(with_sct = TRUE)
  before_meta <- srt@meta.data
  before_tools <- srt@tools
  before_variable_features <- SeuratObject::VariableFeatures(srt)

  g <- SeuratToScopGiotto(
    srt,
    assay = "RNA",
    layer = "counts",
    use_sct = "auto",
    use_official = FALSE,
    verbose = FALSE
  )

  expect_s3_class(g, "scop_giotto")
  expect_false(inherits(g, "Seurat"))
  expect_equal(g$source$assay, "RNA")
  expect_equal(g$source$sct_mode, "auto")
  expect_equal(g$source$cells, colnames(srt))
  expect_equal(srt@meta.data, before_meta)
  expect_equal(srt@tools, before_tools)
  expect_equal(SeuratObject::VariableFeatures(srt), before_variable_features)
})

test_that("SeuratToScopGiotto can use SCT counts as fallback", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_scop_giotto_seurat(with_sct = TRUE)
  SeuratObject::DefaultAssay(srt) <- "SCT"
  srt[["RNA"]] <- NULL

  g <- SeuratToScopGiotto(
    srt,
    layer = "counts",
    use_official = FALSE,
    verbose = FALSE
  )

  expect_s3_class(g, "scop_giotto")
  expect_equal(g$source$assay, "SCT")
  expect_true(any(vapply(g$history, function(x) identical(x$step, "SeuratToScopGiotto"), logical(1))))
})

test_that("AddGiottoToSeurat is the explicit bridge back to Seurat", {
  srt <- make_scop_giotto_seurat()
  g <- make_mock_scop_giotto()

  out <- AddGiottoToSeurat(srt, g, result = "cluster", name = "Giotto_bridge", store_result = FALSE)

  expect_true("Giotto_bridge" %in% colnames(out@meta.data))
  expect_false("Giotto_bridge" %in% colnames(srt@meta.data))
  expect_false("Giotto" %in% names(out@tools))
})

test_that("scop_giotto single-step preprocessing and reduction run", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_scop_giotto_seurat()
  g <- SeuratToScopGiotto(
    srt,
    assay = "RNA",
    layer = "counts",
    use_official = FALSE,
    verbose = FALSE
  )

  g <- GiottoPreprocess(g, verbose = FALSE)
  g <- GiottoReduce(g, reduction = "pca", dims = 1:2, verbose = FALSE)

  expect_s3_class(g, "scop_giotto")
  expect_equal(g$parameters$pca_name, "pca")
})

test_that("RunGiottoWorkflow basic returns a scop_giotto object", {
  testthat::skip_if_not_installed("Giotto")
  srt <- make_scop_giotto_seurat()

  g <- RunGiottoWorkflow(
    srt,
    steps = "basic",
    assay = "RNA",
    layer = "counts",
    use_official = FALSE,
    verbose = FALSE,
    seed = 11
  )

  expect_s3_class(g, "scop_giotto")
  expect_true("cluster" %in% names(g$results))
  expect_true("spatial_network" %in% names(g$results))
  expect_gt(nrow(g$results$spatial_network$table), 0)
  expect_s3_class(GiottoPlot(g, plot_type = "cluster"), "ggplot")
  expect_s3_class(GiottoPlot(g, plot_type = "network"), "ggplot")
  expect_s3_class(GiottoPlot(g, plot_type = "dim"), "ggplot")
})
