test_that("CytoTRACE2 model loading tolerates duplicated CSV row indexes", {
  data_dir <- tempfile()
  dir.create(data_dir)
  saveRDS(list(), file.path(data_dir, "model_parameters.rds"))
  writeLines(
    c(
      ",0",
      "0,A1bg",
      "0,A1cf",
      "1,A4galt"
    ),
    file.path(data_dir, "features_model_training_17.csv")
  )

  expect_identical(
    scop:::load_cytotrace2_data(data_dir, verbose = FALSE)$features,
    c("A1bg", "A1cf", "A4galt")
  )
})

test_that("CytoTRACE2 model loading rejects duplicated feature names", {
  data_dir <- tempfile()
  dir.create(data_dir)
  saveRDS(list(), file.path(data_dir, "model_parameters.rds"))
  writeLines(
    c(
      ",0",
      "0,A1bg",
      "1,A1bg"
    ),
    file.path(data_dir, "features_model_training_17.csv")
  )

  expect_error(
    scop:::load_cytotrace2_data(data_dir, verbose = FALSE),
    "duplicated feature names"
  )
})

test_that("CytoTRACEPlot boxplot keeps angled x-axis labels when wrapped", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("patchwork")

  set.seed(1)
  counts <- Matrix::Matrix(
    stats::rpois(60, lambda = 5),
    nrow = 6,
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

  srt <- SeuratObject::CreateSeuratObject(counts)
  srt[["CytoTRACE2_Potency"]] <- rep(
    c("Multipotent", "Oligopotent"),
    length.out = ncol(srt)
  )
  srt[["CytoTRACE2_Score"]] <- seq(0.1, 0.9, length.out = ncol(srt))
  srt[["CytoTRACE2_Relative"]] <- seq(0.2, 0.8, length.out = ncol(srt))
  srt[["clusters"]] <- rep(
    c("Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"),
    length.out = ncol(srt)
  )

  embedding <- matrix(stats::rnorm(ncol(srt) * 2), ncol = 2)
  rownames(embedding) <- colnames(srt)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(srt)
  )

  plots <- CytoTRACEPlot(
    srt,
    group.by = "clusters",
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    palette = "ChineseSet8",
    theme_use = "theme_blank",
    combine = FALSE,
    verbose = FALSE
  )
  wrapped <- patchwork::wrap_plots(plots[c(4, 1, 2, 3, 5)], ncol = 5)
  x_labels <- ggplot2::ggplot_build(plots[["Boxplot"]])$
    layout$panel_params[[1]]$x$get_labels()

  expect_identical(plots[["Boxplot"]]$theme$axis.text.x$angle, 45)
  expect_false(inherits(plots[["Boxplot"]]$theme$axis.text, "element_blank"))
  expect_false(any(grepl("\n", x_labels, fixed = TRUE)))
  expect_s3_class(patchwork::patchworkGrob(wrapped), "gtable")
})
