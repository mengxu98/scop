test_that("LISI scaling matches scIB-style 0-1 bounds", {
  expect_equal(scale_ilisi(1, 4), 0)
  expect_equal(scale_ilisi(4, 4), 1)
  expect_equal(scale_ilisi(2.5, 4), 0.5)
  expect_equal(scale_clisi(1, 5), 1)
  expect_equal(scale_clisi(5, 5), 0)
  expect_equal(scale_asw_bio(1), 1)
  expect_equal(scale_asw_bio(-1), 0)
  expect_true(is.na(scale_ilisi(2, 1)))
})

test_that("format_integration_metrics adds iLISI and cLISI", {
  srt <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(
      matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:4))),
      sparse = TRUE
    )
  )
  srt$tech <- c("A", "A", "B", "B")
  srt$celltype <- c("alpha", "beta", "alpha", "beta")
  raw <- data.frame(
    metric = c("tech_LISI_mean", "celltype_LISI_mean", "batch_ASW_mixing", "celltype_ASW"),
    value = c(1.8, 1.1, 0.4, 0.5),
    stringsAsFactors = FALSE
  )
  out <- format_integration_metrics(
    raw,
    srt,
    batch_col = "tech",
    celltype_col = "celltype",
    method = "Harmony"
  )
  expect_true(all(c("iLISI", "cLISI") %in% out$metric))
  expect_identical(unique(out$method), "Harmony")
  expect_equal(out$category[out$metric == "iLISI"], "batch")
  expect_equal(out$category[out$metric == "cLISI"], "bio")
  expect_equal(out$scaled[out$metric == "iLISI"], scale_ilisi(1.8, 2))
  expect_equal(out$scaled[out$metric == "cLISI"], scale_clisi(1.1, 2))
})

fake_integration_benchmark_srt <- function() {
  srt <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(
      matrix(1:20, nrow = 2, dimnames = list(c("g1", "g2"), paste0("c", 1:10))),
      sparse = TRUE
    )
  )
  set.seed(1)
  srt$Uncorrected_tech_LISI <- stats::runif(10, 1.0, 1.4)
  srt$Harmony_tech_LISI <- stats::runif(10, 2.0, 2.8)
  srt$Uncorrected_celltype_LISI <- stats::runif(10, 1.0, 1.2)
  srt$Harmony_celltype_LISI <- stats::runif(10, 1.0, 1.3)
  srt@tools$IntegrationBenchmark <- list(
    summary = data.frame(
      method = c("Uncorrected", "Harmony"),
      bio = c(0.80, 0.82),
      batch = c(0.20, 0.70),
      overall = c(0.56, 0.77),
      status = "success",
      stringsAsFactors = FALSE
    ),
    metrics = data.frame(
      method = rep(c("Uncorrected", "Harmony"), each = 2),
      metric = rep(c("iLISI", "cLISI"), 2),
      category = rep(c("batch", "bio"), 2),
      value = c(1.2, 1.1, 2.4, 1.2),
      scaled = c(0.1, 0.9, 0.7, 0.85),
      direction = "higher",
      stringsAsFactors = FALSE
    ),
    runs = data.frame(
      method = c("Uncorrected", "Harmony"),
      status = "success",
      umap = NA_character_,
      stringsAsFactors = FALSE
    ),
    batch = "tech",
    celltype = "celltype"
  )
  srt
}

test_that("IntegrationBenchmarkPlot draws box heatmap and scatter from Seurat", {
  srt <- fake_integration_benchmark_srt()
  expect_s4_class(srt, "Seurat")
  expect_s3_class(IntegrationBenchmarkPlot(srt, plot_type = "box"), "ggplot")
  expect_s3_class(IntegrationBenchmarkPlot(srt, plot_type = "heatmap"), "ggplot")
  expect_s3_class(IntegrationBenchmarkPlot(srt, plot_type = "scatter"), "ggplot")
  plots <- IntegrationBenchmarkPlot(srt, plot_type = "auto")
  expect_true(all(c("box", "heatmap", "scatter") %in% names(plots)))
  expect_false("umap" %in% names(plots))
  scores <- integration_overall_scores(srt@tools$IntegrationBenchmark$metrics)
  expect_equal(nrow(scores), 2)
  expect_true(all(c("bio", "batch", "overall") %in% colnames(scores)))
})

test_that("IntegrationBenchmarkPlot box works from LISI columns alone", {
  srt <- fake_integration_benchmark_srt()
  srt@tools$IntegrationBenchmark <- NULL
  expect_s3_class(IntegrationBenchmarkPlot(srt, plot_type = "box"), "ggplot")
})

test_that("parse_lisi_feature splits method and label", {
  expect_equal(parse_lisi_feature("Harmony_tech_LISI")$method, "Harmony")
  expect_equal(parse_lisi_feature("Harmony_tech_LISI")$label, "tech")
  expect_equal(parse_lisi_feature("Uncorrected_celltype_LISI")$method, "Uncorrected")
  expect_equal(parse_lisi_feature("Uncorrected_celltype_LISI")$label, "celltype")
  expect_null(parse_lisi_feature("not_lisi"))
})
