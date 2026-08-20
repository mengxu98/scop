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

test_that("IntegrationBenchmarkPlot draws overview heatmap and scatter", {
  metrics <- data.frame(
    method = rep(c("Uncorrected", "Harmony"), each = 2),
    metric = rep(c("iLISI", "cLISI"), 2),
    category = rep(c("batch", "bio"), 2),
    value = c(1.2, 1.1, 2.4, 1.2),
    scaled = c(0.1, 0.9, 0.7, 0.85),
    direction = "higher",
    stringsAsFactors = FALSE
  )
  fake <- structure(
    list(
      summary = data.frame(
        method = c("Uncorrected", "Harmony"),
        bio = c(0.80, 0.82),
        batch = c(0.20, 0.70),
        overall = c(0.56, 0.77),
        status = "success",
        stringsAsFactors = FALSE
      ),
      metrics = metrics,
      runs = data.frame(
        method = c("Uncorrected", "Harmony"),
        status = "success",
        umap = NA_character_,
        stringsAsFactors = FALSE
      ),
      srt = NULL,
      batch = "tech",
      celltype = "celltype"
    ),
    class = c("integration_benchmark_result", "list")
  )
  expect_s3_class(IntegrationBenchmarkPlot(fake, plot_type = "overview"), "ggplot")
  expect_s3_class(IntegrationBenchmarkPlot(fake, plot_type = "heatmap"), "ggplot")
  expect_s3_class(IntegrationBenchmarkPlot(fake, plot_type = "scatter"), "ggplot")
  scores <- integration_overall_scores(metrics)
  expect_equal(nrow(scores), 2)
  expect_true(all(c("bio", "batch", "overall") %in% colnames(scores)))
})

test_that("print.integration_benchmark_result returns the object", {
  fake <- structure(
    list(
      summary = data.frame(
        method = "Uncorrected",
        bio = 0.5,
        batch = 0.4,
        overall = 0.46,
        status = "success",
        stringsAsFactors = FALSE
      ),
      metrics = data.frame(
        method = "Uncorrected",
        metric = "iLISI",
        category = "batch",
        value = 1.4,
        scaled = 0.4,
        direction = "higher",
        stringsAsFactors = FALSE
      ),
      runs = data.frame(method = "Uncorrected", status = "success", stringsAsFactors = FALSE),
      srt = NULL
    ),
    class = c("integration_benchmark_result", "list")
  )
  expect_identical(print(fake), fake)
})
