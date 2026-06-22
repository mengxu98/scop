make_gsva_stat_srt <- function() {
  cells <- paste0("Cell", seq_len(8))
  expr <- matrix(seq_len(40), nrow = 5, ncol = 8)
  rownames(expr) <- paste0("Gene", seq_len(5))
  colnames(expr) <- cells
  srt <- suppressWarnings(Seurat::CreateSeuratObject(expr))
  srt$condition <- rep(c("A", "B"), each = 4)
  srt$sample <- rep(c("S1", "S2", "S3", "S4"), each = 2)
  scores <- matrix(
    c(
      2.0, 2.1, 1.8, 2.2, 0.2, 0.1, 0.3, 0.2,
      0.1, 0.2, 0.1, 0.2, 1.5, 1.6, 1.4, 1.7,
      0.4, 0.3, 0.5, 0.4, 0.4, 0.3, 0.5, 0.4
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("PathwayA", "PathwayB", "PathwayC"), cells)
  )
  suppressWarnings({
    srt[["GSVA"]] <- Seurat::CreateAssayObject(data = scores)
  })
  srt
}

test_that("GSVAPlot score mode ranks without placeholder p-values", {
  srt <- make_gsva_stat_srt()

  out <- GSVAPlot(
    srt = srt,
    mode = "score",
    group.by = "condition",
    topTerm = 2,
    return_data = TRUE
  )

  expect_equal(nrow(out), 2)
  expect_false("pvalue" %in% colnames(out))
  expect_false("p.adjust" %in% colnames(out))
  expect_equal(out$Description[[1]], "PathwayA")
})

test_that("GSVAPlot diff mode aggregates by sample and tests pathways", {
  srt <- make_gsva_stat_srt()

  out <- suppressWarnings(GSVAPlot(
    srt = srt,
    mode = "diff",
    group.by = "condition",
    sample.by = "sample",
    group_use = c("A", "B"),
    test.use = "t.test",
    return_data = TRUE
  ))

  expect_named(
    out,
    c(
      "ID", "Description", "Database", "group1", "group2",
      "mean_group1", "mean_group2", "diff", "pvalue", "p.adjust",
      "test.use", "n_group1", "n_group2"
    )
  )
  expect_equal(unique(out$group1), "A")
  expect_equal(unique(out$group2), "B")
  expect_true(all(out$n_group1 == 2))
  expect_true(all(out$n_group2 == 2))
  expect_equal(out$diff[out$Description == "PathwayA"], 1.825, tolerance = 1e-8)
})

test_that("GSVAPlot diff mode accepts sample-level score matrices", {
  srt <- make_gsva_stat_srt()
  sample_scores <- matrix(
    c(
      2.05, 2.00, 0.15, 0.25,
      0.15, 0.15, 1.55, 1.55
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("PathwayA", "PathwayB"), c("S1", "S2", "S3", "S4"))
  )

  out <- suppressWarnings(GSVAPlot(
    srt = srt,
    res = sample_scores,
    mode = "diff",
    group.by = "condition",
    sample.by = "sample",
    group_use = c("A", "B"),
    test.use = "t.test",
    return_data = TRUE
  ))

  expect_equal(out$n_group1[[1]], 2)
  expect_equal(out$n_group2[[1]], 2)
  expect_equal(out$diff[out$Description == "PathwayA"], 1.825, tolerance = 1e-8)
})

test_that("GSVAPlot diff mode validates grouping inputs", {
  srt <- make_gsva_stat_srt()

  expect_error(
    GSVAPlot(srt = srt, mode = "diff", group.by = "condition"),
    "sample.by"
  )
  expect_error(
    GSVAPlot(
      srt = srt,
      mode = "diff",
      group.by = "missing",
      sample.by = "sample"
    ),
    "group.by"
  )
  expect_error(
    GSVAPlot(
      srt = srt,
      mode = "diff",
      group.by = "condition",
      sample.by = "sample",
      group_use = "A"
    ),
    "group_use"
  )
})

test_that("GSVAPlot diff mode rejects group-level RunGSVA results", {
  srt <- make_gsva_stat_srt()
  group_scores <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("PathwayA", "PathwayB"), c("A", "B"))
  )
  res <- list(scores = group_scores, group.by = "condition")

  expect_error(
    GSVAPlot(
      srt = srt,
      res = res,
      mode = "diff",
      group.by = "condition",
      sample.by = "sample",
      group_use = c("A", "B")
    ),
    "cannot produce real p-values"
  )
})
