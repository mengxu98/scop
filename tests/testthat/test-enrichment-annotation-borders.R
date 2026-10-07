enrichment_border_fixture <- function() {
  list(enrichment = data.frame(
    Database = "GO_BP", Groups = "C1", ID = "GO:0000001",
    Description = "test process", pvalue = 0.001, p.adjust = 0.01,
    geneID = "G1/G2", Count = 2, GeneRatio = "2/2", BgRatio = "2/20"
  ))
}

capture_enrichment_borders <- function(...) {
  calls <- list(textboxes = list(), blocks = list(), legends = list())
  original_textbox <- ComplexHeatmap::anno_textbox
  original_block <- ComplexHeatmap::anno_block
  original_legend <- ComplexHeatmap::Legend
  testthat::local_mocked_bindings(
    RunEnrichment = function(...) enrichment_border_fixture(), .package = "scop"
  )
  testthat::local_mocked_bindings(
    Legend = function(...) {
      args <- list(...)
      calls$legends[[length(calls$legends) + 1L]] <<- args
      do.call(original_legend, args)
    },
    anno_textbox = function(...) {
      args <- list(...)
      calls$textboxes[[length(calls$textboxes) + 1L]] <<- args
      do.call(original_textbox, args)
    },
    anno_block = function(...) {
      args <- list(...)
      calls$blocks[[length(calls$blocks) + 1L]] <<- args
      do.call(original_block, args)
    }, .package = "ComplexHeatmap"
  )
  result <- heatmap_enrichment(
    geneID = c("G1", "G2"), geneID_groups = c("C1", "C1"),
    anno_terms = TRUE, db = "GO_BP", padjustCutoff = 0.05, ...
  )
  calls$result <- result
  calls
}

test_that("text-only enrichment terms honor gray border and line width", {
  calls <- capture_enrichment_borders(
    terms_stat = "none", terms_border = TRUE,
    terms_border_palcolor = "grey60", terms_border_size = 2.5,
    feature_annotation_border = FALSE
  )
  expect_length(calls$textboxes, 1L)
  expect_identical(calls$textboxes[[1L]]$background_gp$col, "grey60")
  expect_identical(calls$textboxes[[1L]]$background_gp$lwd, 2.5)
  expect_true(is.na(calls$blocks[[1L]]$gp$col))
  expect_identical(calls$legends[[1L]]$border, "grey60")
  expect_identical(calls$legends[[1L]]$legend_gp$col, "grey60")
  expect_identical(calls$legends[[1L]]$legend_gp$lwd, 2.5)
  expect_identical(calls$result$res, enrichment_border_fixture())
})

test_that("text-only term frames and cluster strips have independent switches", {
  calls <- capture_enrichment_borders(
    terms_stat = "none", terms_border = FALSE,
    feature_annotation_border = TRUE,
    feature_annotation_border_color = "navy", feature_annotation_border_size = 3
  )
  expect_true(is.na(calls$textboxes[[1L]]$background_gp$col))
  expect_true(is.na(calls$legends[[1L]]$border))
  expect_true(is.na(calls$legends[[1L]]$legend_gp$col))
  expect_identical(calls$blocks[[1L]]$gp$col, "navy")
  expect_identical(calls$blocks[[1L]]$gp$lwd, 3)
})

test_that("graphic enrichment cluster strips also honor hidden borders", {
  calls <- capture_enrichment_borders(
    terms_stat = "p.adjust", feature_annotation_border = FALSE
  )
  expect_true(is.na(calls$blocks[[1L]]$gp$col))
})

test_that("enrichment border defaults retain visible black frames", {
  calls <- capture_enrichment_borders(terms_stat = "none")
  expect_identical(calls$textboxes[[1L]]$background_gp$col, "black")
  expect_identical(calls$blocks[[1L]]$gp$col, "black")
})
