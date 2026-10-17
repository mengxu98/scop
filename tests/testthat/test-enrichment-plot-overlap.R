test_that("EnrichmentPlot overlap edges match pairwise shared-gene counts", {
  actual <- enrichment_overlap_edges(
    c("term_a", "term_b", "term_c"),
    list(c("A", "B", "B"), c("B", "C"), "D")
  )
  expected <- data.frame(
    from = "term_a",
    to = "term_b",
    weight = 1,
    stringsAsFactors = FALSE
  )
  expect_equal(actual, expected)
})

test_that("EnrichmentPlot overlap edges omit disjoint and singleton term sets", {
  actual <- enrichment_overlap_edges(c("term_a", "term_b"), list("A", "B"))
  expect_equal(actual, data.frame(
    from = character(), to = character(), weight = numeric()
  ))
})

test_that("enrichment network views execute the selected layouts and clustering", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("shadowtext")
  tab <- data.frame(ID = paste0("t", 1:3), Description = c("Alpha response", "Beta response", "Gamma process"),
    Database = "custom", Groups = "G", pvalue = c(0.001, 0.002, 0.003), p.adjust = c(0.01, 0.02, 0.03),
    Count = 3L, GeneRatio = "3/5", BgRatio = "5/10", geneID = c("a/b/c", "b/c/d", "c/d/e"),
    core_enrichment = c("a/b/c", "b/c/d", "c/d/e"), NES = c(1.3, -1.5, 2), row.names = paste0("t", 1:3))
  res <- list(enrichment = tab, results = list(`G-custom` = methods::new("gseaResult", result = tab)))
  calls <- character()
  original_get <- get_namespace_fun
  testthat::local_mocked_bindings(get_namespace_fun = function(package, fun, ...) {
    backend <- original_get(package, fun, ...)
    if (identical(package, "igraph") && grepl("^(layout_with_|cluster_)", fun)) {
      return(function(graph, ...) {
        calls <<- c(calls, fun)
        backend(graph, ...)
      })
    }
    backend
  }, .package = "scop")
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    for (kind in c("network", "enrichmap")) {
      plots <- suppressWarnings(plotter(res = res, db = "custom", plot_type = kind,
        network_layoutadjust = FALSE, combine = FALSE, verbose = FALSE))
      expect_length(plots, 1L)
      expect_s3_class(plots[[1]], "ggplot")
      expect_gt(length(suppressWarnings(ggplot2::ggplot_build(plots[[1]]))$data), 0L)
    }
  }
  expect_equal(sum(calls == "layout_with_fr"), 4L)
  expect_equal(sum(calls == "cluster_fast_greedy"), 2L)
})
