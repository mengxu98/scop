test_that("GO keyword counts and Fisher results match the native backend", {
  skip_if_not_installed("simplifyEnrichment")
  ids <- c("GO:0006351", "GO:0006355", "GO:0006366", "GO:0006396", "GO:0006412")
  set.seed(11)
  rng <- .Random.seed
  for (query in list(ids, c(ids, ids[[1]], "GO:unknown"), character())) {
    expect_identical(
      enrichment_go_keywords(query),
      simplifyEnrichment::keyword_enrichment_from_GO(query)
    )
  }
  expect_identical(.Random.seed, rng)
})

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

test_that("P value weights preserve positive values and reject invalid probabilities", {
  values <- c(0, .Machine$double.xmin, .Machine$double.xmin / 2, 0.01, 1, NA_real_)
  expected <- c(-log10(.Machine$double.xmin), -log10(values[2:5]), NA_real_)
  expect_equal(enrichment_metric_weight(values), expected)
  for (value in c(-0.1, 1.1, Inf, -Inf)) {
    expect_error(enrichment_metric_weight(value), "between 0 and 1")
  }
  expect_equal(enrichment_metric_weight(c(-2, 1), score_mode = TRUE), c(2, 1))
})

test_that("GSEA maps retain zero P values as finite weights and preserve NES signs", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("shadowtext")
  tab <- data.frame(
    ID = paste0("t", 1:3), Description = c("Alpha response", "Beta response", "Gamma process"),
    Database = "custom", Groups = "G", pvalue = c(0, 0.001, 0.01),
    p.adjust = c(0, 0.01, 0.03), Count = 3L,
    core_enrichment = c("a/b/c", "b/c/d", "c/d/e"), NES = c(1.3, -1.5, 2),
    row.names = paste0("t", 1:3)
  )
  result <- list(enrichment = tab, results = list(`G-custom` = methods::new("gseaResult", result = tab)))
  plot <- suppressWarnings(GSEAPlot(
    res = result, db = "custom", plot_type = "enrichmap", enrichmap_layout = "circle", direction = "both",
    enrichmap_show_keyword = TRUE, combine = FALSE, verbose = FALSE
  ))[[1]]
  nodes <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)[[1L]]$data
  nodes <- nodes[match(tab$ID, nodes$name), ]
  expect_true(all(is.finite(nodes$metric)))
  expect_equal(nodes$metric, c(-log10(.Machine$double.xmin), 2, -log10(0.03)))
  expect_identical(as.character(nodes$Direction), c("Pos", "Neg", "Pos"))
  expect_equal(nodes$NES, tab$NES)
  expect_identical(nodes$p.adjust, tab$p.adjust)
  expect_identical(result$enrichment, tab)
  expect_true(all(vapply(ggplot2::ggplot_build(plot)$data, function(layer) {
    !"size" %in% names(layer) || all(is.finite(layer$size))
  }, logical(1))))
})

test_that("enrichment network views execute the selected layouts and clustering", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("shadowtext")
  tab <- data.frame(
    ID = paste0("t", 1:3), Description = c("Alpha response", "Beta response", "Gamma process"),
    Database = "custom", Groups = "G", pvalue = c(0.001, 0.002, 0.003), p.adjust = c(0.01, 0.02, 0.03),
    Count = 3L, GeneRatio = "3/5", BgRatio = "5/10", geneID = c("a/b/c", "b/c/d", "c/d/e"),
    core_enrichment = c("a/b/c", "b/c/d", "c/d/e"), NES = c(1.3, -1.5, 2), row.names = paste0("t", 1:3)
  )
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
      plots <- suppressWarnings(plotter(
        res = res, db = "custom", plot_type = kind,
        network_layoutadjust = FALSE, combine = FALSE, verbose = FALSE
      ))
      expect_length(plots, 1L)
      expect_s3_class(plots[[1]], "ggplot")
      expect_gt(length(suppressWarnings(ggplot2::ggplot_build(plots[[1]]))$data), 0L)
    }
  }
  expect_equal(sum(calls == "layout_with_fr"), 4L)
  expect_equal(sum(calls == "cluster_fast_greedy"), 2L)
})

test_that("enrichment maps retain weighted labels, counts and GSEA direction", {
  skip_if_not_installed("DOSE")
  tab <- data.frame(
    ID = paste0("t", 1:3), Description = c("Hormone secretion", "Endocrine differentiation", "Vesicle transport"),
    Database = "custom", Groups = "G", pvalue = c(0, 0.001, 0.01), p.adjust = c(0, 0.01, 0.03),
    Count = c(3L, 2L, 2L), geneID = c("a/b/c", "b/c", "a/c"),
    core_enrichment = c("a/b/c", "b/c", "a/c"), NES = c(1.3, -1.5, 2),
    Score = c(-8, 3, 1), row.names = paste0("t", 1:3)
  )
  res <- list(enrichment = tab, results = list(`G-custom` = methods::new("gseaResult", result = tab)))
  original <- res
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    args <- list(res = res, db = "custom", plot_type = "enrichmap", enrichmap_layout = "circle",
      enrichmap_show_keyword = TRUE, enrichmap_nlabel = 2L, words_excluded = character(), combine = FALSE, verbose = FALSE)
    if (identical(plotter, GSEAPlot)) args$direction <- "both"
    for (label in c("term", "feature")) {
      args$enrichmap_label <- label
      p <- do.call(plotter, args)[[1L]]
      nodes <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), p$layers)[[1L]]$data
      nodes <- nodes[match(tab$ID, nodes$name), ]
      expect_equal(nodes$metric, c(-log10(.Machine$double.xmin), 2, -log10(0.03)))
      expect_identical(nodes$Count, tab$Count)
      expect_identical(nodes$keyword1[nodes$name == "t1"], "Hormone Secretion")
      for (cluster in unique(nodes$clusters)) {
        rows <- nodes[nodes$clusters == cluster, ]
        totals <- tapply(rep(rows$metric, lengths(rows$geneID)), unlist(rows$geneID), sum)
        expected <- paste(utils::head(names(sort(totals, decreasing = TRUE)), 2L), collapse = " ")
        expect_identical(unique(rows$keyword2), expected)
      }
      if (identical(plotter, GSEAPlot)) {
        expect_identical(as.character(nodes$Direction), c("Pos", "Neg", "Pos"))
      }
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    }
  }
  p <- EnrichmentPlot(res = res, db = "custom", plot_type = "enrichmap",
    score_col = "Score", enrichmap_layout = "circle", combine = FALSE, verbose = FALSE)[[1L]]
  nodes <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), p$layers)[[1L]]$data
  expect_equal(nodes$metric[match(tab$ID, nodes$name)], abs(tab$Score))
  expect_identical(res, original)
})
