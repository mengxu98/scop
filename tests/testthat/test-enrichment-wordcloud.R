wordcloud_result_fixture <- function() {
  tab <- data.frame(
    ID = paste0("t", 1:3), Description = c("Shared alpha", "Shared beta", "Shared gamma"),
    Database = "custom", Groups = "G", pvalue = c(0.001, 0.01, 0.02),
    p.adjust = c(0.001, 0.01, 0.02), Count = 2L,
    geneID = c("ora/first", "ora/second", "ora/third"),
    core_enrichment = c("a/b", "b/c", "c"), NES = c(1, -2, 3), Score = c(-2, 5, -1),
    row.names = paste0("t", 1:3)
  )
  list(enrichment = tab, results = list(`G-custom` = methods::new("gseaResult", result = tab)))
}

test_that("wordclouds retain method-specific genes, weights and top-word ordering", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  res <- wordcloud_result_fixture()
  original <- res
  ora <- EnrichmentPlot(
    res = res, db = "custom", plot_type = "wordcloud",
    word_type = "feature", topWord = 2L
  )
  gsea <- GSEAPlot(
    res = res, db = "custom", plot_type = "wordcloud",
    word_type = "feature", topWord = 2L, direction = "both"
  )
  expect_identical(ora$data$keyword, c("ora", "first"))
  expect_equal(ora$data$score, c(sum(-log10(res$enrichment$p.adjust)), 3))
  expect_identical(ora$data$count, c(3L, 1L))
  expect_identical(gsea$data$keyword, c("b", "c"))
  expect_equal(gsea$data$score, c(5, -log10(0.01) - log10(0.02)))
  expect_identical(gsea$data$count, c(2L, 2L))
  score <- EnrichmentPlot(
    res = res, db = "custom", plot_type = "wordcloud",
    score_col = "Score", score_label = "Effect", topWord = 2L,
    words_excluded = c("shared", "alpha")
  )
  expect_identical(score$data$keyword, c("beta", "gamma"))
  expect_equal(score$data$score, c(5, 1))
  expect_identical(score$scales$get_scales("colour")$name, "Effect:")
  expect_identical(res, original)
})

test_that("only GO term wordclouds request the GO keyword dependency", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  original_check <- check_r
  testthat::local_mocked_bindings(check_r = function(pkg, ...) {
    if ("simplifyEnrichment" %in% pkg) stop("GO dependency sentinel")
    original_check(pkg, ...)
  }, .package = "scop")
  res <- wordcloud_result_fixture()
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    for (word in c("term", "feature")) {
      expect_s3_class(plotter(res = res, db = "custom", plot_type = "wordcloud", word_type = word), "ggplot")
    }
  }
  res$enrichment$Database <- "GO_BP"
  res$results <- list(`G-GO_BP` = methods::new("gseaResult", result = res$enrichment))
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    expect_s3_class(plotter(res = res, db = "GO_BP", plot_type = "wordcloud", word_type = "feature"), "ggplot")
    expect_error(plotter(res = res, db = "GO_BP", plot_type = "wordcloud", word_type = "term"), "GO dependency sentinel")
  }
  expect_s3_class(EnrichmentPlot(
    res = res, db = "GO_BP", plot_type = "wordcloud",
    word_type = "term", score_col = "Score"
  ), "ggplot")
})

test_that("term and feature wordclouds retain their distinct top-word scopes", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  res <- wordcloud_result_fixture()
  first <- res$enrichment
  second <- first
  second$Groups <- "H"
  second$p.adjust <- second$p.adjust / 2
  res$enrichment <- rbind(first, second)
  terms <- EnrichmentPlot(
    res = res, db = "custom", split_by = "Database", plot_type = "wordcloud",
    word_type = "term", topWord = 2L
  )
  genes <- EnrichmentPlot(
    res = res, db = "custom", split_by = "Database", plot_type = "wordcloud",
    word_type = "feature", topWord = 2L
  )
  expect_equal(nrow(terms$data), 4L)
  expect_equal(as.integer(table(terms$data$Groups)), c(2L, 2L))
  expect_equal(nrow(genes$data), 2L)
  expect_identical(genes$data$keyword, c("ora", "ora"))
})

test_that("an empty GO keyword group retains the following groups", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  skip_if_not_installed("simplifyEnrichment")
  res <- wordcloud_result_fixture()
  first <- res$enrichment
  first$Database <- "GO_BP"
  second <- first
  second$ID <- paste0("h", 1:3)
  second$Groups <- "H"
  res$enrichment <- rbind(first, second)
  local_mocked_bindings(enrichment_go_keywords = function(ids) {
    if (all(startsWith(ids, "t"))) {
      return(data.frame(keyword = character(), padj = numeric(), n_term = integer()))
    }
    data.frame(keyword = "retained", padj = 0.01, n_term = 2L)
  }, .package = "scop")
  plot <- EnrichmentPlot(
    res = res, db = "GO_BP", split_by = "Database", plot_type = "wordcloud"
  )
  expect_identical(plot$data$keyword, "retained")
  expect_equal(as.character(plot$data$Groups), "H")
  expect_equal(plot$data$score, 2)
})
