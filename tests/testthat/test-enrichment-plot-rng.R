enrichment_rng_fixture <- function() {
  tab <- data.frame(
    ID = c("t1", "t2"), Description = c("Alpha response", "Beta process"),
    Database = "custom", Groups = "G", pvalue = c(0.001, 0.002),
    p.adjust = c(0.01, 0.02), Count = c(3L, 2L),
    geneID = c("a/b/c", "b/c"), core_enrichment = c("a/b/c", "b/c"),
    NES = c(1.3, -1.5), row.names = c("t1", "t2")
  )
  list(enrichment = tab, results = list(`G-custom` = methods::new("gseaResult", result = tab)))
}

test_that("enrichment plotting restores caller RNG after success and errors", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  res <- enrichment_rng_fixture()
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    for (has_seed in c(TRUE, FALSE)) {
      for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        expect_rng_restored(
          plotter(res = res, db = "custom", plot_type = "wordcloud", word_type = "feature", seed = 37),
          has_seed, kind
        )
        expect_rng_restored(
          expect_error(plotter(res = res, db = "custom", plot_type = "invalid", seed = 37), "arg"),
          has_seed, kind
        )
      }
    }
  }
})

test_that("enrichment plot seeds reproduce sampled angles and network coordinates", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  withr::local_preserve_seed()
  res <- enrichment_rng_fixture()
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    for (kind in c("wordcloud", "network")) {
      args <- list(
        res = res, db = "custom", plot_type = kind, word_type = "feature",
        seed = 37, network_layoutadjust = FALSE, combine = FALSE
      )
      set.seed(210)
      first <- do.call(plotter, args)[[1L]]
      set.seed(890)
      second <- do.call(plotter, args)[[1L]]
      expect_identical(first$data, second$data)
      expect_equal(ggplot2::ggplot_build(first)$data, ggplot2::ggplot_build(second)$data)
    }
  }
})
