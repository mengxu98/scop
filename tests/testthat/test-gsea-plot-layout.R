gsea_layout_fixture <- function() {
  tabs <- lapply(paste0("Group", 1:4), function(group) {
    data.frame(
      ID = c("t1", "t2"), Description = c("Alpha response", "Beta process"),
      Database = "custom", Groups = group, pvalue = c(0.001, 0.01),
      p.adjust = c(0.01, 0.02), NES = c(2, -1.5), Count = 2L,
      core_enrichment = c("a/b", "b/c"), row.names = c("t1", "t2")
    )
  })
  results <- lapply(tabs, function(tab) methods::new("gseaResult", result = tab))
  names(results) <- paste0("Group", 1:4, "-custom")
  list(enrichment = do.call(rbind, tabs), results = results)
}

test_that("GSEAPlot arranges combined plots in the requested direction", {
  skip_if_not_installed("DOSE")
  skip_if_not_installed("ggwordcloud")
  loadNamespace("DOSE")
  res <- gsea_layout_fixture()
  panels <- GSEAPlot(res = res, db = "custom", plot_type = "wordcloud", combine = FALSE)
  expect_length(panels, 4L)
  for (byrow in c(TRUE, FALSE)) {
    plot <- GSEAPlot(
      res = res, db = "custom", plot_type = "wordcloud",
      nrow = 2, ncol = 2, byrow = byrow
    )
    grob <- patchwork::patchworkGrob(plot)
    layout <- grob$layout[startsWith(grob$layout$name, "panel") & grob$layout$name != "panel-area", ]
    layout <- layout[match(1:4, as.integer(sub(".*-", "", layout$name))), ]
    expect_equal(match(layout$t, sort(unique(layout$t))), if (byrow) c(1L, 1L, 2L, 2L) else c(1L, 2L, 1L, 2L))
    expect_equal(match(layout$l, sort(unique(layout$l))), if (byrow) c(1L, 2L, 1L, 2L) else c(1L, 1L, 2L, 2L))
  }
  columns <- GSEAPlot(res = res, db = "custom", plot_type = "wordcloud", combine = FALSE, byrow = FALSE)
  expect_identical(lapply(columns, function(p) p$data), lapply(panels, function(p) p$data))
})

test_that("enrichment plots retain selected IDs and reject invalid group keys", {
  skip_if_not_installed("DOSE")
  loadNamespace("DOSE")
  res <- gsea_layout_fixture()
  original <- serialize(res, NULL)
  for (plotter in list(EnrichmentPlot, GSEAPlot)) {
    plots <- plotter(
      res = res, db = "custom", plot_type = "bar", topTerm = 1,
      id_use = c("t2", "t1"), combine = FALSE, verbose = FALSE
    )
    expect_setequal(unlist(lapply(plots, function(p) as.character(p$data$ID))), c("t1", "t2"))
    expect_no_warning(invisible(lapply(plots, ggplot2::ggplotGrob)))
    selected <- plotter(
      res = res, db = "custom", plot_type = "bar",
      id_use = list(Group2 = "t2", Group1 = "t1"), combine = FALSE, verbose = FALSE
    )
    actual <- do.call(rbind, lapply(selected, function(p) p$data[, c("Groups", "ID")]))
    expect_setequal(paste(actual$Groups, actual$ID), c("Group2 t2", "Group1 t1"))
    expect_error(
      plotter(res = res, db = "custom", plot_type = "bar", id_use = list("t1"), verbose = FALSE),
      "must be named"
    )
    expect_error(
      plotter(res = res, db = "custom", plot_type = "bar", id_use = list(Unknown = "t1"), verbose = FALSE),
      "Names in 'id_use' is invalid: Unknown",
      fixed = TRUE
    )
  }
  expect_identical(serialize(res, NULL), original)
})

test_that("grouped term selection preserves order, repeated groups and empty tables", {
  tab <- data.frame(
    ID = c("t2", "t1", "t1", "t2"),
    Groups = factor(c("B", "A", "B", "A"), levels = c("B", "A", "unused"))
  )
  ids <- structure(list("t1", c("t2", "t2"), "t1"), names = c("A", "B", "A"))
  selected <- filter_enrichment_plot_ids(tab, ids)
  expect_identical(selected$ID, c("t1", "t2", "t1"))
  expect_identical(as.character(selected$Groups), c("A", "B", "A"))
  expect_identical(levels(selected$Groups), levels(tab$Groups))
  expect_identical(filter_enrichment_plot_ids(tab, list(A = "missing")), tab[FALSE, , drop = FALSE])
})
