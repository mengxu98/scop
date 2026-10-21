annotated_enrichment_fixture <- function() {
  tab <- expand.grid(term = 1:3, Database = c("GO_BP", "GO_CC"), Groups = c("A", "B"), stringsAsFactors = FALSE)
  tab$ID <- paste0("term", tab$term)
  tab$Description <- "Repeated description"
  tab$p.adjust <- c(0.001, 0.02, 0.04)[tab$term]
  tab$pvalue <- tab$p.adjust/2
  tab$Count <- tab$term + 1
  tab$geneID <- c("A/B", "C/D/E", NA_character_)[tab$term]
  list(enrichment = tab)
}

test_that("annotated bars select terms per database and preserve duplicate labels", {
  plots <- EnrichmentPlot(res = annotated_enrichment_fixture(), db = c("GO_CC", "GO_BP"), plot_type = "annotated_bar",
    topTerm = 2, combine = FALSE, color_by = "Database", palcolor = c(GO_BP = "blue", GO_CC = "red"))
  expect_length(plots, 2)
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_equal(as.character(p$data$Database), rep(c("GO_CC", "GO_BP"), each = 2))
    expect_equal(p$data$ID, rep(c("term1", "term2"), 2))
    expect_equal(length(unique(p$data$.y)), 4)
    expect_equal(p$data$.genes, rep(c("A, B", "C, D, E"), 2))
    expect_equal(p$data$.score, -log10(p$data$p.adjust))
    built <- ggplot2::ggplot_build(p)
    expect_equal(unname(built$plot$scales$get_scales("fill")$map(c("GO_CC", "GO_BP"))), c("red",
      "blue"))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_equal(plots[[1]]$scales$get_scales("x")$limits, plots[[2]]$scales$get_scales("x")$limits)
})

test_that("annotated bars use GSEA-style opacity on database colors", {
  plots <- EnrichmentPlot(res = annotated_enrichment_fixture(), db = "GO", plot_type = "annotated_bar",
    combine = FALSE)
  expect_s3_class(plots[[1]]$scales$get_scales("fill"), "ScaleDiscrete")
  scale <- plots[[1]]$scales$get_scales("alpha")
  expect_s3_class(scale, "ScaleContinuous")
  expect_equal(scale$name, "-log10(p.adjust)")
  expect_equal(scale$limits, plots[[2]]$scales$get_scales("alpha")$limits)
  layers <- ggplot2::ggplot_build(plots[[1]])$data
  expect_equal(range(layers[[2]]$alpha), c(0.4, 1))
  expect_equal(layers[[2]]$alpha, layers[[3]]$alpha)
  expect_equal(layers[[2]]$alpha, layers[[5]]$alpha)
  expect_equal(layers[[2]]$fill, layers[[3]]$fill)
  expect_equal(layers[[2]]$fill, layers[[5]]$colour)
  expect_equal(plots[[1]]$scales$get_scales("size")$get_breaks(), c(2, 3, 4))
})

test_that("annotated bar spacing is explicit and validated", {
  p <- EnrichmentPlot(res = annotated_enrichment_fixture(), db = "GO",
    group_use = "A", plot_type = "annotated_bar", row_step = 0.9, db_gap = 0.3)
  expect_equal(diff(p$data$.y), c(-0.9, -0.9, -1.2, -0.9, -0.9))
  expect_error(EnrichmentPlot(res = annotated_enrichment_fixture(),
    plot_type = "annotated_bar", row_step = 0.4), "row_step")
  expect_error(EnrichmentPlot(res = annotated_enrichment_fixture(),
    plot_type = "annotated_bar", db_gap = -1), "db_gap")
})

test_that("annotated bars honor explicit IDs, metrics and panel splitting", {
  res <- annotated_enrichment_fixture()
  plots <- EnrichmentPlot(res = res, db = "GO", plot_type = "annotated_bar", id_use = list(B = "term3"),
    padjustCutoff = NULL, pvalueCutoff = 0.001, split_by = c("Database", "Groups"), combine = FALSE)
  expect_length(plots, 2)
  expect_true(all(vapply(plots, function(p) nrow(p$data) == 1L, logical(1))))
  expect_equal(plots[[1]]$data$ID, "term3")
  expect_equal(plots[[1]]$data$.genes, "")
  expect_equal(plots[[1]]$labels$x, "-log10(pvalue)")
  expect_s3_class(EnrichmentPlot(res = res, db = "GO", plot_type = "annotated_bar"), "patchwork")
})

test_that("annotated bars handle zero, missing and invalid values explicitly", {
  res <- annotated_enrichment_fixture()
  res$enrichment$p.adjust[1:2] <- c(0, NA_real_)
  expect_warning(p <- EnrichmentPlot(res = res, db = "GO_BP", group_use = "A", plot_type = "annotated_bar"),
    "missing p-values")
  expect_true(all(is.finite(p$data$.score)))
  expect_equal(p$data$.score[1], -log10(.Machine$double.xmin))
  res$enrichment$p.adjust[1] <- -1
  expect_error(EnrichmentPlot(res = res, plot_type = "annotated_bar"), "between 0 and 1")
  res$enrichment$p.adjust[] <- 1
  expect_error(EnrichmentPlot(res = res, plot_type = "annotated_bar"), "No term enriched")
  expect_error(EnrichmentPlot(res = annotated_enrichment_fixture(), plot_type = "annotated_bar",
    topTerm = 0), "positive integer")
})

test_that("optional annotations do not require unused columns", {
  res <- annotated_enrichment_fixture()
  res$enrichment$geneID <- res$enrichment$Count <- NULL
  expect_error(EnrichmentPlot(res = res, plot_type = "annotated_bar"), "Missing annotated-bar columns")
  p <- EnrichmentPlot(res = res, plot_type = "annotated_bar", group_use = "A", show_genes = FALSE,
    show_count = FALSE, character_width = 8)
  expect_true(all(nchar(p$data$.term) <= 8))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("term labels sit inside bars and database bands render with rounded corners", {
  p <- EnrichmentPlot(res = annotated_enrichment_fixture(), db = "GO",
    group_use = "A", plot_type = "annotated_bar")
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[4]]$y > built$data[[2]]$ymin))
  expect_true(all(built$data[[4]]$y < built$data[[2]]$ymax))
  expect_true(all(built$data[[5]]$y < built$data[[2]]$ymin))
  grob <- enrichment_round_band$draw_panel(
    built$data[[1]], built$layout$panel_params[[1]], built$plot$coordinates
  )
  expect_true(all(vapply(grob$children, inherits, logical(1), "roundrect")))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("annotated bars read the established Seurat enrichment result key", {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1:12, 3, dimnames = list(c("a",
    "b", "c"), paste0("cell", 1:4))), sparse = TRUE))
  object@tools$Enrichment_celltype_wilcox <- annotated_enrichment_fixture()
  p <- EnrichmentPlot(object, group.by = "celltype", group_use = "B", db = "GO", plot_type = "annotated_bar")
  expect_equal(unique(as.character(p$data$Groups)), "B")
  expect_equal(nrow(p$data), 6)
})
test_that("database colors align across bands, dots and gene labels", {
  plots <- EnrichmentPlot(
    res = annotated_enrichment_fixture(), db = c("GO_CC", "GO_BP"),
    plot_type = "annotated_bar", combine = FALSE,
    palcolor = c(GO_BP = "#228833", GO_CC = "#4477AA")
  )
  for (p in plots) {
    layers <- ggplot2::ggplot_build(p)$data
    expect_equal(layers[[1]]$fill, c("#4477AA", "#228833"))
    expect_equal(layers[[3]]$fill, rep(c("#4477AA", "#228833"), each = 3))
    expect_equal(layers[[5]]$colour, layers[[3]]$colour)
    expect_equal(layers[[2]]$fill, layers[[5]]$colour)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})

