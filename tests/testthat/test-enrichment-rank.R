rank_fixture <- function() {
  list(enrichment = data.frame(
    ID = c("a", "b", "c", "a", "b", "c", "x"),
    Description = paste("Term", c("a", "b", "c", "a", "b", "c", "x")),
    Groups = c(rep("A", 3), rep("B", 3), "C"), Database = "custom",
    p.adjust = c(0, 0.02, 0.2, 0.01, 0.03, 0.04, 0.8),
    FoldEnrichment = c(2, 8, 99, 3, 7, 5, 10)
  ))
}

test_that("rank plots filter before ranking and preserve group/size contracts", {
  res <- rank_fixture()
  original <- res
  p <- EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    group_use = c("B", "A", "C"), topTerm = 2, rank_nlabel = 1,
    rank_pt_size = c(3, 6), combine = FALSE
  )
  expect_named(p, c("B", "A", "C"))
  expect_identical(p$B$data$ID, c("b", "c"))
  expect_identical(p$A$data$ID, c("b", "a"))
  expect_equal(p$A$data$.rank, 1:2)
  expect_equal(sum(p$B$data$.highlight), 1)
  expect_true(all(is.finite(p$A$data$.significance)))
  expect_equal(p$A$data$p.adjust[2], 0)
  expect_identical(res, original)
  built <- ggplot2::ggplot_build(p$A)
  point_layer <- which(vapply(p$A$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))[1]
  expect_true(all(built$data[[point_layer]]$size >= 3 & built$data[[point_layer]]$size <= 6))
  expect_true(all(built$data[[point_layer]]$shape == 16))
  expect_s3_class(p$C, "ggplot")
  expect_silent(ggplot2::ggplotGrob(p$A))
})

test_that("soft rank contours preserve point scores and stay ordered", {
  res <- rank_fixture()
  res$enrichment <- res$enrichment[rep(4, 80), ]
  res$enrichment$ID <- paste0("t", 1:80)
  res$enrichment$FoldEnrichment <- c(seq(12, 11.5, length.out = 40), seq(5, 4.5, length.out = 40))
  p <- EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank", combine = FALSE
  )
  built <- ggplot2::ggplot_build(p$B)
  expect_s3_class(p$B$layers[[1]]$geom, "GeomRibbon")
  expect_true(all(built$data[[1]]$ymax > built$data[[1]]$ymin))
  expect_true(all(is.finite(built$data[[1]]$ymax)))
  expect_equal(p$B$data$.score, res$enrichment$FoldEnrichment)
  ribbons <- Filter(function(layer) inherits(layer$geom, "GeomRibbon"), p$B$layers)
  expect_length(ribbons, 6)
  expect_true(all(vapply(ribbons, function(layer) {
    max(layer$data$.x) <= 40 || min(layer$data$.x) >= 41
  }, logical(1))))
  res <- rank_fixture()
  res$enrichment <- res$enrichment[1, ]
  singleton <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank", combine = FALSE)
  expect_false(any(vapply(
    singleton[[1]]$layers,
    function(layer) inherits(layer$geom, "GeomRibbon"), logical(1)
  )))
})

test_that("rank labels select displayed terms without filtering background", {
  p <- EnrichmentPlot(
    res = rank_fixture(), db = "custom", plot_type = "rank",
    topTerm = Inf, rank_label_terms = list(B = "c"), combine = FALSE
  )
  expect_equal(nrow(p$B$data), 3)
  expect_identical(p$B$data$ID[p$B$data$.highlight], "c")
  expect_false(any(p$A$data$.highlight))
  expect_false(any(grepl("^#", p$B$data$.label)))
})

test_that("rank handles singleton and invalid controls", {
  res <- rank_fixture()
  res$enrichment <- res$enrichment[1, ]
  p <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank", rank_nlabel = 0)
  expect_s3_class(p, "patchwork")
  expect_equal(nrow(attr(p, "enrichment_data")), 1)
  expect_silent(patchwork::patchworkGrob(p))
  expect_error(EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    rank_pt_size = c(5, 2)
  ), "[Rr]ank_pt_size")
  expect_error(EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    rank_nlabel = -1
  ), "[Rr]ank_nlabel")
  expect_error(EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    rank_by = "missing"
  ), "rank_by")
  expect_error(EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    rank_label_terms = list(unknown = "a")
  ), "[Rr]ank_label_terms")
  res$enrichment$FoldEnrichment <- NA_real_
  expect_warning(expect_error(EnrichmentPlot(
    res = res, db = "custom",
    plot_type = "rank"
  ), "No valid terms"), "Omitting")
})

test_that("rank resolves stored results and composite group selectors", {
  res <- rank_fixture()
  res$enrichment$Groups <- paste0(res$enrichment$Groups, "::up")
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(1, 3, 3,
    sparse = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("cell", 1:3))
  ))
  object@tools$Enrichment_custom_wilcox <- res
  p <- EnrichmentPlot(object,
    db = "custom", plot_type = "rank",
    group_use = "A", combine = FALSE
  )
  expect_named(p, "A::up")
  expect_equal(nrow(p[[1]]$data), 2)
})

test_that("rank shares axes and respects colors and raw P filtering", {
  res <- rank_fixture()
  res$enrichment$pvalue <- c(0.01, 0.2, 0.2, 0.01, 0.2, 0.2, 0.2)
  p <- EnrichmentPlot(
    res = res, db = "custom", plot_type = "rank",
    padjustCutoff = NULL, pvalueCutoff = 0.05, rank_shared_y = TRUE,
    palcolor = c(A = "red", B = "blue", C = "grey"), combine = FALSE
  )
  expect_equal(nrow(p$A$data), 1)
  expect_equal(p$A$coordinates$limits$y, p$B$coordinates$limits$y)
  expect_equal(ggplot2::ggplot_build(p$A)$data[[1]]$colour, "red")
  res$enrichment$Database[4:7] <- "other"
  expect_error(EnrichmentPlot(
    res = res, db = c("custom", "other"),
    plot_type = "rank"
  ), "exactly one database")
})

test_that("rank supports grid_major toggle, topTerm, right legend and capital labels", {
  res <- rank_fixture()
  p_default <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank", combine = FALSE)
  expect_true(inherits(p_default$A$theme$panel.grid.major, "element_line"))
  expect_identical(p_default$A$theme$panel.grid.major$linetype, 2)
  expect_true(grepl("^[A-Z]", p_default$A$data$.label[1]))

  # Test grid_major = FALSE
  p_nogrid <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank", grid_major = FALSE, combine = FALSE)
  expect_true(inherits(p_nogrid$A$theme$panel.grid, "element_blank"))

  # Test topTerm
  p_term <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank", topTerm = 1, combine = FALSE)
  expect_equal(nrow(p_term$A$data), 1)

  p_comb <- EnrichmentPlot(res = res, db = "custom", plot_type = "rank")
  expect_identical(p_comb$theme$legend.position, "right")
})
