mirror_fixture <- function() {
  tab <- expand.grid(ID = c("one", "two", "three"), Database = c("GO_BP", "GO_CC"),
    Groups = c("A", "B"), stringsAsFactors = FALSE)
  tab$Description <- paste("Complete pathway", tab$ID)
  tab$p.adjust <- rep(c(0.001, 0.01, 0.03), 4)
  tab$p.adjust[tab$Groups == "B"] <- tab$p.adjust[tab$Groups == "B"]^2
  tab$geneID <- "G1/G2/G3"
  list(enrichment = tab)
}

test_that("mirror results preserve ranking, full labels, colors and numeric geometry", {
  p <- EnrichmentMirrorPlot(res = mirror_fixture(), group_use = c("A", "B"),
    db = c("GO_CC", "GO_BP"), topTerm = 2, palcolor = c(B = "red", A = "blue"))
  expect_equal(nrow(p$data), 8)
  expect_equal(p$data$.term, p$data$Description)
  expect_equal(p$data$.xmax - p$data$.xmin, -log10(p$data$p.adjust))
  expect_true(all(p$data$.xmin[p$data$Groups == "A"] < 0))
  expect_true(all(p$data$.xmax[p$data$Groups == "B"] > 0))
  expect_false(any(vapply(p$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1))))
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[5]]$colour == "black"))
  expect_equal(built$data[[4]]$fill, built$data[[6]]$colour)
  expect_equal(built$data[[4]]$alpha, built$data[[6]]$alpha)
  expect_equal(built$data[[4]]$fill[p$data$Groups == "A"], rep("blue", 4))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("mirror axes and group order adapt without changing result values", {
  p <- EnrichmentMirrorPlot(res = mirror_fixture(), group_use = c("B", "A"), db = "GO")
  expect_true(all(p$data$.xmin[p$data$Groups == "B"] < 0))
  q <- EnrichmentMirrorPlot(res = mirror_fixture(), group_use = c("B", "A"),
    db = "GO", symmetric_axis = TRUE, term_size = 4, font_family = "sans")
  expect_equal(p$data$.score, q$data$.score)
  expect_equal(q$layers[[7]]$data$x, -q$layers[[8]]$data$x)
  expect_no_warning(ggplot2::ggplotGrob(q))
})

test_that("mirror handles zero, missing and unbalanced database results", {
  res <- mirror_fixture()
  res$enrichment$p.adjust[1] <- 0
  res$enrichment$p.adjust[2] <- NA_real_
  res$enrichment <- subset(res$enrichment, !(Groups == "B" & Database == "GO_CC"))
  expect_warning(p <- EnrichmentMirrorPlot(res = res, group_use = c("A", "B"), db = "GO"), "missing adjusted")
  expect_true(all(is.finite(p$data$.score)))
  expect_equal(max(p$data$.score), -log10(.Machine$double.xmin))
  expect_no_warning(ggplot2::ggplotGrob(p))
  expect_error(EnrichmentMirrorPlot(res = res, group_use = "A"), "exactly two")
  expect_error(EnrichmentMirrorPlot(res = mirror_fixture(), group_use = c("A", "C"), db = "GO"), "Both selected")
  expect_error(EnrichmentMirrorPlot(res = mirror_fixture(), group_use = c("A", "B"), band_gap = -1), "band_gap")
})

test_that("mirror reads the saved object key and permits custom result lists", {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1:12, 3,
    dimnames = list(letters[1:3], paste0("c", 1:4))), sparse = TRUE))
  object@tools$Enrichment_custom_wilcox <- mirror_fixture()
  p <- EnrichmentMirrorPlot(object, group_use = c("A", "B"), db = "GO")
  expect_equal(nrow(p$data), 12)
  expect_error(EnrichmentMirrorPlot(object, group.by = "wrong", group_use = c("A", "B")), "Missing saved")
})
