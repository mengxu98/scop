local_edition(3)

cellstate_fixture <- function() {
  counts <- matrix(
    c(2, 0, 0, 0, 4, 0,
      0, 0, 2, 2, 0, 0,
      0, 1, 0, 2, 1, 3,
      1, 1, 1, 1, 1, 1),
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("M", 1:4), paste0("c", 1:6))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- SeuratObject::SetAssayData(
    object, layer = "data", new.data = Matrix::Matrix(counts, sparse = TRUE)
  )
  object$cluster <- rep(c("A", "B", "C"), each = 2)
  object$donor <- c("d1", "d2", "d1", "d2", "d1", "d1")
  object$score1 <- c(3, 1, NA, 0, 2, NA)
  object$score2 <- c(0, 1, NA, 2, 0, NA)
  gsea <- expand.grid(Groups = c("A", "B", "C"), ID = c("H_ONE", "H_TWO"))
  gsea$Database <- "H"
  gsea$NES <- c(1, -1, NA, 2, 0, -2)
  gsea$pvalue <- c(1, .05, NA, .01, 0, .5)
  gsea$p.adjust <- c(1, .1, NA, .02, 0, .8)
  list(object = object, counts = counts, gsea = gsea)
}

cellstate_fixture_plot <- function(fixture, ...) {
  CellStatePlot(
    fixture$object, group.by = "cluster", sample.by = "donor",
    features = c("M1", "M2", "M4"),
    program.scores = c(High = "score1", Low = "score2"),
    gsea = fixture$gsea, return_data = TRUE, ...
  )
}

test_that("CellStatePlot summaries agree with independent expectations", {
  f <- cellstate_fixture()
  original <- f$object
  p <- cellstate_fixture_plot(f)
  m <- p$data$markers
  expect_equal(m$mean[m$feature == "M1"], c(1, 0, 2))
  expect_equal(m$mean_z[m$feature == "M1"], c(0, -1, 1) / sqrt(2 / 3))
  expect_equal(m$proportion[m$feature == "M1"], c(.5, 0, .5))
  expect_equal(m$mean_z[m$feature == "M4"], rep(0, 3))
  expect_equal(unname(p$data$program.assignment),
               c("High", "Ambiguous", "Unassigned", "Low", "High", "Unassigned"))
  expect_equal(as.numeric(tapply(p$data$programs$fraction, p$data$programs$cluster, sum)),
               rep(1, 3))
  expect_equal(as.numeric(tapply(p$data$samples$fraction, p$data$samples$cluster, sum)),
               rep(1, 3))
  c.samples <- subset(p$data$samples, cluster == "C")
  expect_equal(c.samples$fraction[as.character(c.samples$category) == "d1"], 1)
  expect_identical(f$object, original)
  expect_s3_class(p$plot, "grob")
  expect_gt(p$width_mm, 0)
  expect_gt(p$height_mm, 0)
})

test_that("stored RunGSEA, missing results and raw P values retain their meaning", {
  f <- cellstate_fixture()
  f$object@tools$GSEA_cluster_wilcox <- list(enrichment = f$gsea)
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M1",
    program.by = "donor", return_data = TRUE
  )
  g <- p$data$pathways
  expect_equal(g$radius[g$p_value == 1 & !is.na(g$p_value)], 0)
  expect_equal(g$p_value[g$cluster == "B" & g$feature == "H_TWO"], 0)
  expect_equal(g$neg_log_p[g$cluster == "B" & g$feature == "H_TWO"], -log(1e-10))
  expect_equal(g$status[g$cluster == "C" & g$feature == "H_ONE"], "NOT_EVALUABLE")
  expect_equal(sum(is.na(g$neg_log_p)), 1)
})

test_that("cell-ID alignment and counts normalization are preserved", {
  f <- cellstate_fixture()
  a <- cellstate_fixture_plot(f)
  f$object <- f$object[, c(6, 2, 4, 1, 5, 3)]
  b <- cellstate_fixture_plot(f, tree = a$tree)
  expect_equal(a$data$markers, b$data$markers)
  expect_equal(a$data$samples, b$data$samples)
  expect_equal(a$data$programs, b$data$programs)
  expect_equal(a$data$program.assignment[sort(names(a$data$program.assignment))],
               b$data$program.assignment[sort(names(b$data$program.assignment))])
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M1", layer = "counts", input.scale = "counts",
    programs = list(One = "M1", Two = "M2"), min.program.genes = 1,
    gsea = f$gsea, tree = a$tree, return_data = TRUE
  )
  normalized <- log1p(sweep(f$counts, 2, colSums(f$counts) / 1e4, "/"))
  expected <- rowMeans(normalized["M1", c("c1", "c2"), drop = FALSE])
  expect_equal(p$data$markers$mean[p$data$markers$cluster == "A"], unname(expected))
  expect_equal(sort(names(p$data$program.assignment)), colnames(f$counts))
  expect_equal(p$data$program.coverage$status, rep("EVALUABLE", 2))
  standardized <- normalized[c("M1", "M2"), ] -
    rowMeans(normalized[c("M1", "M2"), ])
  standardized <- standardized / sqrt(rowMeans(standardized^2))
  expected.assignment <- ifelse(standardized[1, ] > standardized[2, ], "One", "Two")
  expect_equal(unname(p$data$program.assignment[colnames(f$counts)]),
               unname(expected.assignment))
})

test_that("constant dot domains and zero-height trees draw without side effects", {
  f <- cellstate_fixture()
  tree <- stats::hclust(stats::dist(diag(3)))
  tree$labels <- c("A", "B", "C")
  tree$height[] <- 0
  directory <- tempfile("cell-state-construction-")
  dir.create(directory)
  withr::local_dir(directory)
  devices <- grDevices::dev.list()
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M4", program.by = "donor",
    gsea = f$gsea, tree = tree, return_data = TRUE
  )
  expect_identical(grDevices::dev.list(), devices)
  expect_length(list.files(directory, all.files = TRUE, no.. = TRUE), 0)
  expect_gt(diff(p$plot$encoded$marker.radius(c(.4, .8))), 0)
  grDevices::cairo_pdf(
    tempfile(fileext = ".pdf"), width = p$width_mm / 25.4, height = p$height_mm / 25.4
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_condition(grid::grid.draw(p$plot))
  expect_no_condition(grid::grid.force())
})

test_that("ambiguous annotation and duplicate enrichment fail clearly", {
  f <- cellstate_fixture()
  f$object$annotation <- c("One", "Two", "One", "One", "One", "One")
  expect_snapshot(error = TRUE, {
    cellstate_fixture_plot(f, class.by = "annotation")
  })
  f$gsea <- rbind(f$gsea, f$gsea[1, ])
  expect_snapshot(error = TRUE, {
    cellstate_fixture_plot(f)
  })
})
