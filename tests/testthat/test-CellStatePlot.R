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

test_that("annotation maps take precedence over unused metadata columns", {
  f <- cellstate_fixture()
  map <- c(A = "One", B = "Two", C = "Three")
  p <- cellstate_fixture_plot(
    f, class.by = "old_class", subclass.by = "old_subclass",
    class.map = map, subclass.map = map
  )
  expect_equal(p$data$clusters$class, unname(map[p$data$clusters$cluster]))
  expect_equal(p$data$clusters$subclass, unname(map[p$data$clusters$cluster]))
  expect_error(cellstate_fixture_plot(f, class.by = "old_class"), "Missing metadata")
  expect_error(cellstate_fixture_plot(f, class.map = map[-1]), "cover every cluster")
})

test_that("default PDF export needs no font registration", {
  p <- cellstate_fixture_plot(cellstate_fixture())
  file <- withr::local_tempfile(fileext = ".pdf")
  expect_no_condition(ggplot2::ggsave(
    file, p$plot, width = p$width_mm, height = p$height_mm, units = "mm"
  ))
  expect_gt(file.info(file)$size, 1000)
})

test_that("marker size legends describe the displayed high fractions", {
  f <- cellstate_fixture()
  counts <- rbind(
    M1 = c(rep(1, 18), 0, 0, rep(1, 19), 0, rep(1, 20)),
    M2 = rep(1, 60)
  )
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- SeuratObject::SetAssayData(
    object, layer = "data", new.data = Matrix::Matrix(counts, sparse = TRUE)
  )
  object$cluster <- rep(c("A", "B", "C"), each = 20)
  object$donor <- rep(c("d1", "d2"), 30)
  tree <- stats::hclust(stats::dist(diag(3)))
  tree$labels <- c("A", "B", "C")
  p <- CellStatePlot(
    object, "cluster", "donor", "M1", program.by = "donor",
    gsea = f$gsea, tree = tree, return_data = TRUE
  )
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, family = "sans")
  on.exit(grDevices::dev.off(), add = TRUE)
  p$plot$layout$family <- "sans"
  children <- as.list(grid::grid.force(p$plot)$children[[1]]$children)
  legend <- Filter(function(g) {
    inherits(g, "circle") && identical(unname(g$gp$fill), rep("#D3D3D3", 3))
  }, children)
  labels <- unlist(lapply(Filter(function(g) inherits(g, "text"), children), `[[`, "label"))
  expect_equal(p$data$markers$proportion, c(.9, .95, 1))
  expect_length(legend, 1)
  expect_true(all(diff(as.numeric(legend[[1]]$r)) > 0))
  expect_equal(as.numeric(legend[[1]]$r), p$data$markers$radius)
  expect_true(any(grepl("0\\.9", labels)))
})

test_that("single-row panels reserve space for their vertical titles", {
  f <- cellstate_fixture()
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M1", program.by = "donor",
    gsea = f$gsea, pathways = "H_ONE", family = "sans", return_data = TRUE
  )
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = p$width_mm / 25.4, height = p$height_mm / 25.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- grid::grid.force(p$plot)
  vp <- g$children[[1]]$vp
  grid::pushViewport(vp)
  on.exit(grid::popViewport(), add = TRUE, after = FALSE)
  titles <- Filter(function(g) {
    inherits(g, "text") && all(g$label %in% c("Marker gene", "expression", "Pathway", "enrichment"))
  }, as.list(g$children[[1]]$children))
  expect_length(titles, 2)
  bounds <- lapply(titles, function(title) {
    vapply(seq_along(title$label), function(i) {
      t <- grid::textGrob(
        title$label[i], x = title$x[i], y = title$y[i],
        just = c("left", "bottom"), rot = 90, gp = title$gp
      )
      c(
        north = grid::convertY(grid::grobY(t, "north"), "bigpts", valueOnly = TRUE),
        south = grid::convertY(grid::grobY(t, "south"), "bigpts", valueOnly = TRUE)
      )
    }, numeric(2))
  })
  expect_gt(min(bounds[[1]]["south", ]), max(bounds[[2]]["north", ]))
})

test_that("AnnData objects work with either reticulate conversion mode", {
  skip_if_not_installed("reticulate")
  if (!reticulate::py_available(initialize = FALSE)) {
    skip("Requires an existing initialized Python environment")
  }
  if (!reticulate::py_module_available("anndata") ||
      !reticulate::py_module_available("scipy.sparse")) {
    skip("Existing Python environment needs anndata and scipy")
  }
  f <- cellstate_fixture()
  expected <- cellstate_fixture_plot(f)
  scipy <- reticulate::import("scipy.sparse", convert = FALSE)
  for (convert in c(TRUE, FALSE)) {
    ad <- reticulate::import("anndata", convert = convert)
    for (sparse in c(FALSE, TRUE)) {
      x <- t(f$counts)
      if (sparse) x <- scipy$csr_matrix(x)
      object <- ad$AnnData(
        X = x, obs = f$object[[]], var = data.frame(row.names = rownames(f$counts))
      )
      p <- CellStatePlot(
        object, "cluster", "donor", c("M1", "M2", "M4"),
        program.scores = c(High = "score1", Low = "score2"),
        gsea = f$gsea, return_data = TRUE
      )
      expect_equal(p$data, expected$data)
      expect_equal(p$tree, expected$tree)
    }
  }
})
