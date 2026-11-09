local_edition(3)

cellstate_fixture <- function() {
  counts <- matrix(
    c(
      2, 0, 0, 0, 4, 0,
      0, 0, 2, 2, 0, 0,
      0, 1, 0, 2, 1, 3,
      1, 1, 1, 1, 1, 1
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("M", 1:4), paste0("c", 1:6))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- SeuratObject::SetAssayData(
    object,
    layer = "data", new.data = Matrix::Matrix(counts, sparse = TRUE)
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
    fixture$object,
    group.by = "cluster", sample.by = "donor",
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
  expect_equal(
    unname(p$data$program.assignment),
    c("High", "Ambiguous", "Unassigned", "Low", "High", "Unassigned")
  )
  expect_equal(
    as.numeric(tapply(p$data$programs$fraction, p$data$programs$cluster, sum)),
    rep(1, 3)
  )
  expect_equal(
    as.numeric(tapply(p$data$samples$fraction, p$data$samples$cluster, sum)),
    rep(1, 3)
  )
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
  expect_equal(
    a$data$program.assignment[sort(names(a$data$program.assignment))],
    b$data$program.assignment[sort(names(b$data$program.assignment))]
  )
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M1",
    layer = "counts", input.scale = "counts",
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
  expect_equal(
    unname(p$data$program.assignment[colnames(f$counts)]),
    unname(expected.assignment)
  )
})

test_that("constant dot domains and zero-height trees draw without side effects", {
  f <- cellstate_fixture()
  tree <- stats::hclust(stats::dist(diag(3)))
  tree$labels <- c("A", "B", "C")
  tree$height[] <- 0
  directory <- withr::local_tempdir(pattern = "cell-state-construction-")
  withr::local_dir(directory)
  devices <- grDevices::dev.list()
  p <- CellStatePlot(
    f$object, "cluster", "donor", "M4",
    program.by = "donor",
    gsea = f$gsea, tree = tree, family = "sans", return_data = TRUE
  )
  expect_identical(grDevices::dev.list(), devices)
  expect_length(list.files(directory, all.files = TRUE, no.. = TRUE), 0)
  expect_gt(diff(p$plot$encoded$marker.radius(c(.4, .8))), 0)
  grDevices::cairo_pdf(
    withr::local_tempfile(fileext = ".pdf"),
    width = p$width_mm / 25.4, height = p$height_mm / 25.4
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_condition(grid::grid.draw(p$plot))
  expect_no_condition(grid::grid.force())
})

test_that("ambiguous annotation and duplicate enrichment fail clearly", {
  f <- cellstate_fixture()
  f$object$annotation <- c("One", "Two", "One", "One", "One", "One")
  expect_error(
    cellstate_fixture_plot(f, class.by = "annotation"),
    "Mixed annotations in cluster A; supply an explicit map.",
    fixed = TRUE
  )
  f$gsea <- rbind(f$gsea, f$gsea[1, ])
  expect_error(
    cellstate_fixture_plot(f),
    "Duplicate GSEA cluster/pathway rows; select one result per combination.",
    fixed = TRUE
  )
})

test_that("annotation maps take precedence over unused metadata columns", {
  f <- cellstate_fixture()
  map <- c(A = "One", B = "Two", C = "Three")
  p <- cellstate_fixture_plot(
    f,
    class.by = "old_class", subclass.by = "old_subclass",
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
    file, p$plot,
    width = p$width_mm, height = p$height_mm, units = "mm"
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
    object,
    layer = "data", new.data = Matrix::Matrix(counts, sparse = TRUE)
  )
  object$cluster <- rep(c("A", "B", "C"), each = 20)
  object$donor <- rep(c("d1", "d2"), 30)
  tree <- stats::hclust(stats::dist(diag(3)))
  tree$labels <- c("A", "B", "C")
  p <- CellStatePlot(
    object, "cluster", "donor", "M1",
    program.by = "donor",
    gsea = f$gsea, tree = tree, return_data = TRUE
  )
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, family = "sans")
  on.exit(grDevices::dev.off(), add = TRUE)
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
    f$object, "cluster", "donor", "M1",
    program.by = "donor",
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
        title$label[i],
        x = title$x[i], y = title$y[i],
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

test_that("palette changes preserve summaries and respect named overrides", {
  f <- cellstate_fixture()
  p <- cellstate_fixture_plot(f)
  q <- cellstate_fixture_plot(
    f,
    palette = "Set2", expression_palette = "Blues", enrichment_palette = "Spectral",
    palcolor = list(program = c(High = "black"), expression = c("white", "blue"))
  )
  expect_identical(p$data, q$data)
  expect_identical(p$tree, q$tree)
  expected <- thisplot::palette_colors(
    unique(as.character(p$data$programs$category)),
    palette = "Chinese", type = "discrete"
  )
  expect_identical(p$plot$palette$program, expected)
  expect_identical(unname(q$plot$palette$program["High"]), "black")
  expect_false(identical(p$plot$palette$expression, q$plot$palette$expression))
  expect_false(identical(p$plot$palette$enrichment, q$plot$palette$enrichment))
  expect_error(cellstate_fixture_plot(f, palcolor = list(c("red", "blue"))), "named list")
  expect_error(cellstate_fixture_plot(f, palcolor = list(program = c(High = NA))), "color")
})

test_that("file loading stays outside plotting and leaves the session untouched", {
  f <- cellstate_fixture()
  file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(f$object, file)
  before <- Sys.getenv("RETICULATE_USE_MANAGED_VENV", unset = NA_character_)
  expect_error(CellStatePlot(file, "cluster", "donor", "M1", program.by = "donor"), "readRDS")
  expect_error(CellStatePlot("object.h5ad", "cluster", "donor", "M1", program.by = "donor"), "h5ad_to_srt")
  loaded <- CellStatePlot(
    readRDS(file), "cluster", "donor", c("M1", "M2", "M4"),
    program.scores = c(High = "score1", Low = "score2"),
    gsea = f$gsea, return_data = TRUE
  )
  expect_equal(loaded$data, cellstate_fixture_plot(f)$data)
  expect_identical(Sys.getenv("RETICULATE_USE_MANAGED_VENV", unset = NA_character_), before)
})

test_that("short cluster labels keep the composition and dot tracks together", {
  p <- cellstate_fixture_plot(cellstate_fixture())
  a <- p$plot$layout
  expect_lte(a$program.y - a$cluster.bottom, 9)
  expect_gte(a$bar.height, 14)
  expect_lte(a$bar.height, 18)
  expect_lte(a$sample.y - a$program.y - a$bar.height, 5)
  expect_lte(a$gene.top - a$sample.y - a$bar.height, 13)
  expect_lt(a$W, 244)
  expect_gt(a$H, a$bottom)
})

test_that("long labels and many categories draw without text collisions", {
  f <- cellstate_fixture()
  clusters <- c("Alpha progenitor cluster", "Beta differentiated cluster", "Gamma precursor cluster")
  f$object$cluster <- rep(clusters, each = 2)
  f$gsea$Groups <- clusters[match(f$gsea$Groups, c("A", "B", "C"))]
  categories <- paste("Biological sample", seq_len(ncol(f$object)), "extended label")
  f$object$donor <- categories
  p <- CellStatePlot(
    f$object, "cluster", "donor", c("M1", "M2", "M4"),
    program.by = "donor", gsea = f$gsea,
    subclass.map = stats::setNames(paste("Subclass", seq_along(clusters), "annotation"), clusters),
    class.map = stats::setNames(rep("Broad class", length(clusters)), clusters),
    feature.labels = c(M1 = "Long marker gene label", M2 = "Second marker", M4 = "Third marker"),
    pathway.labels = c(H_ONE = "First pathway\nwith multiple lines", H_TWO = "Second pathway label"),
    return_data = TRUE
  )
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = p$width_mm / 25.4, height = p$height_mm / 25.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- grid::grid.force(p$plot)
  grid::pushViewport(g$children[[1]]$vp)
  on.exit(grid::popViewport(), add = TRUE, after = FALSE)
  labels <- Filter(function(g) {
    inherits(g, "text") && length(g$label) == 1L &&
      g$label %in% c(categories, "Program assignment", "Donors / samples", "Not evaluable")
  }, as.list(g$children[[1]]$children))
  expect_equal(sum(vapply(labels, function(g) g$label %in% categories, logical(1))), 2L * length(categories))
  boxes <- vapply(labels, function(g) {
    c(
      left = grid::convertX(grid::grobX(g, "west"), "bigpts", valueOnly = TRUE),
      right = grid::convertX(grid::grobX(g, "east"), "bigpts", valueOnly = TRUE),
      bottom = grid::convertY(grid::grobY(g, "south"), "bigpts", valueOnly = TRUE),
      top = grid::convertY(grid::grobY(g, "north"), "bigpts", valueOnly = TRUE)
    )
  }, numeric(4))
  pairs <- utils::combn(ncol(boxes), 2)
  overlaps <- apply(pairs, 2, function(pair) {
    a <- boxes[, pair[1]]
    b <- boxes[, pair[2]]
    min(a["right"], b["right"]) > max(a["left"], b["left"]) &&
      min(a["top"], b["top"]) > max(a["bottom"], b["bottom"])
  })
  expect_false(any(overlaps))
  expect_gt(min(boxes["left", ]), 0)
  expect_lt(max(boxes["right", ]), p$plot$layout$W)
})

test_that("large sample legends wrap instead of widening the canvas", {
  f <- cellstate_fixture()
  counts <- f$counts[, rep(seq_len(ncol(f$counts)), 20)]
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object <- SeuratObject::SetAssayData(object, layer = "data", new.data = Matrix::Matrix(counts, sparse = TRUE))
  object$cluster <- rep(c("A", "B", "C"), each = 40)
  widths <- numeric()
  for (n in c(20L, 50L)) {
    object$donor <- paste("Biological sample", rep(seq_len(n), length.out = ncol(object)), "extended label")
    p <- CellStatePlot(
      object, "cluster", "donor", c("M1", "M2", "M4"),
      program.by = "donor", gsea = f$gsea, return_data = TRUE
    )
    expect_setequal(names(p$plot$palette$sample), unique(object$donor))
    expect_equal(nrow(p$data$samples), 3L * n)
    expect_lt(p$plot$layout$W, 500)
    file <- withr::local_tempfile(fileext = ".pdf")
    expect_no_condition(ggplot2::ggsave(
      file, p$plot,
      width = p$width_mm, height = p$height_mm, units = "mm"
    ))
    widths <- c(widths, p$plot$layout$W)
  }
  expect_lte(widths[2], widths[1] * 1.2)
})

test_that("plots omit repeated class headings and keep titles optional", {
  f <- cellstate_fixture()
  map <- c(A = "Progenitor", B = "Progenitor", C = "Differentiated")
  p <- cellstate_fixture_plot(f, class.map = map)
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = p$width_mm / 25.4, height = p$height_mm / 25.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  labels <- function(plot) {
    children <- as.list(grid::grid.force(plot)$children[[1]]$children)
    unlist(lapply(Filter(function(g) inherits(g, "text"), children), `[[`, "label"))
  }
  text <- labels(p$plot)
  expect_null(formals(CellStatePlot)$title)
  expect_false(any(c("Cell states", unname(map)) %in% text))
  expect_true(all(c("A", "B", "C") %in% text))
  expect_equal(p$data$clusters$class, unname(map[p$data$clusters$cluster]))
  q <- cellstate_fixture_plot(f, class.map = map, title = "Custom title\nSecond line")
  expect_true(all(c("Custom title", "Second line") %in% labels(q$plot)))
  expect_identical(p$data, q$data)
  expect_gt(q$height_mm, p$height_mm)
})

test_that("GSEA rows are matched by both cluster and pathway IDs", {
  f <- cellstate_fixture()
  groups <- c("A\rB", "A", "C")
  f$object$cluster <- rep(groups, each = 2)
  f$gsea <- data.frame(
    cluster = c("A\rB", "A", "C", "A\rB"),
    feature = c("C", "B\rC", "C", "B\rC"),
    nes = c(1, -2, 3, 4), p_value = c(.1, .2, NA, .4),
    status = c("EVALUABLE", "EVALUABLE", "NOT_EVALUABLE", "NOT_EVALUABLE")
  )
  original <- f$gsea
  p <- cellstate_fixture_plot(f)
  g <- p$data$pathways
  expect_identical(f$gsea, original)
  expect_equal(g[c("cluster", "feature")], expand.grid(
    cluster = sort(groups), feature = c("C", "B\rC"),
    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE
  ))
  first <- g$cluster == "A\rB" & g$feature == "C"
  second <- g$cluster == "A" & g$feature == "B\rC"
  missing <- g$cluster == "A" & g$feature == "C"
  excluded <- g$cluster == "A\rB" & g$feature == "B\rC"
  expect_equal(g$nes[first], 1)
  expect_equal(g$nes[second], -2)
  expect_equal(g$p_value[first], .1)
  expect_equal(g$p_value[second], .2)
  expect_true(is.na(g$nes[missing]))
  expect_equal(g$status[missing], "NOT_EVALUABLE")
  expect_equal(g$p_value[excluded], .4)
  expect_true(is.na(g$neg_log_p[excluded]))
})

test_that("small enrichment values have distinct measured and drawn size labels", {
  f <- cellstate_fixture()
  f$gsea$pvalue[] <- .99
  p <- cellstate_fixture_plot(f)
  encoded <- p$plot$encoded
  expect_length(unique(encoded$size.labels), 3)
  expect_equal(as.numeric(encoded$size.labels), encoded$size.ticks, tolerance = .01)
  file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = p$width_mm / 25.4, height = p$height_mm / 25.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  children <- as.list(grid::grid.force(p$plot)$children[[1]]$children)
  labels <- Filter(function(g) {
    inherits(g, "text") && identical(g$label, encoded$size.labels)
  }, children)
  expect_length(labels, 1)
  a <- p$plot$layout
  expected.x <- a$key.x[3] + a$measure("-ln(P)") + 6 +
    a$key.offsets(encoded$size.labels) + 4
  expect_equal(as.numeric(labels[[1]]$x), expected.x)
  expect_lte(max(expected.x + a$measure(encoded$size.labels)), a$W - 4)
})

test_that("flip rotates the whole plot and preserves numerical results", {
  f <- cellstate_fixture()
  p <- cellstate_fixture_plot(f)
  q <- cellstate_fixture_plot(f, flip = TRUE)
  expect_false(formals(CellStatePlot)$flip)
  expect_identical(p$data, q$data)
  expect_identical(p$tree, q$tree)
  expect_identical(p$provenance, q$provenance)
  expect_identical(p$plot$palette, q$plot$palette)
  expect_identical(p$plot$colors, q$plot$colors)
  geometry <- names(p$plot$layout)[!vapply(p$plot$layout, is.function, logical(1))]
  expect_identical(p$plot$layout[geometry], q$plot$layout[geometry])
  expect_identical(p$width_mm, q$height_mm)
  expect_identical(p$height_mm, q$width_mm)
  for (result in list(p, q)) {
    file <- withr::local_tempfile(fileext = ".pdf")
    local({
      withr::local_pdf(file, width = result$width_mm / 25.4, height = result$height_mm / 25.4)
      expect_no_condition(grid::grid.draw(result$plot))
      drawn <- grid::grid.force(result$plot)
      vp <- drawn$children[[1]]$vp
      expect_equal(vp$angle, if (result$plot$flip) -90 else 0)
      expect_equal(grid::convertWidth(vp$width, "mm", valueOnly = TRUE), p$width_mm)
      expect_equal(grid::convertHeight(vp$height, "mm", valueOnly = TRUE), p$height_mm)
    })
    expect_gt(file.info(file)$size, 1000)
  }
})

test_that("flip requires a single nonmissing logical value", {
  f <- cellstate_fixture()
  for (value in list(NULL, NA, c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(cellstate_fixture_plot(f, flip = value), "flip must be TRUE or FALSE.", fixed = TRUE)
  }
})

test_that("unmatched marker warnings retain the gene IDs", {
  f <- cellstate_fixture()
  expect_warning(CellStatePlot(
    f$object, "cluster", "donor", c("M1", "absent{marker}"),
    program.by = "donor", gsea = f$gsea
  ), "Unmatched markers retained as missing: absent{marker}", fixed = TRUE)
})
