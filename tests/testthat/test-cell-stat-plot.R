make_cell_stat_srt <- function() {
  srt <- SeuratObject::CreateSeuratObject(
    counts = Matrix::sparseMatrix(
      i = rep(1:2, each = 6),
      j = rep(1:6, times = 2),
      x = rep(1, 12),
      dims = c(2, 6),
      dimnames = list(c("gene1", "gene2"), paste0("cell", 1:6))
    )
  )
  srt$celltype <- c("A", "B", "A", "B", "C", "C")
  srt$group <- c("E1", "E1", "E2", "E2", "E3", "E3")
  srt
}

test_that("CellStatPlot validates chord stat.by length", {
  srt <- make_cell_stat_srt()

  expect_error(
    CellStatPlot(srt, stat.by = "celltype", plot_type = "chord"),
    "exactly two metadata columns"
  )
})

test_that("CellStatPlot exposes direct x-axis text angle control", {
  srt <- make_cell_stat_srt()

  plot <- CellStatPlot(
    srt,
    stat.by = "celltype",
    group.by = "group",
    plot_type = "trend",
    x_text_angle = 0
  )

  expect_identical(plot$theme$axis.text.x$angle, 0)
})

test_that("CellStatPlot validates x_text_angle", {
  srt <- make_cell_stat_srt()

  expect_error(
    CellStatPlot(
      srt,
      stat.by = "celltype",
      group.by = "group",
      plot_type = "trend",
      x_text_angle = NA_real_
    ),
    "finite number"
  )
})

test_that("CellStatPlot forwards major grid settings to StatPlot", {
  srt <- make_cell_stat_srt()

  plot <- CellStatPlot(
    srt,
    stat.by = "celltype",
    group.by = "group",
    plot_type = "trend",
    grid_major = TRUE,
    grid_major_colour = "red",
    grid_major_linetype = 3,
    grid_major_linewidth = 0.7
  )

  expect_s3_class(plot$theme$panel.grid.major, "element_line")
  expect_identical(plot$theme$panel.grid.major$colour, "red")
  expect_identical(plot$theme$panel.grid.major$linetype, 3)
  expect_identical(plot$theme$panel.grid.major$linewidth, 0.7)
})

test_that("CellStatPlot match.arg defaults work without explicit choices", {
  srt <- make_cell_stat_srt()

  expect_no_error(
    CellStatPlot(srt, stat.by = "celltype", group.by = "group")
  )
})

test_that("CellStatPlot resolves the exported theme_scop alias", {
  srt <- make_cell_stat_srt()

  expect_true(is.function(scop::theme_scop))
  expect_equal(scop::theme_scop, thisplot::theme_this)

  plot <- CellStatPlot(
    srt,
    stat.by = "celltype",
    group.by = "group",
    theme_use = "theme_scop"
  )
  expect_s3_class(plot, "ggplot")
})

test_that("sankey selections preserve unweighted rows and explicit weights", {
  data <- data.frame(
    Group = c("A", "A", "B"),
    Type = c("X", "Y", "X"),
    weight = c(1, 2, 3)
  )
  unweighted <- thisplot::make_long(data, Group, Type, value = NULL)
  weighted <- thisplot::make_long(data, Group, Type, value = weight)
  expect_identical(as.character(unweighted$node), c("A", "X", "A", "Y", "B", "X"))
  expect_identical(as.character(unweighted$next_node), c("X", NA, "Y", NA, "X", NA))
  expect_false("value" %in% names(unweighted))
  expect_equal(weighted$value, rep(data$weight, each = 2))
  expect_equal(weighted[names(unweighted)], unweighted)
  expect_equal(thisplot::make_long(data, Group = Group, Type = Type), unweighted)

  object <- make_cell_stat_srt()
  before <- object
  plot <- CellStatPlot(
    object, stat.by = c("celltype", "group"),
    plot_type = "sankey", stat_type = "count"
  )
  expect_s3_class(plot, "patchwork")
  expect_no_warning(patchwork::patchworkGrob(plot))
  expect_identical(object, before)
})

test_that("composition styling retains size units, visible category keys and plot shape", {
  srt <- make_cell_stat_srt()
  args <- list(object = srt, stat.by = "celltype", group.by = "group")
  bar <- do.call(CellStatPlot, args)
  expect_identical(bar$theme$axis.text.x$angle, 0)
  expect_s3_class(bar$theme$panel.grid.major.x, "element_blank")
  expect_s3_class(bar$guides$guides$fill, "GuideLegend")
  custom_grid <- do.call(CellStatPlot, c(args, list(grid_major_colour = "red")))
  expect_identical(custom_grid$theme$panel.grid.major$colour, "red")

  p <- do.call(CellStatPlot, c(args, list(plot_type = "dot")))
  expect_identical(p$theme$axis.text.x$angle, 0)
  expect_identical(p$scales$get_scales("size")$labels(0.25), "25%")
  expect_s3_class(p$guides$guides$fill, "GuideLegend")
  hidden <- do.call(CellStatPlot, c(args, list(plot_type = "dot", legend.position = "none")))
  expect_identical(hidden$theme$legend.position, "none")
  expect_equal(p$scales$get_scales("y")$expand, c(0, 1, 0, 1))

  flipped <- do.call(CellStatPlot, c(args, list(
    plot_type = "dot", flip = TRUE, stat_type = "count", x_text_angle = 30
  )))
  expect_identical(flipped$theme$axis.text.x$angle, 30)
  expect_s3_class(flipped$theme$panel.grid.major.x, "element_line")
  expect_false(is.function(flipped$scales$get_scales("size")$labels))
  expect_no_error(ggplot2::ggplotGrob(flipped))

  custom <- do.call(CellStatPlot, c(args, list(
    plot_type = "dot",
    theme_use = function() ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
  )))
  expect_s3_class(custom$theme$axis.text.y, "element_blank")
  expect_s3_class(custom$guides$guides$fill, "GuideLegend")

  srt$split <- rep(c("left", "right"), 3)
  args$object <- srt
  args$split.by <- "split"
  separate <- do.call(CellStatPlot, c(args, list(plot_type = "dot", combine = FALSE)))
  combined <- do.call(CellStatPlot, c(args, list(plot_type = "dot", nrow = 1)))
  expect_length(separate, 2L)
  expect_s3_class(combined, "patchwork")
  expect_identical(combined[[1]]$data, separate[[1]]$data)
  expect_s3_class(combined[[2]]$guides$guides$fill, "GuideLegend")
  expect_no_error(patchwork::patchworkGrob(combined))
})
