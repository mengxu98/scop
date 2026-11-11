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
