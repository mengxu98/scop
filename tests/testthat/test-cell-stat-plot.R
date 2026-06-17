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
