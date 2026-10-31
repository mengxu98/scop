test_that("GroupHeatmap validates custom legend lists", {
  skip_if_not_installed("ComplexHeatmap")
  for (params in list(NULL, TRUE, list(list()), list(other = list()),
    structure(list(list(), list()), names = c("expression", "expression")))) {
    expect_error(GroupHeatmap(legend_params = params), "legend_params")
  }
  for (params in list(NULL, TRUE, list("horizontal"), list(direcion = "horizontal"),
    structure(list(1, 2), names = c("nrow", "nrow")))) {
    expect_error(GroupHeatmap(legend_params = list(percent = params)), "legend_params")
  }
})

test_that("GroupHeatmap renders custom legends without changing expression data", {
  skip_if_not_installed("ComplexHeatmap")
  counts <- matrix(rep(0:5, 12), nrow = 6,
    dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:12)))
  counts[1:3, 1:6] <- counts[1:3, 1:6] + 5
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  object$celltype <- rep(c("A", "B"), each = 6)
  captured <- list()
  original_legend <- get_namespace_fun("ComplexHeatmap", "Legend")
  capture_legend <- function(...) {
    args <- lapply(as.list(match.call())[-1L], eval, envir = parent.frame())
    legend <- do.call(original_legend, args)
    captured[[args$title]] <<- list(args = args, grob = legend@grob)
    legend
  }
  formals(capture_legend) <- formals(original_legend)
  local_mocked_bindings(Legend = capture_legend, .package = "ComplexHeatmap")
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  args <- list(object = object, features = rownames(object), group.by = "celltype",
    layer = "counts", add_dot = TRUE, feature_split = NULL,
    cluster_rows = FALSE, cluster_columns = FALSE, nlabel = 0, verbose = FALSE,
    legend.position = "bottom")
  default <- do.call(GroupHeatmap, args)
  expect_identical(captured$Percent$args$direction, "vertical")
  default_legends <- captured
  custom <- do.call(GroupHeatmap, c(args, list(legend_params = list(
    expression = list(direction = "horizontal", legend_width = grid::unit(35, "mm"),
      title_position = "topcenter", title = "Expression"),
    percent = list(direction = "horizontal", nrow = 1, title_position = "topcenter")
  ))))
  print(custom$plot)
  expect_equal(custom$matrix_list, default$matrix_list)
  expect_identical(captured$Expression$args$direction, "horizontal")
  expect_equal(captured$Percent$args$nrow, 1)
  expect_equal(captured$Percent$args$size, default_legends$Percent$args$size)
  expect_equal(captured$celltype$args, default_legends$celltype$args)
  for (name in c("Expression", "Percent")) {
    grob <- captured[[name]]$grob
    expect_gt(grid::convertWidth(grid::grobWidth(grob), "mm", valueOnly = TRUE),
      grid::convertHeight(grid::grobHeight(grob), "mm", valueOnly = TRUE))
  }
  args$add_dot <- FALSE
  captured <- list()
  expect_no_error(do.call(GroupHeatmap, c(args,
    list(legend_params = list(percent = list(nrow = 1))))))
  expect_null(captured$Percent)
})
