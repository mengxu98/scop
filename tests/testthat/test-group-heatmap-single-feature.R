test_that("GroupHeatmap single-gene split heatmaps agree with multi-gene heatmaps", {
  counts <- Matrix::Matrix(
    rbind(Rpsa = c(1, 0, 4, 2), Reference = c(9, 10, 6, 8)),
    sparse = TRUE
  )
  colnames(counts) <- paste0("Cell", 1:4)
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$celltype.1 <- factor(c("B", "A", "B", "A"), levels = c("A", "B"))
  object$group <- factor(c("CTRL", "CTRL", "Day1", "Day1"))

  single <- GroupHeatmap(
    object,
    features = "Rpsa",
    group.by = "celltype.1",
    split.by = "group",
    cell_split_palette = "Dark2"
  )
  expected <- matrix(
    (c(0, 2, 1, 4) - 1.75) / sqrt(35 / 12),
    nrow = 1,
    dimnames = list("Rpsa", c("A : CTRL", "A : Day1", "B : CTRL", "B : Day1"))
  )
  expect_identical(dimnames(single$matrix_list$celltype.1), dimnames(expected))
  expect_equal(c(single$matrix_list$celltype.1), c(expected))
  expect_identical(rownames(single$feature_metadata), "Rpsa")

  multiple <- GroupHeatmap(
    object,
    features = c("Rpsa", "Reference"),
    group.by = "celltype.1",
    split.by = "group",
    cell_split_palette = "Dark2"
  )
  expect_equal(
    c(single$matrix_list$celltype.1),
    c(multiple$matrix_list$celltype.1["Rpsa", , drop = FALSE])
  )

  flipped_dots <- GroupHeatmap(
    object,
    features = "Rpsa",
    group.by = "celltype.1",
    split.by = "group",
    cell_split_palette = "Dark2",
    flip = TRUE,
    add_dot = TRUE
  )
  expect_identical(dimnames(flipped_dots$matrix_list$celltype.1), dimnames(expected))
  expect_equal(c(flipped_dots$matrix_list$celltype.1), c(expected))
})
