test_that("CellStatPlot validates chord stat.by length", {
  srt <- SeuratObject::CreateSeuratObject(
    counts = Matrix::sparseMatrix(
      i = c(1, 2),
      j = c(1, 2),
      x = c(1, 1),
      dims = c(2, 2),
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    )
  )
  srt$celltype <- c("A", "B")

  expect_error(
    CellStatPlot(srt, stat.by = "celltype", plot_type = "chord"),
    "exactly two metadata columns"
  )
})
