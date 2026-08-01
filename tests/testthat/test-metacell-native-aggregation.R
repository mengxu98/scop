test_that("metacell native helpers preserve KNN graph and count aggregation", {
  knn_idx <- rbind(
    c(2L, 3L),
    c(1L, 3L),
    c(1L, 2L)
  )
  actual_graph <- metacell_knn_adjacency(knn_idx, n_cells = 3L)
  expected_graph <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 2L, 3L, 3L),
    j = c(2L, 3L, 1L, 3L, 1L, 2L),
    x = 1,
    dims = c(3L, 3L)
  )
  expect_equal(as.matrix(actual_graph), as.matrix(expected_graph))

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 4, 0, 0, 5, 6), nrow = 3L),
    sparse = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3))
  )
  membership <- c("MC2", "MC1", "MC2")
  actual_counts <- metacell_aggregate_counts(counts, membership)
  expected_counts <- cbind(
    MC2 = Matrix::rowSums(counts[, c(1L, 3L), drop = FALSE]),
    MC1 = Matrix::rowSums(counts[, 2L, drop = FALSE])
  )

  expect_equal(as.matrix(actual_counts), expected_counts)
  expect_identical(colnames(actual_counts), c("MC2", "MC1"))
})

test_that("MetaCellPlot supports one or more metacell centroids", {
  cells <- paste0("cell", seq_len(6))
  counts <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("gene1", "gene2"), cells)
  )
  original <- Seurat::CreateSeuratObject(counts = counts)
  original[["CellType"]] <- factor(rep(c("A", "B"), each = 3))
  embedding <- matrix(
    seq_len(12) / 12,
    nrow = 6,
    dimnames = list(cells, c("PC_1", "PC_2"))
  )
  original[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "PC_",
    assay = "RNA"
  )

  metacell_counts <- Matrix::Matrix(
    matrix(seq_len(4), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("gene1", "gene2"), c("mc1", "mc2"))
  )
  metacells <- Seurat::CreateSeuratObject(counts = metacell_counts)
  metacells[["CellType"]] <- factor(c("A", "B"))
  metacells@misc[["original_srt"]] <- original

  metacells@misc[["cell_membership"]] <- stats::setNames(
    rep(c("mc1", "mc2"), each = 3),
    cells
  )
  expect_s3_class(
    MetaCellPlot(
      metacells,
      reduction = "pca",
      group.by = "CellType"
    ),
    "ggplot"
  )

  metacells@misc[["cell_membership"]] <- stats::setNames(
    rep("mc1", length(cells)),
    cells
  )
  expect_s3_class(
    MetaCellPlot(
      metacells,
      reduction = "pca",
      group.by = "CellType"
    ),
    "ggplot"
  )
})
