test_that("GroupHeatmap yields the same library-normalized values from metadata and counts", {
  counts <- Matrix::Matrix(
    rbind(target = c(1, 4, 2, 3), reference = c(3, 0, 6, 7)),
    sparse = TRUE
  )
  colnames(counts) <- paste0("cell", 1:4)
  object <- SeuratObject::CreateSeuratObject(counts)
  object$group <- factor(c("B", "A", "B", "A"), levels = c("A", "B"))
  draw <- function(x) {
    GroupHeatmap(
      x, features = "target", group.by = "group",
      exp_method = "raw", lib_normalize = TRUE, verbose = FALSE
    )$matrix_list$group
  }
  expected <- matrix(c(3.9, 1.5), nrow = 1,
    dimnames = list("target", c("A", "B"))
  )
  expect_equal(draw(object), expected)
  object$nCount_RNA <- NULL
  expect_equal(draw(object), expected)
})
