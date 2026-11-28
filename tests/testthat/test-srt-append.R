make_append_identity_object <- function(cells, labels, levels) {
  counts <- Matrix::Matrix(
    matrix(
      1L,
      nrow = 2L,
      ncol = length(cells),
      dimnames = list(c("g1", "g2"), cells)
    ),
    sparse = TRUE
  )
  srt <- SeuratObject::CreateSeuratObject(counts)
  srt@active.ident <- setNames(factor(labels, levels = levels), cells)
  srt
}

test_that("srt_append aligns overwritten identities to target cell names", {
  target <- make_append_identity_object(c("c1", "c2", "c3"), rep("old", 3L), "old")
  source <- make_append_identity_object(
    c("c3", "c1", "c2"), c("B", "C", "A"), c("C", "B", "A", "unused")
  )
  expect_true(methods::validObject(target, test = TRUE))
  expect_true(methods::validObject(source, test = TRUE))

  out <- srt_append(
    target, source, slots = "active.ident",
    pattern = "^irrelevant$", overwrite = TRUE, verbose = FALSE
  )
  expect_true(methods::validObject(out, test = TRUE))
  expect_identical(
    SeuratObject::Idents(out),
    setNames(
      factor(c("C", "A", "B"), levels = c("C", "B", "A", "unused")),
      colnames(target)
    )
  )
  expect_identical(out[["RNA"]], target[["RNA"]])
  expect_identical(out[[]], target[[]])

  already_aligned <- srt_append(
    target, out, slots = "active.ident", overwrite = TRUE, verbose = FALSE
  )
  expect_identical(SeuratObject::Idents(already_aligned), SeuratObject::Idents(out))
})

test_that("srt_append requires a source identity for each overwritten target cell", {
  target <- make_append_identity_object(c("c1", "c2"), c("X", "Y"), c("X", "Y"))
  extra <- make_append_identity_object(
    c("c3", "c2", "c1"), c("Z", "B", "A"), c("Z", "A", "B")
  )
  covered <- srt_append(
    target, extra, slots = "active.ident", overwrite = TRUE, verbose = FALSE
  )
  expect_true(methods::validObject(covered, test = TRUE))
  expect_identical(SeuratObject::Idents(covered), SeuratObject::Idents(extra)[c("c1", "c2")])

  missing <- make_append_identity_object(c("c1", "c3"), c("A", "B"), c("A", "B"))
  expect_error(
    srt_append(target, missing, slots = "active.ident", overwrite = TRUE, verbose = FALSE),
    "^Cannot overwrite Seurat identities:"
  )
  expect_identical(
    srt_append(target, missing, slots = "active.ident", overwrite = FALSE, verbose = FALSE),
    target
  )
})
