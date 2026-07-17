test_that("ccc_pivot_matrix aggregates numeric interaction tables in place", {
  table <- data.frame(
    sender = factor(
      c("B", "A", "B", "A"),
      levels = c("A", "B", "unused")
    ),
    receiver = c("X", "Y", "X", "Y"),
    score = c(1.5, 2, 3.5, 4),
    stringsAsFactors = FALSE
  )

  out <- ccc_pivot_matrix(table, "sender", "receiver", "score")

  expect_identical(rownames(out), c("A", "B"))
  expect_identical(colnames(out), c("X", "Y"))
  expect_equal(
    out,
    matrix(c(NA, 5, 6, NA), 2, 2, dimnames = list(c("A", "B"), c("X", "Y")))
  )
})

test_that("ccc_pivot_matrix preserves legacy missing-value accumulation", {
  table <- data.frame(
    sender = c("A", "A", "A", "B"),
    receiver = c("X", "X", "X", "Y"),
    score = c(2, NA_real_, 5, 7),
    stringsAsFactors = FALSE
  )

  out <- ccc_pivot_matrix(table, "sender", "receiver", "score")

  expect_equal(out["A", "X"], 5)
  expect_equal(out["B", "Y"], 7)
})

test_that("ccc_ligand_target_matrix preserves ordered finite pair sums", {
  table <- data.frame(
    ligand = c("L2", "L1", "L2", "missing"),
    target = c("T1", "T2", "T1", "T2"),
    weight = c(1, 2, 3, Inf),
    stringsAsFactors = FALSE
  )

  out <- ccc_ligand_target_matrix(
    table,
    ligand_levels = c("L1", "L2"),
    target_levels = c("T1", "T2")
  )

  expect_equal(
    out,
    matrix(c(NA, 4, 2, NA), 2, 2,
      dimnames = list(c("L1", "L2"), c("T1", "T2"))
    )
  )
})
