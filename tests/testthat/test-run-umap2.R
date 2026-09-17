test_that("RunUMAP2 appends spread without shifting existing formals", {
  seurat_args <- names(formals(RunUMAP2.Seurat))
  default_args <- names(formals(RunUMAP2.default))

  expect_identical(
    seurat_args[match("min.dist", seurat_args) + 1L],
    "set.op.mix.ratio"
  )
  expect_identical(
    default_args[match("min.dist", default_args) + 1L],
    "set.op.mix.ratio"
  )
  expect_gt(match("spread", seurat_args), match("seed.use", seurat_args))
  expect_gt(match("spread", default_args), match("seed.use", default_args))
  expect_true("n_threads" %in% seurat_args)
  expect_true("n_threads" %in% default_args)
  expect_true("spread" %in% seurat_args)
  expect_true("spread" %in% default_args)
})
