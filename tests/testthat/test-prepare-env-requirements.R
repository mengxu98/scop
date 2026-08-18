test_that("PrepareEnv module specs keep scanpy and scvelo distinct", {
  specs <- getFromNamespace("env_module_requirements", "scop")()
  scanpy <- specs$scanpy
  scvelo <- specs$scvelo

  expect_true("scanpy" %in% names(scanpy$packages))
  expect_false("scvelo" %in% names(scanpy$packages))
  expect_true("scvelo" %in% names(scvelo$packages))
  expect_false("scanpy" %in% names(scvelo$packages))
})

test_that("PrepareEnv resolves scvelo through modules rather than a dedicated helper", {
  req <- env_requirements(modules = "scvelo")
  expect_true("scanpy" %in% names(req$packages))
  expect_true("scvelo" %in% names(req$packages))
  expect_false(exists("scvelo_python_requirements", envir = asNamespace("scop"), inherits = FALSE))
  expect_false(exists("scanpy_python_requirements", envir = asNamespace("scop"), inherits = FALSE))
})
