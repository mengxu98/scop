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

test_that("PrepareEnv retains Imports without per-method checks", {
  packages <- env_r_packages()
  expect_true(all(c(
    "reticulate", "ComplexHeatmap", "igraph", "Signac",
    "SummarizedExperiment", "S4Vectors", "thisplot"
  ) %in% packages))
})

test_that("scprep modules keep a separate compatible dependency profile", {
  expect_false(any(c("phate", "magic") %in% default_env_modules()))
  expect_identical(env_module_environment("sccoda"), "tensorflow_env")
  expect_identical(env_requirements(modules = "sccoda")$packages[["rpy2"]], "rpy2==3.5.17")
  expect_identical(env_module_environment("scomm"), "tensorflow_env")
  expect_identical(env_module_environment(c("scvelo", "magic")), "scprep_env")
  expect_identical(env_requirements(modules = c("scvelo", "magic"))$packages[["pandas"]], "pandas==2.0.3")
  expect_identical(env_requirements(modules = "phate")$packages[["phate"]], "phate==1.0.11")
  expect_identical(env_requirements(modules = "scrublet")$packages[["pandas"]], "pandas==2.2.0")
})
