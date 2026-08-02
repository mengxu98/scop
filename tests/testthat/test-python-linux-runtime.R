test_that("CUDA runtime discovery is inactive without a PyTorch installation", {
  env_path <- tempfile("scop-conda-env-")
  dir.create(file.path(env_path, "lib"), recursive = TRUE)
  file.create(file.path(env_path, "lib", "libcublas.so.12"))

  expect_identical(find_conda_cuda_libraries(env_path), character())
})

test_that("CUDA runtime discovery finds conda and pip NVIDIA libraries in load order", {
  env_path <- tempfile("scop-conda-env-")
  site_packages <- file.path(
    env_path,
    "lib",
    "python3.10",
    "site-packages"
  )
  cublas_dir <- file.path(site_packages, "nvidia", "cublas", "lib")
  dir.create(file.path(site_packages, "torch"), recursive = TRUE)
  dir.create(cublas_dir, recursive = TRUE)
  target_dir <- file.path(env_path, "targets", "x86_64-linux", "lib")
  dir.create(target_dir, recursive = TRUE)

  file.create(file.path(target_dir, "libnvJitLink.so.12"))
  file.create(file.path(env_path, "lib", "libcudart.so.12"))
  file.create(file.path(cublas_dir, "libcublasLt.so.12"))
  file.create(file.path(cublas_dir, "libcublas.so.12"))
  file.create(file.path(cublas_dir, "libcublas.so.12.9.1.4"))

  libraries <- find_conda_cuda_libraries(env_path)

  expect_identical(
    basename(libraries),
    c(
      "libnvJitLink.so.12",
      "libcudart.so.12",
      "libcublasLt.so.12",
      "libcublas.so.12"
    )
  )
})

test_that("shared-library preloading retries dependencies loaded out of order", {
  # LD_PRELOAD is a Linux-only mechanism; the preload strategy is a
  # no-op elsewhere and the environment variable is not meaningful.
  skip_on_os(c("mac", "windows"))

  first <- tempfile(fileext = ".so")
  second <- tempfile(fileext = ".so")
  file.create(first, second)
  attempts <- new.env(parent = emptyenv())
  attempts$first <- 0L
  old_preload <- Sys.getenv("LD_PRELOAD", unset = NA_character_)
  on.exit(
    {
      if (is.na(old_preload)) {
        Sys.unsetenv("LD_PRELOAD")
      } else {
        Sys.setenv(LD_PRELOAD = old_preload)
      }
    },
    add = TRUE
  )

  testthat::local_mocked_bindings(
    load_shared_library = function(lib) {
      if (identical(lib, first)) {
        attempts$first <- attempts$first + 1L
        if (attempts$first == 1L) {
          stop("dependency not loaded")
        }
      }
      invisible(NULL)
    },
    .package = "scop"
  )

  loaded <- preload_shared_libraries(c(first, second))

  expect_identical(loaded, c(second, first))
  expect_identical(attempts$first, 2L)
  expect_true(all(c(first, second) %in% strsplit(
    Sys.getenv("LD_PRELOAD"),
    .Platform$path.sep,
    fixed = TRUE
  )[[1]]))
})
