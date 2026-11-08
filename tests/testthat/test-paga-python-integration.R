# Run in a fresh R process with RETICULATE_PYTHON pointing to a prepared
# SCOP Scanpy environment. The test never provisions or upgrades Python packages.
test_that("real large-data RunPAGA roundtrips Seurat through AnnData", {
  skip_if_not(nzchar(Sys.getenv("RETICULATE_PYTHON")))
  skip_if_not_installed("callr")

  skip_if_not(reticulate::py_module_available("scanpy"))
  skip_if_not(reticulate::py_module_available("pynndescent"))
  expect_identical(
    reticulate::py_eval(
      "import importlib.metadata; importlib.metadata.version('pynndescent')"
    ),
    "0.6.0"
  )

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  output <- callr::r(function(package_path, libpath) {
    .libPaths(libpath)
    if (file.exists(file.path(package_path, ".Rbuildignore"))) {
      pkgload::load_all(package_path, quiet = TRUE)
    } else {
      library(scop)
    }

    # The Python environment is prepared before this test. Mock only package
    # provisioning so the regression test cannot install or update dependencies.
    testthat::local_mocked_bindings(
      PrepareEnv = function(...) invisible(NULL),
      .package = "scop"
    )

    set.seed(42)
    counts <- Matrix::Matrix(
      matrix(rpois(40L * 8200L, 2), nrow = 40L),
      sparse = TRUE
    )
    dimnames(counts) <- list(
      paste0("gene", seq_len(40)),
      paste0("cell", seq_len(8200))
    )
    object <- Seurat::CreateSeuratObject(counts)
    object$group <- factor(rep(c("A", "B", "C", "D"), each = 2050))
    pcs <- matrix(
      rnorm(8200L * 10L),
      nrow = 8200L,
      dimnames = list(colnames(object), paste0("PC_", seq_len(10)))
    )
    object[["PCA"]] <- SeuratObject::CreateDimReducObject(
      pcs,
      key = "PC_",
      assay = "RNA"
    )
    umap <- pcs[, 1:2]
    colnames(umap) <- c("UMAP_1", "UMAP_2")
    object[["UMAP"]] <- SeuratObject::CreateDimReducObject(
      umap,
      key = "UMAP_",
      assay = "RNA"
    )

    scop::RunPAGA(
      object = object,
      group.by = "group",
      linear_reduction = "PCA",
      nonlinear_reduction = "UMAP",
      n_pcs = 10,
      n_neighbors = 10,
      assay_y = NULL,
      backend = "python",
      show_plot = FALSE,
      save_plot = FALSE,
      dirpath = withr::local_tempdir(),
      verbose = FALSE
    )
  }, args = list(package_path, .libPaths()), user_profile = FALSE, show = TRUE)

  expect_s4_class(output, "Seurat")
  expect_identical(colnames(output), paste0("cell", seq_len(8200)))
  expect_identical(
    as.character(output$group),
    rep(c("A", "B", "C", "D"), each = 2050)
  )
  expect_identical(output@tools$PAGA$backend, "python")
  expect_identical(output@tools$PAGA$groups, "group")
  graph <- as.matrix(output@tools$PAGA$connectivities)
  expect_equal(dim(graph), c(4L, 4L))
  expect_true(all(is.finite(graph)))
  expect_equal(graph, t(graph), ignore_attr = TRUE)
  expect_true(all(graph >= 0 & graph <= 1))
  expect_true(any(graph > 0))
  expect_null(output@misc$paga)
})
