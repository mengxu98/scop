test_that("real large-data RunPAGA roundtrips Seurat through AnnData", {
  skip_if_not_installed("callr")
  config <- reticulate::py_discover_config()
  skip_if(is.null(config), "No existing Python interpreter available")
  python <- config$python

  python_check <- suppressWarnings(system2(
    python,
    c("-c", shQuote(
      "import scanpy, pynndescent; print(pynndescent.__version__)"
    )),
    stdout = TRUE,
    stderr = TRUE
  ))
  python_status <- attr(python_check, "status")
  skip_if(!is.null(python_status) && python_status != 0L)
  expect_identical(tail(python_check, 1L), "0.6.0")

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  output <- callr::r(function(package_path, libpath, python) {
    .libPaths(libpath)
    reticulate::use_python(python, required = TRUE)
    reticulate::py_config()
    if (file.exists(file.path(package_path, ".Rbuildignore"))) {
      pkgload::load_all(package_path, quiet = TRUE, compile = FALSE)
    } else {
      loadNamespace("scop")
    }

    testthat::local_mocked_bindings(
      PrepareEnv = function(...) invisible(NULL),
      .package = "scop"
    )

    getFromNamespace("scop_python_import", "scop")("functions", convert = FALSE)
    reticulate::py_run_string(paste(
      "import sys, types",
      "from unittest.mock import patch",
      "from functions import _check_paga_pynndescent_runtime",
      "with patch.dict(sys.modules, {'pynndescent': None}):",
      "    _check_paga_pynndescent_runtime()",
      "with patch('importlib.metadata.version', return_value='0.6.0'):",
      "    with patch.dict(sys.modules, {'pynndescent': types.SimpleNamespace(__version__='0.6.0')}):",
      "        _check_paga_pynndescent_runtime()",
      "    with patch.dict(sys.modules, {'pynndescent': types.SimpleNamespace(__version__='0.7.0')}):",
      "        try:",
      "            _check_paga_pynndescent_runtime()",
      "        except RuntimeError as error:",
      "            assert 'Restart R' in str(error)",
      "            assert '0.7.0 is still loaded' in str(error)",
      "            assert '0.6.0 is installed' in str(error)",
      "        else:",
      "            raise AssertionError('Stale PyNNDescent was accepted')",
      sep = "\n"
    ))

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
  }, args = list(package_path, .libPaths(), python), user_profile = FALSE, show = TRUE)

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
