make_moscot_test_srt <- function() {
  counts <- Matrix::Diagonal(4, x = rep(1, 4))
  rownames(counts) <- paste0("g", 1:4)
  colnames(counts) <- paste0("c", 1:4)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$cluster <- factor(c("a", "a", "b", "b"))
  srt
}

test_that("moscot is an optional trajectory environment with a Python-version profile", {
  legacy <- env_requirements(version = "3.10-1", modules = "moscot")
  modern <- env_requirements(version = "3.12-1", modules = "moscot")

  expect_identical(legacy$packages[["moscot"]], "moscot==0.4.3")
  expect_identical(modern$packages[["moscot"]], "moscot==0.5.2")
  expect_true(all(c("scanpy", "cellrank", "moscot") %in% names(modern$packages)))
  expect_identical(modern$packages[["anndata"]], "anndata==0.12.19")
  expect_identical(modern$packages[["jax"]], "jax==0.11.1")
  expect_identical(modern$packages[["jaxlib"]], "jaxlib==0.11.1")
  expect_identical(modern$packages[["flax"]], "flax==0.12.9")
})

test_that("RunCellRank forwards the moscot and real-time configuration", {
  srt <- make_moscot_test_srt()
  captured <- new.env(parent = emptyenv())
  adata_out <- srt
  adata_out$cellrank_moscot_time <- factor(c(0, 0, 1, 1), levels = c(0, 1))

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      captured$modules <- modules
      NULL
    },
    check_python = function(package, ...) {
      captured$checked <- c(captured$checked, package)
      TRUE
    },
    py_to_r2 = function(x) x,
    srt_to_adata = function(object, layer_x, ...) {
      captured$layer_x <- layer_x
      list(obs = data.frame(cluster = factor(c("a", "a", "b", "b"))))
    },
    palette_colors = function(...) c(a = "#111111", b = "#222222"),
    scop_python_import = function(...) {
      list(CellRank = function(...) {
        captured$args <- list(...)
        list(
          adata_out,
          "estimator",
          "kernel",
          list(
            temporal = list(
              time_key = "cellrank_moscot_time",
              time_values = c(0, 1)
            ),
            versions = list(cellrank = "2.3.2", moscot = "0.5.2")
          )
        )
      })
    },
    adata_to_srt = function(...) adata_out
  )

  out <- RunCellRank(
    object = srt,
    group.by = "cluster",
    kernel_type = "moscot",
    backend = "python",
    time_field = "day",
    time_values = c(D0 = 0, D1 = 1),
    moscot_args = list(solve = list(epsilon = 0.05, tau_a = 0.95)),
    realtime_args = list(transition = list(self_transitions = "all", conn_weight = 0.2)),
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(captured$layer_x, "data")
  expect_true("moscot" %in% captured$modules)
  expect_true("moscot" %in% captured$checked)
  expect_equal(unlist(captured$args$time_values, use.names = FALSE), c(0, 1))
  expect_equal(captured$args$moscot_args$solve$epsilon, 0.05)
  expect_identical(captured$args$realtime_args$transition$self_transitions, "all")
  expect_identical(out@tools$CellRank$temporal$time_key, "cellrank_moscot_time")
  expect_identical(out@tools$CellRank$versions$moscot, "0.5.2")
})

test_that("moscot refuses the approximate C++ route and backward transitions", {
  srt <- make_moscot_test_srt()
  expect_error(
    RunCellRank(
      object = srt,
      group.by = "cluster",
      kernel_type = "moscot",
      backend = "cpp",
      allow_approximate = TRUE,
      show_plot = FALSE,
      verbose = FALSE
    ),
    "requires.*backend.*python"
  )
})

test_that("CellRankPlot supports stored experimental-time flow and earliest starts", {
  srt <- make_moscot_test_srt()
  srt$cluster <- factor(c("A", "A", "B", "B"))
  srt$cellrank_moscot_time <- factor(c(0, 0, 1, 1), levels = c(0, 1))
  transition <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L),
    j = c(3L, 4L, 3L, 4L),
    x = c(0.8, 0.7, 1, 1),
    dims = c(4, 4),
    dimnames = list(colnames(srt), colnames(srt))
  )
  srt@tools$CellRank$transition_matrix <- transition
  srt@tools$CellRank$temporal <- list(
    time_key = "cellrank_moscot_time",
    time_values = c(0, 1)
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(8), ncol = 2, dimnames = list(colnames(srt), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_",
    assay = "RNA"
  )

  expect_s3_class(
    CellRankPlot(srt, plot_type = "flow", group.by = "cluster"),
    "ggplot"
  )
  expect_s3_class(
    CellRankPlot(
      srt,
      plot_type = "random_walks",
      reduction = "umap",
      start_cells = c("c1", "c2"),
      n_sims = 2L,
      max_iter = 2L
    ),
    "ggplot"
  )
})
