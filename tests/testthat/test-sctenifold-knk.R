make_sctenifold_knk_srt <- function(n_genes = 60L, n_cells = 80L, seed = 1L) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_genes * n_cells, mu = 5, size = 2),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)
  )
}

run_sctenifold_knk <- function(srt, store_networks) {
  set.seed(11)
  RunscTenifoldKnk(
    srt,
    gKO = "g1",
    qc = FALSE,
    nc_nNet = 2,
    nc_nCells = 40,
    td_maxIter = 50,
    store_networks = store_networks,
    store_manifold = TRUE,
    backend = "cpp",
    verbose = FALSE
  )
}

test_that("RunscTenifoldKnk cpp backend builds the network ensemble natively", {
  srt <- make_sctenifold_knk_srt()
  manifold_calls <- 0L
  # Manifold alignment and differential regulation need the optional
  # RSpectra/RhpcBLASctl/MASS packages; stubbing them keeps this guard running
  # where the check installs hard dependencies only.
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(pkg, fun, ...) {
      stop(pkg, "::", fun, "() is not expected on the cpp backend")
    },
    sctenifold_manifold_cpp = function(x, y, d, cores) {
      manifold_calls <<- manifold_calls + 1L
      genes <- rownames(x)
      matrix(
        0,
        nrow = 2L * length(genes),
        ncol = d,
        dimnames = list(
          c(paste0("X_", genes), paste0("Y_", genes)),
          paste0("NLMA ", seq_len(d))
        )
      )
    },
    sctenifold_dregulation = function(manifold_output, gko) {
      genes <- sub("^X_", "", grep("^X_", rownames(manifold_output), value = TRUE))
      data.frame(
        gene = genes,
        distance = rev(seq_along(genes)),
        Z = 0,
        FC = 1,
        p.value = 1,
        p.adj = 1
      )
    }
  )

  out <- run_sctenifold_knk(srt, store_networks = TRUE)
  dr <- out@tools$scTenifoldKnk$diffRegulation

  expect_identical(manifold_calls, 1L)
  expect_identical(dim(out@tools$scTenifoldKnk$result$tensorNetworks$WT), c(60L, 60L))
  expect_s3_class(dr, "data.frame")
  expect_setequal(dr$gene, paste0("g", seq_len(60L)))
})

test_that("RunscTenifoldKnk store_networks only controls what is stored", {
  testthat::skip_if_not_installed("RSpectra")
  testthat::skip_if_not_installed("RhpcBLASctl")
  testthat::skip_if_not_installed("MASS")
  srt <- make_sctenifold_knk_srt()
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) invisible(TRUE)
  )

  stored <- run_sctenifold_knk(srt, store_networks = TRUE)
  dropped <- run_sctenifold_knk(srt, store_networks = FALSE)

  expect_identical(
    stored@tools$scTenifoldKnk$diffRegulation,
    dropped@tools$scTenifoldKnk$diffRegulation
  )
  expect_false(is.null(stored@tools$scTenifoldKnk$result$tensorNetworks))
  expect_null(dropped@tools$scTenifoldKnk$result$tensorNetworks)
  expect_identical(
    stored@tools$scTenifoldKnk$result$manifoldAlignment,
    dropped@tools$scTenifoldKnk$result$manifoldAlignment
  )
})
