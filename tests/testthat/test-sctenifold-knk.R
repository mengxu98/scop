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

# A gene without counts in any cell is dropped by every subsample, so each
# network has to be re-embedded into the full gene space.
make_sctenifold_dropout_counts <- function(
  n_genes = 40L,
  n_cells = 60L,
  seed = 1L
) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_genes * n_cells, mu = 5, size = 2),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  counts["g2", ] <- 0
  counts
}

build_sctenifold_networks <- function(counts, n_cells) {
  set.seed(11)
  sctenifold_make_networks_cpp(
    counts,
    nNet = 2L,
    nCells = as.integer(n_cells),
    nComp = 3L,
    scaleScores = TRUE,
    symmetric = FALSE,
    q = 0.9,
    nCores = 1L
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

test_that("RunscTenifoldKnk re-embeds subsampled networks with their scores", {
  counts <- make_sctenifold_dropout_counts()
  # 30 sampled cells exercise the dual solver, 50 the primal one
  for (n_cells in c(30L, 50L)) {
    nets <- NULL
    expect_no_warning(nets <- build_sctenifold_networks(counts, n_cells))

    expect_length(nets, 2L)
    for (net in nets) {
      # a pattern (ngCMatrix) result silently drops the regression scores and
      # collapses the network to unweighted logical entries
      expect_s4_class(net, "dgCMatrix")
      expect_identical(dim(net), c(nrow(counts), nrow(counts)))
      expect_identical(rownames(net), rownames(counts))
      expect_gt(Matrix::nnzero(net), nrow(counts))
      expect_gt(length(unique(as.vector(net))), 2L)
      expect_true(all(net["g2", ] == 0))
      expect_true(all(net[, "g2"] == 0))
    }
  }
})

test_that("RunscTenifoldKnk cpp ensembles match the upstream construction", {
  testthat::skip_if_not_installed("scTenifoldNet")
  counts <- make_sctenifold_dropout_counts()
  make_networks <- get_namespace_fun("scTenifoldNet", "makeNetworks")

  for (n_cells in c(30L, 50L)) {
    native <- build_sctenifold_networks(counts, n_cells)
    set.seed(11)
    upstream <- suppressMessages(make_networks(
      X = counts,
      nNet = 2L,
      nCells = as.integer(n_cells),
      nComp = 3L,
      scaleScores = TRUE,
      symmetric = FALSE,
      q = 0.9,
      nCores = 1L
    ))

    expect_length(native, length(upstream))
    for (i in seq_along(native)) {
      expect_identical(dim(native[[i]]), dim(upstream[[i]]))
      expect_identical(rownames(native[[i]]), rownames(upstream[[i]]))
      expect_equal(
        as.matrix(native[[i]]),
        as.matrix(upstream[[i]]),
        tolerance = 1e-8
      )
    }
  }
})
