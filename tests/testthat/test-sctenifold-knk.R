make_sctenifold_knk_srt <- function(n_genes = 60L, n_cells = 80L, seed = 1L) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_genes * n_cells, mu = 5, size = 2),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
}

mock_sctenifold_knk_upstream <- function(run, version = "2.0.0") {
  testthat::local_mocked_bindings(
    .package = "scop",
    .env = parent.frame(),
    check_r = function(packages, ...) {
      expect_identical(packages, "scTenifoldKnk")
      invisible(TRUE)
    },
    get_namespace_fun = function(pkg, fun) {
      if (pkg == "utils" && fun == "packageVersion") {
        return(function(...) numeric_version(version))
      }
      if (pkg == "scTenifoldKnk" && fun == "scTenifoldKnk") {
        return(run)
      }
      stop("Unexpected upstream call: ", pkg, "::", fun)
    }
  )
}

test_that("RunscTenifoldKnk passes raw selected counts to the official workflow", {
  srt <- make_sctenifold_knk_srt()
  received <- NULL
  result <- list(
    tensorNetworks = list(WT = Matrix::Diagonal(3)),
    manifoldAlignment = matrix(seq_len(12), 6, 2),
    diffRegulation = data.frame(gene = c("g5", "g3", "g1"), FC = 1, p.adj = 1),
    upstream_field = "retained"
  )
  mock_sctenifold_knk_upstream(function(countMatrix, ...) {
    received <<- c(list(countMatrix = countMatrix), list(...))
    result
  })
  out <- RunscTenifoldKnk(
    srt,
    gKO = c("g1", "g1"), features = c("g5", "g3"),
    qc_min_library_size = 12, qc_min_pct = 0.2,
    qc_remove_outlier_cells = FALSE, qc_mt_threshold = 0.3,
    nc_priorNetwork = data.frame(regulators = "g1", targets = "g3"),
    ma_method = "heat", ma_heatT = 4, dr_empiricalNull = TRUE,
    dr_direction = FALSE, dr_directionT = 2,
    dr_directionRegressLibSize = TRUE, seed = 7, cores = 2,
    tool_name = "knockout", verbose = FALSE
  )
  expect_identical(
    received$countMatrix,
    GetAssayData5(srt, assay = "RNA", layer = "counts")[c("g5", "g3", "g1"), , drop = FALSE]
  )
  expect_identical(received$gKO, "g1")
  expect_identical(received$qc_minLibSize, 12)
  expect_identical(received$qc_minPCT, 0.2)
  expect_identical(received$qc_removeOutlierCells, FALSE)
  expect_identical(received$qc_maxMTratio, 0.3)
  expect_identical(received$nc_priorNetwork, data.frame(regulators = "g1", targets = "g3"))
  expect_identical(received$ma_method, "heat")
  expect_identical(received$ma_heatT, 4)
  expect_identical(received$dr_empiricalNull, TRUE)
  expect_identical(received$dr_direction, FALSE)
  expect_identical(received$dr_directionT, 2)
  expect_identical(received$dr_directionRegressLibSize, TRUE)
  expect_identical(received$nCores, 2)
  expect_identical(received$seed, 7)
  expect_identical(out@tools$knockout$result, result)
  expect_identical(out@tools$knockout$diffRegulation, result$diffRegulation)
  expect_identical(out@tools$knockout$package_version, "2.0.0")
})

test_that("RunscTenifoldKnk storage flags only trim requested upstream fields", {
  srt <- make_sctenifold_knk_srt()
  result <- list(
    tensorNetworks = list(WT = Matrix::Diagonal(3)),
    manifoldAlignment = matrix(seq_len(12), 6, 2),
    diffRegulation = data.frame(gene = "g1", FC = 1, p.adj = 1)
  )
  mock_sctenifold_knk_upstream(function(...) result)
  out <- RunscTenifoldKnk(
    srt,
    gKO = "g1", store_networks = FALSE, store_manifold = FALSE,
    verbose = FALSE
  )
  expect_null(out@tools$scTenifoldKnk$result$tensorNetworks)
  expect_null(out@tools$scTenifoldKnk$result$manifoldAlignment)
  expect_identical(out@tools$scTenifoldKnk$diffRegulation, result$diffRegulation)
})

test_that("RunscTenifoldKnk validates SCOP inputs before invoking upstream", {
  srt <- make_sctenifold_knk_srt()
  mock_sctenifold_knk_upstream(function(...) stop("Unexpected calculation"))
  expect_error(RunscTenifoldKnk(matrix(1), "g1"), "Seurat")
  expect_error(RunscTenifoldKnk(srt), "gKO")
  expect_error(RunscTenifoldKnk(srt, NA_character_), "gene symbols")
  expect_error(RunscTenifoldKnk(srt, "absent"), "absent")
  expect_error(
    RunscTenifoldKnk(srt, transcriptomeWide = TRUE, features = "absent"),
    "No requested"
  )
})

test_that("RunscTenifoldKnk rejects upstream versions before 2.0", {
  srt <- make_sctenifold_knk_srt()
  mock_sctenifold_knk_upstream(function(...) stop("Unexpected calculation"), version = "1.1.5")
  expect_error(RunscTenifoldKnk(srt, "g1"), "scTenifoldKnk.*2.0.0")
})

test_that("RunscTenifoldKnk matches scTenifoldKnk 2.0 across supported modes", {
  skip_if_not_installed("scTenifoldKnk", minimum_version = "2.0.0")
  skip_if_not_installed("scTenifoldNet", minimum_version = "1.4.3")
  srt <- make_sctenifold_knk_srt()
  counts <- GetAssayData5(srt, assay = "RNA", layer = "counts")
  testthat::local_mocked_bindings(check_r = function(...) invisible(TRUE), .package = "scop")
  wrapper_args <- list(
    object = srt, gKO = "g1", qc = FALSE, nc_nNet = 3, nc_nCells = 40,
    td_maxIter = 50, cores = 1, seed = 11, verbose = FALSE
  )
  upstream_args <- list(
    countMatrix = counts, gKO = "g1", qc = FALSE, nc_nNet = 3, nc_nCells = 40,
    td_maxIter = 50, nCores = 1, seed = 11
  )
  cases <- list(
    single = list(wrapper = list(), upstream = list()),
    joint = list(
      wrapper = list(gKO = c("g1", "g2"), dr_direction = FALSE),
      upstream = list(gKO = c("g1", "g2"), dr_direction = FALSE)
    ),
    qc = list(
      wrapper = list(
        qc = TRUE, qc_min_library_size = 0, qc_min_pct = 0.2,
        qc_remove_outlier_cells = FALSE, qc_mt_threshold = 0.4
      ),
      upstream = list(
        qc = TRUE, qc_minLibSize = 0, qc_minPCT = 0.2,
        qc_removeOutlierCells = FALSE, qc_maxMTratio = 0.4
      )
    ),
    heat = list(
      wrapper = list(
        gKO = c("g1", "g2"), ma_method = "heat",
        dr_directionRegressLibSize = TRUE
      ),
      upstream = list(
        gKO = c("g1", "g2"), ma_method = "heat",
        dr_directionRegressLibSize = TRUE
      )
    ),
    screen = list(
      wrapper = list(gKO = c("g1", "g2", "g3"), transcriptomeWide = TRUE),
      upstream = list(gKO = c("g1", "g2", "g3"), transcriptomeWide = TRUE)
    ),
    all_genes = list(
      wrapper = list(gKO = NULL, transcriptomeWide = TRUE, dr_direction = FALSE),
      upstream = list(gKO = NULL, transcriptomeWide = TRUE, dr_direction = FALSE)
    ),
    caller_rng = list(wrapper = list(seed = NULL), upstream = list(seed = NULL))
  )
  for (case in cases) {
    set.seed(123)
    out <- do.call(RunscTenifoldKnk, utils::modifyList(wrapper_args, case$wrapper, keep.null = TRUE))
    wrapper_rng <- .Random.seed
    set.seed(123)
    expected <- suppressMessages(do.call(
      scTenifoldKnk::scTenifoldKnk,
      utils::modifyList(upstream_args, case$upstream, keep.null = TRUE)
    ))
    expect_identical(wrapper_rng, .Random.seed)
    expect_equal(out@tools$scTenifoldKnk$result, expected, tolerance = 1e-12)
    expect_identical(out@tools$scTenifoldKnk$diffRegulation, expected$diffRegulation)
  }
})

test_that("RunscTenifoldKnk preserves feature selection and upstream QC errors", {
  skip_if_not_installed("scTenifoldKnk", minimum_version = "2.0.0")
  srt <- make_sctenifold_knk_srt()
  testthat::local_mocked_bindings(check_r = function(...) invisible(TRUE), .package = "scop")
  out <- RunscTenifoldKnk(
    srt, "g1",
    features = paste0("g", 2:30), qc = FALSE,
    nc_nNet = 3, nc_nCells = 40, td_maxIter = 50, verbose = FALSE
  )
  counts <- GetAssayData5(srt, assay = "RNA", layer = "counts")[c(paste0("g", 2:30), "g1"), ]
  expected <- suppressMessages(scTenifoldKnk::scTenifoldKnk(
    counts, "g1",
    qc = FALSE, nc_nNet = 3, nc_nCells = 40,
    td_maxIter = 50, nCores = 1
  ))
  expect_equal(out@tools$scTenifoldKnk$result, expected, tolerance = 1e-12)
  counts["g1", ] <- 0
  filtered <- Seurat::CreateSeuratObject(counts = counts)
  expect_error(
    RunscTenifoldKnk(filtered, "g1",
      qc_min_library_size = 0,
      qc_remove_outlier_cells = FALSE, verbose = FALSE
    ),
    "after quality control"
  )
})

test_that("scTenifoldKnk 2.0 outputs retain all existing single-knockout plots", {
  skip_if_not_installed("scTenifoldKnk", minimum_version = "2.0.0")
  srt <- make_sctenifold_knk_srt()
  testthat::local_mocked_bindings(check_r = function(...) invisible(TRUE), .package = "scop")
  out <- RunscTenifoldKnk(
    srt, "g1",
    qc = FALSE, nc_nNet = 3, nc_nCells = 40,
    td_maxIter = 50, verbose = FALSE
  )
  expect_true(all(c("direction", "directionScore") %in% names(out@tools$scTenifoldKnk$diffRegulation)))
  expect_named(out@tools$scTenifoldKnk$result$tensorNetworks, "WT")
  for (type in c("qq", "effect", "network", "manifold")) {
    p <- scTenifoldKnkPlot(out, plot_type = type, top_n = 8, label = FALSE)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  out@tools$scTenifoldKnk$diffRegulation$direction <- NULL
  out@tools$scTenifoldKnk$diffRegulation$directionScore <- NULL
  expect_s3_class(scTenifoldKnkPlot(out, plot_type = "effect"), "ggplot")
  out@tools$scTenifoldKnk$parameters$transcriptomeWide <- TRUE
  expect_error(scTenifoldKnkPlot(out), "Transcriptome-wide")
  out@tools$scTenifoldKnk$parameters$transcriptomeWide <- FALSE
  out@tools$scTenifoldKnk$parameters$ma_method <- "heat"
  expect_error(scTenifoldKnkPlot(out, plot_type = "manifold"), "Heat alignment")
  out@tools$scTenifoldKnk$parameters$ma_method <- "manifold"
  out@tools$scTenifoldKnk$result$manifoldAlignment <- NULL
  expect_error(scTenifoldKnkPlot(out, plot_type = "manifold"), "store_manifold")
})
