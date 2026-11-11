rctd_nonzero_shared_features <- function(st_counts, ref_counts) {
  st_counts <- methods::as(st_counts, "dgCMatrix")
  ref_counts <- methods::as(ref_counts, "dgCMatrix")
  quality <- rctd_sparse_quality_cpp(st_counts, ref_counts)
  rownames(st_counts)[quality$keep_features]
}
make_rctd_validation_pair <- function() {
  counts <- matrix(
    c(4, 1, 0, 3, 2, 1),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:2))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  spatial <- SeuratObject::CreateSeuratObject(counts)
  reference <- SeuratObject::CreateSeuratObject(counts)
  reference$celltype <- c("Alpha", "Beta")
  list(spatial = spatial, reference = reference)
}

test_that("RunRCTD validates result names and storage before backend work", {
  pair <- make_rctd_validation_pair()
  backend_touched <- FALSE
  testthat::local_mocked_bindings(
    check_r = function(...) {
      backend_touched <<- TRUE
      stop("backend check should not run")
    },
    rctd_run_spacexr = function(...) {
      backend_touched <<- TRUE
      stop("backend should not run")
    },
    .package = "scop"
  )
  invalid <- list(
    prefix = list("", 1L),
    tool_name = list("", 1L),
    store_results = list(1L, NA)
  )
  base_args <- list(
    object = pair$spatial,
    reference = pair$reference,
    reference_label = "celltype",
    verbose = FALSE
  )
  for (arg in names(invalid)) {
    pattern <- if (identical(arg, "store_results")) {
      "store_results.*TRUE or FALSE"
    } else {
      paste0(arg, ".*single non-empty string")
    }
    for (value in invalid[[arg]]) {
      args <- base_args
      args[[arg]] <- value
      expect_error(do.call(RunRCTD, args), pattern)
    }
  }
  expect_false(backend_touched)
})

test_that("RCTD sparse quality helper matches R sparse summaries", {
  st_dense <- matrix(
    c(
      1, 0, 2,
      0, 0, 0,
      3, 0, 0,
      0, 5, 0,
      0, 0, 0
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("Spot", 1:3))
  )
  ref_dense <- matrix(
    c(
      0, 4, 0, 1,
      2, 0, 0, 0,
      0, 0, 0, 0,
      1, 1, 0, 0,
      0, 0, 7, 0
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:4))
  )
  st <- methods::as(Matrix::Matrix(st_dense, sparse = TRUE), "dgCMatrix")
  ref <- methods::as(Matrix::Matrix(ref_dense, sparse = TRUE), "dgCMatrix")

  out <- rctd_sparse_quality_cpp(st, ref)
  keep <- Matrix::rowSums(st) > 0 & Matrix::rowSums(ref) > 0

  expect_equal(out$keep_features, unname(keep))
  expect_equal(
    out$st_numi,
    unname(Matrix::colSums(st[keep, , drop = FALSE]))
  )
  expect_equal(
    out$ref_numi,
    unname(Matrix::colSums(ref[keep, , drop = FALSE]))
  )
  expect_equal(rctd_nonzero_shared_features(st, ref), rownames(st)[keep])
})

test_that("RCTD normalizes check_r list results without adding a backend helper", {
  requested <- character()
  testthat::local_mocked_bindings(
    check_r = function(pkgs, verbose = FALSE) {
      requested <<- c(requested, pkgs)
      stats::setNames(as.list(rep(TRUE, length(pkgs))), sub("^.*/", "", pkgs))
    },
    .package = "scop"
  )
  expect_invisible(rctd_require_namespaces(c("spacexr", "S4Vectors")))
  expect_identical(requested, c("dmcable/spacexr", "S4Vectors"))
  expect_invisible(rctd_require_namespaces(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors")))

  testthat::local_mocked_bindings(
    check_r = function(pkgs, verbose = FALSE) list(unexpected = TRUE),
    .package = "scop"
  )
  expect_error(
    rctd_require_namespaces("spacexr"),
    "invalid result"
  )
})

test_that("spacexr API variants prefer new and accept complete legacy exports", {
  new_api <- c("createRctd", "runRctd")
  legacy_api <- c("SpatialRNA", "Reference", "create.RCTD", "run.RCTD")
  args <- list(
    st_counts = matrix(1, 2, 2), coords = data.frame(x = 1:2, y = 1:2),
    st_numi = 1:2, ref_counts = matrix(1, 2, 2), ref_labels = c("A", "B"),
    ref_numi = 1:2, rctd_mode = "doublet", max_cores = 1,
    create_rctd_params = list(), run_rctd_params = list()
  )
  testthat::local_mocked_bindings(
    .package = "scop",
    rctd_require_namespaces = function(...) invisible(TRUE),
    rctd_run_spacexr_new = function(...) "new_api",
    rctd_run_spacexr_old = function(...) "old_api",
    rctd_spacexr_exports = function() new_api
  )
  expect_identical(do.call(rctd_run_spacexr, args), "new_api")

  testthat::local_mocked_bindings(
    .package = "scop",
    rctd_require_namespaces = function(...) invisible(TRUE),
    rctd_run_spacexr_new = function(...) stop("new API must not run"),
    rctd_run_spacexr_old = function(...) "old_api",
    rctd_spacexr_exports = function() legacy_api
  )
  expect_identical(do.call(rctd_run_spacexr, args), "old_api")
})
