make_choir_seurat <- function() {
  counts <- matrix(
    c(
      8, 7, 1, 0, 1, 0,
      5, 4, 0, 1, 0, 1,
      0, 1, 7, 8, 1, 0,
      1, 0, 5, 6, 0, 1,
      0, 0, 1, 0, 8, 7,
      1, 0, 0, 1, 6, 5
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$batch <- rep(c("B1", "B2"), each = 3)
  Seurat::NormalizeData(srt, verbose = FALSE)
}

with_mock_choir <- function(code, backend = NULL) {
  received <- new.env(parent = emptyenv())
  if (is.null(backend)) {
    backend <- function(
      object,
      key,
      alpha,
      p_adjust,
      feature_set,
      exclude_features,
      n_iterations,
      n_trees,
      min_accuracy,
      max_clusters,
      normalization_method,
      batch_correction_method,
      batch_labels,
      use_assay,
      use_slot,
      reduction,
      var_features,
      atac,
      n_cores,
      random_seed,
      distance_awareness,
      verbose
    ) {
      received$args <- as.list(environment())
      cluster_col <- paste0("CHOIR_clusters_", alpha)
      object@meta.data[[cluster_col]] <- rep(1:3, each = 2)
      object@misc[[key]] <- list(
        clusters = stats::setNames(list(object@meta.data[[cluster_col]]), cluster_col),
        records = list(mock = TRUE)
      )
      object
    }
  }

  checked <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, dependencies, verbose, ...) {
      checked <<- c(checked, packages)
      expect_identical(dependencies, NA)
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "CHOIR")
      expect_identical(name, "CHOIR")
      backend
    },
    choir_installed_commit = function() {
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    },
    choir_namespace_loaded = function() FALSE,
    choir_mark_namespace_trusted = function() invisible(TRUE),
    choir_validate_platform = function(...) invisible(TRUE)
  )
  force(code)
  list(received = received, checked = checked)
}

test_that("CHOIR rejects unsafe numeric cluster limits", {
  validate <- getFromNamespace("choir_validate_max_clusters", "scop")
  expect_identical(validate("auto"), "auto")
  expect_error(validate(2L), "may not terminate")
})

test_that("RunCHOIR validates inputs before checking the backend", {
  expect_error(
    RunCHOIR(matrix(1, nrow = 2, ncol = 2), verbose = FALSE),
    "Seurat"
  )

  srt <- make_choir_seurat()
  expect_error(
    RunCHOIR(srt, assay = "missing", verbose = FALSE),
    "not present"
  )
  expect_error(
    RunCHOIR(
      srt,
      normalization_method = "SCTransform",
      verbose = FALSE
    ),
    "requires.*counts"
  )
  expect_error(
    RunCHOIR(srt, batch.by = "missing", verbose = FALSE),
    "cell metadata"
  )
  expect_error(
    RunCHOIR(srt, seed = 0, verbose = FALSE),
    "positive integer"
  )
  expect_error(
    RunCHOIR(
      srt,
      layer = "counts",
      normalization_method = "SCTransform",
      verbose = FALSE
    ),
    "does not support.*Assay5"
  )
})

test_that("RunCHOIR calls the optional backend and standardizes results", {
  srt <- make_choir_seurat()
  mocked <- with_mock_choir({
    out <- RunCHOIR(
      srt,
      batch.by = "batch",
      n_iterations = 20,
      n_trees = 10,
      n_cores = 2,
      distance_awareness = 1.8,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_identical(
    mocked$checked,
    "CHOIR"
  )
  expect_identical(mocked$received$args$use_assay, "RNA")
  expect_identical(mocked$received$args$use_slot, "data")
  expect_identical(mocked$received$args$batch_labels, "batch")
  expect_identical(mocked$received$args$batch_correction_method, "Harmony")
  expect_identical(mocked$received$args$n_cores, 2L)
  expect_equal(mocked$received$args$distance_awareness, 1.8)

  expect_true("CHOIR_clusters_0.05" %in% colnames(out@meta.data))
  expect_true("CHOIR_cluster" %in% colnames(out@meta.data))
  expect_s3_class(out$CHOIR_cluster, "factor")
  expect_equal(
    as.character(out$CHOIR_cluster),
    rep(as.character(1:3), each = 2)
  )
  expect_true("CHOIR" %in% names(out@misc))
  expect_true("CHOIR" %in% names(out@tools))
  expect_identical(out@tools$CHOIR$backend$repository, "corceslab/CHOIR")
  expect_identical(out@tools$CHOIR$backend$key, "CHOIR")
  expect_identical(out@tools$CHOIR$parameters$batch.by, "batch")
  expect_equal(rownames(out@tools$CHOIR$clusters), colnames(out))
})

test_that("RunCHOIR uses an existing reduction and variable features", {
  srt <- make_choir_seurat()
  SeuratObject::VariableFeatures(srt) <- rownames(srt)[1:4]
  embeddings <- matrix(
    seq_len(18) / 10,
    nrow = 6,
    dimnames = list(colnames(srt), paste0("PC_", 1:3))
  )
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "PC_",
    assay = "RNA"
  )

  mocked <- with_mock_choir({
    out <- RunCHOIR(
      srt,
      reduction = "pca",
      distance_awareness = 2,
      verbose = FALSE
    )
  })

  expect_equal(
    mocked$received$args$reduction,
    SeuratObject::Embeddings(srt[["pca"]])
  )
  expect_setequal(
    mocked$received$args$var_features,
    SeuratObject::VariableFeatures(srt)
  )
  expect_identical(out@tools$CHOIR$parameters$reduction, "pca")
})

test_that("RunCHOIR rejects unsupported pass-through arguments", {
  srt <- make_choir_seurat()
  expect_error(
    with_mock_choir({
      RunCHOIR(
        srt,
        unsupported_parameter = 1,
        verbose = FALSE
      )
    }),
    "Unsupported CHOIR argument"
  )
})

test_that("RunCHOIR rejects invalid backend results", {
  srt <- make_choir_seurat()
  missing_clusters <- function(object, key, use_assay, use_slot, verbose) {
    object
  }
  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }, backend = missing_clusters),
    "final cluster column"
  )

  changed_cells <- function(object, key, alpha, use_assay, use_slot, verbose) {
    object <- object[, colnames(object)[-1]]
    object@meta.data[[paste0("CHOIR_clusters_", alpha)]] <- 1:ncol(object)
    object
  }
  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }, backend = changed_cells),
    "changed cell identities or cell order"
  )
})

test_that("RunCHOIR does not accept a stale cluster column", {
  srt <- make_choir_seurat()
  srt$CHOIR_clusters_0.05 <- rep(1:3, each = 2)
  unchanged <- function(object, key, alpha, use_assay, use_slot, verbose) {
    object
  }

  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }, backend = unchanged),
    "Existing CHOIR output"
  )
  expect_error(
    with_mock_choir({
      RunCHOIR(srt, overwrite = TRUE, verbose = FALSE)
    }, backend = unchanged),
    "final cluster column"
  )
})

test_that("RunCHOIR protects and explicitly replaces managed output", {
  srt <- make_choir_seurat()
  srt$CHOIR_clusters_0.05 <- "old internal"
  srt$CHOIR_cluster <- "old stable"
  embeddings <- matrix(
    seq_len(12) / 10,
    nrow = 6,
    dimnames = list(colnames(srt), paste0("CHOIRP0_", 1:2))
  )
  srt[["CHOIR_P0_reduction"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "CHOIRP0_",
    assay = "RNA"
  )
  srt@misc$CHOIR <- list(user = "keep")
  srt@tools$CHOIR <- list(user = "keep")
  before <- srt

  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }),
    "Set.*overwrite"
  )
  expect_identical(srt@meta.data, before@meta.data)
  expect_identical(srt@reductions, before@reductions)
  expect_identical(srt@misc, before@misc)
  expect_identical(srt@tools, before@tools)

  mocked <- with_mock_choir({
    out <- RunCHOIR(
      srt,
      store_tool = FALSE,
      overwrite = TRUE,
      verbose = FALSE
    )
  })
  expect_true("CHOIR_cluster" %in% colnames(out@meta.data))
  expect_false("CHOIR_P0_reduction" %in% names(out@reductions))
  expect_identical(out@misc$CHOIR$records, list(mock = TRUE))
  expect_identical(out@tools$CHOIR, list(user = "keep"))
})

test_that("RunCHOIR preserves tools when storage is disabled", {
  srt <- make_choir_seurat()
  srt@tools$CHOIR <- list(user = "keep")

  with_mock_choir({
    out <- RunCHOIR(
      srt,
      store_tool = FALSE,
      verbose = FALSE
    )
  })

  expect_identical(out@tools$CHOIR, list(user = "keep"))
})

test_that("RunCHOIR preserves a stable output column used for batches", {
  srt <- make_choir_seurat()
  srt$CHOIR_cluster <- srt$batch
  backend <- function(
    object,
    key,
    alpha,
    batch_labels,
    use_assay,
    use_slot,
    verbose
  ) {
    expect_true(batch_labels %in% colnames(object@meta.data))
    object@meta.data[[paste0("CHOIR_clusters_", alpha)]] <- rep(1:3, each = 2)
    object@misc[[key]] <- list(
      clusters = list(mock = TRUE),
      records = list(mock = TRUE)
    )
    object
  }

  with_mock_choir({
    out <- RunCHOIR(
      srt,
      batch.by = "CHOIR_cluster",
      cluster_colname = "CHOIR_cluster",
      overwrite = TRUE,
      verbose = FALSE
    )
  }, backend = backend)
  expect_equal(as.character(out$CHOIR_cluster), rep(as.character(1:3), each = 2))

  srt$CHOIR_clusters_0.05 <- srt$batch
  expect_error(
    with_mock_choir({
      RunCHOIR(
        srt,
        batch.by = "CHOIR_clusters_0.05",
        overwrite = TRUE,
        verbose = FALSE
      )
    }),
    "cannot use CHOIR's managed cluster column"
  )
})

test_that("RunCHOIR requires complete upstream records", {
  srt <- make_choir_seurat()
  missing_records <- function(
    object,
    key,
    alpha,
    use_assay,
    use_slot,
    verbose
  ) {
    object@meta.data[[paste0("CHOIR_clusters_", alpha)]] <- seq_len(ncol(object))
    object
  }

  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }, backend = missing_records),
    "complete records"
  )
})

test_that("RunCHOIR rejects invalid cluster labels", {
  srt <- make_choir_seurat()
  invalid_labels <- function(
    object,
    key,
    alpha,
    use_assay,
    use_slot,
    verbose
  ) {
    cluster_col <- paste0("CHOIR_clusters_", alpha)
    object@meta.data[[cluster_col]] <- c(1, 2, 3, 4, 5, Inf)
    object@misc[[key]] <- list(
      clusters = list(object@meta.data[[cluster_col]]),
      records = list(mock = TRUE)
    )
    object
  }

  expect_error(
    with_mock_choir({
      RunCHOIR(srt, verbose = FALSE)
    }, backend = invalid_labels),
    "invalid final cluster labels"
  )
})

test_that("CHOIR rejects unsupported Windows execution", {
  validate <- getFromNamespace("choir_validate_platform", "scop")
  expect_invisible(validate("unix"))
  expect_error(validate("windows"), "not Windows")
})

test_that("CHOIR helper reuses the installed pinned backend", {
  checked <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, dependencies, install, verbose, ...) {
      checked <<- list(
        packages = packages,
        dependencies = dependencies,
        install = install,
        verbose = verbose
      )
      invisible(TRUE)
    },
    choir_installed_commit = function() {
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    },
    choir_namespace_loaded = function() FALSE
  )
  expect_invisible(choir_check_r(verbose = FALSE))
  expect_identical(checked$packages, "CHOIR")
  expect_identical(checked$dependencies, NA)
  expect_false(checked$install)
  expect_identical(checked$verbose, FALSE)
})

test_that("CHOIR helper always skips the upstream telemetry configure", {
  calls <- list()
  fallback_calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      calls[[length(calls) + 1L]] <<- packages
      if (identical(packages, "CHOIR") && length(calls) == 1L) {
        return(list(CHOIR = FALSE))
      }
      if (!identical(packages, "CHOIR")) {
        return(as.list(stats::setNames(rep(TRUE, length(packages)), packages)))
      }
      list(CHOIR = TRUE)
    },
    choir_install_without_configure = function(verbose) {
      fallback_calls <<- fallback_calls + 1L
      expect_false(verbose)
      invisible(TRUE)
    },
    choir_namespace_loaded = function() FALSE,
    choir_installed_commit = function() {
      if (fallback_calls == 0L) {
        return(NULL)
      }
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    }
  )

  expect_invisible(getFromNamespace("choir_check_r", "scop")(verbose = FALSE))
  expect_identical(calls[[1L]], "CHOIR")
  expect_identical(
    calls[[2L]],
    getFromNamespace(".choir_dependencies", "scop")
  )
  expect_identical(calls[[3L]], "CHOIR")
  expect_identical(fallback_calls, 1L)
})

test_that("CHOIR helper stops before source installation when dependencies fail", {
  fallback_calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      if (identical(packages, "CHOIR")) {
        return(list(CHOIR = FALSE))
      }
      as.list(stats::setNames(
        c(FALSE, rep(TRUE, length(packages) - 1L)),
        packages
      ))
    },
    choir_namespace_loaded = function() FALSE,
    choir_installed_commit = function() NULL,
    choir_install_without_configure = function(...) {
      fallback_calls <<- fallback_calls + 1L
    }
  )

  expect_error(
    getFromNamespace("choir_check_r", "scop")(verbose = FALSE),
    "prepare dependencies"
  )
  expect_identical(fallback_calls, 0L)
})

test_that("CHOIR helper refuses to replace a loaded stale namespace", {
  fallback_calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) list(CHOIR = TRUE),
    choir_installed_commit = function() "stale",
    choir_namespace_loaded = function() TRUE,
    choir_loaded_commit = function() "stale",
    choir_install_without_configure = function(...) {
      fallback_calls <<- fallback_calls + 1L
    }
  )

  expect_error(
    getFromNamespace("choir_check_r", "scop")(verbose = FALSE),
    "Restart R"
  )
  expect_identical(fallback_calls, 0L)
})

test_that("CHOIR helper accepts a loaded pinned namespace on the fast path", {
  checked <- FALSE
  marked <- character()
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(...) {
      checked <<- TRUE
      list(CHOIR = TRUE)
    },
    choir_namespace_loaded = function() TRUE,
    choir_loaded_commit = function() {
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    },
    choir_mark_backend_verified = function(commit) {
      marked <<- c(marked, commit)
    },
    choir_mark_namespace_trusted = function() {
      marked <<- c(marked, "namespace")
    }
  )

  expect_invisible(getFromNamespace("choir_check_r", "scop")(verbose = FALSE))
  expect_false(checked)
  expect_identical(
    marked,
    c(
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc",
      "namespace"
    )
  )
})

test_that("CHOIR helper rejects Windows before checking dependencies", {
  checked <- FALSE
  testthat::local_mocked_bindings(
    .package = "scop",
    choir_validate_platform = function(...) {
      stop("windows blocked")
    },
    check_r = function(...) {
      checked <<- TRUE
      TRUE
    }
  )

  expect_error(
    getFromNamespace("choir_check_r", "scop")(verbose = FALSE),
    "windows blocked"
  )
  expect_false(checked)
})
