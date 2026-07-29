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
    }
  )
  force(code)
  list(received = received, checked = checked)
}

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
    RunCHOIR(srt, max_clusters = 2L, verbose = FALSE),
    "may not terminate"
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
    paste0(
      "corceslab/CHOIR@",
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    )
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
    "final cluster column"
  )
})

test_that("CHOIR helper resolves the registered optional repository", {
  checked <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, dependencies, verbose, ...) {
      checked <<- list(
        packages = packages,
        dependencies = dependencies,
        verbose = verbose
      )
      invisible(TRUE)
    },
    choir_installed_commit = function() {
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    }
  )
  expect_invisible(choir_check_r(verbose = FALSE))
  expect_identical(
    checked$packages,
    paste0(
      "corceslab/CHOIR@",
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    )
  )
  expect_identical(checked$dependencies, NA)
  expect_identical(checked$verbose, FALSE)
})

test_that("CHOIR helper skips the upstream telemetry configure on fallback", {
  calls <- character()
  fallback_calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      calls <<- c(calls, packages)
      if (startsWith(packages, "corceslab/CHOIR@")) {
        return(list(CHOIR = FALSE))
      }
      list(CHOIR = TRUE)
    },
    choir_install_without_configure = function(verbose) {
      fallback_calls <<- fallback_calls + 1L
      expect_false(verbose)
      invisible(TRUE)
    },
    choir_installed_commit = function() {
      "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
    }
  )

  expect_invisible(getFromNamespace("choir_check_r", "scop")(verbose = FALSE))
  expect_identical(
    calls,
    c(
      paste0(
        "corceslab/CHOIR@",
        "e9ebfbc9089beeaf4ca088c7b81b18f39758b0bc"
      ),
      "CHOIR"
    )
  )
  expect_identical(fallback_calls, 1L)
})
