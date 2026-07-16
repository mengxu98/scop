make_bayesspace_object <- function() {
  counts <- matrix(
    c(
      3, 1, 0, 2,
      0, 4, 1, 0,
      2, 1, 3, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$col <- c(0, 2, 1, 3)
  srt$row <- c(0, 0, 1, 1)
  srt
}

make_bayesspace_multi_image_object <- function() {
  srt <- make_bayesspace_object()
  assay <- SeuratObject::DefaultAssay(srt)
  slice1 <- data.frame(
    x = c(0, 2), y = c(0, 0),
    row.names = c("Spot1", "Spot2")
  )
  slice2 <- data.frame(
    x = c(1, 3), y = c(1, 1),
    row.names = c("Spot3", "Spot4")
  )
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    slice1, type = "centroids", assay = assay, key = "bs1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    slice2, type = "centroids", assay = assay, key = "bs2_"
  )
  srt
}

with_mock_bayesspace <- function(cluster_fun, code) {
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "BayesSpace")
      if (identical(name, "spatialCluster")) {
        return(cluster_fun)
      }
      stop("unexpected BayesSpace function: ", name)
    },
    .package = "scop"
  )
  force(code)
}

mock_bayesspace_complete <- function(sce, ...) {
  cdata <- as.data.frame(SummarizedExperiment::colData(sce))
  cdata$spatial.cluster <- paste0("domain_", seq_len(nrow(cdata)) %% 2 + 1)
  cdata$cluster.init <- paste0("init_", seq_len(nrow(cdata)) %% 2 + 1)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cdata)
  sce
}

test_that("RunBayesSpace stores raw coordinate provenance and aligned cells", {
  srt <- make_bayesspace_object()
  with_mock_bayesspace(mock_bayesspace_complete, {
    out <- RunBayesSpace(
      srt,
      q = 2,
      preprocess = FALSE,
      store_sce = FALSE,
      verbose = FALSE
    )
  })

  result <- out@tools[["BayesSpace"]]
  expect_identical(result$schema_version, 1L)
  expect_identical(result$source$coordinate_space, "raw")
  expect_identical(result$source$coord.cols, c("col", "row"))
  expect_identical(result$source$image, NA_character_)
  expect_identical(result$source$selection_strategy, "metadata_columns")
  expect_identical(result$cells, colnames(out))
  expect_identical(result$coords$cell_id, colnames(out))
  expect_identical(rownames(result$colData), colnames(out))
  expect_false("sce" %in% names(result))
  expect_false(anyNA(out$BayesSpace_cluster))
})

test_that("RunBayesSpace rejects ambiguous or partial image selection atomically", {
  srt <- make_bayesspace_multi_image_object()
  metadata_before <- srt[[]]
  tools_before <- srt@tools

  with_mock_bayesspace(mock_bayesspace_complete, {
    expect_error(
      RunBayesSpace(
        srt,
        q = 2,
        preprocess = FALSE,
        store_sce = FALSE,
        verbose = FALSE
      ),
      "Multiple spatial images"
    )
    expect_error(
      RunBayesSpace(
        srt,
        q = 2,
        image = "slice1",
        preprocess = FALSE,
        store_sce = FALSE,
        verbose = FALSE
      ),
      "exactly one row for every object\\s+spot"
    )
  })

  expect_identical(srt[[]], metadata_before)
  expect_identical(srt@tools, tools_before)
})

test_that("RunBayesSpace validates backend spot alignment before mutation", {
  srt <- make_bayesspace_object()
  metadata_before <- srt[[]]
  tools_before <- srt@tools
  incomplete_backend <- function(sce, ...) {
    sce <- sce[, -1, drop = FALSE]
    mock_bayesspace_complete(sce)
  }

  with_mock_bayesspace(incomplete_backend, {
    expect_error(
      RunBayesSpace(
        srt,
        q = 2,
        preprocess = FALSE,
        store_sce = FALSE,
        verbose = FALSE
      ),
      "exactly one row for every input spot"
    )
  })

  expect_identical(srt[[]], metadata_before)
  expect_identical(srt@tools, tools_before)
})

test_that("RunBayesSpace realigns reordered backend output by spot ID", {
  srt <- make_bayesspace_object()
  reordered_backend <- function(sce, ...) {
    cdata <- as.data.frame(SummarizedExperiment::colData(sce))
    cdata$spatial.cluster <- paste0("domain_", rownames(cdata))
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cdata)
    sce[, rev(seq_len(ncol(sce))), drop = FALSE]
  }

  with_mock_bayesspace(reordered_backend, {
    out <- RunBayesSpace(
      srt,
      q = 2,
      preprocess = FALSE,
      store_sce = FALSE,
      verbose = FALSE
    )
  })

  expect_identical(
    as.character(out$BayesSpace_cluster),
    paste0("domain_", colnames(out))
  )
  expect_identical(
    rownames(out@tools[["BayesSpace"]]$colData),
    colnames(out)
  )
})

test_that("RunBayesSpace rejects missing cluster labels before mutation", {
  srt <- make_bayesspace_object()
  metadata_before <- srt[[]]
  empty_label_backend <- function(sce, ...) {
    cdata <- as.data.frame(SummarizedExperiment::colData(sce))
    cdata$spatial.cluster <- c("domain_1", NA, "domain_2", "")
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cdata)
    sce
  }

  with_mock_bayesspace(empty_label_backend, {
    expect_error(
      RunBayesSpace(
        srt,
        q = 2,
        preprocess = FALSE,
        store_sce = FALSE,
        verbose = FALSE
      ),
      "missing or empty"
    )
  })
  expect_identical(srt[[]], metadata_before)
})
