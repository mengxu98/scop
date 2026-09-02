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
    slice1,
    type = "centroids", assay = assay, key = "bs1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    slice2,
    type = "centroids", assay = assay, key = "bs2_"
  )
  srt
}

make_bayesspace_visium_v1_object <- function() {
  srt <- make_bayesspace_object()
  scale_factors <- structure(
    list(spot = 1, fiducial = 1, hires = 0.5, lowres = 0.25),
    class = "scalefactors"
  )
  coordinates <- data.frame(
    tissue = 1L,
    row = c(10L, 10L, 11L, 11L),
    col = c(4L, 6L, 5L, 7L),
    imagerow = c(101.3, 99.8, 213.7, 210.2),
    imagecol = c(52.1, 180.4, 61.9, 191.6),
    row.names = colnames(srt)
  )
  srt[["slice"]] <- methods::new(
    "VisiumV1",
    image = array(1, dim = c(120, 140, 3)),
    scale.factors = scale_factors,
    coordinates = coordinates,
    spot.radius = 0.1,
    assay = SeuratObject::DefaultAssay(srt),
    key = "bsv1_"
  )
  list(object = srt, coordinates = coordinates)
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

bayesspace_receipt_plain <- function(messages) {
  cli::ansi_strip(paste(messages, collapse = "\n"))
}

bayesspace_expect_receipt_plot <- function(messages, expected_call, srt) {
  plain <- bayesspace_receipt_plain(messages)
  expect_true(grepl(paste0("Plot `", expected_call, "`"), plain, fixed = TRUE))
  expect_error(parse(text = expected_call), NA)
  plot <- suppressWarnings(eval(
    parse(text = expected_call),
    envir = list(srt = srt),
    enclos = parent.frame()
  ))
  expect_s3_class(plot, "ggplot")
  invisible(plain)
}

test_that("RunBayesSpace stores raw coordinate provenance and aligned cells", {
  skip_if_not_installed("SingleCellExperiment")
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
  expect_identical(result$parameters$coordinate_space, "raw")
  expect_identical(result$parameters$coord.cols, c("col", "row"))
  expect_identical(result$parameters$image, NA_character_)
  expect_identical(result$cells, colnames(out))
  expect_identical(result$coords$cell_id, colnames(out))
  expect_identical(rownames(result$colData), colnames(out))
  expect_false("sce" %in% names(result))
  expect_false(anyNA(out$BayesSpace_cluster))
})

test_that("RunBayesSpace emits an exact receipt and safely quotes its Plot key", {
  skip_if_not_installed("SingleCellExperiment")
  srt <- make_bayesspace_object()
  cluster_colname <- 'Bayes "quoted"\\domain'
  out <- NULL
  messages <- suppressWarnings(testthat::capture_messages(
    out <- with_mock_bayesspace(mock_bayesspace_complete, {
      RunBayesSpace(
        srt,
        q = 2,
        preprocess = FALSE,
        cluster_colname = cluster_colname,
        store_sce = FALSE,
        verbose = TRUE
      )
    })
  ))
  plain <- bayesspace_receipt_plain(messages)
  expect_true(grepl(
    "BayesSpace completed: 4 spots, 2 domains, domain size 2-2",
    plain,
    fixed = TRUE
  ))
  expect_true(grepl("Scope assay", plain, fixed = TRUE))
  expect_true(grepl("raw coordinates", plain, fixed = TRUE))
  expect_true(grepl("Saved metadata column", plain, fixed = TRUE))
  expect_true(grepl(cluster_colname, plain, fixed = TRUE))
  expect_true(grepl("srt@tools[['BayesSpace']]", plain, fixed = TRUE))
  expect_true(cluster_colname %in% colnames(out[[]]))
  plot_call <- paste0(
    "SpatialSpotPlot(srt, group.by = ",
    deparse1(cluster_colname, width.cutoff = 500L),
    ")"
  )
  bayesspace_expect_receipt_plot(messages, plot_call, out)
})

test_that("RunBayesSpace is silent on request and failure has no receipt", {
  skip_if_not_installed("SingleCellExperiment")
  srt <- make_bayesspace_object()
  quiet <- NULL
  quiet_messages <- suppressWarnings(testthat::capture_messages(
    quiet <- with_mock_bayesspace(mock_bayesspace_complete, {
      RunBayesSpace(
        srt,
        q = 2,
        preprocess = FALSE,
        store_sce = FALSE,
        verbose = FALSE
      )
    })
  ))
  expect_length(quiet_messages, 0L)
  expect_s4_class(quiet, "Seurat")

  missing_cluster <- function(sce, ...) sce
  failure_messages <- suppressWarnings(testthat::capture_messages(
    with_mock_bayesspace(missing_cluster, {
      expect_error(
        RunBayesSpace(
          srt,
          q = 2,
          preprocess = FALSE,
          store_sce = FALSE,
          verbose = TRUE
        ),
        "did not return"
      )
    })
  ))
  expect_false(grepl(
    "BayesSpace completed",
    bayesspace_receipt_plain(failure_messages),
    fixed = TRUE
  ))
})

test_that("BayesSpace preserves exact Visium array indices separately from raw pixels", {
  skip_if_not_installed("SingleCellExperiment")
  fixture <- make_bayesspace_visium_v1_object()
  sce <- Seurat::as.SingleCellExperiment(fixture$object)
  input <- bayesspace_add_spatial_coords(
    srt = fixture$object,
    sce = sce,
    image = "slice",
    platform = "Visium"
  )
  cdata <- as.data.frame(SummarizedExperiment::colData(input$sce))
  expect_identical(cdata$array_row, fixture$coordinates$row)
  expect_identical(cdata$array_col, fixture$coordinates$col)
  expect_equal(input$coords$x, fixture$coordinates$imagecol)
  expect_equal(input$coords$y, fixture$coordinates$imagerow)
  expect_identical(input$source$backend_coordinate_space, "array_index")
})

test_that("RunBayesSpace rejects ambiguous or partial image selection atomically", {
  skip_if_not_installed("SpatialExperiment")
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
  skip_if_not_installed("SpatialExperiment")
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
  skip_if_not_installed("SingleCellExperiment")
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
  skip_if_not_installed("SpatialExperiment")
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
