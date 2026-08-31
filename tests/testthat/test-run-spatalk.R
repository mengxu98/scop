make_spatalk_test_object <- function() {
  counts <- Matrix::sparseMatrix(
    i = rep(seq_len(4), 4), j = rep(seq_len(4), each = 4),
    x = rep(c(5, 3, 2, 1), 4), dims = c(4, 4),
    dimnames = list(c("L1", "R1", "TF1", "G1"), paste0("cell", 1:4))
  )
  srt <- SeuratObject::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "A", "B", "B")
  srt$col <- c(0, 1, 2, 3)
  srt$row <- c(0, 1, 1, 2)
  srt
}

mock_spatalk_bindings <- function(dec_seen = NULL, fail_cci = FALSE) {
  testthat::local_mocked_bindings(
    check_r = function(repos, verbose = FALSE, install = TRUE, ...) {
      expect_identical(repos, c("NNLM", "SpaTalk"))
      expect_false(install)
      list(NNLM = TRUE, SpaTalk = TRUE)
    },
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "SpaTalk")
      switch(fun,
        createSpaTalk = function(st_data, st_meta, species, if_st_is_sc, spot_max_cell, celltype = NULL) {
          list(
            data = list(rawdata = st_data), meta = list(rawmeta = st_meta),
            para = list(species = species, if_st_is_sc = if_st_is_sc),
            coef = matrix(), lrpair = data.frame(), tf = data.frame(),
            lr_path = list()
          )
        },
        dec_celltype = function(object, sc_data, sc_celltype, method = 1, dec_result = NULL, ...) {
          if (!is.null(dec_seen)) {
            dec_seen$value <- TRUE
            dec_seen$method <- method
            dec_seen$dec_result <- dec_result
          }
          object$coef <- dec_result %||% matrix(
            rep(c(0.7, 0.3), ncol(object$data$rawdata)),
            ncol = 2,
            byrow = TRUE,
            dimnames = list(colnames(object$data$rawdata), c("A", "B"))
          )
          object
        },
        find_lr_path = function(object, lrpairs, pathways, ...) {
          object$lr_path <- list(lrpairs = lrpairs, pathways = pathways)
          object
        },
        dec_cci_all = function(object, ...) {
          if (isTRUE(fail_cci)) stop("backend failed")
          object$lrpair <- data.frame(
            ligand = "L1", receptor = "R1", species = "Human",
            celltype_sender = "A", celltype_receiver = "B",
            lr_co_exp_num = 4, lr_co_ratio = 0.5,
            lr_co_ratio_pvalue = 0.02, score = 2,
            stringsAsFactors = FALSE
          )
          object$tf <- data.frame(
            celltype_sender = "A", celltype_receiver = "B",
            receptor = "R1", tf = "TF1", n_hop = 2, n_target = 3,
            score = 1.5, stringsAsFactors = FALSE
          )
          object
        },
        stop("unexpected SpaTalk function")
      )
    },
    spatalk_load_data = function(name, species) {
      if (identical(name, "lrpairs")) {
        return(data.frame(ligand = "L1", receptor = "R1", species = species))
      }
      data.frame(
        src = c("R1", "TF1"), dest = c("TF1", "G1"),
        pathway = c("P1", "P1"), species = species,
        src_tf = c("NO", "YES"), dest_tf = c("YES", "NO")
      )
    },
    spatalk_package_version = function() "1.0",
    .package = "scop",
    .env = parent.frame()
  )
}

test_that("SpaTalk preflight accepts installed packages without remote metadata", {
  seen <- new.env(parent = emptyenv())
  seen$calls <- list()
  testthat::local_mocked_bindings(
    check_r = function(repos, verbose = FALSE, install = TRUE, ...) {
      seen$calls[[length(seen$calls) + 1L]] <- list(repos = repos, install = install)
      list(NNLM = TRUE, SpaTalk = TRUE)
    },
    .package = "scop"
  )
  expect_invisible(spatalk_check_r())
  expect_length(seen$calls, 1L)
  expect_identical(seen$calls[[1L]]$repos, c("NNLM", "SpaTalk"))
  expect_false(seen$calls[[1L]]$install)
})

test_that("SpaTalk preflight falls back to official repositories when absent", {
  seen <- new.env(parent = emptyenv())
  seen$calls <- list()
  testthat::local_mocked_bindings(
    check_r = function(repos, verbose = FALSE, install = TRUE, ...) {
      seen$calls[[length(seen$calls) + 1L]] <- list(repos = repos, install = install)
      if (identical(repos, c("NNLM", "SpaTalk"))) {
        return(list(NNLM = FALSE, SpaTalk = FALSE))
      }
      list(NNLM = TRUE, SpaTalk = TRUE)
    },
    .package = "scop"
  )
  expect_invisible(spatalk_check_r())
  expect_length(seen$calls, 2L)
  expect_identical(
    seen$calls[[2L]]$repos,
    c("linxihui/NNLM", "ZJUFanLab/SpaTalk")
  )
})

test_that("RunSpaTalk stores official single-cell results and unified CCC rows", {
  srt <- make_spatalk_test_object()
  mock_spatalk_bindings()
  out <- RunSpaTalk(
    srt,
    group.by = "celltype", mode = "single_cell",
    coord.cols = c("col", "row"), store.object = "minimal",
    backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$SpaTalk$method, "SpaTalk")
  expect_identical(out@tools$SpaTalk$parameters$coordinate_space, "raw")
  expect_identical(out@tools$SpaTalk$long_table$pvalue, 0.02)
  expect_identical(out@tools$SpaTalk$long_table$score_type, "spatalk_score")
  expect_null(out@tools$SpaTalk$results$default$native_object)
  expect_identical(out@tools$SpaTalk$provenance$backend_version, "1.0")
  expect_true(is.list(out@tools$SpaTalk$results$default$spatial_plot))
  expect_identical(
    out@tools$SpaTalk$results$default$spatial_plot$labels,
    stats::setNames(as.character(srt$celltype), colnames(srt))
  )
  expect_true("SpaTalk" %in% out@tools$CCC$methods)
  expect_true(all(out@tools$CCC$long_table$method == "SpaTalk"))
  expect_s3_class(CCCNetworkPlot(out, method = "SpaTalk", plot_type = "spatial"), "ggplot")
  expect_s3_class(SpaTalkPlot(out, plot_type = "tf"), "ggplot")
})

test_that("RunSpaTalk spot mode runs NNLM deconvolution and can retain native output", {
  spatial <- make_spatalk_test_object()
  reference <- make_spatalk_test_object()
  reference$ref_type <- c("A", "A", "B", "B")
  seen <- new.env(parent = emptyenv())
  seen$value <- FALSE
  mock_spatalk_bindings(dec_seen = seen)
  out <- RunSpaTalk(
    spatial,
    group.by = "celltype", mode = "spot",
    reference = reference, reference.group.by = "ref_type",
    coord.cols = c("col", "row"), store.object = "full",
    backend = "r", verbose = FALSE
  )
  expect_true(seen$value)
  expect_identical(seen$method, 1L)
  expect_null(seen$dec_result)
  expect_true(is.list(out@tools$SpaTalk$results$default$native_object))
  expect_true(is.matrix(out@tools$SpaTalk$results$default$deconvolution))
  expect_true(is.matrix(out@tools$SpaTalk$results$default$spatial_plot$composition))
  expect_equal(
    unname(rowSums(out@tools$SpaTalk$results$default$spatial_plot$composition)),
    rep(1, ncol(spatial))
  )
})

test_that("RunSpaTalk reuses stored RCTD weights with SpaTalk method 2", {
  spatial <- make_spatalk_test_object()
  reference <- make_spatalk_test_object()
  reference$ref_type <- c("A", "A", "B", "B")
  weights <- matrix(
    c(0.1, 0.9, 0.2, 0.8, 0.3, 0.7, 0.4, 0.6),
    ncol = 2, byrow = TRUE,
    dimnames = list(rev(colnames(spatial)), c("A", "B"))
  )
  spatial@tools$RCTD <- list(
    weights = weights,
    parameters = list(coordinate_contract_version = 2L)
  )
  seen <- new.env(parent = emptyenv())
  seen$value <- FALSE
  mock_spatalk_bindings(dec_seen = seen)

  out <- RunSpaTalk(
    spatial,
    group.by = "celltype", mode = "spot",
    reference = reference, reference.group.by = "ref_type",
    deconvolution = "RCTD", coord.cols = c("col", "row"),
    backend = "r", verbose = FALSE
  )

  expected <- weights[colnames(spatial), , drop = FALSE]
  expect_true(seen$value)
  expect_identical(seen$method, 2L)
  expect_identical(seen$dec_result, expected)
  expect_identical(
    out@tools$SpaTalk$parameters$deconvolution_source,
    "srt@tools$RCTD$weights"
  )
  expect_identical(out@tools$SpaTalk$results$default$deconvolution, expected)
  expect_null(out@tools$SpaTalk$parameters$backend_args$dec_result)
})

test_that("RunSpaTalk accepts an explicit external deconvolution matrix", {
  spatial <- make_spatalk_test_object()
  reference <- make_spatalk_test_object()
  reference$ref_type <- c("A", "A", "B", "B")
  dec_result <- matrix(
    rep(c(0.6, 0.4), ncol(spatial)),
    ncol = 2, byrow = TRUE,
    dimnames = list(colnames(spatial), c("A", "B"))
  )
  seen <- new.env(parent = emptyenv())
  mock_spatalk_bindings(dec_seen = seen)

  out <- RunSpaTalk(
    spatial,
    group.by = "celltype", mode = "spot",
    reference = reference, reference.group.by = "ref_type",
    deconvolution = "none", dec_result = dec_result,
    coord.cols = c("col", "row"), backend = "r", verbose = FALSE
  )

  expect_identical(seen$method, 2L)
  expect_identical(seen$dec_result, dec_result)
  expect_identical(
    out@tools$SpaTalk$parameters$deconvolution_source,
    "dec_result argument"
  )
})

test_that("RunSpaTalk validates spot inputs before backend execution", {
  spatial <- make_spatalk_test_object()
  expect_error(
    RunSpaTalk(
      spatial,
      group.by = "celltype", mode = "spot",
      coord.cols = c("col", "row"), verbose = FALSE
    ),
    "reference"
  )
  reference <- make_spatalk_test_object()
  reference$ref_type <- reference$celltype
  expect_error(
    RunSpaTalk(
      spatial,
      group.by = "celltype", mode = "spot",
      reference = reference, reference.group.by = "ref_type",
      deconvolution = "none", coord.cols = c("col", "row"), verbose = FALSE
    ),
    "dec_result"
  )
  expect_error(
    RunSpaTalk(
      spatial,
      group.by = "celltype", mode = "spot",
      reference = reference, reference.group.by = "ref_type",
      deconvolution = "RCTD", coord.cols = c("col", "row"), verbose = FALSE
    ),
    "Stored RCTD weights are absent"
  )
  invalid <- matrix(
    c(1, 0, 0, 1, 1, 0),
    ncol = 2, byrow = TRUE,
    dimnames = list(colnames(spatial)[1:3], c("A", "B"))
  )
  spatial@tools$RCTD <- list(weights = invalid)
  expect_error(
    RunSpaTalk(
      spatial,
      group.by = "celltype", mode = "spot",
      reference = reference, reference.group.by = "ref_type",
      deconvolution = "RCTD", coord.cols = c("col", "row"), verbose = FALSE
    ),
    "rerun\\s+RunRCTD"
  )
  spatial@tools$RCTD <- list(
    weights = invalid,
    parameters = list(coordinate_contract_version = 2L)
  )
  expect_error(
    RunSpaTalk(
      spatial,
      group.by = "celltype", mode = "spot",
      reference = reference, reference.group.by = "ref_type",
      deconvolution = "RCTD", coord.cols = c("col", "row"), verbose = FALSE
    ),
    "is missing 1 selected"
  )
})

test_that("RunSpaTalk never silently selects one image", {
  srt <- make_spatalk_test_object()
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 1), row.names = c("cell1", "cell2")),
    type = "centroids", assay = assay, key = "s1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(1, 2), row.names = c("cell3", "cell4")),
    type = "centroids", assay = assay, key = "s2_"
  )
  expect_error(
    RunSpaTalk(srt, group.by = "celltype", verbose = FALSE),
    "Multiple spatial images"
  )
})

test_that("RunSpaTalk accepts one explicitly selected spatial image", {
  srt <- make_spatalk_test_object()
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1, 2, 3), y = c(0, 1, 1, 2), row.names = colnames(srt)),
    type = "centroids", assay = assay, key = "c1_"
  )
  mock_spatalk_bindings()
  out <- RunSpaTalk(
    srt,
    group.by = "celltype", image = "slice1",
    backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$SpaTalk$parameters$image, "slice1")
  expect_identical(out@tools$SpaTalk$cells, colnames(srt))
})

test_that("SpaTalk spatial input rejects invalid IDs and coordinates", {
  srt <- make_spatalk_test_object()
  duplicate <- data.frame(
    cell_id = c("cell1", "cell1"), x = c(0, 1), y = c(0, 1)
  )
  expect_error(
    testthat::with_mocked_bindings(
      spatalk_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = duplicate, source = list()),
      .package = "scop"
    ),
    "valid unique"
  )
  unknown <- data.frame(cell_id = "absent", x = 0, y = 0)
  expect_error(
    testthat::with_mocked_bindings(
      spatalk_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = unknown, source = list()),
      .package = "scop"
    ),
    "valid unique"
  )
  nonfinite <- data.frame(cell_id = "cell1", x = NA_real_, y = 0)
  expect_error(
    testthat::with_mocked_bindings(
      spatalk_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = nonfinite, source = list()),
      .package = "scop"
    ),
    "finite raw"
  )
})

test_that("SpaTalk expression, labels, and coordinates use the same ID order", {
  srt <- make_spatalk_test_object()
  reordered <- data.frame(
    cell_id = rev(colnames(srt)), x = 4:1, y = 1:4,
    stringsAsFactors = FALSE
  )
  input <- testthat::with_mocked_bindings(
    spatalk_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
    SpatialCoordinates = function(...) list(data = reordered, source = list()),
    .package = "scop"
  )
  expect_identical(input$cells, reordered$cell_id)
  expect_identical(colnames(input$expression), reordered$cell_id)
  expect_identical(input$labels, as.character(srt[[]][reordered$cell_id, "celltype"]))
})

test_that("SpaTalk backend failure leaves the input object unchanged", {
  srt <- make_spatalk_test_object()
  before <- list(
    metadata = srt[[]], assays = srt@assays, reductions = srt@reductions,
    graphs = srt@graphs, tools = srt@tools
  )
  mock_spatalk_bindings(fail_cci = TRUE)
  expect_error(
    RunSpaTalk(
      srt,
      group.by = "celltype", coord.cols = c("col", "row"),
      backend = "r", verbose = FALSE
    ),
    "backend failed"
  )
  expect_identical(srt[[]], before$metadata)
  expect_identical(srt@assays, before$assays)
  expect_identical(srt@reductions, before$reductions)
  expect_identical(srt@graphs, before$graphs)
  expect_identical(srt@tools, before$tools)
})

test_that("RunCCC dispatches SpaTalk only when explicitly selected", {
  srt <- make_spatalk_test_object()
  mock_spatalk_bindings()
  out <- RunCCC(
    srt,
    group.by = "celltype", methods = "spatalk",
    method_params = list(SpaTalk = list(coord.cols = c("col", "row"))),
    backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$RunCCC$completed_methods, "SpaTalk")
  expect_true("SpaTalk" %in% out@tools$CCC$methods)
  expect_false("SpaTalk" %in% eval(formals(RunCCC)$methods))
})

test_that("SpaTalkPlot requires an explicit result when several are stored", {
  srt <- make_spatalk_test_object()
  mock_spatalk_bindings()
  out <- RunSpaTalk(
    srt,
    group.by = "celltype", coord.cols = c("col", "row"),
    result.name = "first", backend = "r", verbose = FALSE
  )
  out@tools$SpaTalk$results$first$coordinate_contract_version <- NULL
  out <- RunSpaTalk(
    out,
    group.by = "celltype", coord.cols = c("col", "row"),
    result.name = "second", backend = "r", verbose = FALSE
  )
  expect_error(SpaTalkPlot(out), "select.*result.name")
  expect_error(
    SpaTalkPlot(out, result.name = "first", plot_type = "tf"),
    "rerun\\s+RunSpaTalk"
  )
  expect_s3_class(
    SpaTalkPlot(out, result.name = "second", plot_type = "tf"),
    "ggplot"
  )
})
