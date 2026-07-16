make_spatial_integration_seurat <- function(samples = c("S1", "S2")) {
  counts <- matrix(
    c(
      8, 7, 1, 0, 0, 1,
      0, 1, 7, 8, 2, 0,
      5, 0, 1, 0, 6, 7,
      0, 4, 2, 1, 0, 8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:6))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$sample <- rep(samples, each = 3)
  srt$col <- c(1, 2, 3, 1, 2, 3)
  srt$row <- c(1, 1, 1, 2, 2, 2)
  srt
}

mock_spatial_integration_backend <- function(code) {
  testthat::local_mocked_bindings(
    spatial_integration_run_backend = function(method, input, verbose = TRUE, ...) {
      expect_identical(method, "PRECAST")
      expect_equal(input$samples, c("S1", "S2"))
      expect_s4_class(input$expr, "dgCMatrix")
      cells <- input$cells
      list(
        embedding = matrix(
          seq_len(length(cells) * 2),
          ncol = 2,
          dimnames = list(cells, c("PC1", "PC2"))
        ),
        domains = stats::setNames(rep(c("D1", "D2"), length.out = length(cells)), cells),
        aligned_coords = data.frame(
          x = input$coords$x + ifelse(input$srt@meta.data[cells, "sample"] == "S2", 10, 0),
          y = input$coords$y,
          row.names = cells
        ),
        raw_result = list(mock = TRUE)
      )
    }
  )
  force(code)
}

test_that("RunSpatialIntegration writes standardized results for merged Seurat", {
  srt <- make_spatial_integration_seurat()
  mock_spatial_integration_backend({
    out <- RunSpatialIntegration(
      srt,
      method = "PRECAST",
      sample.by = "sample",
      assay = "RNA",
      layer = "counts",
      coord.cols = c("col", "row"),
      verbose = FALSE
    )
  })

  expect_equal(unname(out$SpatialIntegration_PRECAST_domain), rep(c("D1", "D2"), 3))
  expect_true("SpatialIntegration_PRECAST" %in% SeuratObject::Reductions(out))
  expect_true("SpatialIntegration" %in% names(out@tools))
  expect_equal(out@tools$SpatialIntegration$active_method, "PRECAST")
  expect_equal(out@tools$SpatialIntegration$parameters$sample.by, "sample")
  expect_identical(out@tools$SpatialIntegration$source$coordinate_space, "raw")
  expect_named(out@tools$SpatialIntegration$summary, c("n_cells", "domains", "samples"))
  expect_named(out@tools$SpatialIntegration$summary$domains, c("domain", "count"))
  expect_named(out@tools$SpatialIntegration$summary$samples, c("sample", "count"))
  expect_true(all(c(
    "SpatialIntegration_PRECAST_aligned_x",
    "SpatialIntegration_PRECAST_aligned_y"
  ) %in% colnames(out@meta.data)))
})

test_that("RunSpatialIntegration supports list input and preserves sample labels", {
  srt <- make_spatial_integration_seurat()
  srt_list <- Seurat::SplitObject(srt, split.by = "sample")
  mock_spatial_integration_backend({
    out <- RunSpatialIntegration(
      srt_list,
      method = "PRECAST",
      assay = "RNA",
      layer = "counts",
      coord.cols = c("col", "row"),
      verbose = FALSE
    )
  })

  expect_true(".scop_spatial_sample" %in% colnames(out@meta.data))
  expect_equal(sort(unique(out$.scop_spatial_sample)), c("S1", "S2"))
  expect_true("SpatialIntegration_PRECAST" %in% SeuratObject::Reductions(out))
  expect_true(all(grepl("^S[12]_", colnames(out))))
})

test_that("RunSpatialIntegration validates required inputs before backend work", {
  srt <- make_spatial_integration_seurat()
  testthat::local_mocked_bindings(
    spatial_integration_run_backend = function(...) stop("backend should not run")
  )
  expect_error(
    RunSpatialIntegration(srt, method = "PRECAST", verbose = FALSE),
    "sample.by"
  )
  expect_error(
    RunSpatialIntegration(
      srt,
      method = "PRECAST",
      sample.by = "missing",
      verbose = FALSE
    ),
    "sample.by"
  )
  expect_error(
    RunSpatialIntegration(
      srt,
      method = "PRECAST",
      sample.by = "sample",
      features = "AbsentGene",
      verbose = FALSE
    ),
    "features"
  )
})

test_that("SpatialIntegrationPlot reuses SCOP plot helpers", {
  srt <- make_spatial_integration_seurat()
  mock_spatial_integration_backend({
    out <- RunSpatialIntegration(
      srt,
      method = "PRECAST",
      sample.by = "sample",
      assay = "RNA",
      layer = "counts",
      coord.cols = c("col", "row"),
      verbose = FALSE
    )
  })

  p_spatial <- SpatialIntegrationPlot(
    out,
    plot_type = "spatial",
    coord.cols = c("col", "row"),
    overlay_image = FALSE
  )
  p_embedding <- SpatialIntegrationPlot(out, plot_type = "embedding")
  p_composition <- SpatialIntegrationPlot(out, plot_type = "composition")
  p_alignment <- SpatialIntegrationPlot(out, plot_type = "alignment")

  expect_s3_class(p_spatial, "ggplot")
  expect_s3_class(p_embedding, "ggplot")
  expect_s3_class(p_composition, "ggplot")
  expect_s3_class(p_alignment, "ggplot")
})

test_that("SpatialMNN follows the atlasClustering three-stage API", {
  srt <- make_spatial_integration_seurat()
  observed <- character()
  testthat::local_mocked_bindings(
    spatialmnn_check_r = function() invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "atlasClustering")
      switch(name,
        stage_1 = function(seu_ls, verbose) {
          observed <<- c(observed, "stage_1")
          expect_true(all(vapply(seu_ls, function(x) {
            all(c("coord_x", "coord_y") %in% colnames(x@meta.data)) &&
              "RNA" %in% SeuratObject::Assays(x) &&
              identical(SeuratObject::DefaultAssay(x), "RNA")
          }, logical(1))))
          for (i in seq_along(seu_ls)) {
            seu_ls[[i]]$merged_cluster <- rep(1:2, length.out = ncol(seu_ls[[i]]))
          }
          seu_ls
        },
        stage_2 = function(seu_ls, method, rtn_seurat, verbose) {
          observed <<- c(observed, "stage_2")
          expect_identical(method, "MNN")
          expect_true(rtn_seurat)
          list(cl_df = do.call(rbind, lapply(names(seu_ls), function(sample) {
            data.frame(sample = sample, cluster = 1:2, louvain = c("D1", "D2"))
          })))
        },
        stop("unexpected atlasClustering function")
      )
    }
  )

  out <- RunSpatialIntegration(
    srt,
    method = "SpatialMNN",
    sample.by = "sample",
    assay = "RNA",
    layer = "counts",
    coord.cols = c("col", "row"),
    coordinate_space = "raw",
    verbose = FALSE
  )
  expect_identical(observed, c("stage_1", "stage_2"))
  expect_true(all(!is.na(out$SpatialIntegration_SpatialMNN_domain)))
  expect_identical(out@tools$SpatialIntegration$source$coordinate_space, "raw")
})

test_that("PRECAST receives selected features and its SelectModel object argument", {
  input <- list(
    srt_list = list(S1 = "slice-1", S2 = "slice-2"),
    features = c("Gene1", "Gene2")
  )
  observed <- list()
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "PRECAST")
      switch(name,
        CreatePRECASTObject = function(seuList, project, customGenelist) {
          observed$features <<- customGenelist
          list(step = "created")
        },
        AddAdjList = function(PRECASTObj, ...) PRECASTObj,
        AddParSetting = function(PRECASTObj, ...) PRECASTObj,
        PRECAST = function(PRECASTObj, ...) PRECASTObj,
        SelectModel = function(obj, ...) {
          observed$selected <<- obj
          obj
        },
        stop("unexpected PRECAST function")
      )
    },
    spatial_integration_extract_backend = function(raw_result, method, input) raw_result
  )

  out <- scop:::spatial_integration_run_precast(input, params = list(), verbose = FALSE)
  expect_identical(observed$features, input$features)
  expect_identical(observed$selected, list(step = "created"))
  expect_identical(out, list(step = "created"))
})

test_that("SpatialMNN installation discovery uses the actual package name", {
  calls <- character()
  testthat::local_mocked_bindings(
    check_r = function(package, install = TRUE, verbose = TRUE) {
      calls <<- c(calls, package)
      if (identical(package, "atlasClustering")) {
        return(list(atlasClustering = FALSE))
      }
      list(atlasClustering = TRUE)
    }
  )

  expect_invisible(scop:::spatialmnn_check_r())
  expect_identical(calls, c("atlasClustering", "Pixel-Dream/spatialMNN"))
})

test_that("BASS discovery rejects the unrelated package-name collision", {
  checks <- 0L
  installs <- character()
  testthat::local_mocked_bindings(
    bass_namespace_ready = function() {
      checks <<- checks + 1L
      checks > 1L
    },
    check_r = function(package, force = FALSE, verbose = TRUE) {
      installs <<- c(installs, package)
      expect_true(force)
      invisible(TRUE)
    }
  )

  expect_invisible(scop:::bass_check_r())
  expect_identical(installs, "zhengli09/BASS")
})

test_that("integration backends declare all producer entry points", {
  registry <- scop:::spatial_backend_registry()
  expect_setequal(
    scop:::spatial_backend_required_symbols(registry$precast),
    c("CreatePRECASTObject", "AddAdjList", "AddParSetting", "PRECAST", "SelectModel")
  )
  expect_setequal(
    scop:::spatial_backend_required_symbols(registry$bass),
    c("createBASSObject", "BASS.preprocess", "BASS.run")
  )
  expect_identical(registry$spatialmnn$package, "atlasClustering")
  expect_setequal(
    scop:::spatial_backend_required_symbols(registry$spatialmnn),
    c("stage_1", "stage_2")
  )
})
