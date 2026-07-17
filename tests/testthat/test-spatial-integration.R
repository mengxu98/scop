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

test_that("PRECAST adapter supplies backend row and col coordinates without mutating input", {
  srt <- make_spatial_integration_seurat()
  srt$x <- srt$col * 10
  srt$y <- srt$row * 10
  srt[["row"]] <- NULL
  srt[["col"]] <- NULL
  input <- scop:::spatial_integration_prepare_input(
    object = srt,
    sample.by = "sample",
    assay = "RNA",
    layer = "counts",
    features = NULL,
    image = NULL,
    coord.cols = c("x", "y")
  )
  input$features <- input$features[seq_len(2L)]

  adapted <- scop:::spatial_integration_precast_seurat_list(input)

  expect_false(any(c("row", "col") %in% colnames(input$srt_list[[1L]]@meta.data)))
  expect_identical(rownames(adapted[[1L]][[input$assay]]), input$features)
  expect_equal(
    unname(adapted[[1L]]$row),
    unname(input$coords_list[[1L]][colnames(adapted[[1L]]), "y"])
  )
  expect_equal(
    unname(adapted[[1L]]$col),
    unname(input$coords_list[[1L]][colnames(adapted[[1L]]), "x"])
  )
})

test_that("PRECAST extractor maps its selected latent embedding and domains to cells", {
  skip_if_not_installed("PRECAST")
  data("PRECASTObj", package = "PRECAST")
  selected <- PRECAST::SelectModel(PRECASTObj)
  cells <- unlist(lapply(selected@seulist, colnames), use.names = FALSE)
  input <- list(cells = cells)

  extracted <- scop:::spatial_integration_extract_precast(selected, input)

  expect_equal(rownames(extracted$embedding), cells)
  expect_equal(names(extracted$domains), cells)
  expect_equal(nrow(extracted$embedding), length(cells))
  expect_equal(length(extracted$domains), length(cells))
  expect_match(colnames(extracted$embedding)[1L], "^PRECAST_")
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
