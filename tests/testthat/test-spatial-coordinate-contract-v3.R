coordinate_v3_fixture <- function(images = TRUE) {
  counts <- matrix(seq_len(24), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:8)))
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$sample <- rep(c("S1", "S2"), each = 4)
  srt$type <- rep(c("A", "B"), 4)
  srt$x <- seq_len(8)
  srt$y <- seq_len(8) * 2
  if (images) for (i in 1:2) {
    cells <- colnames(srt)[srt$sample == paste0("S", i)]
    srt[[paste0("slice", i)]] <- SeuratObject::CreateFOV(
      data.frame(x = c(10, 20, 30, 40), y = c(30, 40, 60, 50), row.names = cells),
      type = "centroids", assay = "RNA", key = paste0("v", i, "_"))
  }
  srt
}

test_that("VisiumV2 image provenance survives metadata edits and cell operations", {
  data(visium_human_pancreas_sub, package = "scop")
  srt <- visium_human_pancreas_sub
  expected <- SpatialCoordinates(srt, image = "slice1")$data
  srt$x <- NULL; srt$y <- NULL
  expect_equal(SpatialCoordinates(srt, image = "slice1")$data, expected)
  srt$x <- rev(expected$y); srt$y <- rev(expected$x)
  expect_equal(SpatialCoordinates(srt, image = "slice1")$data, expected)
  cells <- expected$cell_id[c(5, 2, 1)]
  smaller <- suppressWarnings(srt[, cells])
  expect_equal(SpatialCoordinates(smaller)$data[, c("x", "y")], expected[colnames(smaller), c("x", "y")])
  renamed <- SeuratObject::RenameCells(smaller, new.names = paste0("new", 1:3))
  expect_equal(SpatialCoordinates(renamed)$data$x, expected[colnames(smaller), "x"])
  expect_identical(renamed@misc$spatial_image_axes$slice1, "horizontal")
  p <- SpatialSpotPlot(srt, group.by = "coda_label", image = "slice1")
  expect_equal(p$data$x, expected$x * srt[["slice1"]]@scale.factors$lowres)
  srt@misc$spatial_image_axes$slice1 <- "unknown"
  expect_error(SpatialCoordinates(srt), "coords_x_orientation")
})

test_that("numeric text coordinates retain values and true source columns", {
  srt <- coordinate_v3_fixture(FALSE)
  srt$x <- factor(as.character(seq_len(8) * 100))
  srt$y <- factor(as.character(seq_len(8) * 1000))
  coords <- SpatialCoordinates(srt)
  expect_equal(coords$data$x, seq_len(8) * 100)
  expect_equal(coords$data$y, seq_len(8) * 1000)
  expect_identical(coords$source$coord.cols, c("x", "y"))
  expect_identical(coords$transform$raw_x_col, "x")
  srt$x <- c("bad", as.character(2:8))
  expect_error(SpatialCoordinates(srt), "non-finite")
  srt$x <- rep(TRUE, 8)
  expect_error(SpatialCoordinates(srt), "numeric")
})

test_that("explicit image axes migrate legacy objects before subsetting and merging", {
  data(visium_human_pancreas_sub, package = "scop")
  original <- visium_human_pancreas_sub
  original@misc$spatial_image_axes <- NULL
  original <- SetSpatialImageAxes(original, image = "slice1", x_orientation = "horizontal")
  raw <- SpatialCoordinates(original)$data
  a <- suppressWarnings(original[, 1:4])
  b <- suppressWarnings(original[, 5:8])
  parts <- spatial_integration_as_list(list(S1 = a, S2 = b), "sample")
  merged <- spatial_integration_merge_list(parts, "sample")
  coords <- spatial_sample_coords(merged, "sample")$data
  expect_equal(coords$x, raw[1:8, "x"])
  expect_equal(coords$y, raw[1:8, "y"])
  expect_length(merged@misc$spatial_image_axes, 2L)
  expect_error(SetSpatialImageAxes(coordinate_v3_fixture(), image = "slice1"), "VisiumV2")
})

test_that("FOV point pie and network final axes agree", {
  skip_if_not_installed("BiocNeighbors")
  skip_if_not_installed("scatterpie")
  srt <- coordinate_v3_fixture()
  point <- SpatialSpotPlot(srt, image = "slice1", group.by = "type", crop = FALSE)
  pie <- SpatialSpotPlot(srt, image = "slice1", plot_type = "pie",
    values = matrix(1, nrow = 8, ncol = 2, dimnames = list(colnames(srt), c("A", "B"))))
  graph <- RunSpatialNetwork(srt, image = "slice1", k = 1, verbose = FALSE)
  network <- SpatialNetworkPlot(graph)
  point_built <- ggplot2::ggplot_build(point)
  pie_built <- ggplot2::ggplot_build(pie)
  network_built <- ggplot2::ggplot_build(network)
  expect_equal(point_built$data[[1]]$y, network_built$data[[2]]$y)
  expect_identical(pie_built$layout$panel_scales_y[[1]]$get_transformation()$name, "identity")
  manual <- SpatialSpotPlot(srt, image = "slice1", group.by = "type", flip.y = TRUE)
  expect_equal(ggplot2::ggplot_build(manual)$data[[1]]$y, -point_built$data[[1]]$y)
})

test_that("SpatialExperiment preserves raw units through display exports and subsets", {
  skip_if_not_installed("SpatialExperiment")
  data(visium_human_pancreas_sub, package = "scop")
  srt <- suppressWarnings(visium_human_pancreas_sub[, c(5, 2, 1, 4)])
  raw <- SpatialCoordinates(srt)$data
  for (space in c("raw", "legacy_display")) {
    spe <- srt_to_spe(srt, coordinate_space = space)
    expect_identical(S4Vectors::metadata(spe)$scop_spatial_coordinates$source$coordinate_space, space)
    spe <- spe[, c(4, 1)]
    out <- spe_to_srt(spe)
    restored <- SpatialCoordinates(out)$data
    expect_equal(restored$x, raw[colnames(spe), "x"])
    expect_equal(restored$y, raw[colnames(spe), "y"])
    expect_identical(restored$cell_id, colnames(spe))
  }
  default <- srt_to_spe(srt)
  expect_equal(unname(SpatialExperiment::spatialCoords(default)), unname(as.matrix(raw[, c("x", "y")])))
  bad <- srt_to_spe(srt, coordinate_space = "legacy_display")
  S4Vectors::metadata(bad)$scop_spatial_coordinates$transform$scale <- NULL
  expect_error(spe_to_srt(bad), "saved scale")
})

test_that("integration covers all image-backed samples and records their sources", {
  srt <- coordinate_v3_fixture()
  testthat::local_mocked_bindings(spatial_integration_run_backend = function(method, input, ...) {
    expect_setequal(input$samples, c("S1", "S2"))
    expect_length(input$cells, 8)
    expect_equal(input$coords$x, rep(c(10, 20, 30, 40), 2))
    list(domains = stats::setNames(rep("D1", 8), input$cells), raw_result = list())
  })
  for (object in list(srt, Seurat::SplitObject(srt, split.by = "sample"))) {
    out <- RunSpatialIntegration(object, sample.by = "sample", verbose = FALSE)
    sources <- out@tools$SpatialIntegration$parameters$coordinate_sources
    expect_named(sources, c("S1", "S2"))
    expect_identical(sources$S1$image, "slice1")
    expect_identical(sources$S2$image, "slice2")
    expect_false(anyNA(out$SpatialIntegration_PRECAST_domain))
  }
  expect_error(RunSpatialIntegration(srt, sample.by = "sample", image = "slice1", verbose = FALSE), "cover every cell")
  srt[["duplicate"]] <- srt[["slice1"]]
  expect_error(RunSpatialIntegration(srt, sample.by = "sample", verbose = FALSE), "one image covering")
  expect_no_error(RunSpatialIntegration(srt, sample.by = "sample",
    image = c(S1 = "slice1", S2 = "slice2"), verbose = FALSE))
})

test_that("SpatialEcoTyper analysis uses the same selected raw image coordinates", {
  srt <- coordinate_v3_fixture()
  observed <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      if (package != "SpatialEcoTyper") return(get(name, asNamespace(package)))
      if (name %in% c("mostFrequent", "GetSpatialMetacells", "GetPCList")) return(function(...) NULL)
      if (name == "SpatialEcoTyper") return(function(normdata, metadata, ...) {
        observed <<- metadata
        list(metadata = data.frame(SE = rep("SE1", ncol(normdata)), row.names = colnames(normdata)))
      })
      function(data_list, metadata_list, ...) {
        observed <<- metadata_list
        data.frame(CID = unlist(lapply(data_list, colnames), use.names = FALSE),
          Sample = rep(names(data_list), lengths(lapply(data_list, colnames))), InitSE = "I1", SE = "SE1")
      }
    })
  expect_error(RunSpatialEcoTyper(srt, layer = "counts", celltype.by = "type"), "Multiple spatial images")
  out <- RunSpatialEcoTyper(srt, layer = "counts", celltype.by = "type", image = "slice1", verbose = FALSE)
  expect_equal(observed$X, c(10, 20, 30, 40))
  expect_identical(rownames(observed), colnames(srt)[1:4])
  expect_true(all(is.na(out$SpatialEcoTyper_SE[5:8])))
  expect_equal(out@tools$SpatialEcoTyper$coordinates$x, observed$X)
  expect_identical(out@tools$SpatialEcoTyper$source$samples$single$image, "slice1")
  out <- RunSpatialEcoTyper(srt, mode = "multi", sample.by = "sample", layer = "counts",
    celltype.by = "type", outdir = tempfile(), verbose = FALSE)
  expect_named(observed, c("S1", "S2"))
  expect_equal(observed$S1$X, observed$S2$X)
  expect_false(anyNA(out$SpatialEcoTyper_SE))
  expect_identical(out@tools$SpatialEcoTyper$source$samples$S2$image, "slice2")
  expect_error(spatialecotyper_build_metadata(data.frame(X = Inf, Y = 1, ct = "A"), "ct", "X", "Y"), "invalid")
})

test_that("SpotSweeper sample prefixes are restored only with a verified mapping", {
  skip_if_not_installed("SpatialExperiment")
  srt <- coordinate_v3_fixture(FALSE)
  spe <- srt_to_spe(srt)
  output <- spe[, c(8, 2, 5, 1, 6, 3, 7, 4)]
  colnames(output) <- paste(SummarizedExperiment::colData(output)$sample, colnames(output), sep = ".")
  SummarizedExperiment::colData(output)$value_z <- seq_len(8)
  restored <- spot_sweeper_align_local_output(output, spe, "sample", "value")
  expect_identical(colnames(restored), colnames(spe))
  expect_equal(SummarizedExperiment::colData(restored)$value_z, match(seq_len(8), c(8, 2, 5, 1, 6, 3, 7, 4)))
  changed <- output
  SpatialExperiment::spatialCoords(changed)[1, 1] <- 999
  expect_error(spot_sweeper_align_local_output(changed, spe, "sample", "value"), "changed spot")
  colnames(changed) <- paste0("wrong.", colnames(output))
  expect_error(spot_sweeper_align_local_output(changed, spe, "sample", "value"), "changed spot")
  expect_error(spot_sweeper_align_local_output(output[, -1], spe, "sample", "value"), "changed spot")
})

test_that("v2 coordinate results are rejected after the semantic repair", {
  expect_error(spatial_require_coordinate_contract(list(coordinate_contract_version = 2L), "producer"), "rerun")
})
