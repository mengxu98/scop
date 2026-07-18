make_spatial_neighborhood_seurat <- function() {
  counts <- matrix(
    c(
      2, 0, 1, 3, 0, 1,
      0, 4, 0, 1, 2, 0,
      1, 0, 3, 0, 2, 4
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt$CellType <- c("T", "B", "T", "Myeloid", "B", "T")
  srt$condition <- c("A", "A", "B", "B", "A", "B")
  srt$sample <- c("S1", "S1", "S2", "S2", "S1", "S2")
  srt$x <- c(0, 1, 0, 1, 2, 2)
  srt$y <- c(0, 0, 1, 1, 0, 1)
  srt
}

with_mock_spicyr <- function(code) {
  fake_spicy <- function(cells, condition, subject, cellType, imageID, spatialCoords, ...) {
    expect_true(all(c(condition, subject, cellType, imageID, spatialCoords) %in% colnames(cells)))
    expect_identical(condition, ".scop_condition")
    expect_identical(cellType, ".scop_cell_type")
    expect_identical(imageID, ".scop_image_id")
    data.frame(
      cellType1 = c("T", "B", "T"),
      cellType2 = c("B", "T", "Myeloid"),
      coefficient = c(1.2, -0.8, 0.1),
      p.value = c(0.001, 0.02, 0.9),
      stringsAsFactors = FALSE
    )
  }
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_true(packages %in% c("spicyR", "BiocNeighbors"))
      invisible(TRUE)
    },
    get_namespace_fun = function(package, name) {
      if (identical(package, "spicyR")) {
        expect_identical(name, "spicy")
        return(fake_spicy)
      }
      expect_identical(package, "BiocNeighbors")
      getExportedValue(package, name)
    }
  )
  force(code)
}

test_that("RunSpatialNeighborhood stores standardized spicyR results", {
  srt <- make_spatial_neighborhood_seurat()
  with_mock_spicyr({
    out <- RunSpatialNeighborhood(
      srt,
      group.by = "CellType",
      split.by = "condition",
      sample.by = "sample",
      coord.cols = c("x", "y"),
      k = 2,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("SpatialNeighborhood" %in% names(out@tools))
  bundle <- out@tools$SpatialNeighborhood$methods$spicyR
  expect_named(bundle$pair_table[1:10], c(
    "method", "comparison", "condition", "from", "to",
    "estimate", "statistic", "pval", "FDR", "direction"
  ))
  expect_equal(bundle$pair_table$from, c("T", "B", "T"))
  expect_equal(bundle$pair_table$direction[1:2], c("enriched", "depleted"))
  expect_true(all(c("cell", "neighbor", "from", "to", "distance") %in% colnames(bundle$edge_table)))
  expect_named(bundle$summary, c("n_pairs", "n_edges", "top_pairs"))
  expect_equal(out@tools$SpatialNeighborhood$summary$n_pairs, nrow(bundle$pair_table))
  expect_identical(out@tools$SpatialNeighborhood$provenance$backend_id, "spicyr")
  expect_identical(bundle$provenance$backend_id, "spicyr")
})

test_that("RunSpatialNeighborhood defaults to native observed summaries", {
  srt <- make_spatial_neighborhood_seurat()
  out <- RunSpatialNeighborhood(
    srt,
    group.by = "CellType",
    sample.by = "sample",
    coord.cols = c("x", "y"),
    k = 2,
    verbose = FALSE
  )
  bundle <- out@tools$SpatialNeighborhood$methods$observed
  expect_true(nrow(bundle$pair_table) > 0)
  expect_true(all(!is.na(bundle$pair_table$subject)))
  expect_equal(unique(bundle$pair_table$method), "observed")
  expect_identical(out@tools$SpatialNeighborhood$provenance$backend_id, "core")
  expect_identical(out@tools$SpatialNeighborhood$source$coordinate_space, "raw")
  expect_identical(bundle$provenance$backend_id, "core")
  expect_null(bundle$raw)
})

test_that("RunSpatialNeighborhood refuses to impersonate an unrun spicyR backend", {
  srt <- make_spatial_neighborhood_seurat()
  expect_error(
    RunSpatialNeighborhood(
      srt,
      group.by = "CellType",
      method = "spicyR",
      coord.cols = c("x", "y"),
      verbose = FALSE
    ),
    "split.by.*required"
  )
  expect_false("SpatialNeighborhood" %in% names(srt@tools))
})

test_that("spicyR empty and malformed outputs cannot impersonate observed results", {
  observed <- data.frame(
    method = "observed", comparison = "A", condition = "A",
    from = "T", to = "B", estimate = 0.5, statistic = NA_real_,
    pval = NA_real_, FDR = NA_real_, direction = "observed",
    sample = "s1", subject = "s1", count = 1, total = 2,
    fraction = 0.5, stringsAsFactors = FALSE
  )
  expect_error(
    scop:::spatial_neighborhood_standardize_pair_table(
      list(raw = list(), table = data.frame()),
      observed,
      "spicyR"
    ),
    "returned no"
  )
  expect_error(
    scop:::spatial_neighborhood_standardize_pair_table(
      list(raw = list(), table = data.frame(estimate = 1)),
      observed,
      "spicyR"
    ),
    "incompatible result"
  )
  expect_identical(
    scop:::spatial_neighborhood_standardize_pair_table(
      list(raw = NULL, table = NULL),
      observed,
      "observed"
    ),
    observed
  )
})

test_that("SpatialNeighborhoodPlot returns scop-style ggplot objects", {
  srt <- make_spatial_neighborhood_seurat()
  with_mock_spicyr({
    srt <- RunSpatialNeighborhood(
      srt,
      group.by = "CellType",
      method = "spicyR",
      split.by = "condition",
      sample.by = "sample",
      coord.cols = c("x", "y"),
      k = 2,
      verbose = FALSE
    )
  })

  p_heatmap <- SpatialNeighborhoodPlot(srt, plot_type = "heatmap", theme_use = NULL)
  p_stat <- SpatialNeighborhoodPlot(srt, plot_type = "stat", theme_use = NULL)
  p_network <- SpatialNeighborhoodPlot(srt, plot_type = "network", theme_use = NULL)
  p_empty <- SpatialNeighborhoodPlot(
    srt,
    plot_type = "heatmap",
    comparison = "missing",
    theme_use = NULL
  )
  p_spatial <- SpatialNeighborhoodPlot(
    srt,
    plot_type = "spatial",
    pair = c("T", "B"),
    coord.cols = c("x", "y"),
    overlay_image = FALSE,
    theme_use = NULL
  )

  expect_s3_class(p_heatmap, "ggplot")
  expect_s3_class(p_stat, "ggplot")
  expect_s3_class(p_network, "ggplot")
  expect_s3_class(p_empty, "ggplot")
  expect_s3_class(p_spatial, "ggplot")
})

test_that("RunSpatialNeighborhood validates spatial inputs clearly", {
  srt <- make_spatial_neighborhood_seurat()
  expect_error(
    RunSpatialNeighborhood(srt, group.by = "missing", coord.cols = c("x", "y"), verbose = FALSE),
    "Missing metadata"
  )
  expect_error(
    RunSpatialNeighborhood(srt, group.by = "CellType", coord.cols = c("missing_x", "y"), verbose = FALSE),
    "Spatial coordinates"
  )
  expect_error(
    RunSpatialNeighborhood(srt, group.by = "CellType", method = "mistyR", coord.cols = c("x", "y"), verbose = FALSE),
    "should be one of"
  )
})
