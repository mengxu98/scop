test_that("spatial graph core handles sparse edge cases", {
  skip_if_not_installed("BiocNeighbors")
  collinear <- data.frame(cell_id = letters[1:5], x = 1:5, y = 0)
  graph <- scop:::spatial_graph_compute(collinear, k = 2)
  expect_true(all(graph$edges$from < graph$edges$to))
  expect_equal(anyDuplicated(paste(graph$edges$from, graph$edges$to)), 0L)
  expect_true(all(graph$edges$distance > 0))

  duplicate <- data.frame(cell_id = letters[1:4], x = c(0, 0, 1, 2), y = 0)
  duplicate_graph <- scop:::spatial_graph_compute(duplicate, k = 1)
  expect_true(any(duplicate_graph$edges$distance == 0))

  empty <- scop:::spatial_graph_compute(collinear, method = "radius", radius = 0.1)
  expect_equal(nrow(empty$nodes), 5L)
  expect_equal(nrow(empty$edges), 0L)
  expect_error(scop:::spatial_graph_compute(collinear[1, ], k = 1), "At least two")
  expect_error(scop:::spatial_graph_compute(collinear, k = 5), "number of nodes minus one")
})

test_that("spatial graph weights preserve raw distance", {
  skip_if_not_installed("BiocNeighbors")
  coords <- data.frame(cell_id = letters[1:3], x = c(0, 3, 0), y = c(0, 0, 4))
  binary <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE)
  inverse <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE, weight = "inverse_distance")
  gaussian <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE, weight = "gaussian", sigma = 2)
  expect_equal(binary$edges$distance, inverse$edges$distance)
  expect_equal(inverse$edges$weight, 1 / (1 + inverse$edges$distance))
  expect_equal(gaussian$edges$weight, exp(-(gaussian$edges$distance^2) / 8))
})

test_that("raw and display coordinate transforms round trip", {
  raw <- data.frame(cell_id = c("a", "b"), x = c(2, 8), y = c(4, 10))
  transform <- list(scale = 0.5, y_flip = TRUE, image_height = 20)
  display <- scop:::spatial_coords_to_display(raw, transform)
  expect_equal(display$x, c(1, 4))
  expect_equal(display$y, c(18, 15))
  expect_equal(scop:::spatial_coords_to_raw(display, transform), raw)
})

test_that("RunSpatialNetwork uses deterministic slots and collision protection", {
  skip_if_not_installed("BiocNeighbors")
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  out <- RunSpatialNetwork(srt, k = 2, verbose = FALSE)
  expect_identical(out@tools$SpatialNetwork$active_graph, "knn_k2")
  expect_identical(out@tools$SpatialNetwork$graphs$knn_k2$source$coordinate_space, "raw")
  expect_error(RunSpatialNetwork(out, k = 2, verbose = FALSE), "already exists")
  out <- RunSpatialNetwork(out, k = 2, overwrite = TRUE, verbose = FALSE)
  out <- RunSpatialNetwork(out, method = "radius", radius = 1.5, verbose = FALSE)
  expect_identical(out@tools$SpatialNetwork$active_graph, "radius_r1p5")
  expect_true(all(c("knn_k2", "radius_r1p5") %in% names(out@tools$SpatialNetwork$graphs)))
})

test_that("dedicated spatial result plots accept result-only input", {
  misty <- list(results = list(
    improvements = data.frame(target = c("A", "B"), measure = "gain.R2", value = c(0.1, 0.4)),
    contributions = data.frame(target = "A", view = "paraview.3", value = 0.3)
  ))
  expect_s3_class(MistyRPlot(res = misty), "ggplot")
  expect_s3_class(MistyRPlot(res = misty, type = "contributions"), "ggplot")

  statial <- list(table = data.frame(
    imageID = rep(c("s1", "s2"), each = 2), test = rep(c("A:B", "B:A"), 2),
    r = rep(c(10, 20), 2), kontextual = c(-1, 0.5, 0.2, 1)
  ))
  expect_s3_class(StatialKontextualPlot(res = statial), "ggplot")
})

test_that("SpatialCoordinates returns ordered raw and display contracts", {
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(2, 0, 3, 1)
  srt$row <- c(0, 1, 1, 0)
  raw <- SpatialCoordinates(srt)
  expect_identical(raw$data$cell_id, colnames(srt))
  expect_identical(raw$source$coordinate_space, "raw")
  expect_equal(raw$data$x, unname(srt$col))
  display <- SpatialCoordinates(srt, space = "display")
  expect_identical(display$source$coordinate_space, "display")
  expect_equal(display$data[, c("x", "y")], raw$data[, c("x", "y")])

  subsetted <- subset(srt, cells = c("spot4", "spot2"))
  expect_identical(SpatialCoordinates(subsetted)$data$cell_id, colnames(subsetted))
  renamed <- SeuratObject::RenameCells(srt, new.names = paste0("renamed", seq_len(ncol(srt))))
  expect_identical(SpatialCoordinates(renamed)$data$cell_id, colnames(renamed))

  srt$col[[2L]] <- NA_real_
  expect_error(SpatialCoordinates(srt), "non-finite")
})

test_that("analysis coordinates keep one payload contract across spaces", {
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(2, 0, 3, 1)
  srt$row <- c(0, 1, 1, 0)

  legacy <- scop:::spatial_analysis_coords(srt, coordinate_space = "legacy_display")
  raw <- scop:::spatial_analysis_coords(srt, coordinate_space = "raw")
  required <- c("cell_id", "x", "y", "image")
  expect_identical(colnames(legacy$data), required)
  expect_identical(colnames(raw$data), required)
  expect_identical(legacy$data$cell_id, colnames(srt))
  expect_identical(raw$data$cell_id, colnames(srt))
  expect_identical(legacy$data$image, rep(NA_character_, ncol(srt)))
})

test_that("CytoSPACE resolves coordinates before assignment work", {
  data("visium_human_pancreas_sub", package = "scop")
  srt <- visium_human_pancreas_sub
  spot_ids <- colnames(srt)
  legacy <- scop:::cytospace_get_spatial_coords(
    srt,
    spot_ids,
    coordinate_space = "legacy_display"
  )
  raw <- scop:::cytospace_get_spatial_coords(
    srt,
    spot_ids,
    coordinate_space = "raw"
  )
  expect_identical(rownames(legacy), spot_ids)
  expect_identical(rownames(raw), spot_ids)

  bad_counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("bad", 1:4))
  )
  bad <- suppressWarnings(SeuratObject::CreateSeuratObject(bad_counts))
  bad$col <- c(0, 1, NA_real_, 3)
  bad$row <- c(0, 1, 2, 3)
  expect_error(
    scop:::cytospace_get_spatial_coords(bad, colnames(bad)),
    "non-finite"
  )
})

test_that("CytoSPACE composition estimation matches pairwise correlations", {
  st_expr <- rbind(
    c(4, 1, 2, 5),
    c(1, 3, 2, 1),
    c(2, 2, 4, 1),
    c(0, 1, 0, 2),
    c(3, 2, 1, 4)
  )
  ref_expr <- rbind(
    c(4, 1, 2, 5, 3, 1),
    c(1, 3, 2, 1, 2, 4),
    c(2, 2, 4, 1, 3, 2),
    c(0, 1, 0, 2, 1, 0),
    c(3, 2, 1, 4, 2, 3)
  )
  colnames(st_expr) <- paste0("Spot", seq_len(ncol(st_expr)))
  labels <- c("A", "A", "B", "B", "C", "C")
  cell_types <- unique(labels)
  n_cells <- c(2, 3, 2, 1)

  norm_st <- scop:::cytospace_log_cpm_r(st_expr)
  norm_ref <- scop:::cytospace_log_cpm_r(ref_expr)
  profiles <- vapply(cell_types, function(ct) {
    rowMeans(norm_ref[, labels == ct, drop = FALSE])
  }, numeric(nrow(norm_ref)))
  reference_scores <- vapply(seq_len(ncol(norm_st)), function(spot) {
    vapply(seq_along(cell_types), function(ct) {
      value <- suppressWarnings(stats::cor(norm_st[, spot], profiles[, ct]))
      if (is.finite(value)) max(value, 0) else 0
    }, numeric(1))
  }, numeric(length(cell_types)))
  reference_scores <- t(reference_scores)
  row_sums <- rowSums(reference_scores)
  reference_scores[row_sums <= 0, ] <- 1 / ncol(reference_scores)
  reference_scores <- sweep(reference_scores, 1, rowSums(reference_scores), "/")
  expected <- colSums(reference_scores * n_cells) / sum(n_cells)
  names(expected) <- cell_types

  expect_equal(
    scop:::cytospace_estimate_fractions(st_expr, ref_expr, labels, cell_types, n_cells),
    expected,
    tolerance = 1e-12
  )
})

test_that("SpatialCoordinates makes multi-image selection explicit", {
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 1), row.names = c("spot1", "spot2")),
    type = "centroids", assay = assay, key = "s1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(2, 3), row.names = c("spot3", "spot4")),
    type = "centroids", assay = assay, key = "s2_"
  )
  expect_error(SpatialCoordinates(srt), "Multiple spatial images")
  selected <- SpatialCoordinates(srt, image = "slice2")
  expect_identical(selected$data$cell_id, c("spot3", "spot4"))
  expect_identical(selected$source$image, "slice2")
  expect_error(
    SpatialCoordinates(srt, image_policy = "legacy_first"),
    "should be.*strict"
  )
})

test_that("analysis and plotting never silently select the first spatial image", {
  data("visium_mouse_brain_slices_sub", package = "scop")
  srt <- visium_mouse_brain_slices_sub
  expect_identical(length(SeuratObject::Images(srt)), 2L)
  expect_equal(ncol(srt), 2000L)

  expect_error(scop:::spatial_analysis_coords(srt), "Multiple spatial images")
  expect_error(
    SpatialSpotPlot(srt, group.by = "orig.ident"),
    "Multiple spatial images"
  )

  selected <- scop:::spatial_analysis_coords(srt, image = "anterior2")
  expect_identical(nrow(selected$data), 1000L)
  expect_identical(selected$source$image, "anterior2")
  expect_identical(selected$source$coordinate_space, "legacy_display")
  expect_identical(selected$source$coord.cols, c("x", "y"))
  expect_true(all(c("scale", "y_flip", "raw_x_col", "raw_y_col") %in%
    names(selected$transform)))

  plotted <- SpatialSpotPlot(
    srt,
    group.by = "orig.ident",
    image = "anterior2"
  )
  expect_s3_class(plotted, "ggplot")
  expect_identical(nrow(plotted$data), 1000L)
})

test_that("schema-v1 results support custom keys and legacy read-only views", {
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt@tools$custom_misty <- scop:::spatial_result_build(
    bundle = list(results = list(table = data.frame(value = 1))),
    method = "MistyR",
    result_type = "neighborhood",
    provenance = list(producer = "RunMistyR", backend_id = "mistyr")
  )
  info <- SpatialResultInfo(srt)
  expect_identical(info$tool_name, "custom_misty")
  expect_identical(info$schema_version, 1L)
  expect_identical(GetSpatialResult(srt, method = "RunMistyR")$method, "MistyR")
  expect_identical(GetSpatialResult(srt, tool_name = "custom_misty", raw = TRUE), srt@tools$custom_misty)
  srt@tools$custom_misty_2 <- srt@tools$custom_misty
  expect_error(GetSpatialResult(srt, method = "RunMistyR"), "Multiple spatial results")
  expect_identical(GetSpatialResult(srt, tool_name = "custom_misty_2")$method, "MistyR")
  srt@tools$custom_misty_2 <- NULL

  srt@tools$StatialKontextual <- list(
    table = data.frame(value = 1),
    parameters = list(coordinate_space = "raw", image = "slice1")
  )
  before <- srt@tools$StatialKontextual
  normalized <- GetSpatialResult(srt, tool_name = "StatialKontextual")
  expect_identical(normalized$schema_version, 1L)
  expect_identical(normalized$method, "StatialKontextual")
  expect_identical(srt@tools$StatialKontextual, before)
})

test_that("GetSpatialGraph converts without synchronizing Seurat graphs", {
  skip_if_not_installed("BiocNeighbors")
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 0, 1, 2)
  srt$row <- 0
  out <- RunSpatialNetwork(srt, k = 1, verbose = FALSE)
  graph <- GetSpatialGraph(out)
  expect_true(all(c("nodes", "edges", "parameters", "source") %in% names(graph)))
  adjacency <- GetSpatialGraph(out, format = "sparse")
  expect_s4_class(adjacency, "dgCMatrix")
  expect_equal(adjacency, Matrix::t(adjacency))
  expect_length(SeuratObject::Graphs(out), 0L)
  expect_s4_class(GetSpatialGraph(out, format = "seurat"), "Graph")
  expect_error(GetSpatialGraph(out, format = "sparse", value = "distance"), "Zero-distance")
  graph_info <- SpatialResultInfo(out, detail = "graphs")
  expect_identical(graph_info$graph_name, "knn_k1")
  expect_true(graph_info$active)
})

test_that("boundary validator preserves polygon and ring order", {
  boundaries <- data.frame(
    cell_id = rep("cell1", 8),
    polygon_id = "p1",
    ring_id = rep(c("outer", "hole"), each = 4),
    vertex_order = rep(c(3, 1, 4, 2), 2),
    x = c(1, 0, 0, 1, 0.8, 0.2, 0.2, 0.8),
    y = c(1, 0, 1, 0, 0.8, 0.2, 0.8, 0.2)
  )
  valid <- scop:::spatial_boundary_validate(boundaries, image = "slice1")
  expect_identical(unique(valid$image), "slice1")
  expect_true(all(diff(valid$vertex_order[valid$ring_id == "outer"]) >= 0))
  expect_error(scop:::spatial_boundary_validate(boundaries[, setdiff(names(boundaries), "x")]), "required column")
})
