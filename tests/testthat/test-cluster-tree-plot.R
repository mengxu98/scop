make_cluster_tree_seurat <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0, 4, 1,
      0, 3, 0, 4, 1, 0,
      2, 1, 5, 0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:6))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt[["RNA_snn_res.1"]] <- factor(c("0", "1", "1", "2", "3", "4"))
  srt[["RNA_snn_res.0.2"]] <- factor(c("0", "0", "1", "1", "1", "2"))
  srt[["RNA_snn_res.0.8"]] <- factor(c("0", "0", "3", "1", "1", "2"))
  srt[["SCT_snn_res.0.2"]] <- factor(c("A", "A", "B", "B", "B", "C"))
  srt
}

test_that("ClusterTreePlot detects and sorts Seurat resolution columns", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(srt, prefix = "RNA_snn")

  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1"))
  expect_equal(info$resolution, c(0.2, 0.8, 1))
})

test_that("ClusterTreePlot auto-detection keeps one clustering prefix", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(srt)

  expect_equal(unique(info$prefix), "RNA_snn")
  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1"))
})

test_that("ClusterTreePlot auto-detection respects uppercase default SNN prefix", {
  srt <- make_cluster_tree_seurat()
  srt@meta.data[c("RNA_snn_res.0.2", "RNA_snn_res.0.8", "RNA_snn_res.1")] <- NULL
  srt[["RNA_SNN_res.0.2"]] <- factor(c("0", "0", "1", "1", "1", "2"))
  srt[["RNA_SNN_res.0.8"]] <- factor(c("0", "0", "3", "1", "1", "2"))
  srt[["SCT_snn_res.0.8"]] <- factor(c("A", "A", "C", "B", "B", "C"))
  srt[["SCT_snn_res.1"]] <- factor(c("A", "D", "C", "B", "B", "C"))

  info <- clustertree_resolve_cluster_cols(srt)

  expect_equal(unique(info$prefix), "RNA_SNN")
  expect_equal(info$column, c("RNA_SNN_res.0.2", "RNA_SNN_res.0.8"))
})

test_that("ClusterTreePlot respects manual cluster column order", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(
    srt,
    cluster_cols = c("RNA_snn_res.1", "RNA_snn_res.0.2")
  )

  expect_equal(info$column, c("RNA_snn_res.1", "RNA_snn_res.0.2"))
})

test_that("ClusterTreePlot filters prefix and resolutions", {
  srt <- make_cluster_tree_seurat()

  info <- clustertree_resolve_cluster_cols(
    srt,
    prefix = "RNA_snn",
    resolutions = c(0.2, 1)
  )

  expect_equal(info$column, c("RNA_snn_res.0.2", "RNA_snn_res.1"))
})

test_that("ClusterTreePlot builds node and edge statistics", {
  srt <- make_cluster_tree_seurat()
  info <- clustertree_resolve_cluster_cols(srt, prefix = "RNA_snn")
  tree_data <- clustertree_build_data(srt@meta.data, info)

  node <- tree_data$nodes[tree_data$nodes$node_id == "RNA_snn_res.0.2::1", ]
  expect_equal(node$size, 3)

  edge <- tree_data$all_edges[
    tree_data$all_edges$from_node == "RNA_snn_res.0.2::1" &
      tree_data$all_edges$to_node == "RNA_snn_res.0.8::1",
  ]
  expect_equal(edge$count, 2)
  expect_equal(edge$in_prop, 1)
  expect_equal(edge$out_prop, 2 / 3)

  column_bottoms <- tapply(
    tree_data$nodes$y,
    tree_data$nodes$resolution_index,
    min
  )
  expect_equal(as.numeric(column_bottoms), rep(0, length(column_bottoms)))

  column_tops <- tapply(
    tree_data$nodes$y,
    tree_data$nodes$resolution_index,
    max
  )
  expect_gt(tail(column_tops, 1), column_tops[1])

  first_column_y <- tree_data$nodes$y[
    tree_data$nodes$resolution_index == min(tree_data$nodes$resolution_index)
  ]
  expect_gt(min(diff(sort(first_column_y))), 1)
})

test_that("ClusterTreePlot returns ggplot for cluster tree", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    legend.position = "inside"
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(formals(ClusterTreePlot)$legend.position, "auto")
  expect_identical(plot$theme$legend.position, "inside")
  expect_equal(plot$theme$legend.position.inside, c(0.022, 0.978))
  expect_null(plot$labels$title)
  expect_identical(plot$theme$text$family, "Arial")
  expect_identical(plot$theme$text$colour, "black")
  expect_identical(plot$theme$plot.subtitle$colour, "black")
  expect_identical(plot$theme$axis.title$colour, "black")
  expect_identical(plot$theme$axis.text.x$colour, "black")
  expect_identical(plot$theme$axis.line.x$colour, "black")
  expect_identical(plot$theme$axis.ticks.x$colour, "black")
  expect_identical(plot$theme$legend.title$colour, "black")
  expect_identical(plot$theme$legend.text$colour, "black")
  expect_identical(plot$theme$legend.direction, "horizontal")
  expect_identical(plot$theme$legend.box, "vertical")
  text_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomText"),
    plot$layers
  )
  label_colors <- vapply(
    text_layers,
    function(layer) layer$aes_params$colour,
    character(1)
  )
  expect_true(all(label_colors %in% c("black", "white")))
  expect_s3_class(plot$theme$panel.border, "element_blank")
})

test_that("ClusterTreePlot respects explicit legend position", {
  srt <- make_cluster_tree_seurat()

  for (legend_position in c("inside", "bottom", "right", "none")) {
    plot <- ClusterTreePlot(
      srt,
      prefix = "RNA_snn",
      legend.position = legend_position
    )

    expect_s3_class(plot, "ggplot")
    expect_identical(plot$theme$legend.position, legend_position)
  }
})

test_that("ClusterTreePlot supports all tree directions", {
  srt <- make_cluster_tree_seurat()
  directions <- c(
    "left-to-right",
    "right-to-left",
    "top-to-bottom",
    "bottom-to-top"
  )
  expected_anchors <- list(
    "left-to-right" = c(0.022, 0.978),
    "right-to-left" = c(0.978, 0.978),
    "top-to-bottom" = c(0.978, 0.978),
    "bottom-to-top" = c(0.978, 0.022)
  )

  for (tree_direction in directions) {
    plot <- ClusterTreePlot(
      srt,
      prefix = "RNA_snn",
      legend.position = "inside",
      direction = tree_direction
    )

    expect_s3_class(plot, "ggplot")
    expect_equal(
      plot$theme$legend.position.inside,
      expected_anchors[[tree_direction]]
    )
    if (tree_direction %in% c("top-to-bottom", "bottom-to-top")) {
      expect_s3_class(plot$coordinates, "CoordFlip")
      expect_s3_class(plot$theme$axis.line.x, "element_blank")
      expect_s3_class(plot$theme$axis.line.y, "element_line")
    } else {
      expect_false(inherits(plot$coordinates, "CoordFlip"))
      expect_s3_class(plot$theme$axis.line.x, "element_line")
      expect_s3_class(plot$theme$axis.line.y, "element_blank")
    }
  }
})

test_that("ClusterTreePlot moves auto legends away from left-heavy trees", {
  srt <- make_cluster_tree_seurat()

  auto_plot <- ClusterTreePlot(
    srt,
    cluster_cols = c(
      "RNA_snn_res.1",
      "RNA_snn_res.0.8",
      "RNA_snn_res.0.2"
    )
  )
  forced_inside <- ClusterTreePlot(
    srt,
    cluster_cols = c(
      "RNA_snn_res.1",
      "RNA_snn_res.0.8",
      "RNA_snn_res.0.2"
    ),
    legend.position = "inside"
  )

  expect_identical(auto_plot$theme$legend.position, "right")
  expect_identical(forced_inside$theme$legend.position, "inside")
})

test_that("ClusterTreePlot respects explicit titles and font families", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    title = "My cluster tree",
    family = "Helvetica",
    label.fg = "black"
  )
  text_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomText"),
    plot$layers
  )

  expect_identical(plot$labels$title, "My cluster tree")
  expect_identical(plot$theme$text$family, "Helvetica")
  expect_length(text_layers, 1L)
  expect_identical(text_layers[[1L]]$aes_params$family, "Helvetica")
})

test_that("ClusterTreePlot uses compact non-redundant guides", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn")

  expect_s3_class(
    plot$scales$get_scales("linewidth")$guide,
    "GuideLegend"
  )
  expect_s3_class(plot$scales$get_scales("colour")$guide, "GuideColourbar")
  expect_identical(plot$scales$get_scales("colour")$name, "Incoming share")
  expect_equal(plot$scales$get_scales("colour")$breaks, c(0, 0.5, 1))
  expect_null(plot$scales$get_scales("alpha"))
  expect_identical(plot$scales$get_scales("fill")$guide, "none")
  expect_identical(
    rlang::as_label(plot$layers[[1L]]$mapping$linewidth),
    "count"
  )
  expect_identical(
    rlang::as_label(plot$layers[[1L]]$mapping$colour),
    "in_prop"
  )
  expect_identical(plot$layers[[1L]]$aes_params$alpha, 0.9)
  expect_gte(plot$theme$legend.title$size, 9.5)
  expect_gte(plot$theme$legend.text$size, 8.5)
  expect_gte(as.numeric(plot$theme$legend.key.height), 14)
})

test_that("ClusterTreePlot returns ggplot for a single marker overlay", {
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn", features = "Gene1")

  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$title)
  expect_identical(plot$labels$subtitle, "Gene1")
  expect_s3_class(plot$scales$get_scales("fill")$guide, "GuideColourbar")
})

test_that("ClusterTreePlot accepts numeric metadata feature overlays", {
  srt <- make_cluster_tree_seurat()
  srt$G2M_score <- seq_len(ncol(srt))

  plot <- ClusterTreePlot(srt, prefix = "RNA_snn", features = "G2M_score")
  feature_values <- clustertree_feature_values(
    srt,
    features = "G2M_score",
    cluster_cols = c("RNA_snn_res.0.2", "RNA_snn_res.0.8")
  )

  expect_s3_class(plot, "ggplot")
  expect_named(feature_values, "G2M_score")
})

test_that("ClusterTreePlot combines or lists multiple marker overlays", {
  srt <- make_cluster_tree_seurat()

  plots <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = c("Gene1", "Gene2"),
    combine = FALSE
  )
  expect_type(plots, "list")
  expect_named(plots, c("Gene1", "Gene2"))

  testthat::skip_if_not_installed("patchwork")
  plot <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = c("Gene1", "Gene2"),
    combine = TRUE
  )
  expect_s3_class(plot, "patchwork")
})

test_that("inside multi-feature patchwork retains guide titles", {
  testthat::skip_if_not_installed("patchwork")
  srt <- make_cluster_tree_seurat()

  plot <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = c("Gene1", "Gene2"),
    combine = TRUE,
    family = "sans"
  )
  grob <- patchwork::patchworkGrob(plot)
  grob_labels <- function(x) {
    out <- if (!is.null(x$label)) as.character(x$label) else character()
    if (!is.null(x$children)) {
      for (i in seq_along(x$children)) {
        out <- c(out, grob_labels(x$children[[i]]))
      }
    }
    if (!is.null(x$grobs)) {
      for (i in seq_along(x$grobs)) {
        out <- c(out, grob_labels(x$grobs[[i]]))
      }
    }
    out
  }

  expect_true(all(c(
    "Cells moved",
    "Incoming share",
    "Cluster size",
    "Gene1",
    "Gene2"
  ) %in% grob_labels(grob)))
})

test_that("ClusterTreePlot preserves named feature-list groups", {
  srt <- make_cluster_tree_seurat()

  plots <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    features = list(Beta = c("Gene1", "Gene2"), Alpha = "Gene1"),
    combine = FALSE
  )

  expect_type(plots, "list")
  expect_named(plots, c("Beta", "Alpha"))
})

test_that("ClusterTreePlot aligns feature values to metadata cells", {
  srt <- make_cluster_tree_seurat()
  values <- c(Cell3 = 10, Cell1 = 2, Cell5 = 4)

  feature_nodes <- clustertree_aggregate_feature_values(
    values = values,
    meta_data = srt@meta.data,
    cluster_cols = "RNA_snn_res.0.2"
  )

  expect_equal(
    feature_nodes$feature_value[feature_nodes$node_id == "RNA_snn_res.0.2::0"],
    2
  )
  expect_equal(
    feature_nodes$feature_value[feature_nodes$node_id == "RNA_snn_res.0.2::1"],
    7
  )
})

test_that("cluster-tree compact breaks are finite and bounded", {
  breaks <- clustertree_compact_breaks(3)

  expect_equal(breaks(c(5, 5)), 5)
  expect_length(breaks(c(NA_real_, Inf)), 0L)

  values <- breaks(c(1, 4000))
  expect_lte(length(values), 3L)
  expect_true(all(values >= 1 & values <= 4000))
})

test_that("cluster-tree labels choose readable black or white text", {
  expect_identical(
    clustertree_contrast_text(c(
      "#111111",
      "#FFFFFF",
      "#FEE08B",
      "#F46D43",
      "#63BDA6"
    )),
    c("white", "black", "black", "black", "black")
  )

  gradient_colors <- clustertree_gradient_colors(
    values = c(0, 0.5, 1, NA_real_),
    colors = c("#313695", "#FFFFBF", "#A50026")
  )
  expect_length(gradient_colors, 4L)
  expect_identical(gradient_colors[4L], "grey85")
  expect_true(all(
    clustertree_contrast_text(gradient_colors) %in% c("black", "white")
  ))
})

test_that("cluster-tree anchored positions retain size-aware spacing", {
  positions <- clustertree_anchored_positions(
    sizes = rep(100, 20),
    size_limits = c(100, 100)
  )

  expect_equal(min(positions), 0)
  expect_gt(min(diff(sort(positions))), 1)

  heterogeneous <- clustertree_anchored_positions(
    sizes = c(100, 100, 10),
    size_limits = c(10, 100)
  )
  expect_gt(
    heterogeneous[1] - heterogeneous[2],
    heterogeneous[2] - heterogeneous[3]
  )

  large_column <- clustertree_anchored_positions(
    sizes = c(80, 100),
    size_limits = c(10, 100)
  )
  small_column <- clustertree_anchored_positions(
    sizes = c(10, 20),
    size_limits = c(10, 100)
  )
  expect_gt(diff(range(large_column)), diff(range(small_column)))
})

test_that("ClusterTreePlot preserves the return_data contract", {
  srt <- make_cluster_tree_seurat()

  out <- ClusterTreePlot(
    srt,
    prefix = "RNA_snn",
    return_data = TRUE
  )

  expect_named(
    out,
    c("plot", "nodes", "edges", "all_edges", "cluster_info")
  )
  expect_s3_class(out$plot, "ggplot")
  expect_gte(nrow(out$all_edges), nrow(out$edges))
})

test_that("ClusterTreePlot errors without multi-resolution columns", {
  counts <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(paste0("Gene", 1:2), paste0("Cell", 1:3))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    ClusterTreePlot(srt),
    "No multi-resolution Seurat clustering columns"
  )
})
