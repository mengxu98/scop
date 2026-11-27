ccc_layout_fixture <- function() {
  cells <- paste0("cell", seq_len(6))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(matrix(
    seq_len(12), nrow = 2L, dimnames = list(c("G1", "G2"), cells)
  )))
  groups <- c("Ductal", "Ngn3-high-EP", "Pre-endocrine")
  object@tools$MDIC3 <- list(long_table = data.frame(
    sender = factor(groups[c(1, 2, 3, 1)], levels = groups),
    receiver = factor(groups[c(2, 3, 2, 1)], levels = groups),
    ligand = c("L1", "L2", "L3", "L4"),
    receptor = c("R1", "R2", "R3", "R4"),
    interaction_name = c("L1_R1", "L2_R2", "L3_R3", "L4_R4"),
    pathway_name = c("P1", "P1", "P2", "P2"),
    pvalue = 0.01, score = c(.3, .4, .2, .8)
  ))
  object
}

test_that("pathway hierarchy preserves sender/receiver scores and shared positions", {
  x <- ccc_layout_fixture()
  original <- x@tools
  p <- CCCNetworkPlot(
    x, method = "MDIC3", plot_type = "pathway",
    layout = "hierarchy", signaling = c("P1", "P2")
  )
  expect_length(p, 2L)
  expect_s3_class(p[[1]], "ggplot")
  expect_equal(sort(p[[1]]$data$weight), c(.3, .4))
  expect_equal(sort(p[[2]]$data$weight), c(.2, .8))
  expect_equal(p[[2]]$data$weight[p[[2]]$data$sender == p[[2]]$data$receiver], .8)
  expect_identical(p[[1]]$layers[[2]]$data[c("x", "y", "color")],
    p[[2]]$layers[[2]]$data[c("x", "y", "color")])
  expect_identical(p[[1]]$scales$get_scales("linewidth")$limits, c(0, .8))
  expect_identical(p[[2]]$scales$get_scales("linewidth")$limits, c(0, .8))
  expect_false(is.null(p[[1]]$layers[[1]]$geom_params$arrow))

  filtered <- CCCNetworkPlot(
    x, method = "MDIC3", plot_type = "pathway",
    layout = "hierarchy", signaling = "P1", receiver.use = "Pre-endocrine",
    directed = FALSE, link_curvature = 0
  )
  expect_identical(filtered$data$receiver, "Pre-endocrine")
  expect_equal(filtered$data$weight, .4)
  expect_null(filtered$layers[[1]]$geom_params$arrow)
  expect_identical(x@tools, original)
})

test_that("aggregate heatmaps measure the native layout without leaking a device", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  x <- ccc_layout_fixture()
  dev <- grDevices::dev.cur()
  p <- CCCHeatmap(x, method = "MDIC3", title = "Aggregate communication")
  expect_identical(grDevices::dev.cur(), dev)
  size <- attr(p, "size")
  expect_s3_class(p, "ggplot")
  expect_gt(size$width, 0)
  expect_gt(size$height, 0)
  expect_lt(size$height, 5)
  manual <- CCCHeatmap(
    x, method = "MDIC3", width = 2,
    x_text_angle = 45, legend.position = "none", top_anno = NULL, left_anno = NULL
  )
  expect_identical(grDevices::dev.cur(), dev)
  expect_gt(attr(manual, "size")$width, 2)
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)
  ggplot2::ggsave(f, p, width = size$width, height = size$height, units = size$units)
  expect_gt(file.info(f)$size, 4000)
})

test_that("pathway networks default to circle without removing hierarchy", {
  x <- ccc_layout_fixture()
  args <- list(object = x, method = "MDIC3", plot_type = "pathway", signaling = "P1")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grDevices::dev.control(displaylist = "enable")
  expect_s3_class(do.call(CCCNetworkPlot, args), "recordedplot")
  expect_s3_class(do.call(CCCNetworkPlot, c(args, list(layout = "circle"))), "recordedplot")
  expect_s3_class(do.call(CCCNetworkPlot, c(args, list(layout = "hierarchy"))), "ggplot")
})

test_that("aggregate heatmaps keep bars and cell-type strips on separate sides", {
  x <- ccc_layout_fixture()
  local_mocked_bindings(ccc_heatmap_full_plot = function(...) list(...), .package = "scop")
  p <- CCCHeatmap(x, method = "MDIC3")
  expect_identical(p$top_anno, "bar")
  expect_identical(p$left_anno, "bar")
  expect_identical(p$bottom_anno, "cell")
  expect_identical(p$right_anno, "cell")
  expect_identical(p$x_text_angle, 90)
  expect_identical(p$legend.position, "right")
  expect_identical(p$legend.direction, "vertical")
  p <- CCCHeatmap(x, method = "MDIC3", top_anno = c("bar", "cell"),
    right_anno = NULL, x_text_angle = 0, legend.position = "bottom")
  expect_identical(p$top_anno, c("bar", "cell"))
  expect_length(p$right_anno, 0)
  expect_identical(p$x_text_angle, 0)
  expect_identical(p$legend.position, "bottom")
  expect_identical(p$legend.direction, "horizontal")
})

test_that("CellChat bubbles preserve raw values and the significance boundaries", {
  x <- ccc_layout_fixture()
  x@tools$CellChat <- list(results = list(ALL = list()))
  records <- data.frame(
    source = "Ductal",
    target = "Ngn3-high-EP",
    interaction_name = paste0("L", seq_len(6), "-R"),
    interaction_name_2 = paste0("L", seq_len(6), "  - R"),
    prob = c(.2, .2, .1, .1, .05, .05),
    pval = c(0, .01, .010001, .05, .08, NA)
  )
  local_mocked_bindings(subset_cc_table = function(...) records, .package = "scop")
  args <- list(object = x, method = "CellChat", plot_type = "bubble", thresh = .1)
  p <- do.call(CCCHeatmap, args)
  built <- ggplot2::ggplot_build(p)
  layer <- built$data[[1L]]
  expect_identical(p$data[names(records)], records)
  expect_identical(as.integer(p$data$pvalue_band), c(1L, 1L, 2L, 2L, 3L, NA))
  expect_equal(layer$size[c(1L, 3L, 5L, 6L)], c(5.4, 4.1, 2.8, 1.5))
  expect_true(all(layer$shape == 21))
  expect_identical(layer$fill[1L], layer$fill[2L])
  expect_identical(as.character(built$layout$panel_params[[1L]]$x$get_labels()), paste0("L", seq_len(6), " - R"))
  expect_identical(built$layout$panel_params[[1L]]$y$get_labels(), "Ductal \u2192 Ngn3-high-EP")
  expect_equal(ggplot2::calc_element("axis.text.x", p$theme)$angle, 45)

  q <- do.call(CCCHeatmap, c(args, list(
    bubble_color.by = "pval", x_text_angle = 35, legend.position = "bottom",
    bubble_size.range = c(1, 6)
  )))
  q_layer <- ggplot2::ggplot_build(q)$data[[1L]]
  expect_identical(q_layer$fill[1L], q_layer$fill[2L])
  expect_length(unique(q_layer$fill[1:5]), 3L)
  expect_equal(max(q_layer$size), 6)
  expect_equal(ggplot2::calc_element("axis.text.x", q$theme)$angle, 35)
  expect_identical(q$theme$legend.position, "bottom")
  expect_silent(ggplot2::ggplotGrob(q))
  expect_identical(do.call(CCCHeatmap, c(args, list(return.data = TRUE))), records)

  records$pval[] <- 0
  p <- do.call(CCCHeatmap, args)
  expect_silent(ggplot2::ggplotGrob(p))
  expect_true(all(ggplot2::ggplot_build(p)$data[[1L]]$size == 5.4))
})

test_that("CellChat comparison conditions share scales and centered right guides", {
  x <- ccc_layout_fixture()
  x@tools$CellChat <- list(comparisons = list(ConditionA_vs_ConditionB = list(
    object.list = list(ConditionA = list(), ConditionB = list())
  )))
  records <- data.frame(
    source = "Ductal", target = "Ngn3-high-EP",
    interaction_name = c("L1-R1", "L2-R2"),
    prob = c(.2, .1), pval = c(0, .02)
  )
  local_mocked_bindings(
    subset_cc_table = function(object, ..., dataset) {
      transform(records, dataset = dataset, prob = prob * if (dataset == "ConditionA") 1 else .5)
    },
    .package = "scop"
  )
  args <- list(
    object = x, method = "CellChat", plot_type = "bubble",
    condition = "ConditionA_vs_ConditionB", comparison = c(2, 1)
  )
  p <- do.call(CCCHeatmap, args)
  expect_identical(ggplot2::ggplot_build(p)$layout$layout$dataset, factor(
    c("ConditionB", "ConditionA"),
    levels = c("ConditionB", "ConditionA")
  ))
  expect_identical(p$theme$legend.position, "right")
  expect_identical(p$theme$legend.justification.right, c(0, .5))
  expect_identical(p$theme$strip.text.y.right$colour, "black")
  expect_s3_class(p$theme$strip.background, "element_blank")
  g <- ggplot2::ggplotGrob(p)
  guide <- which(g$layout$name == "guide-box-right")
  panels <- grepl("^panel", g$layout$name)
  expect_equal(g$layout$t[guide], min(g$layout$t[panels]))
  expect_equal(g$layout$b[guide], max(g$layout$b[panels]))
  expect_equal(sum(g$grobs[[guide]]$layout$name == "guides"), 2L)
  size <- attr(p, "size", exact = TRUE)
  expect_identical(size$units, "mm")
  expect_true(all(c(size$width, size$height) > 0))
  result <- do.call(CCCHeatmap, c(args, list(return.data = TRUE)))
  expect_identical(names(result), c("ConditionB", "ConditionA"))
  expect_identical(result$ConditionA$prob, records$prob)
  expect_identical(result$ConditionB$prob, records$prob * .5)
})

test_that("dense pathway summaries retain ranks on a shared numeric axis", {
  x <- ccc_layout_fixture()
  records <- x@tools$MDIC3$long_table
  records <- records[rep(1L, 20L), ]
  records$pathway_name <- paste0("Pathway ", seq_len(20L))
  records$score <- seq_len(20L)
  p <- ccc_pathway_summary_plot(records)
  expect_equal(p$data$total_strength, 20:1)
  expect_equal(as.integer(table(p$data$section)), c(7, 7, 6))
  built <- ggplot2::ggplot_build(p)
  ranges <- lapply(built$layout$panel_params, function(panel) panel$x.range)
  expect_identical(ranges, rep(ranges[1L], 3L))
  expect_identical(p$labels$x, "Total pathway communication strength")
  sparse <- ccc_pathway_summary_plot(records, top_n = 8, value = "count")
  expect_false("section" %in% names(sparse$data))
  expect_identical(sparse$labels$x, "Significant interactions")
})

test_that("comparison bar labels have space and legends are shown by default", {
  x <- ccc_layout_fixture()
  records <- x@tools$MDIC3$long_table
  records$dataset <- rep(c("Group1", "Group2"), each = 2)
  x@tools$MDIC3$long_table <- records
  args <- list(object = x, method = "MDIC3", plot_type = "comparison",
    comparison = c("Group1", "Group2"), sample_col = "dataset")
  p <- do.call(CCCStatPlot, args)
  built <- ggplot2::ggplot_build(p)
  expect_identical(p$theme$legend.position, "right")
  expect_gt(max(built$layout$panel_params[[1L]]$y.range), max(built$data[[2L]]$y))
  expect_identical(p$labels$y, "Number of inferred interactions")
  q <- do.call(CCCStatPlot, c(args, list(legend.position = "none")))
  expect_identical(q$theme$legend.position, "none")
  expect_identical(p$data, q$data)
})

test_that("interaction dot plots wrap senders and avoid duplicate score guides", {
  x <- ccc_layout_fixture()
  x@tools$MDIC3$long_table$sender <- paste0("Sender", seq_len(4L))
  args <- list(object = x, method = "MDIC3", plot_type = "dot", display_by = "interaction")
  p <- do.call(CCCHeatmap, args)
  expect_identical(p$facet$params$ncol, 2L)
  expect_setequal(as.character(p$data$source), paste0("Sender", seq_len(4L)))
  expect_identical(p$scales$get_scales("size")$guide, "none")
  expect_identical(p$theme$legend.justification.right, c(0, .5))

  q <- do.call(CCCHeatmap, c(args, list(
    color.by = "pvalue", x_text_angle = 90, legend.position = "bottom", nrows = 1
  )))
  expect_equal(q$facet$params$nrow, 1)
  expect_equal(ggplot2::calc_element("axis.text.x", q$theme)$angle, 90)
  expect_identical(q$theme$legend.direction, "horizontal")
  g <- ggplot2::ggplotGrob(q)
  guide <- g$grobs[[which(g$layout$name == "guide-box-bottom")]]
  expect_equal(sum(guide$layout$name == "guides"), 2L)
})

test_that("sample dot plots align conditions with compact default panels", {
  x <- ccc_layout_fixture()
  df <- x@tools$MDIC3$long_table
  df <- rbind(
    transform(df, dataset = "ConditionA"),
    transform(df, dataset = "ConditionB", score = score / 2)
  )
  x@tools$MDIC3$long_table <- df
  p <- CCCHeatmap(x, method = "MDIC3", plot_type = "sample_dot", sample_key = "dataset", top_n = 3)
  expect_equal(p$facet$params$ncol, 1)
  expect_equal(nrow(p$data), 6L)
  expect_setequal(p$data$sample, c("ConditionA", "ConditionB"))
  expect_equal(ggplot2::calc_element("axis.text.x", p$theme)$angle, 0)
  expect_equal(ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$x$get_labels(),
    c("ConditionA", "ConditionB"))
})

test_that("role and difference heatmaps report rendered export sizes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  x <- ccc_layout_fixture()
  a <- list(scaled = matrix(c(1, .4, NA, 1), 2L,
    dimnames = list(c("WNT", "NOTCH"), c("Ductal", "Endocrine"))
  ))
  b <- a
  b$scaled[2, 1] <- .7
  x@tools$CellChat <- list(
    results = list(ConditionA = list(cellchat_object = a)),
    comparisons = list(ConditionA_vs_ConditionB = list(object.list = list(ConditionA = a, ConditionB = b)))
  )
  local_mocked_bindings(
    ccc_cellchat_role_matrix = function(object, ...) object,
    .package = "scop"
  )
  dev <- grDevices::dev.cur()
  p <- CCCHeatmap(x, method = "CellChat", condition = "ConditionA",
    plot_type = "role_heatmap", signaling = c("WNT", "NOTCH"), width = 1.2
  )
  pair <- CCCHeatmap(x, method = "CellChat", condition = "ConditionA_vs_ConditionB",
    plot_type = "role_heatmap", signaling = c("WNT", "NOTCH"), width = 1.2
  )
  diff <- CCCHeatmap(x, method = "CellChat", condition = "ConditionA_vs_ConditionB",
    plot_type = "diff_heatmap", signaling = c("WNT", "NOTCH"), legend.position = "none"
  )
  expect_identical(grDevices::dev.cur(), dev)
  expect_gt(attr(p, "size")$width, 1.2)
  expect_gt(attr(pair, "size")$width, 2.4)
  expect_gt(attr(pair, "size")$height, 0)
  expect_gt(attr(diff, "size")$height, 0)
})

test_that("CellChat pathway rankings stay in one panel and size to the rows", {
  x <- ccc_layout_fixture()
  paths <- paste("Pathway", seq_len(29L))
  x@tools$CellChat <- list(comparisons = list(ConditionA_vs_ConditionB = list(
    object.list = list(ConditionA = 1, ConditionB = 2)
  )))
  local_mocked_bindings(
    ccc_cellchat_pathway_summary_df = function(object, dataset, signaling, ...) {
      df <- data.frame(dataset = dataset, pathway = paths, value = object * seq_along(paths))
      if (is.null(signaling)) df else df[df$pathway %in% signaling, , drop = FALSE]
    },
    .package = "scop"
  )
  args <- list(object = x, method = "CellChat", plot_type = "ranknet", condition = "ConditionA_vs_ConditionB")
  p <- do.call(CCCStatPlot, args)
  small <- do.call(CCCStatPlot, c(args, list(signaling = paths[seq_len(8L)])))
  built <- ggplot2::ggplot_build(p)
  expect_length(built$layout$panel_params, 1L)
  expect_identical(as.character(rev(built$layout$panel_params[[1L]]$y$get_labels())), rev(paths))
  expect_identical(p$labels$x, "Information flow")
  expect_gt(p$theme$aspect.ratio, 3 * small$theme$aspect.ratio)
  expect_identical(p$theme$axis.text$size, small$theme$axis.text$size)
  raw <- do.call(CCCStatPlot, c(args, list(return.data = TRUE)))
  expect_identical(raw$pathway, rev(paths))
  expect_equal(raw$diff, 29:1)
})

test_that("ligand-receptor contribution plots preserve percentages and a common axis", {
  x <- ccc_layout_fixture()
  p <- x@tools$MDIC3$long_table
  rownames(p) <- NULL
  p$score <- c(.75, .25, .6, .4)
  p$pathway_name <- rep(c("P1", "P2"), each = 2L)
  x@tools$CellChat <- list(results = list(ALL = list(cellchat_object = p)))
  local_mocked_bindings(
    subset_cc_table = function(object, signaling, ...) {
      object[object$pathway_name %in% signaling, , drop = FALSE]
    },
    .package = "scop"
  )
  args <- list(object = x, method = "CellChat", plot_type = "lr_contribution", signaling = c("P1", "P2"))
  raw <- do.call(CCCStatPlot, c(args, list(return.data = TRUE)))
  expect_equal(raw$contribution, c(.75, .25, .6, .4))
  p <- do.call(CCCStatPlot, c(args, list(combine = FALSE)))
  built <- lapply(p, ggplot2::ggplot_build)
  expect_identical(built[[1L]]$layout$panel_params[[1L]]$x.range,
    built[[2L]]$layout$panel_params[[1L]]$x.range)
  expect_identical(built[[1L]]$data[[2L]]$label, c("75.0%", "25.0%"))
  expect_identical(p[[1L]]$labels$x, "Contribution share")
  limited <- do.call(CCCStatPlot, c(args, list(top_n = 1, return.data = TRUE)))
  expect_equal(limited$contribution, c(1, 1))
})

test_that("communication distributions keep their scores and explicit label settings", {
  x <- ccc_layout_fixture()
  args <- list(object = x, method = "MDIC3", plot_type = "box", facet_by = "sender")
  p <- do.call(CCCStatPlot, args)
  expect_equal(sort(p$data$score), sort(x@tools$MDIC3$long_table$score))
  expect_identical(p$labels$y, "Communication score")
  expect_equal(ggplot2::calc_element("axis.text.x", p$theme)$angle, 45)
  q <- do.call(CCCStatPlot, c(args, list(x_text_angle = 0, legend.position = "none")))
  expect_equal(ggplot2::calc_element("axis.text.x", q$theme)$angle, 0)
  expect_identical(q$theme$legend.position, "none")
})

test_that("bipartite cell labels sit outside unchanged edge endpoints", {
  x <- ccc_layout_fixture()
  p <- CCCNetworkPlot(x, method = "MDIC3", plot_type = "bipartite", legend.position = "none")
  expect_s3_class(p, "ggplot")
  edges <- p$layers[[1L]]$data
  expect_true(nrow(edges) > 0L)
  left <- p$layers[[6L]]$data
  right <- p$layers[[8L]]$data
  expect_true(all(left$x_label < left$x))
  expect_true(all(right$x_label > right$x))
  expect_identical(p$theme$legend.position, "none")
})

test_that("dense role-change panels label only the largest points without removing marks", {
  x <- ccc_layout_fixture()
  paths <- paste0("P", seq_len(16L))
  base <- matrix(seq_len(16L), ncol = 1L, dimnames = list(paths, "Ductal"))
  x@tools$CellChat <- list(comparisons = list(ConditionA_vs_ConditionB = list(
    object.list = list(ConditionA = base, ConditionB = 2 * base)
  )))
  local_mocked_bindings(
    ccc_cellchat_role_components = function(object, signaling = NULL) {
      x <- if (is.null(signaling)) object else object[signaling, , drop = FALSE]
      list(outgoing = x, incoming = x)
    },
    .package = "scop"
  )
  args <- list(object = x, method = "CellChat", plot_type = "role_change",
    condition = "ConditionA_vs_ConditionB", idents.use = "Ductal")
  p <- do.call(CCCStatPlot, args)
  expect_identical(p$data, do.call(CCCStatPlot, c(args, list(return.data = TRUE))))
  expect_identical(p$data$signaling, paths)
  expect_identical(p$layers[[4L]]$data$signaling, rev(tail(paths, 8L)))
  expect_s3_class(p$scales$get_scales("fill")$guide, "GuideLegend")
  expect_identical(p$labels$x, "Outgoing change")
  expect_identical(p$labels$subtitle, "ConditionB - ConditionA")
})

test_that("role scatter keeps both keys independently of direct labeling", {
  x <- ccc_layout_fixture()
  args <- list(object = x, method = "MDIC3", plot_type = "role_scatter")
  p <- do.call(CCCStatPlot, args)
  no_labels <- do.call(CCCStatPlot, c(args, list(do.label = FALSE, xlabel = "Custom")))
  hidden <- do.call(CCCStatPlot, c(args, list(do.label = FALSE, show.legend = FALSE)))
  expect_identical(p$data, no_labels$data)
  g <- ggplot2::ggplotGrob(p)
  right <- g$grobs[[which(g$layout$name == "guide-box-right")]]
  expect_equal(sum(right$layout$name == "guides"), 2L)
  expect_s3_class(no_labels$scales$get_scales("fill")$guide, "GuideLegend")
  expect_identical(hidden$theme$legend.position, "none")
  expect_identical(no_labels$labels$x, "Custom")
})

test_that("straight and sigmoid flows put stage headings above the graph", {
  x <- ccc_layout_fixture()
  for (type in c("arrow", "sigmoid")) {
    p <- CCCNetworkPlot(x, method = "MDIC3", plot_type = type, legend.position = "none")
    expect_s3_class(p, "ggplot")
    expect_identical(p$scales$get_scales("x")$position, "top")
    nodes <- p$layers[[3L]]$data
    expect_gt(nrow(nodes), 0L)
  }
})

test_that("communication gene plots forward grouping and stacking exactly once", {
  x <- ccc_layout_fixture()
  x$cell_type <- rep(c("A", "B"), each = 3L)
  x@tools$MDIC3$long_table$ligand <- "G1"
  x@tools$MDIC3$long_table$receptor <- "G2"
  seen <- list()
  local_mocked_bindings(
    FeatureStatPlot = function(object, stat.by, group.by, stack, ...) {
      seen[[length(seen) + 1L]] <<- list(features = stat.by, grouping = group.by, stacked = stack)
      ggplot2::ggplot()
    },
    .package = "scop"
  )
  CCCStatPlot(
    x, method = "MDIC3", plot_type = "gene",
    group.by = "cell_type", stack = FALSE
  )
  expect_true(length(seen) > 0L)
  expect_setequal(unlist(lapply(seen, `[[`, "features")), c("G1", "G2"))
  expect_true(all(vapply(seen, function(x) identical(x$grouping, "cell_type"), logical(1))))
  expect_false(any(vapply(seen, `[[`, logical(1), "stacked")))
})

test_that("Sankey export reserves space for separate stage legends", {
  x <- ccc_layout_fixture()
  p <- CCCStatPlot(x, method = "MDIC3", plot_type = "sankey", edge_value = "count")
  hidden <- CCCStatPlot(
    x, method = "MDIC3", plot_type = "sankey",
    edge_value = "count", legend.position = "none"
  )
  expect_s3_class(p, "patchwork")
  expect_gt(attr(p, "size")$width, attr(hidden, "size")$width)
})
