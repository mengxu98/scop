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

test_that("CellChat bubble labels and guides remain compact without changing records", {
  records <- data.frame(
    source = c("Ductal", "Pre-endocrine"),
    target = c("Ngn3-high-EP", "Ductal"),
    interaction_name = c("MDK-SDC4", "PTN-NCL"),
    prob = c(0.2, 0.1),
    pval = c(0.01, 0.02)
  )
  p <- custom_cc_bubble_plot(records, angle.x = 90)
  expect_equal(p$data$prob, records$prob)
  expect_identical(p$scales$get_scales("size")$guide, "none")
  labels <- ggplot2::ggplot_build(p)$layout$panel_params[[1L]]$x$get_labels()
  expect_identical(gsub("\n", " ", as.character(labels)), levels(p$data$pair))
  q <- custom_cc_bubble_plot(records, angle.x = 35, color.by = "pval", legend.position = "bottom")
  expect_equal(q$data$fill_val, -log10(records$pval))
  expect_false(identical(q$scales$get_scales("size")$guide, "none"))
  expect_equal(ggplot2::calc_element("axis.text.x", q$theme)$angle, 35)
  expect_identical(q$theme$legend.position, "bottom")
})

test_that("CellChat bubbles separate reported zeros from positive p-values", {
  x <- ccc_layout_fixture()
  x@tools$CellChat <- list(results = list(ALL = list()))
  records <- data.frame(
    source = rep("Ductal", 6L),
    target = rep("Ngn3-high-EP", 6L),
    interaction_name = paste0("L", seq_len(6), "-R"),
    prob = c(.2, .2, .1, .1, .05, .05),
    pval = c(0, .01, .02, .03, .04, NA)
  )
  local_mocked_bindings(subset_cc_table = function(...) records, .package = "scop")
  args <- list(object = x, method = "CellChat", plot_type = "bubble")
  p <- do.call(CCCHeatmap, c(args, list(bubble_color.by = "pval")))
  base <- do.call(CCCHeatmap, args)
  p_layer <- ggplot2::ggplot_build(p)$data[[1L]]
  base_layer <- ggplot2::ggplot_build(base)$data[[1L]]
  expect_identical(p$data$pval, records$pval)
  expect_identical(p_layer[c("x", "y", "size")], base_layer[c("x", "y", "size")])
  expect_identical(unname(p_layer$shape), c(23, rep(21, 5)))
  expect_length(unique(p_layer$fill[2:5]), 4)
  expect_equal(p$scales$get_scales("fill")$breaks, -log10(rev(records$pval[2:5])))
  expect_identical(p$scales$get_scales("fill")$labels, c("0.04", "0.03", "0.02", "0.01"))
  expect_identical(base$data$fill_val, records$prob)
  expect_identical(do.call(CCCHeatmap, c(args, list(
    bubble_color.by = "pval", return.data = TRUE
  ))), records)

  records$pval <- rep(0, nrow(records))
  zero <- do.call(CCCHeatmap, c(args, list(bubble_color.by = "pval")))
  expect_silent(ggplot2::ggplotGrob(zero))
  expect_identical(zero$scales$get_scales("fill")$guide, "none")
  expect_true(all(ggplot2::ggplot_build(zero)$data[[1L]]$shape == 23))
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

test_that("comparison bar labels have space and explicit legends are retained", {
  x <- ccc_layout_fixture()
  records <- x@tools$MDIC3$long_table
  records$dataset <- rep(c("Group1", "Group2"), each = 2)
  x@tools$MDIC3$long_table <- records
  args <- list(object = x, method = "MDIC3", plot_type = "comparison",
    comparison = c("Group1", "Group2"), sample_col = "dataset")
  p <- do.call(CCCStatPlot, args)
  built <- ggplot2::ggplot_build(p)
  expect_identical(p$theme$legend.position, "none")
  expect_gt(max(built$layout$panel_params[[1L]]$y.range), max(built$data[[2L]]$y))
  expect_identical(p$labels$y, "Number of inferred interactions")
  q <- do.call(CCCStatPlot, c(args, list(legend.position = "right")))
  expect_identical(q$theme$legend.position, "right")
  expect_identical(p$data, q$data)
})
