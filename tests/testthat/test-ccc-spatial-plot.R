make_ccc_spatial_plot_table <- function() {
  data.frame(
    sender = c("A", "A", "B", "C"),
    receiver = c("B", "B", "A", "A"),
    ligand = c("L1", "L2", "L3", "L4"),
    receptor = c("R1", "R2", "R3", "R4"),
    interaction_name = c("L1_R1", "L2_R2", "L3_R3", "L4_R4"),
    pathway_name = rep("SPP1", 4),
    classification = rep("SPP1", 4),
    score = c(0.8, 0.2, 0.5, 0.3),
    pvalue = c(0.01, 0.02, 0.03, 0.04),
    method = "SpatialCellChat",
    sample = "slice1",
    modality = "spatial",
    spatially_constrained = TRUE,
    stringsAsFactors = FALSE
  )
}

make_ccc_spatial_plot_object <- function(composition = FALSE, two_samples = FALSE) {
  ids <- paste0("spot", seq_len(6))
  coordinates <- data.frame(
    cell_id = ids,
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    x_raw = c(0, 1, 2, 0, 1, 2),
    y_raw = c(0, 0, 0, 1, 1, 1),
    x_display = c(0, 1, 2, 0, 1, 2),
    y_display = c(0, 0, 0, 1, 1, 1),
    image = "slice1",
    stringsAsFactors = FALSE
  )
  if (isTRUE(composition)) {
    composition_matrix <- matrix(
      c(1, 0, 0, 0.5, 0.5, 0, 0, 1, 0, 0, 0.5, 0.5, 0, 0, 1, 0.5, 0, 0.5),
      nrow = length(ids), byrow = TRUE,
      dimnames = list(ids, c("A", "B", "C"))
    )
    payload <- ccc_spatial_plot_payload_from_composition(
      coordinates, composition_matrix,
      source = list(image = "slice1")
    )
  } else {
    payload <- ccc_spatial_plot_payload_from_labels(
      coordinates, c("A", "A", "B", "B", "C", "C"),
      source = list(image = "slice1")
    )
  }
  entry <- spatial_tag_coordinate_contract(list(
    interactions = make_ccc_spatial_plot_table(),
    coordinates = coordinates,
    spatial_plot = payload,
    source = payload$source,
    diagnostics = list(interpretation = "spot/domain-level communication")
  ))
  samples <- if (isTRUE(two_samples)) list(slice1 = entry, slice2 = entry) else list(slice1 = entry)
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(
    matrix(seq_len(2 * length(ids)),
      nrow = 2, ncol = length(ids),
      dimnames = list(c("g1", "g2"), ids)
    )
  ))
  srt@tools$SpatialCellChat <- spatial_tag_coordinate_contract(list(
    method = "SpatialCellChat",
    results = list(default = samples),
    active_result = "default",
    long_table = make_ccc_spatial_plot_table(),
    parameters = list()
  ))
  srt
}

test_that("CCC entrypoints preserve common interaction filters and aggregate values", {
  object <- make_ccc_spatial_plot_object()
  original <- object
  seen <- list()
  testthat::local_mocked_bindings(
    ccc_stat_bar_plot = function(df, pair_df, ...) {
      seen$stat <<- list(long = df, pair = pair_df)
      ggplot2::ggplot()
    },
    ccc_circle_plot = function(pair_df, interaction_df, ...) {
      seen$network <<- list(pair = pair_df, interaction = interaction_df)
      ggplot2::ggplot()
    },
    ccc_interaction_dot_plot = function(df, ...) {
      seen$dot <<- df
      ggplot2::ggplot()
    }, .package = "scop"
  )
  args <- list(
    object = object, method = "SpatialCellChat", sender.use = "A",
    receiver.use = "B", pairLR.use = "L2_R2"
  )
  do.call(CCCStatPlot, c(args, list(plot_type = "bar")))
  do.call(CCCNetworkPlot, c(args, list(plot_type = "circle")))
  do.call(CCCHeatmap, c(args, list(plot_type = "dot", display_by = "interaction")))
  expect_identical(seen$stat$long$ligand, "L2")
  expect_equal(seen$stat$long$score, 0.2)
  expect_true(seen$stat$long$significant)
  expect_equal(seen$stat$pair$sum, 0.2)
  expect_identical(seen$stat$pair, seen$network$pair)
  expect_identical(seen$network$interaction, seen$dot)
  expect_identical(object, original)
})

test_that("CCC heatmap retains focused thresholds and pathway auto-selection", {
  object <- make_ccc_spatial_plot_object()
  object@tools$SpatialCellChat$long_table$pathway_name <- c("P1", "P1", "P2", "P2")
  object@tools$SpatialCellChat$long_table$classification <- c("P1", "P1", "P2", "P2")
  original <- object
  seen <- NULL
  testthat::local_mocked_bindings(
    ccc_heatmap_full_plot = function(df, pair_df, interaction_df, ...) {
      seen <<- list(long = df, pair = pair_df, interaction = interaction_df)
      ggplot2::ggplot()
    }, .package = "scop"
  )
  CCCHeatmap(object, method = "SpatialCellChat", plot_type = "focused_heatmap", min_interaction_threshold = 0.45)
  expect_equal(seen$long$score, c(0.8, 0.5))
  expect_equal(seen$pair$sum, c(0.8, 0.5))
  expect_error(CCCHeatmap(object,
    method = "SpatialCellChat", plot_type = "focused_heatmap",
    min_interaction_threshold = 2
  ), "No communication records remain")
  CCCHeatmap(object, method = "SpatialCellChat", plot_type = "pathway_bubble", top_n = 1L)
  expect_identical(unique(seen$long$pathway_name), "P1")
  expect_equal(seen$long$score, c(0.8, 0.2))
  expect_equal(seen$pair$sum, 1)
  expect_identical(object, original)
})

test_that("CCC plotting aliases preserve explicit choices and reject conflicts", {
  object <- make_ccc_spatial_plot_object()
  for (plotter in list(CCCNetworkPlot, CCCHeatmap, CCCStatPlot)) {
    expect_error(plotter(
      object,
      method = "SpatialCellChat", sender.use = "A", sender_use = "B"
    ), "sender.use.*sender_use|sender_use.*sender.use")
    expect_error(plotter(
      object,
      method = "SpatialCellChat", thresh = 0.01, pvalue_threshold = 0.05
    ), "thresh.*pvalue_threshold|pvalue_threshold.*thresh")
    pairs <- data.frame(ligand = "L1", receptor = "R1")
    swapped <- pairs
    names(swapped) <- c("receptor", "ligand")
    expect_error(plotter(
      object,
      method = "SpatialCellChat", pairLR.use = pairs, pair_lr_use = swapped
    ), "pairLR.use.*pair_lr_use|pair_lr_use.*pairLR.use")
  }
  draw <- function(...) {
    CCCStatPlot(
      object,
      method = "SpatialCellChat", plot_type = "bar", ...
    )
  }
  canonical <- draw(sender.use = "A", thresh = 0.03)
  expect_error(draw(thresh = 0.03, pvalue_threshold = 0.03 + 1e-12), "thresh.*pvalue_threshold|pvalue_threshold.*thresh")
  expect_identical(as.character(canonical$data$group), "A -> B")
  expect_equal(canonical$data$value, 1)
  expect_equal(draw(sender_use = "A", pvalue_threshold = 0.03)$data, canonical$data)
  expect_equal(draw(sender.use = "A", sender_use = "A", thresh = 0.03, pvalue_threshold = 0.03)$data, canonical$data)
  expect_error(draw(sender.use = NULL, sender_use = "B"), "sender.use.*sender_use|sender_use.*sender.use")
  expect_error(draw(sender.use = "A", sender_use = NULL), "sender.use.*sender_use|sender_use.*sender.use")
  expect_error(draw(sample_col = "one", context_col = "two"), "sample_col.*context_col|context_col.*sample_col")
  expect_error(draw(context_col = "one", condition_col = "two"), "context_col.*condition_col|condition_col.*context_col")
  expect_error(draw(context_col = NULL, condition_col = "two"), "context_col.*condition_col|condition_col.*context_col")
  expect_equal(draw(sender_use = NULL)$data, draw()$data)
})

test_that("spatial payloads normalize labels and composition without backend calls", {
  ids <- paste0("spot", 1:2)
  coords <- data.frame(cell_id = ids, x = 1:2, y = 1:2)
  labels <- c(spot2 = "B", spot1 = "A")
  payload <- ccc_spatial_plot_payload_from_labels(coords, labels)
  expect_identical(unname(payload$labels), c("A", "B"))
  expect_identical(payload$group_levels, c("A", "B"))

  composition <- matrix(
    c(2, 1, 1, 3),
    nrow = 2, byrow = TRUE,
    dimnames = list(ids, c("A", "B"))
  )
  payload <- ccc_spatial_plot_payload_from_composition(coords, composition)
  expect_equal(unname(rowSums(payload$composition)), c(1, 1))
  expect_equal(payload$composition["spot1", "A"], 2 / 3)
  expect_equal(payload$source$coordinate_space, "raw")
  expect_equal(payload$source$plot_coordinate_space, "display")
})

test_that("spatial renderer maps abundance, weighted centroids, arrows, and labels", {
  object <- make_ccc_spatial_plot_object()
  stored <- ccc_spatial_stored(object, "SpatialCellChat", "default", "slice1")
  pair_df <- ccc_spatial_network_data(stored, signaling = "SPP1", top_n = 3)
  render <- ccc_spatial_render_data(stored$payload, pair_df)
  expect_equal(render$nodes$abundance, c(2, 2, 2))
  expect_equal(render$nodes$x, c(0.5, 1, 1.5))
  expect_true(all(c("x", "xend", "weight", "sender", "receiver") %in% colnames(render$edges)))

  plot <- CCCNetworkPlot(
    object,
    method = "SpatialCellChat",
    condition = "default",
    sample = "slice1",
    plot_type = "spatial",
    signaling = "SPP1",
    top_n = 3,
    cell_palcolor = c(A = "#d73027", B = "#4575b4", C = "#1a9850")
  )
  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)
  expect_gte(length(built$data), 4L)
  expect_true(any(vapply(built$data, nrow, integer(1)) == 6L))
})

test_that("composition spatial renderer uses pies and proportion-weighted nodes", {
  object <- make_ccc_spatial_plot_object(composition = TRUE)
  stored <- ccc_spatial_stored(object, "SpatialCellChat", "default", "slice1")
  pair_df <- ccc_spatial_network_data(stored, signaling = "SPP1", top_n = 3)
  render <- ccc_spatial_render_data(stored$payload, pair_df)
  expect_equal(render$nodes$abundance, c(2, 2, 2))
  expect_equal(render$nodes$x[render$nodes$group == "A"], 0.75)
  plot <- CCCNetworkPlot(
    object,
    method = "SpatialCellChat",
    plot_type = "spatial",
    signaling = "SPP1",
    composition_display = "pie",
    composition_radius = 0.1
  )
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("spatial CCC selection is strict for samples and unsupported methods", {
  object <- make_ccc_spatial_plot_object(two_samples = TRUE)
  expect_error(
    CCCNetworkPlot(object, method = "SpatialCellChat", plot_type = "spatial", signaling = "SPP1"),
    "select.*sample"
  )
  expect_error(
    CCCNetworkPlot(object, method = "LIANA", plot_type = "spatial", signaling = "SPP1"),
    "not supported"
  )
})


test_that("CCC dot views retain aggregate scores and sample-specific observations", {
  object <- make_ccc_spatial_plot_object(two_samples = TRUE)
  second <- object@tools$SpatialCellChat$long_table
  second$sample <- "slice2"
  second$score <- second$score * 0.5
  object@tools$SpatialCellChat$long_table <- rbind(
    object@tools$SpatialCellChat$long_table, second
  )
  original <- object
  source_target <- CCCHeatmap(
    object,
    method = "SpatialCellChat", plot_type = "dot",
    display_by = "interaction", top_n = 2, verbose = FALSE
  )
  sample <- CCCHeatmap(
    object,
    method = "SpatialCellChat", plot_type = "sample_dot",
    top_n = 2, verbose = FALSE
  )

  expect_equal(source_target$data$score, c(1.2, 0.75))
  expect_equal(source_target$data$plot_metric, source_target$data$score)
  expect_equal(sample$data$score, c(0.8, 0.5, 0.4, 0.25))
  expect_identical(unique(sample$data$sample), c("slice1", "slice2"))
  source_layout <- ggplot2::ggplot_build(source_target)$layout$layout
  sample_layout <- ggplot2::ggplot_build(sample)$layout$layout
  expect_identical(as.character(source_layout$source), c("A", "B"))
  expect_identical(as.character(sample_layout$interaction_label), c("L1 - R1", "L3 - R3"))
  expect_s3_class(ggplot2::ggplotGrob(source_target), "gtable")
  expect_s3_class(ggplot2::ggplotGrob(sample), "gtable")
  expect_identical(object, original)
})
