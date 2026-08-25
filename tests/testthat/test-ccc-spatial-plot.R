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
      coordinates, composition_matrix, source = list(image = "slice1")
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
    matrix(seq_len(2 * length(ids)), nrow = 2, ncol = length(ids),
      dimnames = list(c("g1", "g2"), ids))
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

test_that("spatial payloads normalize labels and composition without backend calls", {
  ids <- paste0("spot", 1:2)
  coords <- data.frame(cell_id = ids, x = 1:2, y = 1:2)
  labels <- c(spot2 = "B", spot1 = "A")
  payload <- ccc_spatial_plot_payload_from_labels(coords, labels)
  expect_identical(unname(payload$labels), c("A", "B"))
  expect_identical(payload$group_levels, c("A", "B"))

  composition <- matrix(
    c(2, 1, 1, 3), nrow = 2, byrow = TRUE,
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
