test_that("DynamicHeatmap cell union matches the legacy per-cell NA scan", {
  cells <- paste0("Cell", 1:8)
  metadata <- data.frame(
    lineage_a = c(0.1, NA, 0.3, NA, NA, 0.6, NA, 0.8),
    lineage_b = c(NA, 0.2, 0.4, NA, NA, NA, 0.7, 0.9),
    other = seq_along(cells),
    row.names = cells
  )
  legacy <- unique(cells[apply(
    metadata[, c("lineage_a", "lineage_b"), drop = FALSE],
    1,
    function(x) !all(is.na(x))
  )])

  expect_identical(
    dynamic_heatmap_cell_union(cells, metadata, c("lineage_a", "lineage_b")),
    legacy
  )
})
