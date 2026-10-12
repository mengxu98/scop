# Exercise Seurat's actual identity operations on outputs of each producer.
expect_spatial_result_lifecycle <- function(object, method, cluster_colname,
                                            explicit_cell = FALSE, integration = FALSE) {
  expected <- as.character(object@meta.data[[cluster_colname]])
  original_tools <- object@tools
  variants <- list(
    full = object,
    subset = subset(object, cells = colnames(object)[c(1, ncol(object))]),
    renamed = SeuratObject::RenameCells(object, new.names = paste0("new_", colnames(object))),
    permuted = SeuratObject::RenameCells(object, new.names = rev(colnames(object)))
  )
  # Rebuild from counts to guarantee a permutation independent of Seurat's
  # version-specific subset ordering; metadata identities travel with cells.
  cells <- rev(colnames(variants$permuted))[c(1, ncol(object))]
  current <- variants$permuted
  variants$reordered <- Seurat::CreateSeuratObject(
    counts = Seurat::GetAssayData(current, layer = "counts")[, cells, drop = FALSE],
    meta.data = current[[]][cells, , drop = FALSE])
  variants$reordered@tools <- current@tools
  permuted_expected <- stats::setNames(expected, rev(colnames(object)))
  expected_values <- list(full = expected,
    subset = expected[c(1, ncol(object))], renamed = expected,
    permuted = expected, reordered = unname(permuted_expected[cells]))
  for (name in names(variants)) {
    variant <- variants[[name]]
    before <- variant@tools
    result <- GetSpatialResult(variant, method)
    labels <- scop:::spatial_result_cluster_labels(result$clusters)
    ids <- if (is.data.frame(result$clusters)) rownames(result$clusters) else names(result$clusters)
    expect_identical(ids, colnames(variant))
    expect_identical(unname(labels), expected_values[[name]])
    expect_identical(unname(labels), unname(as.character(variant@meta.data[[cluster_colname]])))
    if (explicit_cell) expect_identical(result$clusters$cell, colnames(variant))
    if (integration) {
      expect_equal(result$summary$n_cells, ncol(variant))
      expect_equal(sum(result$summary$domains$count), ncol(variant))
      expect_equal(sum(result$summary$samples$count), ncol(variant))
    }
    expect_identical(variant@tools, before)
  }
  expect_identical(object@tools, original_tools)
  invisible(variants)
}
