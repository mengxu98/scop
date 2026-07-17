make_cell_dim_plot_srt <- function() {
  counts <- matrix(rpois(60, lambda = 5), nrow = 6)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- rep(c("A1", "A2", "B1", "B2"), length.out = 10)
  srt$major_type <- ifelse(grepl("^A", srt$celltype), "A", "B")
  embedding <- cbind(
    c(seq(0, 0.4, length.out = 5), seq(2, 2.4, length.out = 5)),
    c(seq(0, 0.4, length.out = 5), seq(2, 2.4, length.out = 5))
  )
  rownames(embedding) <- colnames(srt)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(srt)
  )
  srt
}

test_that("FeatureDimPlot reuses the expression layer for co-expression", {
  srt <- make_cell_dim_plot_srt()
  counts <- GetAssayData5(srt, assay = "RNA", layer = "counts")
  calls <- 0L
  testthat::local_mocked_bindings(
    GetAssayData5 = function(...) {
      calls <<- calls + 1L
      counts
    },
    CheckDataType = function(...) "raw_counts",
    .package = "scop"
  )

  expect_no_error(
    FeatureDimPlot(
      srt,
      features = c("gene1", "gene2"),
      reduction = "umap",
      layer = "counts",
      calculate_coexp = TRUE,
      force = TRUE,
      verbose = FALSE
    )
  )
  expect_identical(calls, 1L)
})

test_that("FeatureDimPlot3D supports sparse co-expression", {
  srt <- make_cell_dim_plot_srt()
  embedding <- cbind(
    srt[["umap"]]@cell.embeddings,
    UMAP_3 = seq_len(ncol(srt)) / ncol(srt)
  )
  srt[["umap3d"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP3D_",
    assay = SeuratObject::DefaultAssay(srt)
  )

  plot <- FeatureDimPlot3D(
    srt,
    features = c("gene1", "gene2"),
    reduction = "umap3d",
    layer = "counts",
    calculate_coexp = TRUE,
    force = TRUE,
    verbose = FALSE
  )
  expect_s3_class(plot, "plotly")
})

test_that("CellDimPlot supports atlas-style grid and marked groups", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    reduction = "umap",
    add_grid = TRUE,
    grid_n = 4,
    add_mark = TRUE,
    mark_alpha = 0,
    mark_linetype = 2,
    mark_linewidth = 0.3,
    mark_border = "grey30",
    label = TRUE,
    label_insitu = TRUE,
    label_repel = TRUE,
    legend.position = "bottom",
    theme_use = "theme_blank",
    force = TRUE
  )

  built <- ggplot2::ggplot_build(plot)
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(built$data[[1]]), 16)
  expect_true(any(vapply(built$data, nrow, integer(1)) == ncol(srt)))
})

test_that("CellDimPlot keeps synthetic split metadata local", {
  srt <- make_cell_dim_plot_srt()
  original_meta <- srt@meta.data

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = FALSE,
    force = TRUE
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(srt@meta.data, original_meta)
  expect_false("All.groups" %in% colnames(srt@meta.data))
})

test_that("CellDimPlot can color subgroups while labeling major groups", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    label.by = "major_type",
    mark.by = "major_type",
    reduction = "umap",
    add_mark = TRUE,
    mark_alpha = 0,
    label = TRUE,
    label_insitu = TRUE,
    legend.position = "bottom",
    force = TRUE
  )

  built <- ggplot2::ggplot_build(plot)
  label_values <- unlist(
    lapply(built$data, function(layer) as.character(layer[["label"]])),
    use.names = FALSE
  )
  expect_true(all(c("A", "B") %in% label_values))
  expect_false(any(c("A1", "A2", "B1", "B2") %in% label_values))
})

test_that("CellDimPlot can color subgroups with major-group nested legend", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    label.by = "major_type",
    mark.by = "major_type",
    legend.by = "major_type",
    reduction = "umap",
    add_grid = FALSE,
    add_mark = TRUE,
    mark_alpha = 0,
    label = TRUE,
    label_insitu = TRUE,
    legend.position = "bottom",
    force = TRUE
  )

  expect_true(inherits(plot, "patchwork") || inherits(plot, "ggplot"))
})

test_that("CellDimPlot nested bottom legend is horizontal by default", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    legend.by = "major_type",
    reduction = "umap",
    add_mark = TRUE,
    mark_alpha = 0,
    legend.position = "bottom",
    force = TRUE
  )

  expect_true(inherits(plot, "patchwork") || inherits(plot, "ggplot"))
})

test_that("CellDimPlot nested right legend stacks grouped subtypes", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    legend.by = "major_type",
    reduction = "umap",
    add_mark = TRUE,
    mark.by = "major_type",
    mark_alpha = 0,
    mark_palcolor = c(A = "#444444", B = "#999999"),
    legend.position = "right",
    force = TRUE
  )

  expect_true(inherits(plot, "patchwork") || inherits(plot, "ggplot"))
})

test_that("CellDimPlot nested right legend wraps dense groups into columns", {
  parents <- paste0("type", seq_len(12))
  children <- paste0("sample", seq_len(8))
  dat <- expand.grid(
    parent = parents,
    child = children,
    stringsAsFactors = FALSE
  )
  dat[["group.by"]] <- factor(dat[["child"]], levels = children)
  dat[["legend.by"]] <- factor(dat[["parent"]], levels = parents)
  colors <- stats::setNames(grDevices::rainbow(length(children)), children)
  parent_colors <- stats::setNames(grDevices::rainbow(length(parents)), parents)

  legend <- cell_dim_nested_legend_grob(
    dat = dat,
    colors = colors,
    parent_colors = parent_colors,
    title = "parent:",
    legend_position = "right"
  )

  expect_s3_class(legend, "gtable")
  expect_true(as.numeric(grid::convertUnit(attr(legend, "legend_space"), "in")) > 1.85)
})

test_that("CellDimPlot nested legend separates parent and subgroup colors", {
  srt <- make_cell_dim_plot_srt()
  dat <- srt@meta.data
  dat[["group.by"]] <- factor(dat[["celltype"]], levels = c("A1", "A2", "B1", "B2"))
  dat[["legend.by"]] <- factor(dat[["major_type"]], levels = c("A", "B"))
  subgroup_colors <- c(A1 = "#111111", A2 = "#222222", B1 = "#333333", B2 = "#444444")
  parent_colors <- c(A = "#AAAAAA", B = "#BBBBBB")

  legend_data <- cell_dim_nested_legend_data(
    dat = dat,
    colors = subgroup_colors,
    parent_colors = parent_colors
  )

  expect_equal(legend_data[["parent_colors"]], parent_colors)
  expect_equal(legend_data[["nested_items"]][["A"]], c("A1", "A2"))
  expect_equal(legend_data[["nested_items"]][["B"]], c("B1", "B2"))
  expect_equal(unname(subgroup_colors[legend_data[["nested_items"]][["A"]]]), c("#111111", "#222222"))
})

test_that("CellDimPlot aligns named colors to group levels", {
  colors <- cell_dim_palette_colors(
    c("B", "A"),
    palcolor = c(A = "#AAAAAA", B = "#BBBBBB")
  )

  expect_equal(names(colors), c("B", "A"))
  expect_equal(unname(colors), c("#BBBBBB", "#AAAAAA"))
})

test_that("CellDimPlot disables default point stroke for small point sizes", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    reduction = "umap",
    pt.size = 0.0001,
    raster = FALSE,
    force = TRUE
  )

  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    plot$layers
  )
  expect_true(length(point_layers) > 0L)
  point_layer <- point_layers[[length(point_layers)]]
  expect_equal(point_layer$aes_params$size, 0.0001)
  expect_equal(point_layer$aes_params$stroke, 0)
})

test_that("CellDimPlot uses readable default point size", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    srt = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = FALSE,
    force = TRUE
  )

  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    plot$layers
  )
  expect_true(length(point_layers) > 0L)
  point_layer <- point_layers[[length(point_layers)]]
  expect_equal(point_layer$aes_params$size, 1)
  expect_equal(point_layer$aes_params$stroke, 0)
})

test_that("CellDimPlot validates atlas grid density", {
  srt <- make_cell_dim_plot_srt()

  expect_error(
    CellDimPlot(
      srt = srt,
      group.by = "celltype",
      reduction = "umap",
      add_grid = TRUE,
      grid_n = 1,
      force = TRUE
    ),
    "grid_n"
  )
})

test_that("FeatureDimPlot subsets graph before dense conversion without changing values", {
  srt <- make_cell_dim_plot_srt()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  graph <- Matrix::Matrix(0, ncol(srt), ncol(srt), sparse = TRUE)
  dimnames(graph) <- list(colnames(srt), colnames(srt))
  graph["cell1", "cell2"] <- 1
  graph["cell2", "cell1"] <- 1
  graph["cell1", "cell6"] <- 2
  graph["cell6", "cell1"] <- 2
  srt@graphs[["test_graph"]] <- methods::as(graph, "Graph")
  cells_use <- paste0("cell", 1:5)

  expect_equal(
    as_matrix(graph[cells_use, cells_use, drop = FALSE]),
    as_matrix(graph)[cells_use, cells_use]
  )

  plot <- FeatureDimPlot(
    srt,
    features = c("gene1", "gene2"),
    reduction = "umap",
    cells = cells_use,
    compare_features = TRUE,
    graph = "test_graph",
    force = TRUE,
    theme_use = theme_blank
  )

  expect_true(inherits(plot, "patchwork") || inherits(plot, "ggplot"))
})

test_that("FeatureDimPlot disables default point stroke for small point sizes", {
  srt <- make_cell_dim_plot_srt()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  plot <- FeatureDimPlot(
    srt,
    features = "gene1",
    reduction = "umap",
    pt.size = 0.0001,
    raster = FALSE,
    force = TRUE,
    theme_use = theme_blank
  )

  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    plot$layers
  )
  expect_true(length(point_layers) > 0L)
  point_layer <- point_layers[[length(point_layers)]]
  expect_equal(point_layer$aes_params$size, 0.0001)
  expect_equal(point_layer$aes_params$stroke, 0)
})

test_that("FeatureDimPlot uses readable default point size", {
  srt <- make_cell_dim_plot_srt()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)

  plot <- FeatureDimPlot(
    srt,
    features = "gene1",
    reduction = "umap",
    raster = FALSE,
    force = TRUE,
    theme_use = theme_blank
  )

  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    plot$layers
  )
  expect_true(length(point_layers) > 0L)
  point_layer <- point_layers[[length(point_layers)]]
  expect_equal(point_layer$aes_params$size, 1)
  expect_equal(point_layer$aes_params$stroke, 0)
})
