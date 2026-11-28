make_cell_dim_plot_srt <- function(n_cells = 10L) {
  counts <- matrix(rpois(6L * n_cells, lambda = 5), nrow = 6)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- rep(c("A1", "A2", "B1", "B2"), length.out = n_cells)
  srt$major_type <- ifelse(grepl("^A", srt$celltype), "A", "B")
  n_left <- ceiling(n_cells / 2)
  n_right <- n_cells - n_left
  embedding <- cbind(
    c(seq(0, 0.4, length.out = n_left), seq(2, 2.4, length.out = n_right)),
    c(seq(0, 0.4, length.out = n_left), seq(2, 2.4, length.out = n_right))
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

test_that("feature dimension plots align partial assays and embeddings by cell ID", {
  object <- make_cell_dim_plot_srt()
  cells <- colnames(object)
  observed <- cells[c(2, 5, 8)]
  counts <- matrix(
    c(2, 3, 4, 5, 6, 7), nrow = 2L,
    dimnames = list(c("adt1", "adt2"), observed)
  )
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(counts[, rev(observed)], sparse = TRUE)
  )
  object$meta_numeric <- seq_along(cells)
  partial <- SeuratObject::Embeddings(object, "umap")[observed, , drop = FALSE]
  colnames(partial) <- c("SUB_1", "SUB_2")
  object[["partial"]] <- SeuratObject::CreateDimReducObject(
    partial, key = "SUB_", assay = "RNA"
  )
  original <- object
  panels <- FeatureDimPlot(
    object, features = c("adt1", "adt2", "meta_numeric", "SUB_1"),
    assay = "ADT", layer = "counts", reduction = "umap",
    calculate_coexp = TRUE, bg_cutoff = -Inf, lower_cutoff = 0, upper_cutoff = 20,
    raster = FALSE, combine = FALSE, verbose = FALSE
  )
  expect_length(panels, 5L)
  values <- lapply(panels, function(panel) panel$data[cells, "value"])
  expected <- rep(NA_real_, length(cells))
  expected[match(observed, cells)] <- counts["adt1", ]
  expect_equal(values[[1L]], expected)
  expected[match(observed, cells)] <- counts["adt2", ]
  expect_equal(values[[2L]], expected)
  expect_equal(values[[3L]], seq_along(cells))
  expected[match(observed, cells)] <- partial[, "SUB_1"]
  expect_equal(values[[4L]], expected)
  expected[match(observed, cells)] <- sqrt(counts["adt1", ] * counts["adt2", ])
  expect_equal(values[[5L]], expected)
  expect_identical(object, original)
})

test_that("3D feature plots align values with the selected partial reduction", {
  skip_if_not_installed("plotly")
  object <- make_cell_dim_plot_srt()
  observed <- colnames(object)[c(2, 5, 8)]
  counts <- matrix(c(2, 4, 6), nrow = 1L, dimnames = list("adt1", observed))
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(counts, sparse = TRUE))
  selected <- rev(observed)
  partial <- cbind(
    SeuratObject::Embeddings(object, "umap")[selected, , drop = FALSE],
    seq_along(selected)
  )
  colnames(partial) <- paste0("SUB_", seq_len(ncol(partial)))
  object[["partial"]] <- SeuratObject::CreateDimReducObject(
    partial, key = "SUB_", assay = "RNA"
  )
  score <- matrix(c(20, 30), ncol = 1L, dimnames = list(observed[-1L], "SCORE_1"))
  object[["score"]] <- SeuratObject::CreateDimReducObject(
    score, key = "SCORE_", assay = "RNA"
  )
  widget <- FeatureDimPlot3D(
    object, features = c("adt1", "SCORE_1"), assay = "ADT", layer = "counts",
    reduction = "partial", verbose = FALSE
  )
  built <- plotly::plotly_build(widget)
  trace <- built$x$data[[1L]]
  expect_equal(as.numeric(trace$x), unname(partial[observed, "SUB_1"]))
  expect_equal(as.numeric(trace$z), unname(partial[observed, "SUB_3"]))
  expect_equal(as.numeric(trace$marker$color), unname(counts["adt1", ]))
  expect_identical(as.vector(trace$text), paste0("Cell:", observed, "<br />Exp:", counts["adt1", ]))
  score_button <- built$x$layout$updatemenus[[2L]]$buttons[[2L]]
  expect_identical(score_button$label, "SCORE_1")
  expect_equal(as.numeric(score_button$args[[1L]]$marker$color), c(NA_real_, 20, 30))
})

test_that("dimension plots preserve caller randomness on success and errors", {
  srt <- make_cell_dim_plot_srt()
  for (has_seed in c(TRUE, FALSE)) {
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      expect_rng_restored(
        CellDimPlot(srt, group.by = "celltype", reduction = "umap", raster = FALSE, verbose = FALSE),
        has_seed, kind
      )
      expect_rng_restored(
        FeatureDimPlot(srt, features = "gene1", layer = "counts", reduction = "umap", raster = FALSE, verbose = FALSE),
        has_seed, kind
      )
      expect_rng_restored(
        expect_error(CellDimPlot(srt, group.by = "missing", reduction = "umap", verbose = FALSE), "missing"),
        has_seed, kind
      )
      expect_rng_restored(
        expect_error(FeatureDimPlot(srt, features = "gene1", show_points = NA, verbose = FALSE), "show_points"),
        has_seed, kind
      )
    }
  }
})

test_that("repelled dimension labels render independently of caller randomness", {
  srt <- make_cell_dim_plot_srt()
  plots <- list(
    CellDimPlot(
      srt,
      group.by = "celltype", reduction = "umap", raster = FALSE,
      label = TRUE, label_repel = TRUE, seed = 29, verbose = FALSE
    ),
    FeatureDimPlot(
      srt,
      features = "gene1", layer = "counts", reduction = "umap", raster = FALSE,
      label = TRUE, label_repel = TRUE, seed = 29, verbose = FALSE
    )
  )
  files <- c(tempfile(fileext = ".png"), tempfile(fileext = ".png"))
  withr::defer(unlink(files))
  for (p in plots) {
    for (i in seq_along(files)) {
      set.seed(i)
      ggplot2::ggsave(files[[i]], p, width = 6, height = 5, dpi = 100, bg = "white")
    }
    expect_identical(unname(tools::md5sum(files[[1L]])), unname(tools::md5sum(files[[2L]])))
  }
})

test_that("CellDimPlot density contours close within shared panel limits", {
  srt <- make_cell_dim_plot_srt(40L)
  for (show_na in c(FALSE, TRUE)) {
    plots <- CellDimPlot(srt,
      group.by = "celltype", split.by = "major_type",
      reduction = "umap", add_density = TRUE, show_na = show_na,
      raster = FALSE, combine = FALSE, verbose = FALSE
    )
    for (p in plots) {
      layer <- which(vapply(p$layers, function(l) inherits(l$stat, "StatDensity2d"), logical(1)))
      paths <- ggplot2::ggplot_build(p)$data[[layer]]
      expect_gt(nrow(paths), 0L)
      closed <- vapply(split(paths, paths$group), function(d) {
        isTRUE(all.equal(as.numeric(d[1L, c("x", "y")]),
          as.numeric(d[nrow(d), c("x", "y")]),
          tolerance = 1e-8
        ))
      }, logical(1))
      expect_true(all(closed))
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
    expect_equal(
      plots[[1L]]$scales$get_scales("x")$limits,
      plots[[2L]]$scales$get_scales("x")$limits
    )
    expect_equal(
      plots[[1L]]$scales$get_scales("y")$limits,
      plots[[2L]]$scales$get_scales("y")$limits
    )
  }
})

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

test_that("CellDimPlot keeps colour scale guides as characters for ggnewscale", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = FALSE,
    force = TRUE
  )

  colour_scale <- plot$scales$get_scales("colour")
  expect_true(is.character(colour_scale$guide))
  expect_false(inherits(colour_scale$guide, "Guide"))
})

test_that("CellDimPlot can overlay lineages without ggplot2 guide comparison errors", {
  srt <- make_cell_dim_plot_srt()
  srt$Lineage1 <- seq(0, 1, length.out = ncol(srt))

  plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    lineages = "Lineage1",
    raster = FALSE,
    force = TRUE
  )

  expect_true(inherits(plot, "ggplot") || inherits(plot, "patchwork"))
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("CellDimPlot supports atlas-style grid and marked groups", {
  srt <- make_cell_dim_plot_srt()

  plot <- CellDimPlot(
    object = srt,
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
    object = srt,
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
    object = srt,
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
    object = srt,
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
    object = srt,
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
    object = srt,
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
    object = srt,
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
    object = srt,
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

test_that("dimension plot defaults use calibrated square-root scaling", {
  cell_counts <- c(
    1000L,
    3000L,
    5000L,
    10000L,
    20000L,
    30000L,
    50000L,
    100000L
  )

  expect_equal(
    vapply(cell_counts, dim_plot_default_pt_size, numeric(1)),
    c(
      1,
      0.6 * sqrt(5000 / 3000),
      0.6,
      0.6 * sqrt(5000 / 10000),
      0.3,
      0.3,
      0.3,
      0.3
    )
  )
})

test_that("dimension plot defaults remain raster-readable at common PNG resolutions", {
  point_size_mm <- dim_plot_default_pt_size(30000L)

  expect_equal(point_size_mm, 0.3)
  expect_gte(point_size_mm / 25.4 * 300, 3.5)
  expect_gte(point_size_mm / 25.4 * 600, 7)
})

test_that("raster point sizes stay continuous and scale with raster resolution", {
  expect_equal(dim_plot_raster_pt_size(0.3, auto = TRUE), 2)
  expect_equal(dim_plot_raster_pt_size(0.3, auto = FALSE), 0.3)
  expect_equal(dim_plot_raster_pt_size(1.2, auto = TRUE), 2)
  expect_equal(dim_plot_raster_pt_size(2.2, auto = TRUE), 2.2)
  expect_equal(
    dim_plot_raster_pt_size(
      0.3,
      auto = FALSE,
      pixels = c(2048, 2048)
    ),
    1.2
  )
  expect_equal(
    dim_plot_raster_pt_size(
      0.3,
      auto = TRUE,
      pixels = c(2048, 2048)
    ),
    8
  )
})

test_that("CellDimPlot applies continuous resolution-aware raster point sizes", {
  skip_if_not_installed("scattermore")
  srt <- make_cell_dim_plot_srt()

  default_plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = TRUE,
    show_stat = FALSE,
    force = TRUE
  )
  explicit_plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    pt.size = 0.3,
    raster = TRUE,
    raster.dpi = c(2048, 2048),
    show_stat = FALSE,
    force = TRUE
  )
  raster_point_sizes <- function(plot) {
    unname(vapply(
      Filter(
        function(layer) inherits(layer$geom, "GeomScattermore"),
        plot$layers
      ),
      function(layer) layer$geom_params$pointsize,
      numeric(1)
    ))
  }

  expect_equal(raster_point_sizes(default_plot), c(2, 2))
  expect_equal(raster_point_sizes(explicit_plot), c(1.2, 1.2))
})

test_that("FeatureDimPlot shares resolution-aware raster point sizing", {
  skip_if_not_installed("scattermore")
  srt <- make_cell_dim_plot_srt()
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  plot <- FeatureDimPlot(
    srt,
    features = "gene1",
    reduction = "umap",
    pt.size = 0.6,
    raster = TRUE,
    raster.dpi = c(2048, 2048),
    force = TRUE,
    theme_use = theme_blank
  )
  point_sizes <- vapply(
    Filter(
      function(layer) inherits(layer$geom, "GeomScattermore"),
      plot$layers
    ),
    function(layer) layer$geom_params$pointsize,
    numeric(1)
  )

  expect_true(length(point_sizes) > 0L)
  expect_true(all(point_sizes == 2.4))
})

test_that("dimension plot raster settings are shared and validated", {
  expect_equal(
    dim_plot_point_config(
      n = 1000L,
      raster = TRUE,
      raster.dpi = c(1024, 2048)
    ),
    list(
      pt.size = 1,
      raster = TRUE,
      raster_pt_size = 2 * sqrt(8),
      raster.dpi = c(1024, 2048)
    )
  )
  expect_error(
    dim_plot_point_config(
      n = 1000L,
      raster = TRUE,
      raster.dpi = c(512, 0)
    ),
    "positive finite"
  )
})

test_that("dimension plots auto-rasterize high-density data", {
  expect_false(dim_plot_auto_raster(100000L))
  expect_true(dim_plot_auto_raster(100001L))
})

test_that("CellDimPlot applies scaled and explicit point sizes", {
  srt <- make_cell_dim_plot_srt(n_cells = 5000L)

  default_plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = FALSE,
    force = TRUE
  )
  explicit_plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    pt.size = 1,
    raster = FALSE,
    force = TRUE
  )
  point_size <- function(plot) {
    point_layers <- Filter(
      function(layer) inherits(layer$geom, "GeomPoint"),
      plot$layers
    )
    point_layers[[length(point_layers)]]$aes_params$size
  }

  expect_equal(point_size(default_plot), 0.6)
  expect_equal(point_size(explicit_plot), 1)
})

test_that("CellDimPlot scales large non-raster plots consistently", {
  srt <- make_cell_dim_plot_srt(n_cells = 50000L)

  plot <- CellDimPlot(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = FALSE,
    force = TRUE
  )
  point_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint"),
    plot$layers
  )

  expect_equal(
    point_layers[[length(point_layers)]]$aes_params$size,
    0.3
  )
})

test_that("CellDimPlot scales raster points with resolution", {
  skip_if_not_installed("scattermore")
  srt <- make_cell_dim_plot_srt()
  plot_args <- list(
    object = srt,
    group.by = "celltype",
    reduction = "umap",
    raster = TRUE,
    raster.dpi = c(2048, 2048),
    show_stat = FALSE,
    force = TRUE
  )
  auto_plot <- do.call(CellDimPlot, plot_args)
  explicit_plot <- do.call(
    CellDimPlot,
    c(plot_args, list(pt.size = 0.3))
  )
  point_sizes <- function(plot) {
    vapply(
      Filter(
        function(layer) inherits(layer$geom, "GeomScattermore"),
        plot$layers
      ),
      function(layer) layer$geom_params$pointsize,
      numeric(1)
    )
  }

  expect_true(all(point_sizes(auto_plot) == 8))
  expect_true(all(point_sizes(explicit_plot) == 1.2))
})

test_that("CellDimPlot validates atlas grid density", {
  srt <- make_cell_dim_plot_srt()

  expect_error(
    CellDimPlot(
      object = srt,
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


test_that("FeatureDimPlot styles feature headings by identity across reordered splits", {
  srt <- make_cell_dim_plot_srt()
  srt$score <- seq_len(ncol(srt))
  plots <- FeatureDimPlot(
    srt,
    features = c("gene2", "score", "gene1"), layer = "counts",
    split.by = "major_type", reduction = "umap", show_stat = FALSE,
    title.color = c(gene1 = "#B14B28", gene2 = "#2A6F97"), combine = FALSE,
    verbose = FALSE
  )
  expect_length(plots, 6L)
  for (p in plots) {
    f <- unique(p$data$features)
    expect_identical(p$theme$strip.text.x$face, if (f == "score") "plain" else "italic")
    expect_identical(
      p$theme$strip.text.x$colour,
      switch(f,
        gene1 = "#B14B28",
        gene2 = "#2A6F97",
        score = NULL
      )
    )
    expect_equal(nrow(p$data), sum(srt$major_type == unique(p$data$major_type)))
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  p <- FeatureDimPlot(srt, "gene1",
    layer = "counts", reduction = "umap",
    title.color = "red", title.face = "bold", combine = FALSE, verbose = FALSE
  )[[1]]
  expect_identical(p$theme$strip.text.x$face, "bold")
  expect_identical(p$theme$strip.text.x$colour, "red")
  expect_error(FeatureDimPlot(srt, "gene1", title.color = c("red", "blue")), "title.color")
  expect_error(FeatureDimPlot(srt, "gene1", title.color = c(gene1 = "red", gene1 = "blue")), "title.color")
  expect_error(FeatureDimPlot(srt, "gene1", title.color = "not-a-color"), "invalid color")
})

test_that("cell density bands render without losing point colours", {
  srt <- make_cell_dim_plot_srt()
  plain <- CellDimPlot(srt, "major_type", reduction = "umap", raster = FALSE, verbose = FALSE)
  filled <- CellDimPlot(srt, "major_type",
    reduction = "umap", raster = FALSE,
    add_density = TRUE, density_filled = TRUE,
    density_color = "tomato", verbose = FALSE
  )
  expect_no_warning(built <- ggplot2::ggplot_build(filled))
  bands <- which(vapply(filled$layers, function(l) inherits(l$geom, "GeomDensity2dFilled"), logical(1)))
  expect_length(bands, 1L)
  expect_gt(nrow(built$data[[bands]]), 0L)
  expect_true(all(is.finite(built$data[[bands]]$x) & is.finite(built$data[[bands]]$y)))
  points <- which(vapply(filled$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  base_points <- which(vapply(plain$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  expect_setequal(built$data[[points]]$colour, ggplot2::ggplot_build(plain)$data[[base_points]]$colour)
  line <- which(vapply(filled$layers, function(l) inherits(l$geom, "GeomDensity2d"), logical(1)))
  expect_true(all(built$data[[line]]$colour == "tomato"))
})


test_that("dimension plot highlights preserve vector and raster point sizes", {
  skip_if_not_installed("scattermore")
  object <- make_cell_dim_plot_srt()
  original <- object
  selected <- c("cell1", "cell3", "cell10")
  for (plotter in list(CellDimPlot, FeatureDimPlot)) {
    for (raster in c(FALSE, TRUE)) {
      args <- list(
        object = object, reduction = "umap", raster = raster,
        cells.highlight = selected, sizes.highlight = 3.2,
        stroke.highlight = 0.7, alpha.highlight = 0.6,
        show_stat = FALSE, force = TRUE, verbose = FALSE
      )
      if (identical(plotter, CellDimPlot)) {
        args$group.by <- "celltype"
      } else {
        args$features <- "gene1"
        args$layer <- "counts"
      }
      plot <- do.call(plotter, args)
      point_layers <- Filter(function(layer) {
        inherits(layer$geom, if (raster) "GeomScattermore" else "GeomPoint")
      }, plot$layers)
      outlined <- tail(point_layers, 2L)
      expect_equal(unname(vapply(outlined, function(layer) nrow(layer$data), integer(1))), c(3L, 3L))
      sizes <- vapply(outlined, function(layer) {
        if (raster) layer$geom_params$pointsize else layer$aes_params$size
      }, numeric(1))
      expect_equal(unname(sizes), if (raster) c(3.7, 3) else c(3.9, 3.2))
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    }
  }
  expect_identical(object, original)
})

test_that("3D dimension plots retain highlighted traces and feature menus", {
  skip_if_not_installed("plotly")
  object <- make_cell_dim_plot_srt()
  embedding <- cbind(SeuratObject::Embeddings(object, "umap"), UMAP_3 = seq_len(ncol(object)))
  object[["umap3d"]] <- SeuratObject::CreateDimReducObject(
    embedding,
    key = "UMAP3D_", assay = "RNA"
  )
  original <- object
  selected <- c("cell1", "cell3")
  categorical <- plotly::plotly_build(CellDimPlot3D(
    object,
    group.by = "celltype", reduction = "umap3d",
    cells.highlight = selected, verbose = FALSE
  ))
  feature <- plotly::plotly_build(FeatureDimPlot3D(
    object,
    features = c("gene1", "gene2"), reduction = "umap3d",
    split.by = "major_type", layer = "counts",
    cells.highlight = selected, verbose = FALSE
  ))
  expect_equal(sum(vapply(categorical$x$data, function(trace) {
    if (identical(trace$name, "highlight")) 0L else length(trace$x)
  }, integer(1))), ncol(object))
  expect_equal(length(tail(categorical$x$data, 1L)[[1L]]$x), length(selected))
  expect_equal(length(feature$x$data), 2L)
  expect_equal(
    as.numeric(feature$x$data[[1L]]$marker$color),
    as.numeric(SeuratObject::LayerData(object, layer = "counts")["gene1", ])
  )
  expect_equal(feature$x$layout$scene$xaxis$range, range(embedding[, 1]))
  expect_identical(vapply(
    feature$x$layout$updatemenus[[2L]]$buttons,
    function(button) button$label, character(1)
  ), c("gene1", "gene2"))
  for (highlight in list(NULL, selected)) {
    view <- if (is.null(highlight)) {
      plotly::plotly_build(FeatureDimPlot3D(
        object,
        features = c("gene1", "gene2"), reduction = "umap3d",
        split.by = "major_type", layer = "counts", verbose = FALSE
      ))
    } else {
      feature
    }
    buttons <- view$x$layout$updatemenus[[1L]]$buttons
    restore <- buttons[[1L]]$args[[1L]]
    expect_identical(buttons[[1L]]$label, "All.cells")
    for (axis in seq_len(3L)) {
      name <- c("x", "y", "z")[[axis]]
      expect_equal(as.numeric(restore[[name]][[1L]]), as.numeric(embedding[, axis]))
      if (!is.null(highlight)) {
        expect_equal(as.numeric(restore[[name]][[2L]]), as.numeric(embedding[highlight, axis]))
      }
      for (button in buttons[-1L]) {
        expected <- as.numeric(embedding[, axis])
        expected[object$major_type != button$label] <- NA_real_
        expect_equal(as.numeric(button$args[[1L]][[name]][[1L]]), expected)
      }
    }
  }
  expect_identical(object, original)
})

test_that("3D HTML saves isolate staging and replace only the requested file", {
  local_dir <- withr::local_tempdir()
  file <- file.path(local_dir, "plot.html")
  other <- file.path(local_dir, "keep.txt")
  stage <- NULL
  writeLines("previous output", file)
  writeLines("keep this", other)

  save_widget <- function(widget, file, selfcontained, libdir) {
    stage <<- dirname(file)
    expect_identical(dirname(libdir), stage)
    expect_true(dir.exists(stage))
    expect_false(identical(stage, local_dir))
    expect_identical(widget, "example")
    expect_true(selfcontained)
    dir.create(libdir)
    writeLines("temporary dependency", file.path(libdir, "dependency.js"))
    writeLines("complete widget", file)
  }

  save_dim_plot_widget("example", file, save_widget)
  expect_identical(readLines(file), "complete widget")
  expect_identical(readLines(other), "keep this")
  expect_false(dir.exists(stage))
  expect_setequal(list.files(local_dir), c("plot.html", "keep.txt"))
})

test_that("3D HTML saves leave existing output intact when rendering fails", {
  local_dir <- withr::local_tempdir()
  file <- file.path(local_dir, "plot.html")
  stage <- NULL
  writeLines("previous output", file)

  save_widget <- function(widget, file, selfcontained, libdir) {
    stage <<- dirname(file)
    dir.create(libdir)
    writeLines("incomplete widget", file)
    stop("render failed")
  }

  expect_error(save_dim_plot_widget("example", file, save_widget), "render failed")
  expect_identical(readLines(file), "previous output")
  expect_false(dir.exists(stage))
})

test_that("2D highlights treat FALSE as no highlight and validate requested cell IDs", {
  withr::local_seed(11L)
  object <- make_cell_dim_plot_srt()
  original <- object
  plotters <- list(
    function(ids) CellDimPlot(
      object, group.by = "celltype", reduction = "umap", raster = FALSE,
      show_stat = FALSE, cells.highlight = ids
    ),
    function(ids) FeatureDimPlot(
      object, features = "gene1", layer = "counts", reduction = "umap",
      raster = FALSE, show_stat = FALSE, bg_cutoff = -Inf, cells.highlight = ids
    )
  )
  for (draw in plotters) {
    without_highlights <- ggplot2::ggplot_build(draw(NULL))$data
    expect_identical(ggplot2::ggplot_build(draw(FALSE))$data, without_highlights)
    expect_identical(ggplot2::ggplot_build(draw(character()))$data, without_highlights)
    expect_error(draw("missing-cell"), "No cells.*found")
    expect_warning(kept <- draw(c("missing-cell", "cell2")), "Some cells.*not found")
    expect_identical(ggplot2::ggplot_build(kept)$data,
                     ggplot2::ggplot_build(draw("cell2"))$data)
    expect_equal(nrow(tail(draw(TRUE)$layers, 1L)[[1L]]$data), ncol(object))
  }
  expect_identical(object, original)
})

test_that("3D highlights preserve traces for the matched IDs only", {
  skip_if_not_installed("plotly")
  withr::local_seed(11L)
  object <- make_cell_dim_plot_srt()
  embedding <- cbind(SeuratObject::Embeddings(object, "umap"), seq_len(ncol(object)))
  colnames(embedding) <- paste0("UMAP3D_", seq_len(ncol(embedding)))
  object[["umap3d"]] <- SeuratObject::CreateDimReducObject(
    embedding, key = "UMAP3D_", assay = "RNA"
  )
  plotters <- list(
    function(ids) CellDimPlot3D(
      object, group.by = "celltype", reduction = "umap3d", cells.highlight = ids
    ),
    function(ids) FeatureDimPlot3D(
      object, features = "gene1", layer = "counts", reduction = "umap3d",
      cells.highlight = ids
    )
  )
  traces <- function(plot) {
    lapply(plotly::plotly_build(plot)$x$data, function(x) x[c("name", "x", "y", "z", "marker")])
  }
  for (draw in plotters) {
    expect_identical(traces(draw(FALSE)), traces(draw(NULL)))
    expect_error(draw("missing-cell"), "No cells.*found")
    expect_warning(kept <- draw(c("missing-cell", "cell3")), "Some cells.*not found")
    expect_identical(traces(kept), traces(draw("cell3")))
    expect_length(tail(traces(draw(TRUE)), 1L)[[1L]]$x, ncol(object))
  }
})

test_that("correlation and expression statistic plots normalize highlight inputs", {
  withr::local_seed(11L)
  object <- make_cell_dim_plot_srt()
  meta <- object[[]]
  meta$cells <- rownames(meta)
  plotters <- list(
    function(ids) FeatureCorPlot(
      object, features = c("gene1", "gene2"), group.by = "celltype",
      layer = "counts", raster = FALSE, cells.highlight = ids
    ),
    function(ids) ExpressionStatPlot(
      meta.data = meta, stat.by = "nCount_RNA", group.by = "celltype",
      plot_type = "box", add_point = TRUE, cells.highlight = ids
    ),
    function(ids) FeatureStatPlot(
      object, stat.by = "nCount_RNA", group.by = "celltype", layer = "counts",
      plot_type = "box", add_point = TRUE, cells.highlight = ids
    )
  )
  for (draw in plotters) {
    expect_no_error(draw(FALSE))
    expect_error(draw("missing-cell"), "No cells.*found")
    expect_warning(draw(c("missing-cell", "cell3")), "Some cells.*not found")
  }
})
