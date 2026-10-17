make_feature_density_srt <- function() {
  xy <- expand.grid(x = seq(-1, 1, length.out = 5L), y = seq(-1, 1, length.out = 5L))
  counts <- rbind(gene1 = rep(c(0, 0, 0, 0, 10), 5L), gene2 = rep(2, 25L))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  embedding <- as.matrix(xy)
  dimnames(embedding) <- list(colnames(srt), c("UMAP_1", "UMAP_2"))
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "UMAP_", assay = "RNA"
  )
  srt$group <- factor(rep(c("A", "B"), length.out = ncol(srt)))
  srt$zero <- 0
  srt
}

density_path_layer <- function(plot) {
  Filter(function(l) inherits(l$geom, "GeomPath"), plot$layers)[[1L]]
}

test_that("expression density includes zero values and preserves closed contours", {
  srt <- make_feature_density_srt()
  p <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                     add_density = TRUE, show_points = FALSE,
                     theme_use = "theme_blank", cells.highlight = TRUE, verbose = FALSE)
  paths <- density_path_layer(p)$data
  expect_true(all(is.finite(paths$value)))
  expect_true(all(paths$value >= 0 & paths$value <= 10))
  expect_true(any(paths$value < 9))
  expect_setequal(unique(paths$width), c(0.32, 0.75))
  for (d in split(paths, paths$path)) {
    expect_equal(as.numeric(d[1L, c("x", "y")]),
                 as.numeric(d[nrow(d), c("x", "y")]), tolerance = 1e-8)
  }
  types <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_false(any(types %in% c("GeomPoint", "GeomScattermore", "GeomHex")))
  expect_true("GeomCustomAnn" %in% types)
  built <- ggplot2::ggplot_build(p)
  line <- built$data[[which(types == "GeomPath")]]
  expect_true(all(is.finite(line$x) & is.finite(line$y)))
  expect_true(length(unique(line$colour)) > 1L)
  scale <- built$plot$scales$get_scales("colour")
  expect_identical(line$colour, scale$map(paths$value))
  expect_false(any(line$colour == "grey80"))
  p2 <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                      add_density = TRUE, show_points = FALSE,
                      bg_cutoff = 100, verbose = FALSE)
  expect_equal(density_path_layer(p2)$data$value, paths$value)
})

test_that("expression density handles zero, constant and invalid panels", {
  srt <- make_feature_density_srt()
  for (f in c("gene2", "zero")) {
    p <- FeatureDimPlot(srt, f, layer = "counts", reduction = "umap",
                       add_density = TRUE, show_points = FALSE,
                       verbose = FALSE)
    expected <- if (f == "zero") 0 else 2
    expect_equal(density_path_layer(p)$data$value,
                 rep(expected, nrow(density_path_layer(p)$data)), tolerance = 1e-12)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  expect_warning(
    p <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                        cells = colnames(srt)[1:2], add_density = TRUE,
                        show_points = FALSE, verbose = FALSE),
    "three finite cells"
  )
  expect_equal(nrow(density_path_layer(p)$data), 0L)
  expect_no_error(ggplot2::ggplot_build(p))
  dat <- data.frame(x = rep(c(-1, 0, 1), 400L),
                    y = rep(c(-1, 0, 0, 1), 300L), value = 3)
  expect_equal(feature_dim_expression_density(dat)$value,
               rep(3, nrow(feature_dim_expression_density(dat))), tolerance = 1e-12)
  dat$value <- 0
  expect_true(all(feature_dim_expression_density(dat)$value == 0))
})

test_that("expression density respects feature splits, selected cells and shared scales", {
  srt <- make_feature_density_srt()
  cells <- colnames(srt)[1:20]
  plots <- FeatureDimPlot(srt, c("gene1", "gene2"), layer = "counts", reduction = "umap",
                         split.by = "group", cells = cells, add_density = TRUE,
                         show_points = FALSE, keep_scale = "feature",
                         combine = FALSE, verbose = FALSE)
  expect_length(plots, 4L)
  for (p in plots) {
    paths <- density_path_layer(p)$data
    f <- unique(p$data$features)
    group <- as.character(unique(p$data$group))
    expect_true(all(rownames(p$data) %in% cells))
    expect_identical(unique(paths$features), f)
    expect_identical(unique(paths$group), group)
    dat <- p$data
    expected <- feature_dim_expression_density(data.frame(
      x = dat$x, y = dat$y, value = as.numeric(GetAssayData5(srt, layer = "counts")[f, rownames(dat)])
    ))
    expect_equal(paths$value, expected$value)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  expect_equal(plots[[1]]$scales$get_scales("colour")$limits,
               plots[[2]]$scales$get_scales("colour")$limits)
  shared <- FeatureDimPlot(srt, c("gene1", "gene2"), layer = "counts", reduction = "umap",
                          split.by = "group", add_density = TRUE,
                          show_points = FALSE, keep_scale = "all", combine = FALSE, verbose = FALSE)
  expect_true(all(vapply(shared, function(p) {
    identical(p$scales$get_scales("colour")$limits, shared[[1]]$scales$get_scales("colour")$limits)
  }, logical(1))))
  independent <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                               split.by = "group", add_density = TRUE,
                               show_points = FALSE, keep_scale = NULL, combine = FALSE, verbose = FALSE)
  for (p in independent) {
    values <- density_path_layer(p)$data$value
    expect_equal(p$scales$get_scales("colour")$limits,
                 unname(c(min(values), stats::quantile(values, 0.99) + 0.001)))
  }
})

test_that("feature density defaults to expression and CellDimPlot keeps cell density", {
  srt <- make_feature_density_srt()
  p <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                     add_density = TRUE, raster = FALSE, verbose = FALSE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))))
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPath"), logical(1))))
  expect_false(any(vapply(p$layers, function(l) inherits(l$stat, "StatDensity2d"), logical(1))))
  cells <- CellDimPlot(srt, group.by = "group", reduction = "umap", add_density = TRUE,
                       raster = FALSE, verbose = FALSE)
  expect_true(any(vapply(cells$layers, function(l) inherits(l$stat, "StatDensity2d"), logical(1))))
  plain <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                         raster = FALSE, verbose = FALSE)
  expect_true(any(vapply(plain$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))))
  expect_false(any(vapply(plain$layers, function(l) inherits(l$geom, "GeomPath"), logical(1))))
  expect_error(FeatureDimPlot(srt, "gene1", show_points = NA), "show_points")
})

test_that("filled feature density displays local mean expression", {
  srt <- make_feature_density_srt()
  plots <- FeatureDimPlot(srt, c("gene1", "gene2"), layer = "counts", reduction = "umap",
                         split.by = "group", add_density = TRUE, density_filled = TRUE,
                         show_points = FALSE, combine = FALSE, verbose = FALSE)
  for (p in plots) {
    raster <- Filter(function(l) inherits(l$geom, "GeomRaster"), p$layers)[[1L]]
    values <- raster$data$value[is.finite(raster$data$value)]
    expect_gt(length(values), 0L)
    expect_true(anyNA(raster$data$value))
    expect_true(isTRUE(raster$geom_params$interpolate))
    expect_true(all(raster$data$alpha >= 0 & raster$data$alpha <= 1))
    expect_true(any(raster$data$alpha > 0 & raster$data$alpha < 1))
    expect_identical(unique(raster$data$features), unique(p$data$features))
    expect_identical(unique(raster$data$group), as.character(unique(p$data$group)))
    if (unique(p$data$features) == "gene2") {
      expect_equal(values, rep(2, length(values)), tolerance = 1e-12)
    }
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})

test_that("density overlays preserve the original point values, colours and scale", {
  srt <- make_feature_density_srt()
  srt$ramp <- rep(c(0, 1, 2, 5, 10), 5L)
  base <- FeatureDimPlot(srt, "ramp", reduction = "umap", bg_cutoff = -Inf,
                         raster = FALSE, verbose = FALSE)
  for (filled in c(FALSE, TRUE)) {
    p <- FeatureDimPlot(srt, "ramp", reduction = "umap", bg_cutoff = -Inf,
                       add_density = TRUE, density_filled = filled,
                       raster = FALSE, verbose = FALSE)
    expect_identical(p$data$value, base$data$value)
    expect_equal(p$scales$get_scales("colour")$limits,
                 base$scales$get_scales("colour")$limits)
    point <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
    base_point <- which(vapply(base$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
    expect_identical(ggplot2::ggplot_build(p)$data[[point]]$colour,
                     ggplot2::ggplot_build(base)$data[[base_point]]$colour)
  }
})

test_that("legacy blended and fixed-colour density calls remain supported", {
  srt <- make_feature_density_srt()
  for (filled in c(FALSE, TRUE)) {
    draw <- function() {
      p <- FeatureDimPlot(srt, c("gene1", "gene2"), layer = "counts", reduction = "umap",
                         compare_features = TRUE, add_density = TRUE,
                         density_filled = filled, bg_cutoff = -Inf,
                         raster = FALSE, verbose = FALSE)
      ggplot2::ggplot_build(p)
    }
    expect_no_warning(draw())
  }
  fixed <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                          add_density = TRUE, density_color = "tomato", verbose = FALSE)
  layers <- Filter(function(l) inherits(l$stat, "StatDensity2d"), fixed$layers)
  expect_length(layers, 1L)
  if (length(layers)) expect_identical(layers[[1L]]$aes_params$colour, "tomato")
})

test_that("zero-expression splits retain shared point scales", {
  srt <- make_feature_density_srt()
  srt$zero_split <- factor(ifelse(SeuratObject::Embeddings(srt, "umap")[, 1] <= 0, "A", "B"))
  for (scope in c("feature", "all")) {
    base <- FeatureDimPlot(srt, "gene1", layer = "counts", split.by = "zero_split",
                           keep_scale = scope, combine = FALSE, verbose = FALSE)
    expected <- base[[2L]]$scales$get_scales("colour")$limits
    for (filled in c(FALSE, TRUE)) {
      plots <- FeatureDimPlot(srt, "gene1", layer = "counts", split.by = "zero_split",
                              add_density = TRUE, density_filled = filled,
                              keep_scale = scope, combine = FALSE, verbose = FALSE)
      for (p in plots) {
        expect_equal(p$scales$get_scales("colour")$limits, expected)
        expect_no_warning(ggplot2::ggplotGrob(p))
      }
    }
  }
})

test_that("expression fills inherit the feature palette unless explicitly overridden", {
  srt <- make_feature_density_srt()
  p <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                     add_density = TRUE, density_filled = TRUE, show_points = FALSE,
                     palcolor = c("#112233", "#DD6677"), verbose = FALSE)
  fill <- Filter(function(s) any(grepl("^fill", s$aesthetics)), p$scales$scales)[[1L]]
  expect_identical(fill$palette(c(0, 0.5, 1)), p$scales$get_scales("colour")$palette(c(0, 0.5, 1)))
  override <- FeatureDimPlot(srt, "gene1", layer = "counts", reduction = "umap",
                            add_density = TRUE, density_filled = TRUE, show_points = FALSE,
                            palcolor = c("#112233", "#DD6677"),
                            density_filled_palette = "Blues", verbose = FALSE)
  override_fill <- Filter(function(s) any(grepl("^fill", s$aesthetics)), override$scales$scales)[[1L]]
  expect_false(identical(override_fill$palette(c(0, 0.5, 1)), fill$palette(c(0, 0.5, 1))))
})
