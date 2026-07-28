non_white_fraction <- function(filename, tolerance = 0.98) {
  image <- png::readPNG(filename)
  rgb <- image[, , seq_len(3), drop = FALSE]
  mean(apply(rgb, c(1, 2), min) < tolerance)
}

save_point_plot <- function(plot, filename, width, height = width, dpi = 120) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
}

raster_point_size <- function(plot) {
  layers <- Filter(
    function(layer) inherits(layer$geom, "GeomScattermore"),
    plot$layers
  )
  unique(vapply(
    layers,
    function(layer) layer$geom_params$pointsize,
    numeric(1)
  ))
}

run_cell_dim_point_size_real_case <- function(
  source_dir = ".",
  output_dir = file.path("artifacts", "cell-dim-point-size")
) {
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  options(scop_env_init = FALSE, log_message.verbose = FALSE)
  devtools::load_all(source_dir, quiet = TRUE, compile = FALSE)
  data("pancreas_sub", package = "scop", envir = environment())

  set.seed(277)
  pancreas_case <- Seurat::NormalizeData(pancreas_sub, verbose = FALSE)
  pancreas_case <- Seurat::FindVariableFeatures(
    pancreas_case,
    nfeatures = 500,
    verbose = FALSE
  )
  pancreas_case <- Seurat::ScaleData(pancreas_case, verbose = FALSE)
  pancreas_case <- Seurat::RunPCA(
    pancreas_case,
    npcs = 15,
    verbose = FALSE
  )
  pancreas_case <- Seurat::RunUMAP(
    pancreas_case,
    reduction = "pca",
    dims = 1:10,
    reduction.name = "pointsize_umap",
    reduction.key = "POINTSIZEUMAP_",
    seed.use = 277,
    verbose = FALSE
  )

  base_args <- list(
    srt = pancreas_case,
    group.by = "CellType",
    reduction = "pointsize_umap",
    show_stat = FALSE,
    label = FALSE,
    legend.position = "none",
    theme_use = "theme_blank",
    force = TRUE
  )
  make_plot <- function(...) {
    do.call(CellDimPlot, c(base_args, list(...)))
  }

  vector_default <- make_plot(raster = FALSE)
  vector_large_adjusted <- make_plot(raster = FALSE, pt.size = 3)
  canvas_files <- file.path(
    output_dir,
    c(
      "canvas-4in-default.png",
      "canvas-12in-default.png",
      "canvas-12in-adjusted.png"
    )
  )
  save_point_plot(vector_default, canvas_files[[1]], width = 4)
  save_point_plot(vector_default, canvas_files[[2]], width = 12)
  save_point_plot(vector_large_adjusted, canvas_files[[3]], width = 12)
  canvas_labels <- c(
    "4 in, auto pt.size = 1",
    "12 in, auto pt.size = 1",
    "12 in, explicit pt.size = 3"
  )
  canvas_images <- Map(function(filename, label) {
    image <- magick::image_read(filename)
    image <- magick::image_resize(image, "480x480!")
    label_image <- magick::image_blank(480, 42, color = "white")
    label_image <- magick::image_annotate(
      label_image,
      label,
      gravity = "center",
      size = 20,
      color = "black"
    )
    magick::image_append(c(label_image, image), stack = TRUE)
  }, canvas_files, canvas_labels)
  magick::image_write(
    magick::image_append(do.call(c, canvas_images)),
    file.path(output_dir, "canvas-size-comparison.png")
  )

  raster_reference <- make_plot(
    raster = TRUE,
    raster.dpi = c(512, 512)
  )
  raster_legacy_large <- make_plot(
    raster = TRUE,
    raster.dpi = c(2048, 2048),
    pt.size = 0.5
  )
  raster_fixed_large <- make_plot(
    raster = TRUE,
    raster.dpi = c(2048, 2048)
  )
  raster_files <- file.path(
    output_dir,
    c(
      "raster-512-auto.png",
      "raster-2048-legacy-equivalent.png",
      "raster-2048-fixed-auto.png"
    )
  )
  Map(
    function(plot, filename) save_point_plot(plot, filename, width = 4),
    list(raster_reference, raster_legacy_large, raster_fixed_large),
    raster_files
  )
  raster_comparison <- patchwork::wrap_plots(
    raster_reference + ggplot2::ggtitle("512 px auto (radius 2)"),
    raster_legacy_large + ggplot2::ggtitle("2048 px legacy-equivalent (radius 2)"),
    raster_fixed_large + ggplot2::ggtitle("2048 px fixed auto (radius 8)"),
    nrow = 1
  )
  save_point_plot(
    raster_comparison,
    file.path(output_dir, "raster-resolution-comparison.png"),
    width = 12,
    height = 4
  )

  explicit_sizes <- c(0.3, 0.6, 0.9)
  explicit_plots <- lapply(explicit_sizes, function(size) {
    plot <- make_plot(
      raster = TRUE,
      raster.dpi = c(2048, 2048),
      pt.size = size
    )
    plot + ggplot2::ggtitle(sprintf(
      "pt.size %.1f (radius %.1f px)",
      size,
      raster_point_size(plot)
    ))
  })
  explicit_comparison <- patchwork::wrap_plots(explicit_plots, nrow = 1)
  save_point_plot(
    explicit_comparison,
    file.path(output_dir, "raster-continuous-control.png"),
    width = 12,
    height = 4
  )

  metrics <- data.frame(
    case = c(
      "vector_4in_default",
      "vector_12in_default",
      "vector_12in_adjusted",
      "raster_512_auto",
      "raster_2048_legacy_equivalent",
      "raster_2048_fixed_auto"
    ),
    pt.size = c(1, 1, 3, 1, 0.5, 1),
    raster_radius_px = c(
      NA,
      NA,
      NA,
      raster_point_size(raster_reference),
      raster_point_size(raster_legacy_large),
      raster_point_size(raster_fixed_large)
    ),
    colored_pixel_fraction = c(
      vapply(canvas_files, non_white_fraction, numeric(1)),
      vapply(raster_files, non_white_fraction, numeric(1))
    )
  )
  utils::write.csv(
    metrics,
    file.path(output_dir, "metrics.csv"),
    row.names = FALSE
  )

  invisible(list(object = pancreas_case, metrics = metrics))
}

if (sys.nframe() == 0L) {
  run_cell_dim_point_size_real_case()
}
