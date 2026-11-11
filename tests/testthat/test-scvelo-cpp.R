make_velocity_plot_object <- function(n_cells = 6L) {
  counts <- matrix(
    1,
    nrow = 2,
    ncol = n_cells,
    dimnames = list(paste0("gene", 1:2), paste0("cell", seq_len(n_cells)))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  embedding <- cbind(
    UMAP_1 = seq(-1, 1, length.out = n_cells),
    UMAP_2 = stats::approx(
      seq_len(6), c(-1, -0.4, 0, 0.3, 0.7, 1),
      xout = seq(1, 6, length.out = n_cells)
    )$y
  )
  rownames(embedding) <- colnames(srt)
  velocity <- cbind(
    velocity_1 = rep(0.25, n_cells),
    velocity_2 = rep(0.1, n_cells)
  )
  rownames(velocity) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  srt[["stochastic_umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = velocity,
    key = "velocity_",
    assay = "RNA"
  )
  srt
}

test_that("vector projections preserve raw and grid arrow styles", {
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  object$group <- factor(rep(c("A", "B"), each = 6L))
  original <- object
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  default_segment <- ggplot2::ggplot_build(
    ggplot2::ggplot(data.frame(x = 0, y = 0)) +
      ggplot2::geom_segment(ggplot2::aes(x, y, xend = x + 1, yend = y + 1))
  )$data[[1L]]
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    for (type in c("raw", "grid")) {
      for (group in list(NULL, "group")) {
        args <- list(
          object = object, reduction = "umap", plot_type = type,
          group.by = group, density = if (type == "raw") 1 else 0.15,
          n_neighbors = 3, min_mass = 0, return_layer = TRUE, verbose = FALSE
        )
        if (identical(plotter, PseudotimeProjectionPlot)) {
          args$time_key <- "time"
          args$method <- "gradient"
          args$label <- FALSE
        }
        layers <- do.call(plotter, args)
        expect_length(layers$velocity_layer, if (type == "raw" && !is.null(group)) 2L else 1L)
        arrow <- layers$velocity_layer[[1L]]$geom_params$arrow
        expect_identical(grid::unitType(arrow$length), "mm")
        expect_equal(as.numeric(arrow$length), 1.5)
        built <- ggplot2::ggplot_build(ggplot2::ggplot() + layers$velocity_layer)$data[[1L]]
        if (identical(plotter, VelocityPlot)) {
          expect_true(all(built$linewidth == if (type == "raw") 0.25 else 0.3))
          expect_true(all(built$alpha == if (type == "raw") 0.75 else 0.8))
          if (type == "raw") {
            expect_equal(built$xend, unname(embedding[, 1L]) + 0.25)
            expect_equal(built$yend, unname(embedding[, 2L]) + 0.1)
          }
        } else {
          expect_true(all(built$linewidth == default_segment$linewidth))
          expect_true(all(is.na(built$alpha)))
          if (type == "grid") {
            expect_identical(rownames(layers$cell_layer[[1L]]$data), colnames(object))
          }
        }
        if (type == "grid" || is.null(group)) {
          expect_true(all(built$colour == "black"))
        } else {
          colors <- thisplot::palette_colors(object$group, palette = "Chinese")
          expect_equal(built$colour, unname(colors[as.character(object$group)]))
        }
      }
    }
  }
  expect_identical(object, original)
})


test_that("vector projection calls restore caller RNG", {
  skip_if_not_installed("withr")
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  original <- object
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    args <- list(
      object = object, reduction = "umap", plot_type = "raw", density = 0.5,
      n_neighbors = 3, seed = 37, verbose = FALSE
    )
    if (identical(plotter, PseudotimeProjectionPlot)) {
      args$time_key <- "time"
      args$method <- "gradient"
      args$label <- FALSE
    }
    for (has_seed in c(TRUE, FALSE)) {
      for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        for (return_layer in c(FALSE, TRUE)) {
          expect_rng_restored(
            do.call(plotter, c(args, list(return_layer = return_layer))),
            has_seed, kind
          )
        }
        expect_rng_restored(
          expect_error(do.call(plotter, utils::modifyList(args,
            list(plot_type = "invalid")
          )), "arg"),
          has_seed, kind
        )
        expect_rng_restored(
          expect_error(do.call(plotter, c(args,
            list(theme_use = function() stop("theme failure"))
          )), "theme failure"),
          has_seed, kind
        )
      }
    }
  }
  expect_identical(object, original)
})

test_that("velocity and pseudotime stream layers remain composable", {
  skip_if_not_installed("metR")
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  original <- object
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    for (color in list(NULL, "black")) {
      args <- list(
        object = object, reduction = "umap", plot_type = "stream",
        density = 0.2, n_neighbors = 3, streamline_n = 4,
        streamline_L = 0.6, streamline_minL = 0.1,
        streamline_color = color, return_layer = TRUE, verbose = FALSE
      )
      if (identical(plotter, PseudotimeProjectionPlot)) {
        args$time_key <- "time"
        args$method <- "gradient"
        args$show_cells <- FALSE
      }
      layers <- do.call(plotter, args)
      expect_length(layers$velocity_layer, if (is.null(color)) 5L else 3L)
      expect_s3_class(
        ggplot2::ggplotGrob(ggplot2::ggplot() + layers), "gtable"
      )
    }
  }
  expect_identical(object, original)
})
