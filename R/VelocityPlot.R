#' @title Velocity Plot
#'
#' @description
#' This function creates a velocity plot for a given Seurat object.
#' The plot shows the velocity vectors of the cells in a specified reduction space.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams GraphPlot
#' @param velocity Name of the velocity to use for plotting.
#' Default is `"stochastic"`.
#' @param plot_type Type of plot to create.
#' Can be `"raw"`, `"grid"`, or `"stream"`.
#' @param group_palette Name of the palette to use for coloring the groups.
#' Defaults is `"Chinese"`.
#' @param group_palcolor Colors to use for coloring the groups.
#' Defaults is `NULL`.
#' @param n_neighbors Number of neighbors to include for the density estimation.
#' Defaults is `ceiling(ncol(srt@assays[[1]]) / 50)`.
#' @param density Proportion of cells to plot.
#' Defaults is `1` (plot all cells).
#' @param smooth Smoothing parameter for density estimation.
#' Defaults is `0.5`.
#' @param scale Scaling factor for the velocity vectors.
#' Defaults is `1`.
#' @param min_mass Minimum mass value for the density-based cutoff.
#' Defaults is `1`.
#' @param cutoff_perc Percentile value for the density-based cutoff.
#' Defaults is `5`.
#' @param arrow_angle Angle of the arrowheads.
#' Defaults is `20`.
#' @param arrow_color Color of the arrowheads.
#' Defaults is `"black"`.
#' @param streamline_L Length of the streamlines.
#' Defaults is `5`.
#' @param streamline_minL Minimum length of the streamlines.
#' Defaults is `1`.
#' @param streamline_res Resolution of the streamlines.
#' Defaults is `1`.
#' @param streamline_n Number of streamlines to plot.
#' Defaults is `15`.
#' @param streamline_width Width of the streamlines.
#' Defaults is `c(0, 0.8)`.
#' @param streamline_alpha Alpha transparency of the streamlines.
#' Defaults is `1`.
#' @param streamline_color Color of the streamlines.
#' Defaults is `NULL`.
#' @param streamline_palette Name of the palette to use for coloring the streamlines.
#' Defaults is `"RdYlBu"`.
#' @param streamline_palcolor Colors to use for coloring the streamlines.
#' Defaults is `NULL`.
#' @param streamline_bg_color Background color of the streamlines.
#' Defaults is `"white"`.
#' @param streamline_bg_stroke Stroke width of the streamlines background.
#' Defaults is `0.5`.
#' @param title The text for the title.
#' Defaults is `"Cell velocity"`.
#'
#' @seealso
#' [RunSCVELO], [CellDimPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCVELO(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   linear_reduction = "pca",
#'   nonlinear_reduction = "umap",
#'   return_seurat = TRUE
#' )
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP"
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   group.by = "SubCellType"
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "grid"
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "stream"
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "stream",
#'   streamline_color = "black"
#' )
#'
#' VelocityPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   plot_type = "stream",
#'   streamline_color = "black",
#'   arrow_color = "red"
#' )
#' }
VelocityPlot <- function(
    srt,
    reduction,
    dims = c(1, 2),
    cells = NULL,
    velocity = "stochastic",
    plot_type = c("raw", "grid", "stream"),
    group.by = NULL,
    group_palette = "Chinese",
    group_palcolor = NULL,
    n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50),
    density = 1,
    smooth = 0.5,
    scale = 1,
    min_mass = 1,
    cutoff_perc = 5,
    arrow_angle = 20,
    arrow_color = "black",
    streamline_L = 5,
    streamline_minL = 1,
    streamline_res = 1,
    streamline_n = 15,
    streamline_width = c(0, 0.8),
    streamline_alpha = 1,
    streamline_color = NULL,
    streamline_palette = "RdYlBu",
    streamline_palcolor = NULL,
    streamline_bg_color = "white",
    streamline_bg_stroke = 0.5,
    aspect.ratio = 1,
    title = "Cell velocity",
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    return_layer = FALSE,
    seed = 11) {
  set.seed(seed)

  plot_type <- match.arg(plot_type)
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }

  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "{.val {reduction}} is not in the srt reduction names",
      message_type = "error"
    )
  }
  v_reduction <- paste0(velocity, "_", reduction)
  if (!v_reduction %in% names(srt@reductions)) {
    log_message(
      "Cannot find the velocity embedding {.val {v_reduction}}",
      message_type = "error"
    )
  }
  x_emb <- srt@reductions[[reduction]]@cell.embeddings[, dims]
  v_emb <- srt@reductions[[v_reduction]]@cell.embeddings[, dims]
  if (!is.null(cells)) {
    x_emb <- x_emb[intersect(rownames(x_emb), cells), , drop = FALSE]
    v_emb <- v_emb[intersect(rownames(v_emb), cells), , drop = FALSE]
  }

  reduction_key <- srt@reductions[[reduction]]@key
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  if (plot_type == "raw") {
    if (!is.null(density) && (density > 0 && density < 1)) {
      s <- ceiling(density * nrow(x_emb))
      ix_choice <- sample(seq_len(nrow(x_emb)), size = s, replace = FALSE)
      x_emb <- x_emb[ix_choice, ]
      v_emb <- v_emb[ix_choice, ]
    }
    if (!is.null(scale)) {
      v_emb <- v_emb * scale
    }
    df_field <- cbind.data.frame(x_emb, v_emb)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(
      max(df_field[["x"]], na.rm = TRUE)^2 +
        max(df_field[["y"]], na.rm = TRUE)^2
    )
    df_field[["length_perc"]] <- df_field[["length"]] / global_size

    arrow_length <- grid::unit(
      mean(df_field[["length_perc"]], na.rm = TRUE), "npc"
    )

    if (!is.null(group.by)) {
      df_field[["group.by"]] <- srt@meta.data[
        rownames(df_field),
        group.by,
        drop = TRUE
      ]
      velocity_layer <- list(
        geom_segment(
          data = df_field,
          aes(x = x, y = y, xend = x + u, yend = y + v, color = group.by),
          arrow = grid::arrow(
            length = arrow_length,
            type = "closed",
            angle = arrow_angle
          ),
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        scale_color_manual(
          name = group.by,
          values = palette_colors(
            df_field[["group.by"]],
            palette = group_palette,
            palcolor = group_palcolor
          ),
          guide = guide_legend(
            title.hjust = 0,
            order = 1,
            override.aes = list(linewidth = 2, alpha = 1)
          )
        )
      )
    } else {
      velocity_layer <- list(
        geom_segment(
          data = df_field,
          aes(x = x, y = y, xend = x + u, yend = y + v),
          color = arrow_color,
          arrow = grid::arrow(
            length = arrow_length,
            type = "closed",
            angle = arrow_angle
          ),
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        )
      )
    }
  }
  if (plot_type == "grid") {
    res <- compute_velocity_on_grid(
      x_emb,
      v_emb,
      density = density,
      smooth = smooth,
      n_neighbors = n_neighbors,
      min_mass = min_mass,
      scale = scale
    )
    x_grid <- res$x_grid
    v_grid <- res$v_grid

    df_field <- cbind.data.frame(x_grid, v_grid)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(
      max(df_field[["x"]], na.rm = TRUE)^2 +
        max(df_field[["y"]], na.rm = TRUE)^2
    )
    df_field[["length_perc"]] <- df_field[["length"]] / global_size

    arrow_length <- grid::unit(
      mean(df_field[["length_perc"]], na.rm = TRUE), "npc"
    )
    velocity_layer <- list(
      geom_segment(
        data = df_field,
        aes(x = x, y = y, xend = x + u, yend = y + v),
        color = arrow_color,
        arrow = grid::arrow(
          length = arrow_length,
          type = "closed",
          angle = arrow_angle
        ),
        lineend = "round",
        linejoin = "mitre",
        inherit.aes = FALSE
      )
    )
  }
  if (plot_type == "stream") {
    check_r("metR", verbose = FALSE)
    res <- compute_velocity_on_grid(
      x_emb,
      v_emb,
      density = density,
      smooth = smooth,
      n_neighbors = n_neighbors,
      min_mass = min_mass,
      scale = 1,
      cutoff_perc = cutoff_perc,
      adjust_for_stream = TRUE
    )
    x_grid <- res$x_grid
    v_grid <- res$v_grid

    df_field <- expand.grid(x_grid[1, ], x_grid[2, ])
    colnames(df_field) <- c("x", "y")
    u <- reshape2::melt(Matrix::t(v_grid[1, , ]))
    v <- reshape2::melt(Matrix::t(v_grid[2, , ]))
    df_field[, "u"] <- u$value
    df_field[, "v"] <- v$value
    df_field[is.na(df_field)] <- 0

    if (!is.null(streamline_color)) {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field,
          aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          linewidth = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke,
          color = streamline_bg_color,
          alpha = streamline_alpha,
          arrow.type = "closed",
          arrow.angle = arrow_angle,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field,
          aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          linewidth = max(streamline_width, na.rm = TRUE),
          color = streamline_color,
          alpha = streamline_alpha,
          arrow.type = "closed",
          arrow.angle = arrow_angle,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field,
          aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          linetype = 0,
          color = arrow_color,
          arrow.type = "closed",
          arrow.angle = arrow_angle,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        )
      )
    } else {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field,
          aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          linewidth = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke,
          color = streamline_bg_color,
          alpha = streamline_alpha,
          arrow.type = "closed",
          arrow.angle = arrow_angle,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field,
          aes(
            x = x,
            y = y,
            dx = u,
            dy = v,
            linewidth = after_stat(step),
            color = sqrt(after_stat(dx)^2 + after_stat(dy)^2)
          ),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          alpha = streamline_alpha,
          arrow = NULL,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field,
          aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L,
          min.L = streamline_minL,
          res = streamline_res,
          n = streamline_n,
          linetype = 0,
          color = arrow_color,
          arrow.type = "closed",
          arrow.angle = arrow_angle,
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        scale_color_gradientn(
          name = "Velocity",
          colors = palette_colors(
            palette = streamline_palette,
            palcolor = streamline_palcolor
          ),
          guide = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            title.hjust = 0,
            order = 1
          )
        ),
        scale_linewidth(range = range(streamline_width), guide = "none")
      )
    }
  }

  lab_layer <- list(
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    )
  )
  theme_layer <- list(
    do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  )

  if (isTRUE(return_layer)) {
    return(
      list(
        velocity_layer = velocity_layer,
        lab_layer = lab_layer,
        theme_layer = theme_layer
      )
    )
  } else {
    return(
      ggplot() +
        velocity_layer +
        lab_layer +
        lab_layer +
        theme_layer
    )
  }
}
