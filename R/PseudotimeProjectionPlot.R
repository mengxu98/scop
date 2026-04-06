#' @title Pseudotime Projection Plot
#'
#' @description
#' This function creates a projection plot similar to [VelocityPlot],
#' but uses pseudotime data instead of RNA velocity analysis results.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams VelocityPlot
#' @param time_key Name of the column in the Seurat object metadata containing pseudotime values.
#' @param method Method to compute velocity vectors from pseudotime.
#' Can be `"gradient"` or `"knn"`. Default is `"knn"`.
#' @param k Number of nearest neighbors to use when `method = "knn"`.
#' Default is `30`.
#' @param graph_name Name of the KNN graph in the Seurat object to use.
#' If `NULL`, a new graph will be computed. Default is `NULL`.
#' @param density Scale for the streamline grid (number of points per axis is \code{ceiling(50 * density)}).
#' Default is `2` (CellRank/scvelo-style).
#' @param title The text for the title.
#' Defaults is `"Pseudotime projection"`.
#' @param palette Deprecated alias of `group_palette`.
#' @param palcolor Deprecated alias of `group_palcolor`.
#' @param show_cells Whether to show cell points on the plot.
#' Defaults is `TRUE`.
#' @param pt.size Size of cell points.
#' Defaults is `2` (CellRank-style overlapping patches).
#' @param pt.alpha The transparency of the data points.
#' Default is `0.3` (CellRank/scvelo-style).
#' @param label Whether to label the cell groups.
#' Defaults is `TRUE` when `group.by` is specified.
#' @param label.fg Foreground color of labels.
#' Defaults is `"black"`.
#' @param label.bg Background color of labels.
#' Defaults is `"white"`.
#'
#' @seealso
#' [VelocityPlot], [CellDimPlot]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   group.by = "SubCellType"
#' )
#'
#' PseudotimeProjectionPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   time_key = "Lineage1",
#'   group.by = "SubCellType",
#'   plot_type = "stream",
#'   show_cells = TRUE,
#'   label = TRUE
#' )
#'
#' PseudotimeProjectionPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   time_key = "Lineage2",
#'   plot_type = "grid"
#' )
#'
#' PseudotimeProjectionPlot(
#'   pancreas_sub,
#'   reduction = "UMAP",
#'   time_key = "Lineage1",
#'   method = "gradient",
#'   plot_type = "raw"
#' )
PseudotimeProjectionPlot <- function(
  srt,
  reduction,
  time_key,
  dims = c(1, 2),
  cells = NULL,
  method = c("knn", "gradient"),
  k = 30,
  graph_name = NULL,
  plot_type = c("raw", "grid", "stream"),
  group.by = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50),
  density = 2,
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
  title = "Pseudotime projection",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  palette = NULL,
  palcolor = NULL,
  show_cells = TRUE,
  pt.size = 2,
  pt.alpha = 0.3,
  label = NULL,
  label.size = 4,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  seed = 11
) {
  set.seed(seed)

  plot_type <- match.arg(plot_type)
  method <- match.arg(method)

  # Preserve backward compatibility with the old palette/palcolor arguments.
  if (!is.null(palette)) {
    group_palette <- palette
  }
  if (!is.null(palcolor)) {
    group_palcolor <- palcolor
  }

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

  if (!time_key %in% colnames(srt@meta.data)) {
    log_message(
      "Cannot find the pseudotime column {.val {time_key}} in srt metadata",
      message_type = "error"
    )
  }

  x_emb <- srt@reductions[[reduction]]@cell.embeddings[, dims, drop = FALSE]
  pseudotime <- srt@meta.data[[time_key]]

  if (!is.null(cells)) {
    cell_idx <- intersect(rownames(x_emb), cells)
    x_emb <- x_emb[cell_idx, , drop = FALSE]
    pseudotime <- pseudotime[cell_idx]
  }

  if (any(is.na(pseudotime))) {
    valid_idx <- !is.na(pseudotime)
    x_emb <- x_emb[valid_idx, , drop = FALSE]
    pseudotime <- pseudotime[valid_idx]
    log_message(
      "Removed {sum(!valid_idx)} cells with NA pseudotime values",
      message_type = "warning"
    )
  }

  reduction_key <- srt@reductions[[reduction]]@key
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  if (is.null(label)) {
    label <- !is.null(group.by)
  }

  if (method == "knn") {
    v_emb <- compute_pseudotime_on_knn(
      srt = srt,
      x_emb = x_emb,
      pseudotime = pseudotime,
      k = k,
      graph_name = graph_name,
      reduction = reduction
    )
  } else if (method == "gradient") {
    v_emb <- compute_pseudotime_on_gradient(
      x_emb = x_emb,
      pseudotime = pseudotime,
      smooth = smooth
    )
  }

  if (plot_type == "raw") {
    if (!is.null(density) && (density > 0 && density < 1)) {
      s <- ceiling(density * nrow(x_emb))
      ix_choice <- sample(seq_len(nrow(x_emb)), size = s, replace = FALSE)
      x_emb <- x_emb[ix_choice, , drop = FALSE]
      v_emb <- v_emb[ix_choice, , drop = FALSE]
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
      mean(df_field[["length_perc"]], na.rm = TRUE),
      "npc"
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
      mean(df_field[["length_perc"]], na.rm = TRUE),
      "npc"
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
          linewidth = max(streamline_width, na.rm = TRUE) +
            streamline_bg_stroke,
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
          linewidth = max(streamline_width, na.rm = TRUE) +
            streamline_bg_stroke,
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

  velocity_uses_group_color <- plot_type == "raw" && !is.null(group.by)
  velocity_uses_continuous_color <- plot_type == "stream" &&
    is.null(streamline_color)

  cell_layer <- list()
  cell_scale_layer <- list()
  color_reset_layer <- list()
  if (isTRUE(show_cells)) {
    cell_df <- data.frame(
      x = x_emb[, 1],
      y = x_emb[, 2],
      row.names = rownames(x_emb)
    )
    if (!is.null(group.by)) {
      cell_df[["group.by"]] <- srt@meta.data[
        rownames(cell_df),
        group.by,
        drop = TRUE
      ]
      cell_layer <- list(
        geom_point(
          data = cell_df,
          aes(x = x, y = y, color = group.by),
          size = pt.size,
          alpha = pt.alpha,
          inherit.aes = FALSE
        )
      )
      if (!velocity_uses_group_color) {
        cell_scale_layer <- list(
          scale_color_manual(
            name = group.by,
            values = palette_colors(
              cell_df[["group.by"]],
              palette = group_palette,
              palcolor = group_palcolor
            ),
            guide = if (isTRUE(label)) {
              "none"
            } else {
              guide_legend(
                title.hjust = 0,
                order = 1,
                override.aes = list(size = 3, alpha = 1)
              )
            }
          )
        )
      }
      if (velocity_uses_continuous_color) {
        check_r("ggnewscale", verbose = FALSE)
        color_reset_layer <- list(ggnewscale::new_scale_color())
      }
    } else {
      cell_layer <- list(
        geom_point(
          data = cell_df,
          aes(x = x, y = y),
          size = pt.size,
          alpha = pt.alpha,
          color = "gray70",
          inherit.aes = FALSE
        )
      )
    }
  }

  label_layer <- list()
  if (isTRUE(label) && !is.null(group.by)) {
    cell_df_label <- data.frame(
      x = x_emb[, 1],
      y = x_emb[, 2],
      group.by = srt@meta.data[rownames(x_emb), group.by, drop = TRUE],
      row.names = rownames(x_emb)
    )
    label_df <- stats::aggregate(
      cbind(x, y) ~ group.by,
      data = cell_df_label,
      FUN = function(x) stats::median(x, na.rm = TRUE)
    )
    colnames(label_df)[1] <- "label"

    label_layer <- list(
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = x, y = y, label = label),
        fontface = "bold",
        min.segment.length = 0,
        segment.color = "gray50",
        point.size = NA,
        max.overlaps = 100,
        force = 5,
        color = label.fg,
        bg.color = label.bg,
        bg.r = label.bg.r,
        size = label.size,
        inherit.aes = FALSE
      )
    )
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
        cell_layer = cell_layer,
        cell_scale_layer = cell_scale_layer,
        color_reset_layer = color_reset_layer,
        velocity_layer = velocity_layer,
        label_layer = label_layer,
        lab_layer = lab_layer,
        theme_layer = theme_layer
      )
    )
  } else {
    return(
      ggplot() +
        cell_layer +
        cell_scale_layer +
        color_reset_layer +
        velocity_layer +
        label_layer +
        lab_layer +
        theme_layer
    )
  }
}

compute_pseudotime_on_knn <- function(
  srt,
  x_emb,
  pseudotime,
  k = 30,
  graph_name = NULL,
  reduction = NULL
) {
  cell_names <- rownames(x_emb)

  if (!is.null(graph_name) && graph_name %in% names(srt@graphs)) {
    knn_graph <- srt@graphs[[graph_name]]
    if (!inherits(knn_graph, "Graph")) {
      log_message(
        "Graph {.val {graph_name}} is not a Graph object",
        message_type = "error"
      )
    }
    knn_graph_sparse <- SeuratObject::as.sparse(knn_graph)
    neighbors_list <- lapply(seq_len(ncol(knn_graph_sparse)), function(i) {
      neighbors <- which(knn_graph_sparse[, i] > 0)
      if (length(neighbors) > 0) {
        dists <- knn_graph_sparse[neighbors, i]
        neighbors <- neighbors[order(dists, decreasing = FALSE)]
        neighbors <- neighbors[1:min(k, length(neighbors))]
      }
      neighbors
    })
  } else {
    log_message(
      "Computing KNN graph from embedding..."
    )
    d <- proxyC::dist(
      x = SeuratObject::as.sparse(x_emb),
      y = SeuratObject::as.sparse(x_emb),
      method = "euclidean",
      use_nan = TRUE
    )
    neighbors_list <- lapply(seq_len(ncol(d)), function(i) {
      neighbors <- order(d[, i], decreasing = FALSE)[
        seq_len(min(k, ncol(d) - 1)) + 1
      ]
      neighbors
    })
  }

  n_cells <- nrow(x_emb)
  n_dims <- ncol(x_emb)
  v_emb <- matrix(0, nrow = n_cells, ncol = n_dims)

  for (i in seq_len(n_cells)) {
    neighbors <- neighbors_list[[i]]
    if (length(neighbors) == 0) {
      next
    }

    neighbor_pos <- x_emb[neighbors, , drop = FALSE]
    neighbor_time <- pseudotime[neighbors]
    self_pos <- x_emb[i, , drop = FALSE]
    self_time <- pseudotime[i]

    pos_diff <- neighbor_pos -
      matrix(
        self_pos,
        nrow = length(neighbors),
        ncol = n_dims,
        byrow = TRUE
      )
    time_diff <- neighbor_time - self_time

    dists <- sqrt(rowSums(pos_diff^2))
    dists[dists == 0] <- 1e-10

    weights <- time_diff / dists
    weights[is.na(weights)] <- 0
    weights[is.infinite(weights)] <- 0

    velocity <- colSums(pos_diff * weights) / length(neighbors)
    if (any(is.na(velocity)) || any(is.infinite(velocity))) {
      velocity <- rep(0, n_dims)
    }

    v_emb[i, ] <- velocity
  }

  v_norm <- sqrt(rowSums(v_emb^2))
  v_norm[v_norm == 0] <- 1
  v_emb <- v_emb / v_norm * mean(v_norm)

  rownames(v_emb) <- cell_names
  return(v_emb)
}

compute_pseudotime_on_gradient <- function(
  x_emb,
  pseudotime,
  smooth = 0.5
) {
  n_cells <- nrow(x_emb)
  n_dims <- ncol(x_emb)
  v_emb <- matrix(0, nrow = n_cells, ncol = n_dims)

  d <- proxyC::dist(
    x = SeuratObject::as.sparse(x_emb),
    y = SeuratObject::as.sparse(x_emb),
    method = "euclidean",
    use_nan = TRUE
  )

  k_local <- max(10, ceiling(n_cells * smooth / 100))

  for (i in seq_len(n_cells)) {
    neighbors <- order(d[, i], decreasing = FALSE)[2:(k_local + 1)]

    if (length(neighbors) == 0) {
      next
    }

    neighbor_pos <- x_emb[neighbors, , drop = FALSE]
    neighbor_time <- pseudotime[neighbors]
    self_pos <- x_emb[i, , drop = FALSE]
    self_time <- pseudotime[i]

    pos_diff <- neighbor_pos -
      matrix(
        self_pos,
        nrow = length(neighbors),
        ncol = n_dims,
        byrow = TRUE
      )
    time_diff <- neighbor_time - self_time

    dists <- sqrt(rowSums(pos_diff^2))
    dists[dists == 0] <- 1e-10

    sigma <- mean(dists) * smooth
    weights <- exp(-dists^2 / (2 * sigma^2))
    weights <- weights * time_diff

    gradient <- colSums(pos_diff * weights) / sum(abs(weights) + 1e-10)
    if (any(is.na(gradient)) || any(is.infinite(gradient))) {
      gradient <- rep(0, n_dims)
    }

    v_emb[i, ] <- gradient
  }

  v_norm <- sqrt(rowSums(v_emb^2))
  v_norm[v_norm == 0] <- 1
  v_emb <- v_emb / v_norm * mean(v_norm)

  rownames(v_emb) <- rownames(x_emb)
  return(v_emb)
}
