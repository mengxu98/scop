#' Velocity Plot
#'
#' This function creates a velocity plot for a given Seurat object. The plot shows the velocity vectors of the cells in a specified reduction space.
#'
#' @param srt A Seurat object.
#' @param reduction Name of the reduction in the Seurat object to use for plotting.
#' @param dims Indices of the dimensions to use for plotting.
#' @param cells Cells to include in the plot. If NULL, all cells will be included.
#' @param velocity Name of the velocity to use for plotting.
#' @param plot_type Type of plot to create. Can be "raw", "grid", or "stream".
#' @param group_by Name of the column in the Seurat object metadata to group the cells by. Defaults to NULL.
#' @param group_palette Name of the palette to use for coloring the groups. Defaults to "Paired".
#' @param group_palcolor Colors to use for coloring the groups. Defaults to NULL.
#' @param n_neighbors Number of neighbors to include for the density estimation. Defaults to ceiling(ncol(srt@assays[[1]]) / 50).
#' @param density Propotion of cells to plot. Defaults to 1 (plot all cells).
#' @param smooth Smoothing parameter for density estimation. Defaults to 0.5.
#' @param scale Scaling factor for the velocity vectors. Defaults to 1.
#' @param min_mass Minimum mass value for the density-based cutoff. Defaults to 1.
#' @param cutoff_perc Percentile value for the density-based cutoff. Defaults to 5.
#' @param arrow_angle Angle of the arrowheads. Defaults to 20.
#' @param arrow_color Color of the arrowheads. Defaults to "black".
#' @param streamline_L Length of the streamlines. Defaults to 5.
#' @param streamline_minL Minimum length of the streamlines. Defaults to 1.
#' @param streamline_res Resolution of the streamlines. Defaults to 1.
#' @param streamline_n Number of streamlines to plot. Defaults to 15.
#' @param streamline_width Width of the streamlines. Defaults to c(0, 0.8).
#' @param streamline_alpha Alpha transparency of the streamlines. Defaults to 1.
#' @param streamline_color Color of the streamlines. Defaults to NULL.
#' @param streamline_palette Name of the palette to use for coloring the streamlines. Defaults to "RdYlBu".
#' @param streamline_palcolor Colors to use for coloring the streamlines. Defaults to NULL.
#' @param streamline_bg_color Background color of the streamlines. Defaults to "white".
#' @param streamline_bg_stroke Stroke width of the streamlines background. Defaults to 0.5.
#' @param aspect.ratio Aspect ratio of the plot. Defaults to 1.
#' @param title Title of the plot. Defaults to "Cell velocity".
#' @param subtitle Subtitle of the plot. Defaults to NULL.
#' @param xlab x-axis label. Defaults to NULL.
#' @param ylab y-axis label. Defaults to NULL.
#' @param legend.position Position of the legend. Defaults to "right".
#' @param legend.direction Direction of the legend. Defaults to "vertical".
#' @param theme_use Name of the theme to use for plotting. Defaults to "theme_scop".
#' @param theme_args List of theme arguments for customization. Defaults to list().
#' @param return_layer Whether to return the plot layers as a list. Defaults to FALSE.
#' @param seed Random seed for reproducibility. Defaults to 11.
#'
#' @seealso \code{\link{RunSCVELO}} \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' VelocityPlot(pancreas_sub, reduction = "UMAP")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", group_by = "SubCellType")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "grid")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream", streamline_color = "black")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream", streamline_color = "black", arrow_color = "red")
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom ggplot2 ggplot aes geom_segment scale_size scale_alpha scale_color_gradientn scale_color_manual guide_legend
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @export
VelocityPlot <- function(
    srt,
    reduction,
    dims = c(1, 2),
    cells = NULL,
    velocity = "stochastic",
    plot_type = c("raw", "grid", "stream"),
    group_by = NULL,
    group_palette = "Paired",
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
  check_r("metR")

  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  V_reduction <- paste0(velocity, "_", reduction)
  if (!V_reduction %in% names(srt@reductions)) {
    stop("Cannot find the velocity embedding ", V_reduction, ".")
  }
  X_emb <- srt@reductions[[reduction]]@cell.embeddings[, dims]
  V_emb <- srt@reductions[[V_reduction]]@cell.embeddings[, dims]
  if (!is.null(cells)) {
    X_emb <- X_emb[intersect(rownames(X_emb), cells), , drop = FALSE]
    V_emb <- V_emb[intersect(rownames(V_emb), cells), , drop = FALSE]
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
      s <- ceiling(density * nrow(X_emb))
      ix_choice <- sample(seq_len(nrow(X_emb)), size = s, replace = FALSE)
      X_emb <- X_emb[ix_choice, ]
      V_emb <- V_emb[ix_choice, ]
    }
    if (!is.null(scale)) {
      V_emb <- V_emb * scale
    }
    df_field <- cbind.data.frame(X_emb, V_emb)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(
      max(df_field[["x"]], na.rm = TRUE)^2 +
        max(df_field[["y"]], na.rm = TRUE)^2
    )
    df_field[["length_perc"]] <- df_field[["length"]] / global_size

    if (!is.null(group_by)) {
      df_field[["group_by"]] <- srt@meta.data[
        rownames(df_field),
        group_by,
        drop = TRUE
      ]
      velocity_layer <- list(
        geom_segment(
          data = df_field,
          aes(x = x, y = y, xend = x + u, yend = y + v, color = group_by),
          arrow = arrow(
            length = unit(df_field[["length_perc"]], "npc"),
            type = "closed",
            angle = arrow_angle
          ),
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ),
        scale_color_manual(
          name = group_by,
          values = palette_scop(
            df_field[["group_by"]],
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
          arrow = arrow(
            length = unit(df_field[["length_perc"]], "npc"),
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
      X_emb,
      V_emb,
      density = density,
      smooth = smooth,
      n_neighbors = n_neighbors,
      min_mass = min_mass,
      scale = scale
    )
    X_grid <- res$X_grid
    V_grid <- res$V_grid

    df_field <- cbind.data.frame(X_grid, V_grid)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(
      max(df_field[["x"]], na.rm = TRUE)^2 +
        max(df_field[["y"]], na.rm = TRUE)^2
    )
    df_field[["length_perc"]] <- df_field[["length"]] / global_size
    velocity_layer <- list(
      geom_segment(
        data = df_field,
        aes(x = x, y = y, xend = x + u, yend = y + v),
        color = arrow_color,
        arrow = arrow(
          length = unit(df_field[["length_perc"]], "npc"),
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
    check_r("metR")
    res <- compute_velocity_on_grid(
      X_emb,
      V_emb,
      density = density,
      smooth = smooth,
      n_neighbors = n_neighbors,
      min_mass = min_mass,
      scale = 1,
      cutoff_perc = cutoff_perc,
      adjust_for_stream = TRUE
    )
    X_grid <- res$X_grid
    V_grid <- res$V_grid

    # if (!is.null(density) && (density > 0 & density < 1)) {
    #   s <- ceiling(density * ncol(X_grid))
    #   ix_choice <- sample(1:ncol(X_grid), size = s, replace = FALSE)
    #   X_grid <- X_grid[, ix_choice]
    #   V_grid <- V_grid[, ix_choice, ix_choice]
    # }

    df_field <- expand.grid(X_grid[1, ], X_grid[2, ])
    colnames(df_field) <- c("x", "y")
    u <- melt(t(V_grid[1, , ]))
    v <- melt(t(V_grid[2, , ]))
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
          size = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke,
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
          size = max(streamline_width, na.rm = TRUE),
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
          size = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke,
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
            size = after_stat(step),
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
          colors = palette_scop(
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
        scale_size(range = range(streamline_width), guide = "none")
      )
    }
  }

  lab_layer <- list(labs(
    title = title,
    subtitle = subtitle,
    x = xlab,
    y = ylab
  ))
  theme_layer <- list(
    do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  )

  if (isTRUE(return_layer)) {
    return(list(
      velocity_layer = velocity_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
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

#' Compute velocity on grid
#' The original python code is on https://github.com/theislab/scvelo/blob/master/scvelo/plotting/velocity_embedding_grid.py
#'
#' @param X_emb A matrix of dimension n_obs x n_dim specifying the embedding coordinates of the cells.
#' @param V_emb A matrix of dimension n_obs x n_dim specifying the velocity vectors of the cells.
#' @param density An optional numeric value specifying the density of the grid points along each dimension. Default is 1.
#' @param smooth An optional numeric value specifying the smoothing factor for the velocity vectors. Default is 0.5.
#' @param n_neighbors An optional numeric value specifying the number of nearest neighbors for each grid point. Default is ceiling(n_obs / 50).
#' @param min_mass An optional numeric value specifying the minimum mass required for a grid point to be considered. Default is 1.
#' @param scale An optional numeric value specifying the scaling factor for the velocity vectors. Default is 1.
#' @param adjust_for_stream A logical value indicating whether to adjust the velocity vectors for streamlines. Default is FALSE.
#' @param cutoff_perc An optional numeric value specifying the percentile cutoff for removing low-density grid points. Default is 5.
#'
#' @importFrom SeuratObject as.sparse
#' @importFrom proxyC dist
#' @export
compute_velocity_on_grid <- function(
    X_emb,
    V_emb,
    density = NULL,
    smooth = NULL,
    n_neighbors = NULL,
    min_mass = NULL,
    scale = 1,
    adjust_for_stream = FALSE,
    cutoff_perc = NULL) {
  n_obs <- nrow(X_emb)
  n_dim <- ncol(X_emb)

  density <- density %||% 1
  smooth <- smooth %||% 0.5
  n_neighbors <- n_neighbors %||% ceiling(n_obs / 50)
  min_mass <- min_mass %||% 1
  cutoff_perc <- cutoff_perc %||% 5

  grs <- list()
  for (dim_i in 1:n_dim) {
    m <- min(X_emb[, dim_i], na.rm = TRUE)
    M <- max(X_emb[, dim_i], na.rm = TRUE)
    # m <- m - 0.01 * abs(M - m)
    # M <- M + 0.01 * abs(M - m)
    gr <- seq(m, M, length.out = ceiling(50 * density))
    grs <- c(grs, list(gr))
  }
  X_grid <- Matrix::as.matrix(expand.grid(grs))

  d <- dist(
    x = as.sparse(X_emb),
    y = as.sparse(X_grid),
    method = "euclidean",
    use_nan = TRUE
  )
  neighbors <- t(Matrix::as.matrix(apply(
    d,
    2,
    function(x) order(x, decreasing = FALSE)[1:n_neighbors]
  )))
  dists <- t(Matrix::as.matrix(apply(
    d,
    2,
    function(x) x[order(x, decreasing = FALSE)[1:n_neighbors]]
  )))

  # ggplot() +
  #   annotate(geom = "point", x = X_grid[, 1], y = X_grid[, 2]) +
  #   annotate(geom = "point", x = X_grid[1, 1], y = X_grid[1, 2], color = "blue") +
  #   annotate(geom = "point", x = X_grid[neighbors[1, ], 1], y = X_grid[neighbors[1, ], 2], color = "red")

  weight <- dnorm(
    dists,
    sd = mean(sapply(grs, function(g) g[2] - g[1])) * smooth
  )
  p_mass <- p_mass_V <- rowSums(weight)
  p_mass_V[p_mass_V < 1] <- 1

  # qplot(dists[,1],weight[,1])
  # qplot(py$dists[,1],py$weight[,1])

  neighbors_emb <- array(
    V_emb[neighbors, seq_len(ncol(V_emb))],
    dim = c(dim(neighbors), dim(V_emb)[2])
  )
  V_grid <- apply((neighbors_emb * c(weight)), c(1, 3), sum)
  V_grid <- V_grid / p_mass_V

  # qplot(V_grid[,1],V_grid[,2])
  # qplot(py$V_grid[,1],py$V_grid[,2])

  if (isTRUE(adjust_for_stream)) {
    X_grid <- matrix(
      c(unique(X_grid[, 1]), unique(X_grid[, 2])),
      nrow = 2,
      byrow = TRUE
    )
    ns <- floor(sqrt(length(V_grid[, 1])))
    V_grid <- reticulate::array_reshape(t(V_grid), c(2, ns, ns))

    mass <- sqrt(apply(V_grid**2, c(2, 3), sum))
    min_mass <- 10**(min_mass - 6) # default min_mass = 1e-5
    min_mass[min_mass > max(mass, na.rm = TRUE) * 0.9] <- max(
      mass,
      na.rm = TRUE
    ) *
      0.9
    cutoff <- reticulate::array_reshape(mass, dim = c(ns, ns)) < min_mass

    length <- t(apply(apply(abs(neighbors_emb), c(1, 3), mean), 1, sum))
    length <- reticulate::array_reshape(length, dim = c(ns, ns))
    cutoff <- cutoff | length < quantile(length, cutoff_perc / 100)
    V_grid[1, , ][cutoff] <- NA
  } else {
    min_mass <- min_mass * quantile(p_mass, 0.99) / 100
    X_grid <- X_grid[p_mass > min_mass, ]
    V_grid <- V_grid[p_mass > min_mass, ]
    if (!is.null(scale)) {
      V_grid <- V_grid * scale
    }
  }
  return(list(X_grid = X_grid, V_grid = V_grid))
}
