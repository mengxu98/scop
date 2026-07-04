#' @title Run VECTOR developmental direction inference
#'
#' @md
#' @inheritParams RunFitDevo
#' @param reduction Embedding reduction used for the direction field.
#' @param pca.reduction PCA-like reduction used to score local polarization.
#' @param dims Dimensions from `reduction`.
#' @param pca.dims Dimensions from `pca.reduction`.
#' @param grid.n Number of bins per embedding axis.
#' @param arrow.p Distance decay parameter used by the original VECTOR arrow
#' weighting. Larger values keep more distant grid centers influential.
#' @param arrow.ol Arrow vector length multiplier relative to grid spacing.
#' @param score.name Metadata column for the VECTOR score.
#'
#' @return A modified `Seurat` object.
#'
#' @references
#' Zhang F, Li X, Tian W. Unsupervised inference of developmental directions
#' for single cells using VECTOR. Cell Reports, 2020.
#' doi:10.1016/j.celrep.2020.108069.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunVECTOR(
#'   pancreas_sub,
#'   reduction = "umap",
#'   pca.reduction = "pca",
#'   verbose = FALSE
#' )
#' FeatureDimPlot(pancreas_sub, features = "VECTOR_Score")
#'
#' VECTORPlot(pancreas_sub, plot_type = "grid")
#'
#' VECTORPlot(
#'   pancreas_sub,
#'   plot_type = "raw",
#'   group.by = "SubCellType",
#'   background = "none"
#' )
RunVECTOR <- function(
  object,
  reduction = NULL,
  pca.reduction = "pca",
  dims = 1:2,
  pca.dims = 1:30,
  grid.n = 30,
  arrow.p = 0.9,
  arrow.ol = 1.5,
  score.name = "VECTOR_Score",
  tool_name = "VECTOR",
  verbose = TRUE
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object.", message_type = "error")
  }
  reduction <- reduction %||% DefaultReduction(object, min_dim = 2, verbose = FALSE)
  reduction <- resolve_reduction_name(object, reduction)
  pca.reduction <- resolve_reduction_name(object, pca.reduction)
  if (is.null(reduction)) {
    log_message("{.arg reduction} is not available in the object.", message_type = "error")
  }
  if (is.null(pca.reduction)) {
    log_message("{.arg pca.reduction} is not available in the object.", message_type = "error")
  }
  emb <- SeuratObject::Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  pca <- SeuratObject::Embeddings(object, reduction = pca.reduction)
  pca.dims <- intersect(pca.dims, seq_len(ncol(pca)))
  pca <- pca[, pca.dims, drop = FALSE]
  log_message(
    "Run VECTOR field on {.val {nrow(emb)}} cells using {.val {ncol(pca)}} PC dimensions",
    verbose = verbose
  )
  bundle <- vector_field(emb = emb, pca = pca, grid.n = grid.n, arrow.p = arrow.p, arrow.ol = arrow.ol)
  object <- Seurat::AddMetaData(
    object,
    data.frame(VECTOR_Score = bundle$score[colnames(object)], row.names = colnames(object))
  )
  if (!identical(score.name, "VECTOR_Score")) {
    object[[score.name]] <- object[["VECTOR_Score", drop = TRUE]]
  }
  bundle$parameters <- list(
    reduction = reduction,
    pca.reduction = pca.reduction,
    dims = dims,
    pca.dims = pca.dims,
    grid.n = grid.n,
    arrow.p = arrow.p,
    arrow.ol = arrow.ol
  )
  object@tools[[tool_name]] <- bundle
  object <- suppressWarnings(Seurat::LogSeuratCommand(object))
  object
}

#' @title Plot VECTOR results
#'
#' @description
#' Visualize VECTOR grid-level scores or direction arrows on the embedding used
#' by [RunVECTOR()].
#'
#' @md
#' @param object A `Seurat` object processed by [RunVECTOR()].
#' @param plot_type Plot type. `"grid"` colors occupied grid centers by grid
#' score, and `"raw"` draws VelocityPlot-style grid arrows.
#' @param tool_name Name used in `srt@tools`.
#' @param score.name Metadata column containing VECTOR scores.
#' @param group.by Optional metadata column used as the background cell color
#' for direction-field plots.
#' @param background Background for plots. `"score"` uses
#' [FeatureDimPlot()], `"group"` uses [CellDimPlot()] and requires `group.by`,
#' and `"none"` draws the flow field without cell points.
#' @param point.size Cell point size. If `NULL`, uses the same default as
#' [FeatureDimPlot()].
#' @param point.alpha Cell point alpha.
#' @param grid.size Grid-center point size.
#' @param arrow.linewidth Direction arrow line width.
#' @param arrow.length Arrow head length passed to [grid::arrow()].
#' @param arrow.angle Arrow head angle passed to [grid::arrow()].
#' @param arrow.color Direction arrow color.
#' @param title,subtitle,xlab,ylab Plot labels. By default no title is shown.
#' @param aspect.ratio Fixed aspect ratio.
#' @param legend.position,legend.direction Legend position and direction.
#' @param theme_use,theme_args Theme function and arguments.
#' @param ... Additional arguments passed to [FeatureDimPlot()] or
#' [CellDimPlot()] when a cell background is requested.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunVECTOR(pancreas_sub, verbose = FALSE)
#' VECTORPlot(pancreas_sub, plot_type = "grid")
#' VECTORPlot(pancreas_sub, plot_type = "raw", group.by = "SubCellType")
VECTORPlot <- function(
  object,
  plot_type = c("grid", "raw"),
  tool_name = "VECTOR",
  score.name = "VECTOR_Score",
  group.by = NULL,
  background = c("auto", "score", "group", "none"),
  point.size = NULL,
  point.alpha = 0.7,
  grid.size = 2,
  arrow.linewidth = 0.5,
  arrow.length = grid::unit(0.035, "inches"),
  arrow.angle = 20,
  arrow.color = "grey20",
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  aspect.ratio = 1,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  plot_type <- match.arg(plot_type)
  background <- match.arg(background)
  draw_raw <- identical(plot_type, "raw")
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object.", message_type = "error")
  }
  bundle <- object@tools[[tool_name]]
  if (is.null(bundle)) {
    log_message("Cannot find VECTOR results in {.code object@tools[[{tool_name}]]}.", message_type = "error")
  }
  emb <- as.data.frame(bundle$embedding)
  colnames(emb) <- c("x", "y")
  emb$cell <- rownames(emb)
  grid_df <- bundle$grid
  arrows <- bundle$arrows
  reduction <- bundle$parameters$reduction %||% NULL
  dims <- bundle$parameters$dims %||% c(1, 2)
  reduction_key <- object@reductions[[reduction]]@key %||% paste0(reduction, "_")
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])

  if (identical(background, "auto")) {
    background <- if (draw_raw && !is.null(group.by)) {
      "group"
    } else {
      "none"
    }
  }

  if (identical(background, "score")) {
    p <- FeatureDimPlot(
      object,
      features = score.name,
      reduction = reduction,
      dims = dims,
      pt.size = point.size,
      pt.alpha = point.alpha,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = FALSE,
      ...
    )[[1]]
  } else if (identical(background, "group")) {
    if (is.null(group.by)) {
      log_message("{.arg group.by} is required when {.arg background = 'group'}.", message_type = "error")
    }
    p <- CellDimPlot(
      object,
      group.by = group.by,
      reduction = reduction,
      dims = dims,
      pt.size = point.size,
      pt.alpha = point.alpha,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = FALSE,
      ...
    )[[1]]
  } else {
    emb_df <- as.data.frame(bundle$embedding)
    colnames(emb_df) <- c("x", "y")
    p <- ggplot2::ggplot(emb_df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_blank() +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = xlab,
        y = ylab
      ) +
      do.call(theme_use, theme_args) +
      ggplot2::theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  }

  if (draw_raw && !is.null(arrows) && nrow(arrows) > 0L) {
    raw_df <- arrows
    if (!is.null(group.by)) {
      cell_grid <- bundle$cell_grid
      cell_group <- object@meta.data[names(cell_grid), group.by, drop = TRUE]
      grid_group <- vapply(unique(cell_grid), function(g) {
        vals <- cell_group[cell_grid == g]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0L) {
          return(NA_character_)
        }
        names(sort(table(vals), decreasing = TRUE))[1L]
      }, character(1))
      raw_df$group <- unname(grid_group[raw_df$grid])
      p <- p +
        ggnewscale::new_scale_color() +
        ggplot2::geom_segment(
          data = raw_df,
          ggplot2::aes(
            x = .data$x,
            y = .data$y,
            xend = .data$xend,
            yend = .data$yend,
            color = .data$group
          ),
          linewidth = arrow.linewidth,
          arrow = grid::arrow(length = arrow.length, type = "closed", angle = arrow.angle),
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        ) +
        ggplot2::scale_color_manual(
          name = group.by,
          values = palette_colors(raw_df$group),
          na.value = arrow.color,
          guide = ggplot2::guide_legend(
            title.hjust = 0,
            order = 1,
            override.aes = list(linewidth = 2, alpha = 1)
          )
        )
    } else {
      p <- p +
        ggplot2::geom_segment(
          data = raw_df,
          ggplot2::aes(
            x = .data$x,
            y = .data$y,
            xend = .data$xend,
            yend = .data$yend
          ),
          linewidth = arrow.linewidth,
          color = arrow.color,
          arrow = grid::arrow(length = arrow.length, type = "closed", angle = arrow.angle),
          lineend = "round",
          linejoin = "mitre",
          inherit.aes = FALSE
        )
    }
  }

  if (plot_type == "grid") {
    p <- p +
      ggplot2::geom_point(
        data = grid_df,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$score),
        inherit.aes = FALSE,
        shape = 21,
        color = "grey20",
        stroke = 0.15,
        size = grid.size
      ) +
      ggplot2::scale_fill_gradientn(
        colors = palette_colors(palette = "RdYlBu"),
        name = "Grid score"
      )
  }

  p
}
