#' GraphPlot
#'
#' A function to plot a graph with nodes and edges.
#'
#' @param node A data frame representing the nodes of the graph.
#' @param edge A matrix representing the edges of the graph.
#' @param transition A matrix representing the transitions between nodes.
#' @param node_coord A character vector specifying the names of the columns in \code{node} that represent the x and y coordinates.
#' @param node_group A character vector specifying the name of the column in \code{node} that represents the grouping of the nodes.
#' @param node_palette A character vector specifying the name of the color palette for node groups.
#' @param node_palcolor A character vector specifying the names of the colors for each node group.
#' @param node_size A numeric value or column name of \code{node} specifying the size of the nodes.
#' @param node_alpha A numeric value or column name of \code{node} specifying the transparency of the nodes.
#' @param node_highlight A character vector specifying the names of nodes to highlight.
#' @param node_highlight_color A character vector specifying the color for highlighting nodes.
#' @param label A logical value indicating whether to show labels for the nodes.
#' @param label.size A numeric value specifying the size of the labels.
#' @param label.fg A character vector specifying the foreground color of the labels.
#' @param label.bg A character vector specifying the background color of the labels.
#' @param label.bg.r A numeric value specifying the background color transparency of the labels.
#' @param label_insitu A logical value indicating whether to display the node group labels in situ or as numeric values.
#' @param label_repel A logical value indicating whether to use force-directed label repulsion.
#' @param label_repulsion A numeric value specifying the repulsion force for labels.
#' @param label_point_size A numeric value specifying the size of the label points.
#' @param label_point_color A character vector specifying the color of the label points.
#' @param label_segment_color A character vector specifying the color for the label segments.
#' @param edge_threshold A numeric value specifying the threshold for removing edges.
#' @param use_triangular A character vector specifying which part of the edge matrix to use (upper, lower, both).
#' @param edge_line A character vector specifying the type of line for edges (straight, curved).
#' @param edge_line_curvature A numeric value specifying the curvature of curved edges.
#' @param edge_line_angle A numeric value specifying the angle of curved edges.
#' @param edge_color A character vector specifying the color of the edges.
#' @param edge_size A numeric vector specifying the range of edge sizes.
#' @param edge_alpha A numeric value specifying the transparency of the edges.
#' @param edge_shorten A numeric value specifying the length of the edge shorten.
#' @param edge_offset A numeric value specifying the length of the edge offset.
#' @param edge_highlight A character vector specifying the names of edges to highlight.
#' @param edge_highlight_color A character vector specifying the color for highlighting edges.
#' @param transition_threshold A numeric value specifying the threshold for removing transitions.
#' @param transition_line A character vector specifying the type of line for transitions (straight, curved).
#' @param transition_line_curvature A numeric value specifying the curvature of curved transitions.
#' @param transition_line_angle A numeric value specifying the angle of curved transitions.
#' @param transition_color A character vector specifying the color of the transitions.
#' @param transition_size A numeric vector specifying the range of transition sizes.
#' @param transition_alpha A numeric value specifying the transparency of the transitions.
#' @param transition_arrow_type A character vector specifying the type of arrow for transitions (closed, open).
#' @param transition_arrow_angle A numeric value specifying the angle of the transition arrow.
#' @param transition_arrow_length A numeric value specifying the length of the transition arrow.
#' @param transition_shorten A numeric value specifying the length of the transition shorten.
#' @param transition_offset A numeric value specifying the length of the transition offset.
#' @param transition_highlight A character vector specifying the names of transitions to highlight.
#' @param transition_highlight_color A character vector specifying the color for highlighting transitions.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' @param title A character value specifying the title of the plot.
#' @param subtitle A character value specifying the subtitle of the plot.
#' @param xlab A character value specifying the label for the x-axis.
#' @param ylab A character value specifying the label for the y-axis.
#' @param legend.position A character value specifying the position of the legend.
#' @param legend.direction A character value specifying the direction of the legend.
#' @param theme_use A character value specifying the theme to use.
#' @param theme_args A list of arguments to be passed to the theme.
#' @param return_layer A logical value indicating whether to return the layers of the plot instead of the plot itself.
#'
#' @seealso \code{\link{CellDimPlot}}
#'
#' @importFrom ggplot2 scale_size_identity scale_size_continuous scale_size_discrete scale_alpha_identity scale_alpha_continuous scale_alpha_discrete geom_curve geom_segment geom_point scale_color_manual guide_legend guides labs aes scale_linewidth_continuous
#' @importFrom ggnewscale new_scale
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow
#' @importFrom reshape2 melt
#' @export
GraphPlot <- function(
    node,
    edge,
    transition = NULL,
    node_coord = c("x", "y"),
    node_group = NULL,
    node_palette = "Paired",
    node_palcolor = NULL,
    node_size = 4,
    node_alpha = 1,
    node_highlight = NULL,
    node_highlight_color = "red",
    label = FALSE,
    label.size = 3.5,
    label.fg = "white",
    label.bg = "black",
    label.bg.r = 0.1,
    label_insitu = FALSE,
    label_repel = FALSE,
    label_repulsion = 20,
    label_point_size = 1,
    label_point_color = "black",
    label_segment_color = "black",
    edge_threshold = 0.01,
    use_triangular = c("upper", "lower", "both"),
    edge_line = c("straight", "curved"),
    edge_line_curvature = 0.3,
    edge_line_angle = 90,
    edge_color = "grey40",
    edge_size = c(0.2, 1),
    edge_alpha = 0.5,
    edge_shorten = 0,
    edge_offset = 0,
    edge_highlight = NULL,
    edge_highlight_color = "red",
    transition_threshold = 0.01,
    transition_line = c("straight", "curved"),
    transition_line_curvature = 0.3,
    transition_line_angle = 90,
    transition_color = "black",
    transition_size = c(0.2, 1),
    transition_alpha = 1,
    transition_arrow_type = "closed",
    transition_arrow_angle = 20,
    transition_arrow_length = unit(0.02, "npc"),
    transition_shorten = 0.05,
    transition_offset = 0,
    transition_highlight = NULL,
    transition_highlight_color = "red",
    aspect.ratio = 1,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    return_layer = FALSE) {
  use_triangular <- match.arg(use_triangular)
  edge_line <- match.arg(edge_line)
  transition_line <- match.arg(transition_line)
  if (!is.data.frame(node)) {
    stop("'node' must be a data.frame object.")
  }
  if (!is.matrix(edge)) {
    stop("'edge' must be a matrix object.")
  }
  if (!identical(nrow(edge), ncol(edge))) {
    stop("nrow and ncol is not identical in edge matrix")
  }
  if (!identical(nrow(edge), nrow(node))) {
    stop("nrow is not identical between edge and node.")
  }
  if (!identical(rownames(edge), rownames(node))) {
    warning(
      "rownames of node is not identical with edge matrix. They will correspond according to the order.",
      immediate. = TRUE
    )
    colnames(edge) <- rownames(edge) <- rownames(node) <- rownames(node) %||%
      colnames(edge) %||%
      rownames(edge)
  }
  if (!all(node_coord %in% colnames(node))) {
    stop(
      "Cannot find the node_coord ",
      paste(node_coord[!node_coord %in% colnames(node)], collapse = ","),
      " in the node column"
    )
  }
  if (!is.null(transition)) {
    if (!identical(nrow(transition), nrow(node))) {
      stop("nrow is not identical between transition and node.")
    }
    if (!identical(rownames(transition), rownames(node))) {
      warning(
        "rownames of node is not identical with transition matrix. They will correspond according to the order.",
        immediate. = TRUE
      )
      colnames(transition) <- rownames(transition) <- rownames(
        node
      ) <- rownames(node) %||% colnames(transition) %||% rownames(transition)
    }
  }
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  node <- as.data.frame(node)
  node[["x"]] <- node[[node_coord[1]]]
  node[["y"]] <- node[[node_coord[2]]]
  node[["node_name"]] <- rownames(node)
  node_group <- node_group %||% "node_name"
  node_size <- node_size %||% 5
  node_alpha <- node_alpha %||% 1
  if (!node_group %in% colnames(node)) {
    node[["node_group"]] <- node_group
  } else {
    node[["node_group"]] <- node[[node_group]]
  }
  if (!is.factor(node[["node_group"]])) {
    node[["node_group"]] <- factor(
      node[["node_group"]],
      levels = unique(node[["node_group"]])
    )
  }
  if (!node_size %in% colnames(node)) {
    if (!is.numeric(node_size)) {
      node_size <- 5
    }
    node[["node_size"]] <- node_size
    scale_size <- scale_size_identity()
  } else {
    node[["node_size"]] <- node[[node_size]]
    if (is.numeric(node[[node_size]])) {
      scale_size <- scale_size_continuous(name = node_size)
    } else {
      scale_size <- scale_size_discrete()
    }
  }
  if (!node_alpha %in% colnames(node)) {
    if (!is.numeric(node_alpha)) {
      node_alpha <- 1
    }
    node[["node_alpha"]] <- node_alpha
    scale_alpha <- scale_alpha_identity()
  } else {
    node[["node_alpha"]] <- node[[node_alpha]]
    if (is.numeric(node[[node_alpha]])) {
      scale_alpha <- scale_alpha_continuous()
    } else {
      scale_alpha <- scale_alpha_discrete()
    }
  }

  if (isTRUE(label) && !isTRUE(label_insitu)) {
    label_use <- paste0(
      1:nlevels(node[["node_group"]]),
      ": ",
      levels(node[["node_group"]])
    )
  } else {
    label_use <- levels(node[["node_group"]])
  }
  global_size <- sqrt(
    max(node[["x"]], na.rm = TRUE)^2 + max(node[["y"]], na.rm = TRUE)^2
  )

  edge[edge <= edge_threshold] <- NA
  if (use_triangular == "upper") {
    edge[lower.tri(edge)] <- NA
  } else if (use_triangular == "lower") {
    edge[upper.tri(edge)] <- NA
  }
  edge_df <- melt(edge, na.rm = TRUE, stringsAsFactors = FALSE)
  if (nrow(edge_df) == 0) {
    edge_layer <- NULL
  } else {
    colnames(edge_df) <- c("from", "to", "size")
    edge_df[, "from"] <- as.character(edge_df[, "from"])
    edge_df[, "to"] <- as.character(edge_df[, "to"])
    edge_df[, "size"] <- as.numeric(edge_df[, "size"])
    edge_df[, "x"] <- node[edge_df[, "from"], "x"]
    edge_df[, "y"] <- node[edge_df[, "from"], "y"]
    edge_df[, "xend"] <- node[edge_df[, "to"], "x"]
    edge_df[, "yend"] <- node[edge_df[, "to"], "y"]
    rownames(edge_df) <- edge_df[, "edge_name"] <- paste0(
      edge_df[, "from"],
      "-",
      edge_df[, "to"]
    )
    edge_df <- segementsDf(
      edge_df,
      global_size * edge_shorten,
      global_size * edge_shorten,
      global_size * edge_offset
    )

    linetype <- ifelse(is.null(transition), 1, 2)
    if (edge_line == "straight") {
      edge_layer <- list(
        geom_segment(
          data = edge_df,
          mapping = aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend,
            linewidth = size
          ),
          lineend = "round",
          linejoin = "mitre",
          linetype = linetype,
          color = edge_color,
          alpha = edge_alpha,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[
          edge_df[["edge_name"]] %in% edge_highlight, ,
          drop = FALSE
        ]
        edge_layer <- c(
          edge_layer,
          list(
            geom_segment(
              data = edge_df_highlight,
              mapping = aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend,
                linewidth = size
              ),
              lineend = "round",
              linejoin = "mitre",
              linetype = linetype,
              color = edge_highlight_color,
              alpha = 1,
              inherit.aes = FALSE,
              show.legend = FALSE
            )
          )
        )
      }
    } else {
      edge_layer <- list(
        geom_curve(
          data = edge_df,
          mapping = aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend,
            linewidth = size
          ),
          curvature = edge_line_curvature,
          angle = edge_line_angle,
          lineend = "round",
          linetype = linetype,
          color = edge_color,
          alpha = edge_alpha,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[
          edge_df[["edge_name"]] %in% edge_highlight, ,
          drop = FALSE
        ]
        edge_layer <- c(
          edge_layer,
          list(
            geom_curve(
              data = edge_df_highlight,
              mapping = aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend,
                linewidth = size
              ),
              curvature = edge_line_curvature,
              angle = edge_line_angle,
              lineend = "round",
              linetype = linetype,
              color = edge_highlight_color,
              alpha = 1,
              inherit.aes = FALSE,
              show.legend = FALSE
            )
          )
        )
      }
    }
    edge_layer <- c(
      edge_layer,
      list(
        scale_linewidth_continuous(range = range(edge_size), guide = "none"),
        new_scale("linewidth")
      )
    )
  }

  if (!is.null(transition)) {
    trans2 <- trans1 <- Matrix::as.matrix(transition)
    trans1[lower.tri(trans1)] <- 0
    trans2[upper.tri(trans2)] <- 0
    trans <- t(trans1) - trans2
    trans[abs(trans) <= transition_threshold] <- NA
    trans_df <- melt(trans, na.rm = TRUE, stringsAsFactors = FALSE)
    if (nrow(trans_df) == 0) {
      trans_layer <- NULL
    } else {
      trans_df <- as.data.frame(t(apply(trans_df, 1, function(x) {
        if (as.numeric(x[3]) < 0) {
          return(c(x[c(2, 1)], -as.numeric(x[3])))
        } else {
          return(x)
        }
      })))
      colnames(trans_df) <- c("from", "to", "size")
      trans_df[, "from"] <- as.character(trans_df[, "from"])
      trans_df[, "to"] <- as.character(trans_df[, "to"])
      trans_df[, "size"] <- as.numeric(trans_df[, "size"])
      trans_df[, "x"] <- node[trans_df[, "from"], "x"]
      trans_df[, "y"] <- node[trans_df[, "from"], "y"]
      trans_df[, "xend"] <- node[trans_df[, "to"], "x"]
      trans_df[, "yend"] <- node[trans_df[, "to"], "y"]
      rownames(trans_df) <- trans_df[, "trans_name"] <- paste0(
        trans_df[, "from"],
        "-",
        trans_df[, "to"]
      )
      trans_df <- segementsDf(
        trans_df,
        global_size * transition_shorten,
        global_size * transition_shorten,
        global_size * transition_offset
      )

      if (transition_line == "straight") {
        trans_layer <- list(
          geom_segment(
            data = trans_df,
            mapping = aes(
              x = x,
              y = y,
              xend = xend,
              yend = yend,
              linewidth = size
            ),
            arrow = arrow(
              angle = transition_arrow_angle,
              type = transition_arrow_type,
              length = transition_arrow_length
            ),
            lineend = "round",
            linejoin = "mitre",
            color = transition_color,
            alpha = transition_alpha,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
        )
        if (!is.null(transition_highlight)) {
          trans_df_highlight <- trans_df[
            trans_df[["trans_name"]] %in% transition_highlight, ,
            drop = FALSE
          ]
          trans_layer <- c(
            trans_layer,
            list(
              geom_segment(
                data = trans_df_highlight,
                mapping = aes(
                  x = x,
                  y = y,
                  xend = xend,
                  yend = yend,
                  linewidth = size
                ),
                arrow = arrow(
                  angle = transition_arrow_angle,
                  type = transition_arrow_type,
                  length = transition_arrow_length
                ),
                lineend = "round",
                linejoin = "mitre",
                color = transition_highlight_color,
                alpha = 1,
                inherit.aes = FALSE,
                show.legend = FALSE
              )
            )
          )
        }
      } else {
        trans_layer <- list(
          geom_curve(
            data = trans_df,
            mapping = aes(
              x = x,
              y = y,
              xend = xend,
              yend = yend,
              linewidth = size
            ),
            arrow = arrow(
              angle = transition_arrow_angle,
              type = transition_arrow_type,
              length = transition_arrow_length
            ),
            curvature = transition_line_curvature,
            angle = transition_line_angle,
            lineend = "round",
            color = transition_color,
            alpha = transition_alpha,
            inherit.aes = FALSE,
            show.legend = FALSE
          )
        )
        if (!is.null(edge_highlight)) {
          trans_df_highlight <- trans_df[
            trans_df[["trans_name"]] %in% transition_highlight, ,
            drop = FALSE
          ]
          trans_layer <- c(
            trans_layer,
            list(
              geom_curve(
                data = trans_df_highlight,
                mapping = aes(
                  x = x,
                  y = y,
                  xend = xend,
                  yend = yend,
                  linewidth = size
                ),
                arrow = arrow(
                  angle = transition_arrow_angle,
                  type = transition_arrow_type,
                  length = transition_arrow_length
                ),
                curvature = transition_line_curvature,
                angle = transition_line_angle,
                lineend = "round",
                color = transition_highlight_color,
                alpha = 1,
                inherit.aes = FALSE,
                show.legend = FALSE
              )
            )
          )
        }
      }
      trans_layer <- c(
        trans_layer,
        list(
          scale_linewidth_continuous(
            range = range(transition_size),
            guide = "none"
          ),
          new_scale("linewidth")
        )
      )
    }
  } else {
    trans_layer <- NULL
  }

  node_layer <- list(
    geom_point(
      data = node,
      aes(x = x, y = y, size = node_size * 1.2),
      color = "black",
      show.legend = FALSE,
      inherit.aes = FALSE
    ),
    geom_point(
      data = node,
      aes(
        x = x,
        y = y,
        size = node_size,
        color = node_group,
        alpha = node_alpha
      ),
      inherit.aes = FALSE
    )
  )
  if (!is.null(node_highlight)) {
    node_highlight <- node[
      node[["node_name"]] %in% node_highlight, ,
      drop = FALSE
    ]
    node_layer <- c(
      node_layer,
      list(
        geom_point(
          data = node_highlight,
          aes(x = x, y = y, size = node_size * 1.3),
          color = node_highlight_color,
          show.legend = FALSE,
          inherit.aes = FALSE
        ),
        geom_point(
          data = node_highlight,
          aes(
            x = x,
            y = y,
            size = node_size,
            color = node_group,
            alpha = node_alpha
          ),
          inherit.aes = FALSE
        )
      )
    )
  }
  node_layer <- c(
    node_layer,
    list(
      scale_color_manual(
        name = node_group,
        values = palette_scop(
          node[["node_group"]],
          palette = node_palette,
          palcolor = node_palcolor
        ),
        labels = label_use,
        guide = guide_legend(
          title.hjust = 0,
          order = 1,
          override.aes = list(size = 4, alpha = 1)
        )
      ),
      scale_size,
      scale_alpha
    )
  )

  if (isTRUE(label)) {
    if (isTRUE(label_insitu)) {
      node[, "label"] <- node[["node_group"]]
    } else {
      node[, "label"] <- as.numeric(node[["node_group"]])
    }
    if (isTRUE(label_repel)) {
      label_layer <- list(
        geom_point(
          data = node,
          mapping = aes(x = .data[["x"]], y = .data[["y"]]),
          color = label_point_color,
          size = label_point_size,
          inherit.aes = FALSE
        ),
        geom_text_repel(
          data = node,
          aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
          fontface = "bold",
          min.segment.length = 0,
          segment.color = label_segment_color,
          point.size = label_point_size,
          max.overlaps = 100,
          force = label_repulsion,
          color = label.fg,
          bg.color = label.bg,
          bg.r = label.bg.r,
          size = label.size,
          inherit.aes = FALSE
        )
      )
    } else {
      label_layer <- list(geom_text_repel(
        data = node,
        aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
        fontface = "bold",
        min.segment.length = 0,
        segment.color = label_segment_color,
        point.size = NA,
        max.overlaps = 100,
        force = 0,
        color = label.fg,
        bg.color = label.bg,
        bg.r = label.bg.r,
        size = label.size,
        inherit.aes = FALSE
      ))
    }
  } else {
    label_layer <- NULL
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
      edge_layer = edge_layer,
      trans_layer = trans_layer,
      node_layer = node_layer,
      label_layer = label_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  } else {
    return(
      ggplot() +
        edge_layer +
        trans_layer +
        node_layer +
        label_layer +
        lab_layer +
        theme_layer
    )
  }
}

#' Shorten and offset the segment
#'
#' This function takes a data frame representing segments in a plot and shortens
#' and offsets them based on the provided arguments.
#'
#' @param data A data frame containing the segments. It should have columns
#'             'x', 'y', 'xend', and 'yend' representing the start and end points
#'             of each segment.
#' @param shorten_start The amount to shorten the start of each segment by.
#' @param shorten_end The amount to shorten the end of each segment by.
#' @param offset The amount to offset each segment by.
#'
#' @return The modified data frame with the shortened and offset segments.
#'
#' @examples
#' library(ggplot2)
#' tempNodes <- data.frame(
#'   "x" = c(10, 40),
#'   "y" = c(10, 30)
#' )
#' data <- data.frame(
#'   "x" = c(10, 40),
#'   "y" = c(10, 30),
#'   "xend" = c(40, 10),
#'   "yend" = c(30, 10)
#' )
#' ggplot(tempNodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(
#'     data = data,
#'     aes(x = x, xend = xend, y = y, yend = yend)
#'   )
#'
#' ggplot(tempNodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(
#'     data = segementsDf(
#'       data,
#'       shorten_start = 2,
#'       shorten_end = 3,
#'       offset = 1
#'     ),
#'     aes(x = x, xend = xend, y = y, yend = yend)
#'   )
#' @export
segementsDf <- function(data, shorten_start, shorten_end, offset) {
  data$dx <- data$xend - data$x
  data$dy <- data$yend - data$y
  data$dist <- sqrt(data$dx^2 + data$dy^2)
  data$px <- data$dx / data$dist
  data$py <- data$dy / data$dist

  data$x <- data$x + data$px * shorten_start
  data$y <- data$y + data$py * shorten_start
  data$xend <- data$xend - data$px * shorten_end
  data$yend <- data$yend - data$py * shorten_end
  data$x <- data$x - data$py * offset
  data$xend <- data$xend - data$py * offset
  data$y <- data$y + data$px * offset
  data$yend <- data$yend + data$px * offset

  return(data)
}

fc_matrix <- function(matrix) {
  matrix / rowMeans(matrix)
}

zscore_matrix <- function(matrix, ...) {
  t(scale(t(matrix), ...))
}

log2fc_matrix <- function(matrix) {
  log2(matrix / rowMeans(matrix))
}

log1p_matrix <- function(matrix) {
  log1p(matrix)
}

matrix_process <- function(
    matrix,
    method = c("raw", "zscore", "fc", "log2fc", "log1p"),
    ...) {
  if (is.function(method)) {
    matrix_processed <- method(matrix, ...)
  } else if (method == "raw") {
    matrix_processed <- matrix
  } else if (method == "fc") {
    matrix_processed <- fc_matrix(matrix)
  } else if (method == "zscore") {
    matrix_processed <- zscore_matrix(matrix, ...)
  } else if (method == "log2fc") {
    matrix_processed <- log2fc_matrix(matrix)
  } else if (method == "log1p") {
    matrix_processed <- log1p_matrix(matrix)
  }
  if (!identical(dim(matrix_processed), dim(matrix))) {
    stop("The dimensions of the matrix are changed after processing")
  }
  return(matrix_processed)
}

extractgrobs <- function(vlnplots, x_nm, y_nm, x, y) {
  grobs <- vlnplots[paste0(x_nm[x], ":", y_nm[y])]
  if (length(grobs) == 1) {
    grobs <- grobs[[1]]
  }
  return(grobs)
}

#' @importFrom grid viewport grid.draw is.grob
grid_draw <- function(groblist, x, y, width, height) {
  if (is.grob(groblist)) {
    groblist <- list(groblist)
  }
  for (i in seq_along(groblist)) {
    groblist[[i]]$vp <- viewport(
      x = x[i],
      y = y[i],
      width = width[i],
      height = height[i]
    )
    grid.draw(groblist[[i]])
  }
}

#' @importFrom stats hclust as.dendrogram order.dendrogram
#' @importFrom proxyC dist
#' @importFrom ComplexHeatmap merge_dendrogram
cluster_within_group2 <- function(mat, factor) {
  check_r("dendextend")
  if (!is.factor(factor)) {
    factor <- factor(factor, levels = unique(factor))
  }
  dend_list <- list()
  order_list <- list()
  for (le in unique(levels(factor))) {
    m <- mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] <- which(factor == le)
      dend_list[[le]] <- structure(
        which(factor == le),
        class = "dendrogram",
        leaf = TRUE, # height = 0,
        label = 1,
        members = 1
      )
    } else if (ncol(m) > 1) {
      hc1 <- hclust(as.dist(dist(t(m))))
      dend_list[[le]] <- as.dendrogram(hc1)
      order_list[[le]] <- which(factor == le)[order.dendrogram(dend_list[[le]])]
      dendextend::order.dendrogram(dend_list[[le]]) <- order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") <- le
  }
  parent <- as.dendrogram(hclust(as.dist(dist(t(sapply(
    order_list,
    function(x) rowMeans(mat[, x, drop = FALSE])
  ))))))
  dend_list <- lapply(dend_list, function(dend) {
    dendrapply(
      dend,
      function(node) {
        if (is.null(attr(node, "height"))) {
          attr(node, "height") <- 0
        }
        node
      }
    )
  })
  # print(sapply(dend_list, function(x) attr(x, "height")))
  dend <- merge_dendrogram(parent, dend_list)
  order.dendrogram(dend) <- unlist(order_list[order.dendrogram(parent)])
  return(dend)
}

#' @importFrom ComplexHeatmap HeatmapAnnotation anno_empty anno_block anno_textbox
#' @importFrom grid gpar unit
#' @importFrom dplyr %>% filter group_by arrange desc across reframe mutate distinct n .data "%>%"
heatmap_enrichment <- function(
    geneID,
    geneID_groups,
    feature_split_palette = "simspec",
    feature_split_palcolor = NULL,
    ha_right = NULL,
    flip = FALSE,
    anno_terms = FALSE,
    anno_keys = FALSE,
    anno_features = FALSE,
    terms_width = unit(4, "in"),
    terms_fontsize = 8,
    keys_width = unit(2, "in"),
    keys_fontsize = c(6, 10),
    features_width = unit(2, "in"),
    features_fontsize = c(6, 10),
    IDtype = "symbol",
    species = "Homo_sapiens",
    db_update = FALSE,
    db_combine = FALSE,
    db_version = "latest",
    convert_species = FALSE,
    Ensembl_version = 103,
    mirror = NULL,
    db = "GO_BP",
    TERM2GENE = NULL,
    TERM2NAME = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    GO_simplify = FALSE,
    GO_simplify_cutoff = "p.adjust < 0.05",
    simplify_method = "Wang",
    simplify_similarityCutoff = 0.7,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = 5,
    show_termid = FALSE,
    topWord = 20,
    words_excluded = NULL) {
  res <- NULL
  words_excluded <- words_excluded %||% scop::words_excluded

  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    if (isTRUE(flip)) {
      stop(
        "anno_keys, anno_features and anno_terms can only be used when flip is FALSE."
      )
    }
    if (all(is.na(geneID_groups))) {
      geneID_groups <- rep(1, length(geneID))
    }
    if (!is.factor(geneID_groups)) {
      geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
    }
    fill_split <- palette_scop(
      levels(geneID_groups),
      type = "discrete",
      palette = feature_split_palette,
      palcolor = feature_split_palcolor
    )[levels(geneID_groups) %in% geneID_groups]
    res <- RunEnrichment(
      geneID = geneID,
      geneID_groups = geneID_groups,
      IDtype = IDtype,
      species = species,
      db_update = db_update,
      db_version = db_version,
      db_combine = db_combine,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      db = db,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      GO_simplify = GO_simplify,
      GO_simplify_cutoff = GO_simplify_cutoff,
      simplify_method = simplify_method,
      simplify_similarityCutoff = simplify_similarityCutoff
    )
    if (isTRUE(db_combine)) {
      db <- "Combined"
    }
    if (isTRUE(GO_simplify) && any(db %in% c("GO_BP", "GO_CC", "GO_MF"))) {
      db[db %in% c("GO_BP", "GO_CC", "GO_MF")] <- paste0(
        db[db %in% c("GO_BP", "GO_CC", "GO_MF")],
        "_sim"
      )
    }
    if (nrow(res$enrichment) == 0) {
      warning("No enrichment result found.", immediate. = TRUE)
    } else {
      metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
      metric_value <- ifelse(
        is.null(padjustCutoff),
        pvalueCutoff,
        padjustCutoff
      )
      pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
      padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

      df <- res$enrichment
      df <- df[df[["Database"]] %in% db, , drop = FALSE]
      df <- df[df[[metric]] < metric_value, , drop = FALSE]
      df <- df[order(df[[metric]]), , drop = FALSE]
      if (nrow(df) == 0) {
        warning(
          "No term enriched using the threshold: ",
          paste0("pvalueCutoff = ", pvalueCutoff),
          "; ",
          paste0("padjustCutoff = ", padjustCutoff),
          immediate. = TRUE
        )
      } else {
        df_list <- split.data.frame(df, ~ Database + Groups)
        df_list <- df_list[lapply(df_list, nrow) > 0]

        for (enrich in db) {
          nm <- strsplit(names(df_list), "\\.")
          subdf_list <- df_list[
            unlist(lapply(nm, function(x) x[[1]])) %in% enrich
          ]
          if (length(subdf_list) == 0) {
            warning(
              "No ",
              enrich,
              " term enriched using the threshold: ",
              paste0("pvalueCutoff = ", pvalueCutoff),
              "; ",
              paste0("padjustCutoff = ", padjustCutoff),
              immediate. = TRUE
            )
            next
          }
          nm <- strsplit(names(subdf_list), "\\.")

          ha_terms <- NULL
          if (isTRUE(anno_terms)) {
            terms_list <- lapply(subdf_list, function(df) {
              if (isTRUE(show_termid)) {
                terms <- paste(
                  head(df$ID, topTerm),
                  head(df$Description, topTerm)
                )
              } else {
                terms <- head(df$Description, topTerm)
                terms <- capitalize(terms)
              }
              df_out <- data.frame(keyword = terms)
              df_out[["col"]] <- palette_scop(
                -log10(head(df[, metric], topTerm)),
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df_out[["col"]] <- sapply(
                df_out[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df_out[["fontsize"]] <- rep(terms_fontsize, nrow(df_out))
              return(df_out)
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(terms_list))) > 0) {
              ha_terms <- HeatmapAnnotation(
                "terms_empty" = anno_empty(
                  width = unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "terms_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "terms" = anno_textbox(
                  align_to = geneID_groups,
                  text = terms_list,
                  max_width = terms_width,
                  word_wrap = TRUE,
                  add_new_line = TRUE,
                  background_gp = gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = unit(0, "points")
              )
              names(ha_terms) <- paste0(names(ha_terms), "_", enrich)
            }
          }

          ha_keys <- NULL
          if (isTRUE(anno_keys)) {
            check_r("jokergoo/simplifyEnrichment")
            keys_list <- lapply(subdf_list, function(df) {
              if (all(df$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
                df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df[[
                  "ID"
                ]])
                if (nrow(df0) > 0) {
                  df <- df0 %>%
                    reframe(
                      keyword = .data[["keyword"]],
                      score = -(log10(.data[["padj"]])),
                      count = .data[["n_term"]],
                      Database = df[["Database"]][1],
                      Groups = df[["Groups"]][1]
                    ) %>%
                    filter(
                      !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                    ) %>%
                    filter(nchar(.data[["keyword"]]) >= 1) %>%
                    filter(
                      !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                    ) %>%
                    distinct() %>%
                    mutate(
                      angle = 90 *
                        sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
                    ) %>%
                    as.data.frame()
                  df <- df[
                    head(order(df[["score"]], decreasing = TRUE), topWord), ,
                    drop = FALSE
                  ]
                } else {
                  df <- NULL
                }
              } else {
                df <- df %>%
                  mutate(
                    keyword = strsplit(
                      tolower(as.character(.data[["Description"]])),
                      " "
                    )
                  ) %>%
                  unnest(cols = "keyword") %>%
                  group_by(.data[["keyword"]], Database, Groups) %>%
                  reframe(
                    keyword = .data[["keyword"]],
                    score = sum(-(log10(.data[[metric]]))),
                    count = n(),
                    Database = .data[["Database"]],
                    Groups = .data[["Groups"]],
                    .groups = "keep"
                  ) %>%
                  filter(
                    !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                  ) %>%
                  filter(nchar(.data[["keyword"]]) >= 1) %>%
                  filter(
                    !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                  ) %>%
                  distinct() %>%
                  mutate(
                    angle = 90 *
                      sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
                  ) %>%
                  as.data.frame()
                df <- df[
                  head(order(df[["score"]], decreasing = TRUE), topWord), ,
                  drop = FALSE
                ]
              }
              if (isTRUE(nrow(df) > 0)) {
                df[["col"]] <- palette_scop(
                  df[, "score"],
                  type = "continuous",
                  palette = "Spectral",
                  matched = TRUE
                )
                df[["col"]] <- sapply(
                  df[["col"]],
                  function(x) blendcolors(c(x, "black"))
                )
                df[["fontsize"]] <- rescale(df[, "count"], to = keys_fontsize)
                return(df)
              } else {
                return(NULL)
              }
            })
            names(keys_list) <- unlist(lapply(nm, function(x) x[[2]]))
            keys_list <- keys_list[lapply(keys_list, length) > 0]
            if (length(intersect(geneID_groups, names(keys_list))) > 0) {
              ha_keys <- HeatmapAnnotation(
                "keys_empty" = anno_empty(
                  width = unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "keys_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "keys" = anno_textbox(
                  align_to = geneID_groups,
                  text = keys_list,
                  max_width = keys_width,
                  background_gp = gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = unit(0, "points")
              )
              names(ha_keys) <- paste0(names(ha_keys), "_", enrich)
            }
          }

          ha_features <- NULL
          if (isTRUE(anno_features)) {
            features_list <- lapply(subdf_list, function(df) {
              df <- df %>%
                mutate(
                  keyword = strsplit(as.character(.data[["geneID"]]), "/")
                ) %>%
                unnest(cols = "keyword") %>%
                group_by(.data[["keyword"]], Database, Groups) %>%
                reframe(
                  keyword = .data[["keyword"]],
                  score = sum(-(log10(.data[[metric]]))),
                  count = n(),
                  Database = .data[["Database"]],
                  Groups = .data[["Groups"]],
                  .groups = "keep"
                ) %>%
                distinct() %>%
                mutate(
                  angle = 90 *
                    sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))
                ) %>%
                as.data.frame()
              df <- df[
                head(order(df[["score"]], decreasing = TRUE), topWord), ,
                drop = FALSE
              ]
              df[["col"]] <- palette_scop(
                df[, "score"],
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df[["col"]] <- sapply(
                df[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df[["fontsize"]] <- rescale(df[, "count"], to = features_fontsize)
              return(df)
            })
            names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(features_list))) > 0) {
              ha_features <- HeatmapAnnotation(
                "features_empty" = anno_empty(
                  width = unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "features_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "features" = anno_textbox(
                  align_to = geneID_groups,
                  text = features_list,
                  max_width = features_width,
                  background_gp = gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = unit(0, "points")
              )
              names(ha_features) <- paste0(names(ha_features), "_", enrich)
            }
          }

          ha_enrichment <- list(ha_terms, ha_keys, ha_features)
          ha_enrichment <- ha_enrichment[sapply(ha_enrichment, length) > 0]
          ha_enrichment <- do.call(c, ha_enrichment)

          if (is.null(ha_right)) {
            ha_right <- ha_enrichment
          } else {
            ha_right <- c(ha_right, ha_enrichment)
          }
        }
      }
    }
  }
  return(list(ha_right = ha_right, res = res))
}

#' @importFrom grid convertWidth convertHeight unit
#' @importFrom ComplexHeatmap width.HeatmapAnnotation height.HeatmapAnnotation width.Legends
heatmap_rendersize <- function(
    width,
    height,
    units,
    ha_top_list,
    ha_left,
    ha_right,
    ht_list,
    legend_list,
    flip) {
  width_annotation <- height_annotation <- 0
  if (isTRUE(flip)) {
    width_sum <- width[1] %||%
      convertWidth(unit(1, "in"), units, valueOnly = TRUE)
    height_sum <- sum(
      height %||% convertHeight(unit(1, "in"), units, valueOnly = TRUE)
    )
    if (length(ha_top_list) > 0) {
      width_annotation <- convertWidth(
        unit(width_annotation, units) +
          width.HeatmapAnnotation(ha_top_list[[1]]),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_left)) {
      height_annotation <- convertHeight(
        unit(height_annotation, units) + height.HeatmapAnnotation(ha_left),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_right)) {
      height_annotation <- convertHeight(
        unit(height_annotation, units) + height.HeatmapAnnotation(ha_right),
        units,
        valueOnly = TRUE
      )
    }
  } else {
    width_sum <- sum(
      width %||% convertWidth(unit(1, "in"), units, valueOnly = TRUE)
    )
    height_sum <- height[1] %||%
      convertHeight(unit(1, "in"), units, valueOnly = TRUE)
    if (length(ha_top_list) > 0) {
      height_annotation <- convertHeight(
        unit(height_annotation, units) +
          height.HeatmapAnnotation(ha_top_list[[1]]),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_left)) {
      width_annotation <- convertWidth(
        unit(width_annotation, units) + width.HeatmapAnnotation(ha_left),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_right)) {
      width_annotation <- convertWidth(
        unit(width_annotation, units) + width.HeatmapAnnotation(ha_right),
        units,
        valueOnly = TRUE
      )
    }
  }
  dend_width <- name_width <- NULL
  dend_height <- name_height <- NULL
  if (inherits(ht_list, "HeatmapList")) {
    for (nm in names(ht_list@ht_list)) {
      ht <- ht_list@ht_list[[nm]]
      dend_width <- max(ht@row_dend_param$width, dend_width)
      dend_height <- max(ht@column_dend_param$height, dend_height)
      name_width <- max(ht@row_names_param$max_width, name_width)
      name_height <- max(ht@column_names_param$max_height, name_height)
    }
  } else if (inherits(ht_list, "Heatmap")) {
    ht <- ht_list
    dend_width <- max(ht@row_dend_param$width, dend_width)
    dend_height <- max(ht@column_dend_param$height, dend_height)
    name_width <- max(ht@row_names_param$max_width, name_width)
    name_height <- max(ht@column_names_param$max_height, name_height)
  } else {
    stop("ht_list is not a class of HeatmapList or Heatmap.")
  }

  lgd_width <- grid::convertWidth(
    grid::unit(
      unlist(lapply(legend_list, width.Legends)),
      grid::unitType(width.Legends(legend_list[[1]]))
    ),
    unitTo = units,
    valueOnly = TRUE
  )
  width_sum <- grid::convertWidth(
    grid::unit(width_sum, units) +
      grid::unit(width_annotation, units) +
      dend_width +
      name_width,
    units,
    valueOnly = TRUE
  ) +
    sum(lgd_width)
  height_sum <- max(
    convertHeight(
      unit(height_sum, units) +
        unit(height_annotation, units) +
        dend_height +
        name_height,
      units,
      valueOnly = TRUE
    ),
    convertHeight(unit(0.95, "npc"), units, valueOnly = TRUE)
  )
  return(list(width_sum = width_sum, height_sum = height_sum))
}

#' @importFrom stats sd
standardise <- function(data) {
  data[] <- t(apply(data, 1, scale))
  return(data)
}

mestimate <- function(data) {
  N <- nrow(data)
  D <- ncol(data)
  m.sj <- 1 +
    (1418 / N + 22.05) * D^(-2) +
    (12.33 / N + 0.243) *
      D^(-0.0406 * log(N) - 0.1134)
  return(m.sj)
}

adjustlayout <- function(
    graph,
    layout,
    width,
    height = 2,
    scale = 100,
    iter = 100) {
  w <- width / 2
  layout[, 1] <- layout[, 1] / diff(range(layout[, 1])) * scale
  layout[, 2] <- layout[, 2] / diff(range(layout[, 2])) * scale

  adjusted <- c()
  for (v in order(igraph::degree(graph), decreasing = TRUE)) {
    adjusted <- c(adjusted, v)
    neighbors <- as.numeric(neighbors(graph, V(graph)[v]))
    neighbors <- setdiff(neighbors, adjusted)
    x <- layout[v, 1]
    y <- layout[v, 2]
    r <- w[v]
    for (neighbor in neighbors) {
      nx <- layout[neighbor, 1]
      ny <- layout[neighbor, 2]
      ndist <- sqrt((nx - x)^2 + (ny - y)^2)
      nr <- w[neighbor]
      expect <- r + nr
      if (ndist < expect) {
        dx <- (x - nx) * (expect - ndist) / ndist
        dy <- (y - ny) * (expect - ndist) / ndist
        layout[neighbor, 1] <- nx - dx
        layout[neighbor, 2] <- ny - dy
        adjusted <- c(adjusted, neighbor)
      }
    }
  }

  for (i in seq_len(iter)) {
    dist_matrix <- Matrix::as.matrix(dist(layout))
    nearest_neighbors <- apply(
      dist_matrix,
      2,
      function(x) which(x == min(x[x > 0])),
      simplify = FALSE
    )

    for (v in sample(seq_len(nrow(layout)))) {
      neighbors <- unique(nearest_neighbors[[v]])
      x <- layout[v, 1]
      y <- layout[v, 2]
      r <- w[v]
      for (neighbor in neighbors) {
        nx <- layout[neighbor, 1]
        ny <- layout[neighbor, 2]
        nr <- w[neighbor]
        if (abs(nx - x) < (r + nr) && abs(ny - y) < height) {
          dx <- r + nr - (nx - x)
          dy <- height - (ny - y)
          if (sample(c(1, 0), 1) == 1) {
            dx <- 0
          } else {
            dy <- 0
          }
          layout[neighbor, 1] <- nx - dx
          layout[neighbor, 2] <- ny - dy
        }
      }
    }
  }
  return(layout)
}
