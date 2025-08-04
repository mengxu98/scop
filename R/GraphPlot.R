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
#' @seealso \link{CellDimPlot}
#'
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
    transition_arrow_length = grid::unit(0.02, "npc"),
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
    log_message(
      "'node' must be a data.frame object.",
      message_type = "error"
    )
  }
  if (!is.matrix(edge)) {
    log_message(
      "'edge' must be a matrix object.",
      message_type = "error"
    )
  }
  if (!identical(nrow(edge), ncol(edge))) {
    log_message(
      "nrow and ncol is not identical in edge matrix",
      message_type = "error"
    )
  }
  if (!identical(nrow(edge), nrow(node))) {
    log_message(
      "nrow is not identical between edge and node.",
      message_type = "error"
    )
  }
  if (!identical(rownames(edge), rownames(node))) {
    log_message(
      "rownames of node is not identical with edge matrix. They will correspond according to the order.",
      message_type = "warning"
    )
    colnames(edge) <- rownames(edge) <- rownames(node) <- rownames(node) %||%
      colnames(edge) %||%
      rownames(edge)
  }
  if (!all(node_coord %in% colnames(node))) {
    log_message(
      "Cannot find the node_coord ",
      paste(node_coord[!node_coord %in% colnames(node)], collapse = ","),
      " in the node column",
      message_type = "error"
    )
  }
  if (!is.null(transition)) {
    if (!identical(nrow(transition), nrow(node))) {
      log_message(
        "nrow is not identical between transition and node.",
        message_type = "error"
      )
    }
    if (!identical(rownames(transition), rownames(node))) {
      log_message(
        "rownames of node is not identical with transition matrix. They will correspond according to the order.",
        message_type = "warning"
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

  if (isTRUE(label) && isFALSE(label_insitu)) {
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
  edge_df <- reshape2::melt(
    edge,
    na.rm = TRUE,
    stringsAsFactors = FALSE
  )
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
    edge_df <- segements_df(
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
        scale_linewidth_continuous(
          range = range(edge_size),
          guide = "none"
        ),
        ggnewscale::new_scale("linewidth")
      )
    )
  }

  if (!is.null(transition)) {
    trans2 <- trans1 <- Matrix::as.matrix(transition)
    trans1[lower.tri(trans1)] <- 0
    trans2[upper.tri(trans2)] <- 0
    trans <- Matrix::t(trans1) - trans2
    trans[abs(trans) <= transition_threshold] <- NA
    trans_df <- reshape2::melt(
      trans,
      na.rm = TRUE,
      stringsAsFactors = FALSE
    )
    if (nrow(trans_df) == 0) {
      trans_layer <- NULL
    } else {
      trans_df <- as.data.frame(
        Matrix::t(
          apply(trans_df, 1, function(x) {
            if (as.numeric(x[3]) < 0) {
              return(c(x[c(2, 1)], -as.numeric(x[3])))
            } else {
              return(x)
            }
          })
        )
      )
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
      trans_df <- segements_df(
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
            arrow = grid::arrow(
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
                arrow = grid::arrow(
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
            arrow = grid::arrow(
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
                arrow = grid::arrow(
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
          ggnewscale::new_scale("linewidth")
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
        ggrepel::geom_text_repel(
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
      label_layer <- list(
        ggrepel::geom_text_repel(
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
        )
      )
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
