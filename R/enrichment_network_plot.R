enrichment_gene_network <- function(
  df, df_unnest, network_layout, network_layoutadjust, network_adjscale,
  network_adjiter, palette, palcolor, network_blendmode, network_labelsize
) {
  nodes <- rbind(
    data.frame(
      "ID" = df[["Description"]],
      class = "term",
      metric = df[["metric"]]
    ),
    data.frame("ID" = unique(df_unnest$geneID), class = "gene", metric = 0)
  )
  nodes$Database <- df$Database[1]
  nodes$Groups <- df$Groups[1]
  edges <- as.data.frame(
    df_unnest[, c("Description", "geneID")]
  )
  colnames(edges) <- c("from", "to")
  edges[["weight"]] <- 1
  graph <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  if (network_layout %in% c("circle", "tree", "grid")) {
    layout <- switch(network_layout,
      "circle" = igraph::layout_in_circle(graph),
      "tree" = igraph::layout_as_tree(graph),
      "grid" = igraph::layout_on_grid(graph)
    )
  } else {
    layout_fun <- get_namespace_fun(
      "igraph", paste0("layout_with_", network_layout)
    )
    layout <- layout_fun(graph)
  }
  df_graph <- igraph::as_data_frame(graph, what = "both")

  df_nodes <- df_graph$vertices
  if (isTRUE(network_layoutadjust)) {
    width <- nchar(df_nodes$name)
    width[df_nodes$class == "term"] <- 8
    layout <- adjustlayout(
      graph = graph,
      layout = layout,
      width = width,
      height = 2,
      scale = network_adjscale,
      iter = network_adjiter
    )
  }
  df_nodes[["dim1"]] <- layout[, 1]
  df_nodes[["dim2"]] <- layout[, 2]

  df_edges <- df_graph$edges
  df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
  df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
  df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
  df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

  colors <- palette_colors(
    levels(df[["Description"]]),
    palette = palette,
    palcolor = palcolor
  )
  df_edges[["color"]] <- colors[df_edges$from]
  node_colors <- stats::aggregate(
    df_unnest$Description,
    by = list(df_unnest$geneID),
    FUN = function(x) {
      blendcolors(colors = colors[x], mode = network_blendmode)
    }
  )
  colors <- c(colors, stats::setNames(node_colors[, 2], node_colors[, 1]))
  label_colors <- ifelse(
    Matrix::colSums(grDevices::col2rgb(colors)) > 255 * 2,
    "black",
    "white"
  )
  df_nodes[["color"]] <- colors[df_nodes$name]
  df_nodes[["label_color"]] <- label_colors[df_nodes$name]
  df_nodes[["label"]] <- NA
  df_nodes[levels(df[["Description"]]), "label"] <- seq_len(
    nlevels(
      df[[
        "Description"
      ]]
    )
  )

  check_r("shadowtext")
  shadowtextGrob <- get_namespace_fun("shadowtext", "shadowtextGrob")
  draw_key_cust <- function(data, params, size) {
    label <- which(
      levels(df[["Description"]]) %in%
        names(colors)[colors == data$fill]
    )
    grid::grobTree(
      draw_key_point(data, list(color = "white", shape = 21)),
      shadowtextGrob(
        label = label,
        bg.colour = "black",
        bg.r = 0.1,
        gp = grid::gpar(col = "white", fontface = "bold")
      )
    )
  }

  p <- ggplot() +
    geom_segment(
      data = df_edges,
      aes(
        x = from_dim1,
        y = from_dim2,
        xend = to_dim1,
        yend = to_dim2,
        color = color
      ),
      alpha = 1,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_label(
      data = df_nodes[df_nodes$class == "gene", ],
      aes(
        x = dim1,
        y = dim2,
        label = name,
        fill = color,
        color = label_color
      ),
      size = 3,
      show.legend = FALSE
    ) +
    geom_point(
      data = df_nodes[df_nodes$class == "term", ],
      aes(x = dim1, y = dim2),
      size = 8,
      color = "black",
      fill = "black",
      stroke = 1,
      shape = 21,
      show.legend = FALSE
    ) +
    geom_point(
      data = df_nodes[df_nodes$class == "term", ],
      aes(x = dim1, y = dim2, fill = color),
      size = 7,
      color = "white",
      stroke = 1,
      shape = 21,
      key_glyph = draw_key_cust
    ) +
    ggrepel::geom_text_repel(
      data = df_nodes[df_nodes$class == "term", ],
      aes(x = dim1, y = dim2, label = label),
      fontface = "bold",
      min.segment.length = 0,
      segment.color = "black",
      point.size = NA,
      max.overlaps = 100,
      force = 0,
      color = "white",
      bg.color = "black",
      bg.r = 0.1,
      size = network_labelsize
    ) +
    scale_color_identity(guide = "none") +
    scale_fill_identity(
      name = "Term:",
      guide = "legend",
      labels = levels(df[["Description"]]),
      breaks = colors[levels(df[["Description"]])]
    )
  p
}
