enrichment_map_plot <- function(
  df, enrichmap_layout, enrichmap_cluster, enrichmap_show_keyword,
  enrichmap_label, enrichmap_nlabel, character_width, words_excluded,
  enrichmap_mark, enrichmap_labelsize, enrichmap_expand, palette, palcolor
) {
  nodes <- df
  edges <- enrichment_overlap_edges(nodes$ID, df[["geneID"]])
  graph <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  if (enrichmap_layout %in% c("circle", "tree", "grid")) {
    layout <- switch(enrichmap_layout,
      "circle" = igraph::layout_in_circle(graph),
      "tree" = igraph::layout_as_tree(graph),
      "grid" = igraph::layout_on_grid(graph)
    )
  } else {
    layout_fun <- get_namespace_fun(
      "igraph", paste0("layout_with_", enrichmap_layout)
    )
    layout <- layout_fun(graph)
  }
  cluster_fun <- get_namespace_fun(
    "igraph", paste0("cluster_", enrichmap_cluster)
  )
  clusters <- cluster_fun(graph)
  df_graph <- igraph::as_data_frame(graph, what = "both")

  df_nodes <- df_graph$vertices
  df_nodes[["dim1"]] <- layout[, 1]
  df_nodes[["dim2"]] <- layout[, 2]
  df_nodes[["clusters"]] <- factor(
    paste0("C", clusters$membership),
    paste0("C", unique(sort(clusters$membership)))
  )

  if (isTRUE(enrichmap_show_keyword)) {
    df_keyword1 <- df_nodes |>
      dplyr::mutate(
        keyword = strsplit(
          tolower(as.character(.data[["Description"]])),
          "\\s|\\n",
          perl = TRUE
        )
      ) |>
      thisutils::unnest_fun(cols = "keyword") |>
      dplyr::group_by(.data[["keyword"]], Database, Groups, clusters) |>
      dplyr::reframe(
        keyword = capitalize(.data[["keyword"]]),
        score = sum(.data[["metric"]]),
        count = dplyr::n(),
        Database = .data[["Database"]],
        Groups = .data[["Groups"]],
        .groups = "keep"
      ) |>
      dplyr::filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) |>
      dplyr::filter(nchar(.data[["keyword"]]) >= 1) |>
      dplyr::filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) |>
      dplyr::distinct() |>
      dplyr::group_by(Database, Groups, clusters) |>
      dplyr::arrange(dplyr::desc(score)) |>
      dplyr::slice_head(n = enrichmap_nlabel) |>
      dplyr::reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) |>
      as.data.frame()
    rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
    df_keyword1[["keyword"]] <- stringr::str_wrap(
      df_keyword1[["keyword"]],
      width = character_width
    )
    df_keyword1[["label"]] <- paste0(
      df_keyword1[["clusters"]],
      ":\n",
      df_keyword1[["keyword"]]
    )
  } else {
    if (enrichmap_label == "term") {
      df_nodes[["Description"]] <- stringr::str_wrap(
        df_nodes[["Description"]],
        width = character_width
      )
    }
    df_keyword1 <- df_nodes |>
      dplyr::group_by(Database, Groups, clusters) |>
      dplyr::arrange(dplyr::desc(metric)) |>
      dplyr::reframe(keyword = Description) |>
      dplyr::distinct() |>
      dplyr::group_by(Database, Groups, clusters) |>
      dplyr::slice_head(n = enrichmap_nlabel) |>
      dplyr::reframe(keyword = paste0(.data[["keyword"]], collapse = "\n")) |>
      as.data.frame()
    rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
    df_keyword1[["label"]] <- paste0(
      df_keyword1[["clusters"]],
      ":\n",
      df_keyword1[["keyword"]]
    )
  }

  df_keyword2 <- df_nodes |>
    dplyr::mutate(keyword = .data[["geneID"]]) |>
    thisutils::unnest_fun(cols = "keyword") |>
    dplyr::group_by(.data[["keyword"]], Database, Groups, clusters) |>
    dplyr::reframe(
      keyword = .data[["keyword"]],
      score = sum(.data[["metric"]]),
      count = dplyr::n(),
      Database = .data[["Database"]],
      Groups = .data[["Groups"]],
      .groups = "keep"
    ) |>
    dplyr::distinct() |>
    dplyr::group_by(Database, Groups, clusters) |>
    dplyr::arrange(dplyr::desc(score)) |>
    dplyr::slice_head(n = enrichmap_nlabel) |>
    dplyr::reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) |>
    as.data.frame()
  rownames(df_keyword2) <- as.character(df_keyword2[["clusters"]])
  df_keyword2[["keyword"]] <- stringr::str_wrap(
    df_keyword2[["keyword"]],
    width = character_width
  )
  df_keyword2[["label"]] <- paste0(
    df_keyword2[["clusters"]],
    ":\n",
    df_keyword2[["keyword"]]
  )

  df_nodes[["keyword1"]] <- df_keyword1[
    as.character(df_nodes$clusters),
    "keyword"
  ]
  df_nodes[["keyword2"]] <- df_keyword2[
    as.character(df_nodes$clusters),
    "keyword"
  ]

  df_edges <- df_graph$edges
  df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
  df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
  df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
  df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

  check_r("ggforce", verbose = FALSE)
  if (enrichmap_mark == "hull") {
    check_r("concaveman", verbose = FALSE)
  }
  mark_layer <- do.call(
    switch(enrichmap_mark,
      "ellipse" = ggforce::geom_mark_ellipse,
      "hull" = ggforce::geom_mark_hull
    ),
    list(
      data = df_nodes,
      aes(
        x = dim1,
        y = dim2,
        color = clusters,
        fill = clusters,
        label = clusters,
        description = if (enrichmap_label == "term") keyword1 else keyword2
      ),
      expand = grid::unit(3, "mm"),
      alpha = 0.1,
      label.margin = margin(1, 1, 1, 1, "mm"),
      label.fontsize = enrichmap_labelsize * 2,
      label.fill = "grey95",
      label.minwidth = grid::unit(character_width, "in"),
      label.buffer = grid::unit(0, "mm"),
      con.size = 1,
      con.cap = 0
    )
  )

  p <- ggplot() +
    mark_layer +
    geom_segment(
      data = df_edges,
      aes(
        x = from_dim1,
        y = from_dim2,
        xend = to_dim1,
        yend = to_dim2,
        linewidth = weight
      ),
      alpha = 0.1,
      lineend = "round"
    ) +
    geom_point(
      data = df_nodes,
      aes(x = dim1, y = dim2, size = Count, fill = clusters),
      color = "black",
      shape = 21
    ) +
    labs(x = "", y = "") +
    scale_size(
      name = "Count",
      range = c(2, 6),
      scales::breaks_extended(n = 4)
    ) +
    guides(
      size = guide_legend(
        override.aes = list(fill = "grey30", shape = 21),
        order = 1
      )
    ) +
    scale_linewidth(
      name = "Intersection",
      range = c(0.3, 3),
      scales::breaks_extended(n = 4)
    ) +
    guides(
      linewidth = guide_legend(
        override.aes = list(alpha = 1, color = "grey"),
        order = 2
      )
    ) +
    scale_fill_manual(
      name = switch(enrichmap_label,
        "term" = "Feature:",
        "feature" = "Term:"
      ),
      values = palette_colors(
        levels(df_nodes[["clusters"]]),
        palette = palette,
        palcolor = palcolor
      ),
      labels = if (enrichmap_label == "term") {
        df_keyword2[levels(df_nodes[["clusters"]]), "label"]
      } else {
        df_keyword1[levels(df_nodes[["clusters"]]), "label"]
      },
      na.value = "grey80",
      aesthetics = c("colour", "fill")
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(alpha = 1, color = "black", shape = NA),
        byrow = TRUE,
        order = 3
      )
    ) +
    guides(color = guide_none()) +
    scale_x_continuous(
      expand = expansion(c(enrichmap_expand[1], enrichmap_expand[1]), 0)
    ) +
    scale_y_continuous(
      expand = expansion(c(enrichmap_expand[2], enrichmap_expand[2]), 0)
    )
  p
}
