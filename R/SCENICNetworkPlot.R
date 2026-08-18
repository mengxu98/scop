scenic_plot_network_graph <- function(
  srt,
  tool_name,
  features = NULL,
  highlight_tf = NULL,
  max_targets = 30,
  max_edges = Inf,
  network_layout = "auto",
  label_nodes = "auto",
  network_label_top_n = 60,
  network_include_regions = TRUE,
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL,
  rank_table = NULL
) {
  network_layout <- scenic_network_resolve_layout(network_layout, default = "kk")
  edge_data <- scenic_network_edges(
    srt = srt,
    tool_name = tool_name,
    tfs = scenic_network_tf_candidates(features),
    max_targets = max_targets,
    max_edges = max_edges,
    include_regions = network_include_regions,
    filter_tfs = FALSE
  )
  network_data <- scenic_network_plot_data(edge_data = edge_data, layout = network_layout)
  label_data <- scenic_network_label_data(
    node_data = network_data[["nodes"]],
    label_nodes = label_nodes,
    highlight_tf = highlight_tf,
    top_n = network_label_top_n,
    default = "tfs"
  )
  plot <- scenic_network_ggplot(
    node_data = network_data[["nodes"]],
    edge_plot = network_data[["edge_plot"]],
    edge_data = network_data[["edges"]],
    label_data = label_data,
    label_nodes = label_nodes,
    layout = network_layout,
    palette = palette,
    palcolor = palcolor,
    title = title %||% "SCENIC TF–target network",
    curvature = if (network_layout %in% c("star", "hub", "kk", "fr")) 0.08 else 0.12,
    rank_table = rank_table
  )
  list(
    plot = plot,
    plots = list(plot),
    data = list(edges = network_data[["edges"]], nodes = network_data[["nodes"]])
  )
}

scenic_plot_network <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  network_tf = NULL,
  highlight_tf = NULL,
  max_targets = 30,
  network_layout = "auto",
  label_nodes = "auto",
  network_include_regions = TRUE,
  palette = "RdYlBu",
  palcolor = NULL,
  combine = TRUE,
  ncol = 3,
  title = NULL,
  rank_table = NULL
) {
  tf_candidates <- scenic_network_tf_candidates(
    network_tf %||% features %||% highlight_tf %||% unique(top_table[["TF"]])
  )
  if (length(tf_candidates) == 0L) {
    log_message(
      "No TFs were available for the SCENIC network plot",
      message_type = "error"
    )
  }
  layout_use <- scenic_network_resolve_layout(network_layout, default = "star")
  split_panels <- identical(layout_use, "star") && length(tf_candidates) > 1L
  panel_palcolor <- palcolor
  tf_groups <- scenic_network_tf_groups(rank_table %||% top_table, tf_candidates)
  if (isTRUE(split_panels)) {
    panel_palcolor <- palette_colors(
      tf_candidates,
      palette = scenic_network_palette(palette),
      palcolor = palcolor
    )
  }

  if (isTRUE(split_panels)) {
    plots <- list()
    edge_list <- list()
    node_list <- list()
    for (tf in tf_candidates) {
      grp <- tf_groups[[tf]]
      panel_title <- if (!is.na(grp) && nzchar(grp)) {
        paste0(tf, " (", grp, ")")
      } else {
        tf
      }
      one <- scenic_plot_one_network(
        srt = srt,
        tool_name = tool_name,
        tfs = tf,
        max_targets = max_targets,
        max_edges = Inf,
        network_layout = "star",
        label_nodes = label_nodes,
        highlight_tf = highlight_tf,
        include_regions = network_include_regions,
        filter_tfs = TRUE,
        palette = palette,
        palcolor = panel_palcolor,
        title = panel_title,
        curvature = 0,
        rank_table = rank_table %||% top_table
      )
      plots[[tf]] <- one[["plot"]]
      edge_list[[tf]] <- one[["data"]][["edges"]]
      node_list[[tf]] <- one[["data"]][["nodes"]]
    }
    plot <- scenic_combine_plots(
      plots,
      combine = combine,
      ncol = ncol,
      title = title %||% "SCENIC regulon networks"
    )
    edges <- scenic_bind_named_dfs(edge_list, "TF")
    nodes <- scenic_bind_named_dfs(node_list, "TF")
    return(list(
      plot = plot,
      plots = plots,
      data = list(edges = edges, nodes = nodes)
    ))
  }

  scenic_plot_one_network(
    srt = srt,
    tool_name = tool_name,
    tfs = tf_candidates,
    max_targets = max_targets,
    max_edges = Inf,
    network_layout = layout_use,
    label_nodes = label_nodes,
    highlight_tf = highlight_tf,
    include_regions = network_include_regions,
    filter_tfs = TRUE,
    palette = palette,
    palcolor = palcolor,
    title = if (!is.null(title)) {
      title
    } else if (length(tf_candidates) == 1L) {
      grp <- tf_groups[[tf_candidates[[1]]]]
      if (!is.na(grp) && nzchar(grp)) {
        paste0(tf_candidates, " (", grp, ")")
      } else {
        tf_candidates
      }
    } else {
      "SCENIC TF–target network"
    },
    curvature = if (layout_use %in% c("star", "hub")) 0 else 0.16,
    rank_table = rank_table %||% top_table
  )
}

scenic_plot_egrn <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  network_tf = NULL,
  highlight_tf = NULL,
  max_targets = 20,
  network_layout = "auto",
  label_nodes = "auto",
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL,
  rank_table = NULL
) {
  tf_candidates <- scenic_network_tf_candidates(
    network_tf %||% features %||% highlight_tf %||% unique(top_table[["TF"]])
  )
  layout_use <- if (is.null(network_layout) || identical(network_layout, "auto")) {
    if (length(tf_candidates) <= 1L) "star" else "kk"
  } else {
    scenic_network_resolve_layout(network_layout, default = "kk")
  }
  if (!layout_use %in% c("tripartite", "star", "hub", "bipartite", "fr", "kk", "circle", "nicely", "lgl", "drl")) {
    layout_use <- if (length(tf_candidates) <= 1L) "star" else "kk"
  }
  triplets <- scenic_get_triplets(srt, tool_name, required = FALSE)
  if (is.null(triplets)) {
    log_message(
      "{.val egrn} plots need TF-region-gene triplets; drawing a TF-gene network instead",
      message_type = "warning"
    )
  }
  plot_title <- title
  if (is.null(plot_title) && length(tf_candidates) == 1L) {
    grp <- scenic_network_tf_groups(rank_table %||% top_table, tf_candidates)[[1]]
    plot_title <- if (!is.na(grp) && nzchar(grp)) {
      paste0(tf_candidates, " (", grp, ")")
    } else {
      tf_candidates
    }
  }
  scenic_plot_one_network(
    srt = srt,
    tool_name = tool_name,
    tfs = tf_candidates,
    max_targets = max_targets,
    max_edges = Inf,
    network_layout = layout_use,
    label_nodes = label_nodes,
    highlight_tf = highlight_tf,
    include_regions = TRUE,
    filter_tfs = TRUE,
    palette = palette,
    palcolor = palcolor,
    title = plot_title %||% "SCENIC+ enhancer-driven GRN",
    curvature = if (layout_use %in% c("star")) 0 else 0.08,
    prefer_triplets = TRUE,
    rank_table = rank_table %||% top_table
  )
}

scenic_plot_overlap <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  top_n = 12,
  heatmap_palette = NULL,
  heatmap_palcolor = NULL,
  heatmap_group_palette = "Chinese",
  heatmap_group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list(),
  title = NULL
) {
  regulon_list <- scenic_get_regulon_list(srt, tool_name)
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = names(regulon_list),
    top_table = top_table,
    max_features = if (is.null(features)) top_n else NULL
  )
  mat <- scenic_regulon_jaccard(regulon_list[regulons])
  heatmap_result <- scenic_plot_feature_heatmap_from_matrix(
    mat = mat,
    group_names = colnames(mat),
    features = rownames(mat),
    legend_title = "Jaccard",
    title = title %||% "SCENIC eRegulon overlap",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = "right",
    column_names_side = "top",
    row_names_rot = 0,
    column_names_rot = 45,
    border = TRUE,
    heatmap_palette = heatmap_palette %||% "YlOrRd",
    heatmap_palcolor = heatmap_palcolor,
    group_palette = heatmap_group_palette,
    group_palcolor = heatmap_group_palcolor,
    heatmap_limits = heatmap_limits %||% c(0, 1),
    heatmap_args = heatmap_args
  )
  plot <- heatmap_result[["plot"]]
  plot_data <- scenic_matrix_to_long(mat, "regulon_1", "regulon_2", "jaccard")
  list(plot = plot, plots = list(plot), data = plot_data, heatmap = heatmap_result)
}

scenic_plot_one_network <- function(
  srt,
  tool_name,
  tfs,
  max_targets,
  max_edges,
  network_layout,
  label_nodes,
  highlight_tf,
  include_regions,
  filter_tfs,
  palette,
  palcolor,
  title,
  curvature,
  prefer_triplets = FALSE,
  rank_table = NULL
) {
  edge_data <- scenic_network_edges(
    srt = srt,
    tool_name = tool_name,
    tfs = tfs,
    max_targets = max_targets,
    max_edges = max_edges,
    include_regions = include_regions,
    filter_tfs = filter_tfs,
    prefer_triplets = prefer_triplets
  )
  network_data <- scenic_network_plot_data(edge_data = edge_data, layout = network_layout)
  default_labels <- if (
    network_layout %in% c("star", "tripartite", "bipartite") ||
      isTRUE(prefer_triplets)
  ) {
    "all"
  } else {
    "tfs"
  }
  label_data <- scenic_network_label_data(
    node_data = network_data[["nodes"]],
    label_nodes = label_nodes,
    highlight_tf = highlight_tf,
    top_n = Inf,
    default = default_labels
  )
  plot <- scenic_network_ggplot(
    node_data = network_data[["nodes"]],
    edge_plot = network_data[["edge_plot"]],
    edge_data = network_data[["edges"]],
    label_data = label_data,
    label_nodes = label_nodes,
    layout = network_layout,
    palette = palette,
    palcolor = palcolor,
    title = title,
    curvature = curvature,
    rank_table = rank_table
  )
  list(
    plot = plot,
    plots = list(plot),
    data = list(edges = network_data[["edges"]], nodes = network_data[["nodes"]])
  )
}

scenic_network_tf_candidates <- function(...) {
  values <- unique(unlist(lapply(list(...), function(x) {
    if (is.null(x)) {
      return(character())
    }
    scenic_tf_from_regulon(as.character(x))
  })))
  values[nzchar(values)]
}

scenic_get_triplets <- function(srt, tool_name, required = FALSE) {
  result <- srt@tools[[tool_name]]
  triplets <- result[["triplets"]]
  if (is.null(triplets) || nrow(triplets) == 0L) {
    if (isTRUE(required)) {
      log_message(
        "Cannot find TF–region–gene triplets in tools slot {.val {tool_name}}",
        message_type = "error"
      )
    }
    return(NULL)
  }
  triplets <- as.data.frame(triplets, stringsAsFactors = FALSE)
  if (!all(c("TF", "region", "gene") %in% colnames(triplets))) {
    if (isTRUE(required)) {
      log_message(
        "SCENIC+ triplets must contain {.val TF}, {.val region}, and {.val gene} columns",
        message_type = "error"
      )
    }
    return(NULL)
  }
  triplets
}

scenic_network_edges <- function(
  srt,
  tool_name,
  tfs = NULL,
  max_targets = 30,
  max_edges = Inf,
  include_regions = TRUE,
  filter_tfs = FALSE,
  prefer_triplets = FALSE
) {
  triplets <- if (isTRUE(include_regions)) {
    scenic_get_triplets(srt, tool_name, required = FALSE)
  } else {
    NULL
  }
  use_triplets <- !is.null(triplets)
  if (isTRUE(use_triplets) && !is.null(tfs)) {
    use_triplets <- any(triplets[["TF"]] %in% tfs)
  }
  if (isTRUE(prefer_triplets) && is.null(triplets)) {
    use_triplets <- FALSE
  }
  if (isTRUE(use_triplets)) {
    if (isTRUE(filter_tfs) && !is.null(tfs)) {
      triplets <- triplets[triplets[["TF"]] %in% tfs, , drop = FALSE]
    }
    if (nrow(triplets) == 0L) {
      use_triplets <- FALSE
    } else {
      return(scenic_triplet_edges(triplets, max_targets = max_targets))
    }
  }

  adjacency <- scenic_get_adjacency(srt, tool_name)
  cols <- scenic_adjacency_columns(adjacency)
  if (isTRUE(filter_tfs)) {
    if (is.null(tfs) || length(tfs) == 0L) {
      log_message(
        "No SCENIC adjacency edges found for requested TFs",
        message_type = "error"
      )
    }
    adjacency <- adjacency[adjacency[[cols[["tf"]]]] %in% tfs, , drop = FALSE]
  } else if (!is.null(tfs) && length(tfs) > 0L) {
    adjacency <- adjacency[
      adjacency[[cols[["tf"]]]] %in% tfs |
        adjacency[[cols[["target"]]]] %in% tfs,
      ,
      drop = FALSE
    ]
  }
  if (nrow(adjacency) == 0L) {
    log_message(
      "No SCENIC adjacency edges found for the requested network",
      message_type = "error"
    )
  }
  adjacency <- scenic_top_edges(
    adjacency,
    tf_col = cols[["tf"]],
    weight_col = cols[["weight"]],
    max_targets = max_targets
  )
  adjacency <- scenic_limit_edges(
    adjacency,
    weight_col = cols[["weight"]],
    max_edges = max_edges
  )
  scenic_adjacency_edges(adjacency, cols)
}

scenic_adjacency_edges <- function(adjacency, cols) {
  weight <- if (!is.null(cols[["weight"]])) {
    as.numeric(adjacency[[cols[["weight"]]]])
  } else {
    rep(1, nrow(adjacency))
  }
  weight[!is.finite(weight)] <- 1
  edge_data <- data.frame(
    from = as.character(adjacency[[cols[["tf"]]]]),
    to = as.character(adjacency[[cols[["target"]]]]),
    weight = weight,
    edge_type = "tf_gene",
    edge_sign = ifelse(weight < 0, "negative", "positive"),
    stringsAsFactors = FALSE
  )
  scenic_clean_edges(edge_data)
}

scenic_triplet_edges <- function(triplets, max_targets = 20) {
  score_col <- scenic_first(
    intersect(c("score", "importance", "weight"), colnames(triplets)),
    NULL
  )
  if (!is.null(score_col)) {
    triplets <- triplets[order(abs(triplets[[score_col]]), decreasing = TRUE), , drop = FALSE]
  }
  kept <- do.call(
    rbind,
    lapply(split(triplets, triplets[["TF"]]), function(df) {
      genes <- unique(as.character(df[["gene"]]))
      genes <- utils::head(genes, max_targets)
      df <- df[df[["gene"]] %in% genes, , drop = FALSE]
      regions <- unique(as.character(df[["region"]]))
      if (length(regions) > max_targets) {
        df <- df[df[["region"]] %in% utils::head(regions, max_targets), , drop = FALSE]
      }
      df
    })
  )
  rownames(kept) <- NULL
  if (is.null(kept) || nrow(kept) == 0L) {
    log_message(
      "No SCENIC+ triplets remain after target filtering",
      message_type = "error"
    )
  }
  weight <- if (!is.null(score_col)) as.numeric(kept[[score_col]]) else rep(1, nrow(kept))
  weight[!is.finite(weight)] <- 1
  sign_col <- scenic_first(
    intersect(c("r2g_sign", "tf_sign"), colnames(kept)),
    NULL
  )
  edge_sign <- if (!is.null(sign_col)) {
    ifelse(kept[[sign_col]] %in% c("-", "negative", "-1", -1), "negative", "positive")
  } else {
    ifelse(weight < 0, "negative", "positive")
  }
  tf_region <- data.frame(
    from = as.character(kept[["TF"]]),
    to = as.character(kept[["region"]]),
    weight = abs(weight),
    edge_type = "tf_region",
    edge_sign = edge_sign,
    stringsAsFactors = FALSE
  )
  region_gene <- data.frame(
    from = as.character(kept[["region"]]),
    to = as.character(kept[["gene"]]),
    weight = abs(weight),
    edge_type = "region_gene",
    edge_sign = edge_sign,
    stringsAsFactors = FALSE
  )
  scenic_aggregate_edges(rbind(tf_region, region_gene))
}

scenic_aggregate_edges <- function(edge_data) {
  edge_data <- scenic_clean_edges(edge_data)
  if (nrow(edge_data) == 0L) {
    return(edge_data)
  }
  key <- paste(
    edge_data[["from"]],
    edge_data[["to"]],
    edge_data[["edge_type"]],
    sep = "\r"
  )
  split_edges <- split(edge_data, key)
  out <- do.call(
    rbind,
    lapply(split_edges, function(df) {
      idx <- which.max(abs(df[["weight"]]))
      df[idx, , drop = FALSE]
    })
  )
  rownames(out) <- NULL
  out
}

scenic_clean_edges <- function(edge_data) {
  edge_data <- edge_data[edge_data[["from"]] != edge_data[["to"]], , drop = FALSE]
  edge_data <- edge_data[
    stats::complete.cases(edge_data[, c("from", "to"), drop = FALSE]),
    ,
    drop = FALSE
  ]
  if (nrow(edge_data) == 0L) {
    log_message(
      "No valid TF-target edges remain for SCENIC network plotting",
      message_type = "error"
    )
  }
  edge_data
}

scenic_network_plot_data <- function(edge_data, layout = "star") {
  graph <- igraph::graph_from_data_frame(edge_data, directed = TRUE)
  node_names <- igraph::V(graph)$name
  node_type <- scenic_network_node_types(node_names, edge_data)
  xy <- scenic_network_layout(graph, layout = layout, node_type = node_type)
  node_data <- data.frame(
    name = node_names,
    x = xy[node_names, 1],
    y = xy[node_names, 2],
    node_type = unname(node_type[node_names]),
    degree = as.numeric(igraph::degree(graph, mode = "all")),
    out_degree = as.numeric(igraph::degree(graph, mode = "out")),
    in_degree = as.numeric(igraph::degree(graph, mode = "in")),
    stringsAsFactors = FALSE
  )
  node_data[["label"]] <- scenic_network_node_labels(
    node_data[["name"]],
    node_data[["node_type"]]
  )
  node_data[["node_type"]] <- factor(
    node_data[["node_type"]],
    levels = intersect(c("TF", "region", "gene"), unique(node_data[["node_type"]]))
  )
  edge_plot <- merge(
    edge_data,
    node_data[, c("name", "x", "y"), drop = FALSE],
    by.x = "from",
    by.y = "name"
  )
  edge_plot <- merge(
    edge_plot,
    node_data[, c("name", "x", "y"), drop = FALSE],
    by.x = "to",
    by.y = "name",
    suffixes = c("", "_end")
  )
  list(edges = edge_data, nodes = node_data, edge_plot = edge_plot, graph = graph)
}

scenic_network_node_types <- function(node_names, edge_data) {
  tfs <- unique(edge_data[["from"]][edge_data[["edge_type"]] %in% c("tf_gene", "tf_region")])
  regions <- unique(c(
    edge_data[["to"]][edge_data[["edge_type"]] == "tf_region"],
    edge_data[["from"]][edge_data[["edge_type"]] == "region_gene"]
  ))
  node_type <- stats::setNames(rep("gene", length(node_names)), node_names)
  node_type[node_names %in% regions] <- "region"
  node_type[node_names %in% tfs] <- "TF"
  node_type
}

scenic_network_node_labels <- function(names, node_type) {
  vapply(
    seq_along(names),
    function(idx) {
      if (!identical(node_type[[idx]], "region")) {
        return(names[[idx]])
      }
      parsed <- scenic_parse_genomic_region(names[[idx]])
      if (nrow(parsed) == 0L || is.na(parsed[["start"]][[1]])) {
        return(names[[idx]])
      }
      pos <- parsed[["start"]][[1]]
      pos_lab <- if (pos >= 1e6) {
        sprintf("%.2fMb", pos / 1e6)
      } else if (pos >= 1e3) {
        sprintf("%.0fkb", pos / 1e3)
      } else {
        as.character(as.integer(pos))
      }
      paste0(parsed[["seqnames"]][[1]], ":", pos_lab)
    },
    character(1)
  )
}

scenic_network_resolve_layout <- function(layout, default = "star") {
  if (is.null(layout) || identical(layout, "auto")) {
    return(default)
  }
  layout
}

scenic_network_layout <- function(graph, layout = "fr", node_type = NULL) {
  names <- igraph::V(graph)$name
  if (is.null(node_type)) {
    node_type <- stats::setNames(
      ifelse(igraph::degree(graph, mode = "out") > 0, "TF", "gene"),
      names
    )
  }
  node_type <- node_type[names]
  xy <- switch(
    layout,
    star = scenic_layout_star(graph, node_type),
    hub = scenic_layout_hubs(graph, node_type),
    tripartite = scenic_layout_tripartite(graph, node_type),
    bipartite = scenic_layout_bipartite(graph, node_type),
    circle = scenic_layout_named(igraph::layout_in_circle(graph), names),
    fr = scenic_layout_force(graph, method = "fr"),
    nicely = scenic_layout_named(igraph::layout_nicely(graph), names),
    kk = scenic_layout_force(graph, method = "kk"),
    lgl = scenic_layout_named(igraph::layout_with_lgl(graph), names),
    drl = scenic_layout_named(igraph::layout_with_drl(graph), names),
    scenic_layout_force(graph, method = "fr")
  )
  if (is.null(dim(xy))) {
    xy <- matrix(xy, ncol = 2)
  }
  if (ncol(xy) < 2) {
    xy <- cbind(xy[, 1], 0)
  }
  if (is.null(rownames(xy))) {
    rownames(xy) <- names
  }
  xy[names, , drop = FALSE]
}

scenic_layout_named <- function(xy, names) {
  if (is.null(rownames(xy))) {
    rownames(xy) <- names
  }
  xy
}

scenic_layout_force <- function(graph, method = "fr") {
  names <- igraph::V(graph)$name
  layout_weights <- NULL
  if ("weight" %in% igraph::edge_attr_names(graph)) {
    layout_weights <- abs(igraph::E(graph)$weight)
    layout_weights[!is.finite(layout_weights) | layout_weights <= 0] <- 1
  }
  start <- igraph::layout_in_circle(graph)
  xy <- if (identical(method, "kk")) {
    igraph::layout_with_kk(graph, coords = start, weights = layout_weights)
  } else {
    igraph::layout_with_fr(
      graph,
      coords = start,
      niter = 500,
      weights = layout_weights
    )
  }
  scenic_layout_named(xy, names)
}

scenic_circle_coords <- function(n, radius = 1, start = pi / 2) {
  if (n <= 0) {
    return(matrix(numeric(), ncol = 2))
  }
  if (n == 1L) {
    return(matrix(c(0, radius), ncol = 2))
  }
  angles <- start - 2 * pi * (seq_len(n) - 1) / n
  cbind(radius * cos(angles), radius * sin(angles))
}

scenic_layout_star <- function(graph, node_type) {
  names <- igraph::V(graph)$name
  tfs <- names[node_type[names] == "TF"]
  if (length(tfs) != 1L) {
    return(scenic_layout_hubs(graph, node_type))
  }
  xy <- matrix(0, nrow = length(names), ncol = 2, dimnames = list(names, c("x", "y")))
  xy[tfs, ] <- c(0, 0)
  regions <- names[node_type[names] == "region"]
  genes <- names[node_type[names] == "gene"]
  edges <- igraph::as_data_frame(graph, what = "edges")
  if (length(regions) > 0L) {
    region_genes <- vapply(regions, function(region) {
      hits <- as.character(edges[["to"]][edges[["from"]] == region])
      scenic_first(hits[hits %in% genes], region)
    }, character(1))
    regions <- regions[order(region_genes, regions)]
    xy[regions, ] <- scenic_circle_coords(length(regions), radius = 1)
    if (length(genes) > 0L) {
      for (gene in genes) {
        regs <- unique(as.character(edges[["from"]][edges[["to"]] == gene & edges[["from"]] %in% regions]))
        if (length(regs) == 0L) {
          next
        }
        pos <- colMeans(xy[regs, , drop = FALSE])
        r <- sqrt(sum(pos^2))
        if (!is.finite(r) || r < 1e-8) {
          xy[gene, ] <- c(0, 1.85)
        } else {
          xy[gene, ] <- pos / r * 1.85
        }
      }
      placed <- genes[rowSums(abs(xy[genes, , drop = FALSE])) > 0]
      missing <- setdiff(genes, placed)
      if (length(missing) > 0L) {
        xy[missing, ] <- scenic_circle_coords(length(missing), radius = 1.85)
      }
      dup <- duplicated(round(xy[genes, , drop = FALSE], 6))
      if (any(dup)) {
        n <- length(genes)
        xy[genes, ] <- scenic_circle_coords(n, radius = 1.85)
        if (length(regions) == n) {
          xy[regions, ] <- scenic_circle_coords(n, radius = 1)
        }
      }
    }
  } else if (length(genes) > 0L) {
    xy[genes, ] <- scenic_circle_coords(length(genes), radius = 1.35)
  }
  xy
}

scenic_layout_hubs <- function(graph, node_type) {
  names <- igraph::V(graph)$name
  tfs <- names[node_type[names] == "TF"]
  xy <- matrix(0, nrow = length(names), ncol = 2, dimnames = list(names, c("x", "y")))
  n_tf <- length(tfs)
  if (n_tf == 0L) {
    return(scenic_layout_named(igraph::layout_in_circle(graph), names))
  }
  tf_radius <- if (n_tf == 1L) 0 else max(2.4, 0.58 * n_tf)
  xy[tfs, ] <- scenic_circle_coords(n_tf, radius = tf_radius)
  parents <- scenic_network_parent_tfs(graph, node_type)
  local_r <- list(region = 0.55, gene = 1.05)
  for (tf in tfs) {
    for (type_use in c("region", "gene")) {
      members <- names[node_type[names] == type_use]
      unique_members <- members[vapply(members, function(node) {
        p <- parents[[node]]
        length(p) == 1L && identical(p[[1]], tf)
      }, logical(1))]
      if (length(unique_members) == 0L) {
        next
      }
      local <- scenic_circle_coords(length(unique_members), radius = local_r[[type_use]])
      xy[unique_members, 1] <- xy[tf, 1] + local[, 1]
      xy[unique_members, 2] <- xy[tf, 2] + local[, 2]
    }
  }
  shared <- setdiff(names, tfs)
  shared <- shared[vapply(shared, function(node) length(parents[[node]]) > 1L, logical(1))]
  for (node in shared) {
    tf_xy <- xy[parents[[node]], , drop = FALSE]
    center <- colMeans(tf_xy)
    if (!any(is.finite(center)) || identical(as.numeric(center), c(0, 0))) {
      center <- c(0.35, 0)
    } else {
      center <- 0.55 * center
    }
    xy[node, ] <- center
  }
  at_origin <- abs(xy[, 1]) < 1e-12 & abs(xy[, 2]) < 1e-12
  orphans <- setdiff(names[at_origin], tfs)
  if (length(orphans) > 0L) {
    xy[orphans, ] <- scenic_circle_coords(length(orphans), radius = tf_radius + 1.6)
  }
  xy
}

scenic_layout_tripartite <- function(graph, node_type) {
  names <- igraph::V(graph)$name
  tfs <- names[node_type[names] == "TF"]
  regions <- names[node_type[names] == "region"]
  genes <- names[node_type[names] == "gene"]
  if (length(regions) == 0L) {
    return(scenic_layout_bipartite(graph, node_type))
  }
  parents <- scenic_network_parent_tfs(graph, node_type)
  tfs <- scenic_stable_sort(tfs)
  regions <- scenic_order_by_barycenter(regions, parents, tfs)
  genes <- scenic_order_by_barycenter(genes, parents, tfs)
  xy <- matrix(0, nrow = length(names), ncol = 2, dimnames = list(names, c("x", "y")))
  scenic_place_column <- function(nodes, x) {
    n <- length(nodes)
    if (n == 0L) {
      return(invisible())
    }
    y <- if (n == 1L) 0 else seq(1, -1, length.out = n)
    xy[nodes, 1] <<- x
    xy[nodes, 2] <<- y
  }
  scenic_place_column(tfs, 0)
  scenic_place_column(regions, 1.55)
  scenic_place_column(genes, 3.1)
  xy
}

scenic_layout_bipartite <- function(graph, node_type) {
  names <- igraph::V(graph)$name
  tfs <- scenic_stable_sort(names[node_type[names] == "TF"])
  others <- names[node_type[names] != "TF"]
  parents <- scenic_network_parent_tfs(graph, node_type)
  others <- scenic_order_by_barycenter(others, parents, tfs)
  xy <- matrix(0, nrow = length(names), ncol = 2, dimnames = list(names, c("x", "y")))
  if (length(tfs) > 0L) {
    xy[tfs, 1] <- 0
    xy[tfs, 2] <- if (length(tfs) == 1L) 0 else seq(1, -1, length.out = length(tfs))
  }
  if (length(others) > 0L) {
    xy[others, 1] <- 1
    xy[others, 2] <- if (length(others) == 1L) 0 else seq(1, -1, length.out = length(others))
  }
  xy
}

scenic_network_parent_tfs <- function(graph, node_type) {
  names <- igraph::V(graph)$name
  tfs <- names[node_type[names] == "TF"]
  edges <- igraph::as_data_frame(graph, what = "edges")
  parents <- stats::setNames(vector("list", length(names)), names)
  for (node in names) {
    if (node %in% tfs) {
      parents[[node]] <- node
      next
    }
    frontier <- node
    found <- character()
    seen <- node
    while (length(frontier) > 0L) {
      incoming <- unique(edges[["from"]][edges[["to"]] %in% frontier])
      incoming <- setdiff(incoming, seen)
      found <- unique(c(found, intersect(incoming, tfs)))
      frontier <- setdiff(incoming, tfs)
      seen <- unique(c(seen, incoming))
    }
    parents[[node]] <- found
  }
  parents
}

scenic_order_by_barycenter <- function(nodes, parents, parent_order) {
  if (length(nodes) <= 1L) {
    return(nodes)
  }
  parent_idx <- stats::setNames(seq_along(parent_order), parent_order)
  scores <- vapply(
    nodes,
    function(node) {
      p <- intersect(parents[[node]], parent_order)
      if (length(p) == 0L) {
        return(Inf)
      }
      mean(parent_idx[p])
    },
    numeric(1)
  )
  nodes[order(scores, nodes)]
}

scenic_stable_sort <- function(x) {
  sort(unique(as.character(x)))
}

scenic_network_label_data <- function(
  node_data,
  label_nodes = c("auto", "tfs", "all", "none"),
  highlight_tf = NULL,
  top_n = 60,
  default = "tfs"
) {
  if (identical(label_nodes, "auto")) {
    label_nodes <- default
  }
  label_nodes <- match.arg(label_nodes, c("auto", "tfs", "all", "none"))
  if (identical(label_nodes, "auto")) {
    label_nodes <- default
  }
  label_data <- switch(
    label_nodes,
    all = node_data,
    tfs = {
      tf_nodes <- node_data[as.character(node_data[["node_type"]]) == "TF", , drop = FALSE]
      tf_nodes <- tf_nodes[order(tf_nodes[["degree"]], decreasing = TRUE), , drop = FALSE]
      utils::head(tf_nodes, top_n)
    },
    none = node_data[FALSE, , drop = FALSE],
    node_data
  )
  if (!is.null(highlight_tf)) {
    highlight_tf <- scenic_tf_from_regulon(highlight_tf)
    highlight_data <- node_data[node_data[["name"]] %in% highlight_tf, , drop = FALSE]
    label_data <- unique(rbind(label_data, highlight_data))
  }
  label_data
}

scenic_network_palette <- function(palette) {
  if (identical(palette, "RdYlBu")) {
    return("Chinese")
  }
  palette
}

scenic_network_label_contrast <- function(color) {
  rgb <- grDevices::col2rgb(color)[, 1]
  if (sum(rgb) > 255 * 2) {
    "black"
  } else {
    "white"
  }
}

scenic_network_edge_width <- function(weight, width_range = c(0.32, 1.05)) {
  weight <- abs(as.numeric(weight))
  weight[!is.finite(weight)] <- 1
  if (length(weight) == 0L) {
    return(numeric())
  }
  if (length(unique(weight)) == 1L) {
    return(rep(mean(width_range), length(weight)))
  }
  scales::rescale(weight, to = width_range)
}

scenic_network_region_to_tf <- function(edge_data) {
  tf_region <- edge_data[
    edge_data[["edge_type"]] == "tf_region",
    c("from", "to"),
    drop = FALSE
  ]
  if (nrow(tf_region) == 0L) {
    return(stats::setNames(character(), character()))
  }
  stats::setNames(tf_region[["from"]], tf_region[["to"]])
}

scenic_network_edge_source_tf <- function(edge_data) {
  region_to_tf <- scenic_network_region_to_tf(edge_data)
  vapply(
    seq_len(nrow(edge_data)),
    function(idx) {
      edge_type <- edge_data[["edge_type"]][[idx]]
      if (edge_type %in% c("tf_gene", "tf_region")) {
        return(as.character(edge_data[["from"]][[idx]]))
      }
      if (identical(edge_type, "region_gene")) {
        region <- as.character(edge_data[["from"]][[idx]])
        tf <- region_to_tf[[region]]
        if (!is.null(tf) && nzchar(tf)) {
          return(tf)
        }
      }
      as.character(edge_data[["from"]][[idx]])
    },
    character(1)
  )
}

scenic_network_gene_sources <- function(gene, edge_data, region_to_tf) {
  direct <- unique(as.character(edge_data[["from"]][
    edge_data[["to"]] == gene & edge_data[["edge_type"]] == "tf_gene"
  ]))
  if (length(direct) > 0L) {
    return(direct)
  }
  regions <- unique(as.character(edge_data[["from"]][
    edge_data[["to"]] == gene & edge_data[["edge_type"]] == "region_gene"
  ]))
  unique(na.omit(unname(region_to_tf[regions])))
}

scenic_network_style_data <- function(
  node_data,
  edge_plot,
  edge_data,
  palette = "Chinese",
  palcolor = NULL,
  network_blendmode = "average"
) {
  palette <- scenic_network_palette(palette)
  tfs <- sort(unique(as.character(node_data[["name"]][
    as.character(node_data[["node_type"]]) == "TF"
  ])))
  if (length(tfs) == 0L) {
    tfs <- sort(unique(scenic_network_edge_source_tf(edge_data)))
  }
  tf_cols <- palette_colors(tfs, palette = palette, palcolor = palcolor)
  region_to_tf <- scenic_network_region_to_tf(edge_data)

  fill_color <- character(nrow(node_data))
  border_color <- character(nrow(node_data))
  node_size <- numeric(nrow(node_data))
  for (idx in seq_len(nrow(node_data))) {
    node_type <- as.character(node_data[["node_type"]][[idx]])
    node_name <- as.character(node_data[["name"]][[idx]])
    if (identical(node_type, "TF")) {
      fill_color[[idx]] <- unname(tf_cols[[node_name]])
      border_color[[idx]] <- "#1A1A1A"
      node_size[[idx]] <- 8.2
    } else if (identical(node_type, "region")) {
      fill_color[[idx]] <- "#2E2E2E"
      border_color[[idx]] <- "#2E2E2E"
      node_size[[idx]] <- 3.2
    } else {
      fill_color[[idx]] <- "#E0E0E0"
      border_color[[idx]] <- "#3F3F3F"
      node_size[[idx]] <- 3.8
    }
  }
  node_data[["node_color"]] <- fill_color
  node_data[["border_color"]] <- border_color
  node_data[["node_size"]] <- node_size
  node_data[["label_color"]] <- vapply(fill_color, scenic_network_label_contrast, character(1))

  edge_plot[["source_tf"]] <- if ("edge_type" %in% colnames(edge_plot)) {
    scenic_network_edge_source_tf(edge_plot)
  } else {
    as.character(edge_plot[["from"]])
  }
  edge_plot[["edge_color"]] <- unname(tf_cols[edge_plot[["source_tf"]]])
  edge_plot[["edge_color"]][is.na(edge_plot[["edge_color"]]) | !nzchar(edge_plot[["edge_color"]])] <- "#B3B3B3"
  edge_plot[["linewidth_scaled"]] <- scenic_network_edge_width(edge_plot[["weight"]])
  if ("edge_type" %in% colnames(edge_plot)) {
    r2g <- edge_plot[["edge_type"]] == "region_gene"
    edge_plot[["linewidth_scaled"]][r2g] <- edge_plot[["linewidth_scaled"]][r2g] * 0.82
  }

  list(
    nodes = node_data,
    edges = edge_plot,
    tf_cols = tf_cols,
    blendmode = network_blendmode,
    region_to_tf = region_to_tf
  )
}

scenic_is_region_coordinate <- function(name) {
  grepl("^(chr|CHR)[^[:space:]]*[:_-][0-9]", as.character(name))
}

scenic_network_tf_groups <- function(rank_table, tfs) {
  tfs <- unique(as.character(tfs))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  out <- stats::setNames(rep(NA_character_, length(tfs)), tfs)
  if (length(tfs) == 0L || is.null(rank_table) || nrow(rank_table) == 0L) {
    return(out)
  }
  df <- rank_table
  tf_id <- if ("TF" %in% colnames(df)) {
    as.character(df[["TF"]])
  } else if ("regulon" %in% colnames(df)) {
    scenic_tf_from_regulon(df[["regulon"]])
  } else {
    return(out)
  }
  df <- df[tf_id %in% tfs, , drop = FALSE]
  tf_id <- tf_id[tf_id %in% tfs]
  if (nrow(df) == 0L || !"group" %in% colnames(df)) {
    return(out)
  }
  df[["tf_id"]] <- tf_id
  ord <- seq_len(nrow(df))
  if ("specificity_score" %in% colnames(df)) {
    score <- as.numeric(df[["specificity_score"]])
    score[!is.finite(score)] <- -Inf
    ord <- order(df[["tf_id"]], -score, df[["tf_id"]])
  } else if ("rank" %in% colnames(df)) {
    rank <- as.numeric(df[["rank"]])
    rank[!is.finite(rank)] <- Inf
    ord <- order(df[["tf_id"]], rank)
  }
  df <- df[ord, , drop = FALSE]
  df <- df[!duplicated(df[["tf_id"]]), , drop = FALSE]
  out[df[["tf_id"]]] <- as.character(df[["group"]])
  out
}

scenic_network_legend_labels <- function(tf_cols, tf_groups) {
  tfs <- names(tf_cols)
  groups <- unname(tf_groups[tfs])
  ifelse(
    is.na(groups) | !nzchar(groups),
    tfs,
    paste0(tfs, " (", groups, ")")
  )
}

scenic_network_ggplot <- function(
  node_data,
  edge_plot,
  edge_data,
  label_data,
  label_nodes = "auto",
  title,
  layout = "star",
  palette = "RdYlBu",
  palcolor = NULL,
  curvature = 0,
  edge_width_range = c(0.32, 1.05),
  network_blendmode = "average",
  rank_table = NULL
) {
  styled <- scenic_network_style_data(
    node_data = node_data,
    edge_plot = edge_plot,
    edge_data = edge_data,
    palette = palette,
    palcolor = palcolor,
    network_blendmode = network_blendmode
  )
  node_data <- styled[["nodes"]]
  edge_plot <- styled[["edges"]]
  tf_cols <- styled[["tf_cols"]]
  if (length(edge_width_range) == 2L) {
    edge_plot[["linewidth_scaled"]] <- scenic_network_edge_width(
      edge_plot[["weight"]],
      width_range = edge_width_range
    )
  }

  tf_nodes <- node_data[as.character(node_data[["node_type"]]) == "TF", , drop = FALSE]
  region_nodes <- node_data[as.character(node_data[["node_type"]]) == "region", , drop = FALSE]
  gene_nodes <- node_data[as.character(node_data[["node_type"]]) == "gene", , drop = FALSE]
  use_radial <- identical(layout, "star") && nrow(tf_nodes) == 1L
  show_gene_text <- !identical(label_nodes, "none") && !identical(label_nodes, "tfs")
  if (identical(label_nodes, "tfs")) {
    show_gene_text <- FALSE
  }
  if (identical(label_nodes, "auto") && layout %in% c("kk", "fr", "nicely", "lgl", "drl") && nrow(gene_nodes) > 48L) {
    show_gene_text <- FALSE
  }

  p <- ggplot2::ggplot()
  if (abs(curvature) < 1e-8) {
    p <- p + ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["x_end"]],
        yend = .data[["y_end"]],
        color = .data[["edge_color"]],
        linewidth = .data[["linewidth_scaled"]]
      ),
      alpha = 0.62,
      lineend = "round",
      show.legend = FALSE
    )
  } else {
    p <- p + ggplot2::geom_curve(
      data = edge_plot,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["x_end"]],
        yend = .data[["y_end"]],
        color = .data[["edge_color"]],
        linewidth = .data[["linewidth_scaled"]]
      ),
      curvature = curvature,
      alpha = 0.55,
      lineend = "round",
      show.legend = FALSE
    )
  }

  if (nrow(region_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = region_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_color"]],
        color = .data[["border_color"]],
        size = .data[["node_size"]]
      ),
      shape = 23,
      stroke = 0.35,
      show.legend = FALSE
    )
  }
  if (nrow(gene_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = gene_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_color"]],
        color = .data[["border_color"]],
        size = .data[["node_size"]]
      ),
      shape = 21,
      stroke = 0.45,
      show.legend = FALSE
    )
  }
  if (nrow(tf_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = tf_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_color"]],
        size = .data[["node_size"]]
      ),
      shape = 21,
      color = "#1A1A1A",
      stroke = 0.9,
      show.legend = FALSE
    )
  }

  if (nrow(tf_nodes) > 0L && !identical(label_nodes, "none")) {
    p <- p + ggplot2::geom_text(
      data = tf_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        label = .data[["label"]],
        color = .data[["label_color"]]
      ),
      fontface = "bold",
      size = 2.6,
      show.legend = FALSE
    )
  }

  gene_labels <- gene_nodes
  if (!isTRUE(show_gene_text)) {
    gene_labels <- gene_nodes[FALSE, , drop = FALSE]
  }
  region_labels <- region_nodes
  if (identical(label_nodes, "all")) {
    region_labels <- region_nodes[!scenic_is_region_coordinate(region_nodes[["name"]]), , drop = FALSE]
  } else {
    region_labels <- region_nodes[FALSE, , drop = FALSE]
  }
  other_labels <- rbind(gene_labels, region_labels)
  if (nrow(other_labels) > 0L) {
    if (isTRUE(use_radial)) {
      origin <- c(tf_nodes[["x"]][[1]], tf_nodes[["y"]][[1]])
      other_labels <- scenic_radial_label_coords(other_labels, origin = origin, expand = 1.22)
      p <- p + ggplot2::geom_text(
        data = other_labels,
        ggplot2::aes(
          x = .data[["label_x"]],
          y = .data[["label_y"]],
          label = .data[["label"]],
          hjust = .data[["hjust"]]
        ),
        size = 2.85,
        color = "grey15",
        show.legend = FALSE
      )
    } else {
      p <- p + ggrepel::geom_text_repel(
        data = other_labels,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
        size = 2.7,
        color = "grey20",
        max.overlaps = Inf,
        box.padding = 0.18,
        point.padding = 0.22,
        segment.color = "grey75",
        segment.size = 0.15,
        min.segment.length = 0.08,
        show.legend = FALSE
      )
    }
  }

  if (length(tf_cols) > 1L) {
    tf_groups <- scenic_network_tf_groups(rank_table, names(tf_cols))
    legend_labels <- scenic_network_legend_labels(tf_cols, tf_groups)
    legend_df <- data.frame(
      x = mean(node_data[["x"]]),
      y = mean(node_data[["y"]]),
      label = legend_labels,
      stringsAsFactors = FALSE
    )
    legend_name <- if (any(!is.na(unname(tf_groups)) & nzchar(unname(tf_groups)))) {
      "TF (top cell type)"
    } else {
      "TF"
    }
    p <- p +
      ggplot2::scale_color_identity(guide = "none") +
      ggplot2::scale_fill_identity(guide = "none") +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(
        data = legend_df,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["label"]]),
        shape = 21,
        size = 4,
        alpha = 0,
        color = "#1A1A1A",
        stroke = 0.5,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        name = legend_name,
        values = stats::setNames(unname(tf_cols), legend_labels),
        guide = ggplot2::guide_legend(
          override.aes = list(alpha = 1, size = 4, shape = 21, color = "#1A1A1A")
        )
      )
  } else {
    p <- p + ggplot2::scale_color_identity(guide = "none") + ggplot2::scale_fill_identity(guide = "none")
  }

  extra <- NULL
  if (isTRUE(use_radial) && nrow(other_labels) > 0L && "label_x" %in% colnames(other_labels)) {
    extra <- other_labels
  }
  limits <- scenic_network_limits(node_data, extra %||% other_labels, use_radial = isTRUE(use_radial))
  p +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_identity(guide = "none") +
    ggplot2::coord_equal(xlim = limits[["x"]], ylim = limits[["y"]], clip = "off") +
    theme_scop() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      plot.margin = ggplot2::margin(8, 16, 8, 16)
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}

scenic_radial_label_coords <- function(label_data, origin = c(0, 0), expand = 1.18) {
  dx <- label_data[["x"]] - origin[[1]]
  dy <- label_data[["y"]] - origin[[2]]
  r <- sqrt(dx^2 + dy^2)
  center <- r < 1e-8
  r[center] <- 1
  label_data[["label_x"]] <- origin[[1]] + dx / r * (r * expand)
  label_data[["label_y"]] <- origin[[2]] + dy / r * (r * expand)
  label_data[["label_x"]][center] <- origin[[1]]
  label_data[["label_y"]][center] <- origin[[2]]
  label_data[["hjust"]] <- ifelse(label_data[["label_x"]] >= origin[[1]], 0, 1)
  label_data[["hjust"]][center] <- 0.5
  label_data[["vjust"]] <- 0.5
  label_data
}

scenic_network_limits <- function(node_data, label_data, use_radial = FALSE) {
  xs <- node_data[["x"]]
  ys <- node_data[["y"]]
  if (isTRUE(use_radial) && nrow(label_data) > 0L && "label_x" %in% colnames(label_data)) {
    xs <- c(xs, label_data[["label_x"]])
    ys <- c(ys, label_data[["label_y"]])
  }
  pad <- 0.32 * max(diff(range(xs)), diff(range(ys)), 1)
  list(
    x = range(xs) + c(-pad, pad),
    y = range(ys) + c(-pad, pad)
  )
}

scenic_regulon_jaccard <- function(regulon_list) {
  regulons <- names(regulon_list)
  n <- length(regulons)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(regulons, regulons))
  sets <- lapply(regulon_list, function(x) unique(as.character(x)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      a <- sets[[i]]
      b <- sets[[j]]
      union_n <- length(union(a, b))
      value <- if (union_n == 0L) 0 else length(intersect(a, b)) / union_n
      mat[i, j] <- value
      mat[j, i] <- value
    }
  }
  mat
}

scenic_bind_named_dfs <- function(df_list, name_col) {
  if (length(df_list) == 0L) {
    return(NULL)
  }
  out <- do.call(
    rbind,
    lapply(names(df_list), function(nm) {
      df <- df_list[[nm]]
      df[[name_col]] <- nm
      df
    })
  )
  rownames(out) <- NULL
  out
}
