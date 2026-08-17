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
  title = NULL
) {
  network_layout <- scenic_network_resolve_layout(network_layout, default = "hub")
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
    label_data = label_data,
    layout = network_layout,
    palette = palette,
    palcolor = palcolor,
    title = title %||% "SCENIC TF–target network",
    curvature = if (network_layout %in% c("star", "hub")) 0.08 else 0.16
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
  title = NULL
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

  if (isTRUE(split_panels)) {
    plots <- list()
    edge_list <- list()
    node_list <- list()
    for (tf in tf_candidates) {
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
        palcolor = palcolor,
        title = tf,
        curvature = 0
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
    title = title %||% "SCENIC TF–target network",
    curvature = if (layout_use %in% c("star", "hub")) 0 else 0.16
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
  title = NULL
) {
  tf_candidates <- scenic_network_tf_candidates(
    network_tf %||% features %||% highlight_tf %||% unique(top_table[["TF"]])
  )
  layout_use <- scenic_network_resolve_layout(network_layout, default = "tripartite")
  if (!layout_use %in% c("tripartite", "star", "hub", "bipartite", "fr", "kk", "circle")) {
    layout_use <- "tripartite"
  }
  triplets <- scenic_get_triplets(srt, tool_name, required = FALSE)
  if (is.null(triplets)) {
    log_message(
      "{.val egrn} plots need TF-region-gene triplets; drawing a TF-gene network instead",
      message_type = "warning"
    )
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
    title = title %||% "SCENIC+ enhancer-driven GRN",
    curvature = if (layout_use %in% c("star", "hub", "tripartite", "bipartite")) 0.12 else 0.16,
    prefer_triplets = TRUE
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
  prefer_triplets = FALSE
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
  default_labels <- if (network_layout %in% c("star", "tripartite", "bipartite")) {
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
    label_data = label_data,
    layout = network_layout,
    palette = palette,
    palcolor = palcolor,
    title = title,
    curvature = curvature
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
  if (length(regions) > 0L) {
    xy[regions, ] <- scenic_circle_coords(length(regions), radius = 1)
    if (length(genes) > 0L) {
      xy[genes, ] <- scenic_circle_coords(length(genes), radius = 2)
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
  scenic_place_column(regions, 1)
  scenic_place_column(genes, 2)
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

scenic_network_fill_colors <- function(node_types, palette = "RdYlBu", palcolor = NULL) {
  present <- intersect(c("TF", "region", "gene"), as.character(node_types))
  defaults <- c(TF = "#E69F00", region = "#BDBDBD", gene = "#0072B2")
  if (!is.null(palcolor) && length(palcolor) >= length(present)) {
    cols <- stats::setNames(as.character(palcolor)[seq_along(present)], present)
    return(cols)
  }
  if (!identical(palette, "RdYlBu")) {
    return(palette_colors(present, palette = palette, palcolor = palcolor))
  }
  defaults[present]
}

scenic_network_ggplot <- function(
  node_data,
  edge_plot,
  label_data,
  title,
  layout = "star",
  palette = "RdYlBu",
  palcolor = NULL,
  curvature = 0,
  edge_width_range = c(0.25, 1.4)
) {
  node_data[["node_size"]] <- scenic_network_node_size(node_data)
  fill_cols <- scenic_network_fill_colors(node_data[["node_type"]], palette, palcolor)
  edge_plot[["edge_group"]] <- scenic_network_edge_group(edge_plot)
  edge_cols <- c(
    `TF–gene` = "#636363",
    `TF–region` = "#BDBDBD",
    `region–gene` = "#6BAED6",
    positive = "#4D4D4D",
    negative = "#D55E00"
  )
  edge_cols <- edge_cols[intersect(names(edge_cols), unique(edge_plot[["edge_group"]]))]
  tf_nodes <- node_data[as.character(node_data[["node_type"]]) == "TF", , drop = FALSE]
  region_nodes <- node_data[as.character(node_data[["node_type"]]) == "region", , drop = FALSE]
  gene_nodes <- node_data[as.character(node_data[["node_type"]]) == "gene", , drop = FALSE]
  use_radial <- layout %in% c("star") && nrow(tf_nodes) == 1L
  if (isTRUE(use_radial) && nrow(label_data) > 0L) {
    label_data <- scenic_radial_label_coords(label_data, origin = c(tf_nodes[["x"]][[1]], tf_nodes[["y"]][[1]]))
  }

  p <- ggplot2::ggplot()
  arrow_use <- if (layout %in% c("star", "hub") && abs(curvature) < 1e-8) {
    NULL
  } else {
    grid::arrow(length = grid::unit(0.012, "npc"), type = "closed")
  }
  if (abs(curvature) < 1e-8) {
    p <- p + ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["x_end"]],
        yend = .data[["y_end"]],
        color = .data[["edge_group"]],
        linewidth = abs(.data[["weight"]])
      ),
      alpha = 0.7,
      lineend = "round",
      arrow = arrow_use
    )
  } else {
    p <- p + ggplot2::geom_curve(
      data = edge_plot,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["x_end"]],
        yend = .data[["y_end"]],
        color = .data[["edge_group"]],
        linewidth = abs(.data[["weight"]])
      ),
      curvature = curvature,
      alpha = 0.65,
      lineend = "round",
      arrow = arrow_use
    )
  }
  if (nrow(region_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = region_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_type"]],
        size = .data[["node_size"]]
      ),
      shape = 24,
      color = "grey25",
      stroke = 0.35
    )
  }
  if (nrow(gene_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = gene_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_type"]],
        size = .data[["node_size"]]
      ),
      shape = 21,
      color = "white",
      stroke = 0.4
    )
  }
  if (nrow(tf_nodes) > 0L) {
    p <- p + ggplot2::geom_point(
      data = tf_nodes,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_type"]],
        size = .data[["node_size"]]
      ),
      shape = 23,
      color = "grey10",
      stroke = 0.7
    )
  }
  present_nodes <- names(fill_cols)
  p <- p +
    ggplot2::scale_fill_manual(
      values = fill_cols,
      breaks = present_nodes,
      drop = FALSE,
      name = "Node"
    ) +
    ggplot2::scale_color_manual(
      values = edge_cols,
      breaks = names(edge_cols),
      drop = FALSE,
      name = "Edge"
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth(range = edge_width_range, guide = "none") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = unname(c(TF = 23, region = 24, gene = 21)[present_nodes]),
          size = 4,
          color = "black",
          stroke = 0.5
        ),
        order = 1
      ),
      color = ggplot2::guide_legend(order = 2)
    )
  p <- scenic_network_add_labels(p, label_data, tf_nodes, use_radial = use_radial)
  limits <- scenic_network_limits(node_data, label_data, use_radial = use_radial)
  p +
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
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(8, 18, 8, 18)
    ) +
    ggplot2::labs(title = title)
}

scenic_network_node_size <- function(node_data) {
  base <- c(TF = 9.5, region = 3.2, gene = 4.6)
  type <- as.character(node_data[["node_type"]])
  size <- unname(base[type])
  size[is.na(size)] <- 4
  size + 0.35 * log1p(pmax(node_data[["degree"]] - 1, 0))
}

scenic_network_edge_group <- function(edge_plot) {
  if ("edge_type" %in% colnames(edge_plot) &&
    any(edge_plot[["edge_type"]] %in% c("tf_region", "region_gene"))) {
    out <- edge_plot[["edge_type"]]
    out[out == "tf_gene"] <- "TF–gene"
    out[out == "tf_region"] <- "TF–region"
    out[out == "region_gene"] <- "region–gene"
    return(out)
  }
  if ("edge_sign" %in% colnames(edge_plot) && any(edge_plot[["edge_sign"]] == "negative")) {
    return(edge_plot[["edge_sign"]])
  }
  rep("TF–gene", nrow(edge_plot))
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

scenic_network_add_labels <- function(p, label_data, tf_nodes, use_radial = FALSE) {
  if (nrow(label_data) == 0L) {
    return(p)
  }
  tf_names <- tf_nodes[["name"]]
  tf_labels <- label_data[label_data[["name"]] %in% tf_names, , drop = FALSE]
  other_labels <- label_data[!label_data[["name"]] %in% tf_names, , drop = FALSE]
  if (nrow(tf_labels) > 0L) {
    p <- p + ggplot2::geom_text(
      data = tf_labels,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
      fontface = "bold",
      size = 3.6,
      color = "grey10",
      vjust = -1.35,
      show.legend = FALSE
    )
  }
  if (nrow(other_labels) == 0L) {
    return(p)
  }
  if (isTRUE(use_radial) && all(c("label_x", "label_y") %in% colnames(other_labels))) {
    p <- p + ggplot2::geom_text(
      data = other_labels,
      ggplot2::aes(
        x = .data[["label_x"]],
        y = .data[["label_y"]],
        label = .data[["label"]],
        hjust = .data[["hjust"]]
      ),
      size = 3,
      color = "grey15",
      show.legend = FALSE
    )
  } else {
    p <- p + ggrepel::geom_text_repel(
      data = other_labels,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.22,
      point.padding = 0.18,
      segment.color = "grey70",
      segment.size = 0.2,
      min.segment.length = 0.05,
      show.legend = FALSE
    )
  }
  p
}

scenic_network_limits <- function(node_data, label_data, use_radial = FALSE) {
  xs <- node_data[["x"]]
  ys <- node_data[["y"]]
  if (isTRUE(use_radial) && nrow(label_data) > 0L && "label_x" %in% colnames(label_data)) {
    xs <- c(xs, label_data[["label_x"]])
    ys <- c(ys, label_data[["label_y"]])
  }
  pad <- 0.28 * max(diff(range(xs)), diff(range(ys)), 1)
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
