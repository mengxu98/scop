#' @title Plot a gene regulatory network
#'
#' @md
#' @inheritParams GetGRN
#' @param type `"network"` for edges, `"layered"` for native chains, or
#'   `"overlap"` for shared edges within the common identifier scope.
#' @param top_n Maximum edges or chains to show per result. `Inf` shows all.
#' @param layout `"auto"`, `"star"`, `"hub"`, `"tripartite"`, `"kk"`, or `"fr"`.
#' @param label_nodes Show node labels.
#' @param compare Another GRN result or a named list of results for `"overlap"`.
#' @param compare_group One group per comparison, or a single group for all.
#'   Defaults to `group`.
#' @param title Plot title.
#' @return A ggplot object.
#' @seealso [RunGRN()], [GetGRN()]
#' @export
GRNPlot <- function(
  object,
  group = NULL,
  layer = "tf_gene",
  type = c("network", "layered", "overlap"),
  top_n = 100,
  layout = "auto",
  label_nodes = TRUE,
  compare = NULL,
  compare_group = NULL,
  tool_name = NULL,
  title = NULL
) {
  type <- match.arg(type)
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1 ||
    (!is.infinite(top_n) && (top_n != floor(top_n) || top_n > .Machine$integer.max))) {
    stop("`top_n` must be a positive integer or Inf.", call. = FALSE)
  }
  if (!is.logical(label_nodes) || length(label_nodes) != 1L || is.na(label_nodes)) {
    stop("`label_nodes` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(layout) || length(layout) != 1L ||
    !layout %in% c("auto", "star", "hub", "tripartite", "kk", "fr")) {
    stop("`layout` must be auto, star, hub, tripartite, kk or fr.", call. = FALSE)
  }
  result <- GetGRN(object, what = "result", tool_name = tool_name)
  group <- grn_select_group(names(result$groups), group)
  if (identical(type, "overlap")) {
    if (is.null(compare)) stop("`compare` must supply at least one GRN result for an overlap plot.", call. = FALSE)
    comparisons <- if (inherits(compare, c("GRNResult", "Seurat"))) list(compare) else compare
    if (!is.list(comparisons) || !length(comparisons)) stop("`compare` must be a GRN result or a list of results.", call. = FALSE)
    comparison_results <- lapply(comparisons, function(comparison) {
      comparison_tool <- if (inherits(comparison, "Seurat") &&
        inherits(comparison@tools[[tool_name %||% ""]], "GRNResult")) {
        tool_name
      } else {
        NULL
      }
      GetGRN(comparison, what = "result", tool_name = comparison_tool)
    })
    results <- c(list(result), comparison_results)
    labels <- c(result$method, names(comparisons))
    method_names <- vapply(results, `[[`, character(1), "method")
    if (length(labels) != length(results)) labels <- method_names
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- method_names[missing]
    labels <- make.unique(labels)
    if (!is.null(compare_group) && (!is.character(compare_group) || !length(compare_group) %in% c(1L, length(comparisons)))) {
      stop("`compare_group` must name one group or one group per comparison result.", call. = FALSE)
    }
    groups <- c(group, rep_len(compare_group %||% group, length(comparisons)))
    tables <- lapply(seq_along(results), function(i) GetGRN(results[[i]], group = groups[[i]], layer = layer, what = "edges"))
    endpoints <- lapply(tables, function(table) {
      list(source = grn_plot_node_key(table$source, table$source_type),
        target = grn_plot_node_key(table$target, table$target_type))
    })
    source_universe <- Reduce(intersect, Map(function(table, ids) {
      unique(c(ids$source, ids$target[!table$directed]))
    }, tables, endpoints))
    target_universe <- Reduce(intersect, Map(function(table, ids) {
      unique(c(ids$target, ids$source[!table$directed]))
    }, tables, endpoints))
    counts <- data.frame(method = labels, total = integer(length(tables)), in_scope = integer(length(tables)), shown = integer(length(tables)))
    keys <- lapply(seq_along(tables), function(i) {
      table <- tables[[i]]
      source <- endpoints[[i]]$source
      target <- endpoints[[i]]$target
      edge_keys <- grn_plot_edge_key(source, target, table$directed)
      counts$total[[i]] <<- length(unique(edge_keys))
      in_scope <- source %in% source_universe & target %in% target_universe
      in_scope <- in_scope | (!table$directed & target %in% source_universe & source %in% target_universe)
      edge_keys <- edge_keys[in_scope]
      counts$in_scope[[i]] <<- length(unique(edge_keys))
      rank_by <- results[[i]]$groups[[groups[[i]]]]$layer_specs[[layer]]$rank_by %||% "absolute_score"
      table <- table[in_scope, , drop = FALSE]
      rows <- grn_plot_order(table, top_n, "raw_score", rank_by, edge_keys)
      counts$shown[[i]] <<- length(rows)
      edge_keys[rows]
    })
    overlap <- expand.grid(method_1 = labels, method_2 = labels, stringsAsFactors = FALSE)
    overlap$shared_edges <- vapply(seq_len(nrow(overlap)), function(i) {
      length(intersect(keys[[match(overlap$method_1[[i]], labels)]], keys[[match(overlap$method_2[[i]], labels)]]))
    }, integer(1))
    overlap$method_1 <- factor(overlap$method_1, levels = labels)
    overlap$method_2 <- factor(overlap$method_2, levels = rev(labels))
    scope_label <- sprintf("Exact identifiers; %d common sources and %d common targets\nTop %s per result within this scope", length(source_universe), length(target_universe), format(top_n, trim = TRUE))
    plot <- ggplot2::ggplot(overlap, ggplot2::aes(x = .data$method_1, y = .data$method_2, fill = .data$shared_edges)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = .data$shared_edges), size = 4) +
      ggplot2::scale_fill_gradient(low = "#F5F8FB", high = "#2177B8", name = "Shared edges") +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = title %||% paste(layer, "edge overlap"), x = NULL, y = NULL,
        subtitle = scope_label,
        caption = paste(sprintf("%s: %d of %d in-scope edges shown (%d total)", counts$method, counts$shown, counts$in_scope, counts$total), collapse = "\n")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    attr(plot, "grn_data") <- list(overlap = overlap, counts = counts, source_universe = source_universe, target_universe = target_universe)
    return(plot)
  }
  if (!is.null(compare)) stop("`compare` is used only with type = 'overlap'.", call. = FALSE)
  chains <- NULL
  if (identical(type, "layered")) {
    chains <- GetGRN(result, group = group, what = "chains")
    descriptor <- result$groups[[group]]$chains
    source_col <- intersect(descriptor$source %||% c("regulator", "TF", "tf", "source"), colnames(chains))
    target_col <- intersect(descriptor$target %||% c("target", "gene"), colnames(chains))
    region_col <- descriptor$region %||% "region"
    if (!length(source_col) || !length(target_col) || !region_col %in% colnames(chains)) {
      stop("A layered plot requires native chains with regulator, region and target identifiers.", call. = FALSE)
    }
    score_col <- intersect(descriptor$score %||% c("weight", "chain_delta_bic", "score", "estimate"), colnames(chains))
    all_chains <- chains
    if (length(score_col)) {
      rows <- grn_plot_order(chains, top_n, score_col[[1L]], descriptor$rank_by %||% "absolute_score")
      chains <- chains[rows, , drop = FALSE]
    } else if (nrow(chains) > top_n) {
      stop("Native chains need a score column to select `top_n`; use top_n = Inf to display all chains.", call. = FALSE)
    }
    layers <- lapply(c("tf_region", "region_gene"), function(name) GetGRN(result, group = group, layer = name, what = "edges"))
    source <- as.character(chains[[source_col[[1L]]]])
    target <- as.character(chains[[target_col[[1L]]]])
    region <- as.character(chains[[region_col]])
    tf_keys <- grn_plot_pair_key(layers[[1L]]$source, layers[[1L]]$target)
    gene_keys <- grn_plot_pair_key(layers[[2L]]$source, layers[[2L]]$target)
    selected_tf <- grn_plot_pair_key(source, region)
    selected_gene <- grn_plot_pair_key(region, target)
    if (!all(selected_tf %in% tf_keys) || !all(selected_gene %in% gene_keys)) {
      stop("Selected native chains reference edges absent from their native layers.", call. = FALSE)
    }
    keep_tf <- tf_keys %in% selected_tf
    keep_gene <- gene_keys %in% selected_gene
    common_columns <- Reduce(intersect, lapply(layers, colnames))
    edges <- rbind(layers[[1L]][keep_tf, common_columns, drop = FALSE], layers[[2L]][keep_gene, common_columns, drop = FALSE])
    subtitle <- sprintf("%s / %s: %d of %d native chains; %d layer edges", result$method, group, nrow(chains), nrow(all_chains), nrow(edges))
  } else {
    edges <- GetGRN(result, group = group, layer = layer, what = "edges")
    rank_by <- result$groups[[group]]$layer_specs[[layer]]$rank_by %||% "absolute_score"
    keys <- NULL
    if (any(!edges$directed)) {
      source <- grn_plot_node_key(edges$source, edges$source_type)
      target <- grn_plot_node_key(edges$target, edges$target_type)
      keys <- grn_plot_edge_key(source, target, edges$directed)
      keys[edges$directed] <- paste0("d", which(edges$directed), recycle0 = TRUE)
    }
    total <- if (is.null(keys)) nrow(edges) else length(unique(keys))
    rows <- grn_plot_order(edges, top_n, "raw_score", rank_by, keys)
    edges <- edges[rows, , drop = FALSE]
    subtitle <- sprintf("%s / %s: %d of %d %s edges", result$method, group, nrow(edges), total, layer)
  }
  plot <- grn_plot_network(edges, identical(type, "layered"), layout, label_nodes, title %||% "Gene regulatory network", subtitle)
  view <- attr(plot, "grn_data")
  view$chains <- chains
  view$group <- group
  view$method <- result$method
  view$top_n <- top_n
  attr(plot, "grn_data") <- view
  plot
}

grn_plot_pair_key <- function(first, second) {
  first <- enc2utf8(first)
  paste0(nchar(first, type = "bytes"), ":", first, enc2utf8(second), recycle0 = TRUE)
}

grn_plot_node_key <- function(ids, types, layered = FALSE) {
  types[tolower(types) == "tf"] <- if (layered) "TF" else "gene"
  grn_plot_pair_key(types, ids)
}

grn_plot_edge_key <- function(source, target, directed) {
  if (!length(source)) return(character())
  ids <- sort(unique(c(source, target)), method = "radix")
  swap <- !directed & match(source, ids) > match(target, ids)
  paste0(ifelse(directed, "d", "u"),
    grn_plot_pair_key(ifelse(swap, target, source), ifelse(swap, source, target)))
}

grn_plot_order <- function(table, top_n, score, rank_by, keys = NULL) {
  if (!nrow(table) || (is.infinite(top_n) && !anyDuplicated(keys))) {
    return(seq_len(nrow(table)))
  }
  values <- table[[score]]
  if (identical(rank_by, "absolute_score")) values <- abs(values)
  if (!identical(rank_by, "ascending_score")) values <- -values
  rows <- order(values, na.last = TRUE, method = "radix")
  if (!is.null(keys)) rows <- rows[!duplicated(keys[rows])]
  utils::head(rows, top_n)
}

grn_plot_network <- function(edges, layered, layout, label_nodes, title, subtitle) {
  if (!nrow(edges)) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No edges to display", colour = "#828B91") +
      ggplot2::theme_void() +
      ggplot2::labs(title = title, subtitle = subtitle)
    attr(plot, "grn_data") <- list(edges = edges, nodes = data.frame())
    return(plot)
  }
  node_type <- function(x) {
    x[tolower(x) == "tf"] <- "TF"
    x
  }
  source_types <- node_type(as.character(edges$source_type))
  target_types <- node_type(as.character(edges$target_type))
  edges$from <- grn_plot_node_key(edges$source, source_types, layered)
  edges$to <- grn_plot_node_key(edges$target, target_types, layered)
  nodes <- rbind(
    data.frame(name = edges$from, label = edges$source, node_type = source_types),
    data.frame(name = edges$to, label = edges$target, node_type = target_types)
  )
  nodes <- nodes[order(nodes$node_type != "TF"), , drop = FALSE]
  nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]
  rownames(nodes) <- NULL
  graph_edges <- data.frame(from = edges$from, to = edges$to, weight = abs(edges$raw_score))
  graph_edges$weight[!is.finite(graph_edges$weight) | graph_edges$weight <= 0] <- 1
  graph <- igraph::graph_from_data_frame(graph_edges, directed = any(edges$directed), vertices = nodes)
  layout <- if (identical(layout, "auto")) if (layered) "tripartite" else "kk" else layout
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(
    {
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    },
    add = TRUE
  )
  set.seed(1234)
  xy <- scenic_network_layout(graph, layout = layout, node_type = stats::setNames(nodes$node_type, nodes$name))
  nodes$x <- xy[nodes$name, 1L]
  nodes$y <- xy[nodes$name, 2L]
  edges$x <- nodes$x[match(edges$from, nodes$name)]
  edges$y <- nodes$y[match(edges$from, nodes$name)]
  edges$xend <- nodes$x[match(edges$to, nodes$name)]
  edges$yend <- nodes$y[match(edges$to, nodes$name)]
  distance <- sqrt((edges$xend - edges$x)^2 + (edges$yend - edges$y)^2)
  trim <- pmin(0.012 * max(diff(range(nodes$x)), diff(range(nodes$y))), distance / 3)
  fraction <- ifelse(distance > 0, trim / distance, 0)
  dx <- edges$xend - edges$x
  dy <- edges$yend - edges$y
  edges$xend <- edges$xend - dx * fraction
  edges$yend <- edges$yend - dy * fraction
  edges$x <- edges$x + dx * fraction
  edges$y <- edges$y + dy * fraction
  edges$evidence <- ifelse(
    edges$score_type == "membership", "Membership",
    ifelse(is.na(edges$sign), "Unknown", ifelse(edges$sign == 0, "Zero", ifelse(edges$sign > 0, "Positive", "Negative")))
  )
  edges$evidence <- factor(edges$evidence, levels = c("Positive", "Negative", "Zero", "Unknown", "Membership"))
  edges$edge_width <- abs(edges$raw_score)
  edges$edge_width[!is.finite(edges$edge_width)] <- 0
  if (length(unique(edges$edge_width)) == 1L) edges$edge_width <- 1
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges[!edges$directed, , drop = FALSE],
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = .data$evidence, linetype = .data$layer, linewidth = .data$edge_width),
      alpha = 0.7
    ) +
    ggplot2::geom_segment(
      data = edges[edges$directed, , drop = FALSE],
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = .data$evidence, linetype = .data$layer, linewidth = .data$edge_width),
      arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed"), alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = nodes, ggplot2::aes(x = x, y = y, fill = node_type, shape = node_type),
      size = 3.2, colour = "#38434D", stroke = 0.4
    ) +
    ggplot2::scale_colour_manual(
      values = c(Positive = "#B34740", Negative = "#2177B8", Zero = "#A7ADB1", Unknown = "#828B91", Membership = "#9570A6"),
      name = "Sign or evidence", drop = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c(TF = "#E7B97B", gene = "#D6E6F2", region = "#C8DDBB"), name = "Node type") +
    ggplot2::scale_shape_manual(values = c(TF = 21, gene = 22, region = 23), name = "Node type") +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.1), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, linetype = "Layer") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = grid::unit(c(8, 14, 8, 14), "mm"), legend.position = "right")
  if (isTRUE(label_nodes)) {
    nodes$label <- scenic_network_node_labels(nodes$label, nodes$node_type)
    plot <- plot + ggrepel::geom_text_repel(
      data = nodes, ggplot2::aes(x = x, y = y, label = label),
      size = 3, seed = 1234, max.overlaps = Inf, min.segment.length = 0
    )
  }
  attr(plot, "grn_data") <- list(edges = edges, nodes = nodes)
  plot
}
