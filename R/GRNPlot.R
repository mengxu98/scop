#' @title Plot networks from a unified GRN result
#'
#' @description Read network views through \code{\link{GetGRN}} and draw nodes, edges and native chains while retaining direction, sign and membership evidence.
#'
#' @param object A \code{GRNResult}, or a Seurat object containing one in \code{tools}.
#' @param group One result group. Required when the result contains multiple groups.
#' @param layer Network layer to display or compare.
#' @param type A network, a layered view of native chains, or a pairwise edge-overlap heatmap.
#' @param top_n Maximum edges in a network view, native chains in a layered view, or edges per method within the shared comparison scope. Use \code{Inf} for all entries.
#' @param layout One of \code{"auto"}, \code{"star"}, \code{"hub"}, \code{"tripartite"}, \code{"kk"} or \code{"fr"}. The automatic layout uses tripartite columns for layered views and Kamada--Kawai for other networks.
#' @param label_nodes Whether to label displayed nodes.
#' @param compare For \code{type = "overlap"}, another result or a named list of results. List names provide plot labels.
#' @param compare_group Group in each comparison result. Defaults to \code{group}. A vector supplies one group per comparison result.
#' @param tool_name Name of the Seurat \code{tools} entry. If omitted, the only
#'   stored GRN result is used. Specify a name when several results are stored.
#'   An overlap comparison uses the same name if present; otherwise it selects a
#'   single stored result. Use \code{GetGRN(what = "result")} to select a differently
#'   named result from a comparison with multiple networks.
#' @param title Optional plot title.
#' @return A ggplot object. Its \code{grn_data} attribute contains the displayed edges and nodes, or overlap counts and the comparison scope. The original result is unchanged.
#'
#' @details Signs are labelled positive, negative, zero or unknown. Nonnegative importance scores do not acquire an activation label. Membership edges have a separate evidence style. Arrowheads are drawn only for directed edges. Node shapes distinguish TFs, genes and regions.
#'
#'   A layered view selects existing native chains and then selects the corresponding TF--region and region--gene edges from their native layers. Each segment retains its own sign. Shared nodes can connect several selected chains; the displayed chain table is retained in the plot attribute.
#'
#'   Ranking and \code{top_n} apply only to the displayed view. Captions report the selected size. An overlap heatmap restricts each result to exactly matched source and target identifiers present in all compared results, then selects the top edges independently within each method. TF and gene endpoints share a gene identifier space. Directed and undirected edges remain distinct. Reversed endpoints of an undirected edge denote the same edge. Declared \code{input_summary$entity_namespace} values must agree. No identifier conversion, score aggregation or consensus inference is performed.
#'
#'   The network layout reuses the existing SCENIC graph-layout implementation. Layout and label placement are reproducible and preserve the caller's random-number state. Use existing \code{CellScoring}, \code{FeatureDimPlot} and \code{GroupHeatmap} functions to display extracted gene-set scores.
#' @seealso \code{\link{RunGRN}}, \code{\link{GetGRN}}
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
  if (is.null(group)) {
    if (length(result$groups) != 1L) stop("Choose `group` for a result with multiple groups.", call. = FALSE)
    group <- names(result$groups)[[1L]]
  }
  if (!is.character(group) || length(group) != 1L || !group %in% names(result$groups)) {
    stop("`group` must name one result group.", call. = FALSE)
  }
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
    if (length(labels) != length(results)) labels <- vapply(results, function(x) x$method, character(1))
    labels[is.na(labels) | !nzchar(labels)] <- vapply(results, function(x) x$method, character(1))[is.na(labels) | !nzchar(labels)]
    labels <- make.unique(labels)
    if (!is.null(compare_group) && (!is.character(compare_group) || !length(compare_group) %in% c(1L, length(comparisons)))) {
      stop("`compare_group` must name one group or one group per comparison result.", call. = FALSE)
    }
    groups <- c(group, rep_len(compare_group %||% group, length(comparisons)))
    tables <- lapply(seq_along(results), function(i) GetGRN(results[[i]], group = groups[[i]], layer = layer, what = "edges"))
    namespaces <- lapply(results, function(x) x$input_summary$entity_namespace)
    declared <- namespaces[!vapply(namespaces, is.null, logical(1))]
    if (length(declared) > 1L && !all(vapply(declared, identical, logical(1), declared[[1L]]))) {
      stop("Compared results declare different entity namespaces; align identifiers before plotting.", call. = FALSE)
    }
    entity_key <- function(ids, types) paste(ifelse(tolower(types) == "tf", "gene", types), ids, sep = "\r")
    source_universe <- Reduce(intersect, lapply(tables, function(x) unique(c(entity_key(x$source, x$source_type), entity_key(x$target[!x$directed], x$target_type[!x$directed])))))
    target_universe <- Reduce(intersect, lapply(tables, function(x) unique(c(entity_key(x$target, x$target_type), entity_key(x$source[!x$directed], x$source_type[!x$directed])))))
    counts <- data.frame(method = labels, total = vapply(tables, nrow, integer(1)), in_scope = integer(length(tables)), shown = integer(length(tables)))
    keys <- lapply(seq_along(tables), function(i) {
      table <- tables[[i]]
      source <- entity_key(table$source, table$source_type)
      target <- entity_key(table$target, table$target_type)
      in_scope <- source %in% source_universe & target %in% target_universe
      in_scope <- in_scope | (!table$directed & target %in% source_universe & source %in% target_universe)
      table <- table[in_scope, , drop = FALSE]
      counts$in_scope[[i]] <<- nrow(table)
      rank_by <- results[[i]]$groups[[groups[[i]]]]$layer_specs[[layer]]$rank_by %||% "absolute_score"
      table <- grn_plot_rank(table, top_n, "raw_score", rank_by)
      counts$shown[[i]] <<- nrow(table)
      source <- entity_key(table$source, table$source_type)
      target <- entity_key(table$target, table$target_type)
      swap <- !table$directed & source > target
      from <- ifelse(swap, target, source)
      to <- ifelse(swap, source, target)
      unique(paste(from, to, table$directed, sep = "\r"))
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
      chains <- grn_plot_rank(chains, top_n, score_col[[1L]], descriptor$rank_by %||% "absolute_score")
    } else if (nrow(chains) > top_n) {
      stop("Native chains need a score column to select `top_n`; use top_n = Inf to display all chains.", call. = FALSE)
    }
    layers <- lapply(c("tf_region", "region_gene"), function(name) GetGRN(result, group = group, layer = name, what = "edges"))
    source <- as.character(chains[[source_col[[1L]]]])
    target <- as.character(chains[[target_col[[1L]]]])
    region <- as.character(chains[[region_col]])
    tf_keys <- paste(layers[[1L]]$source, layers[[1L]]$target, sep = "\r")
    gene_keys <- paste(layers[[2L]]$source, layers[[2L]]$target, sep = "\r")
    selected_tf <- paste(source, region, sep = "\r")
    selected_gene <- paste(region, target, sep = "\r")
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
    total <- nrow(edges)
    rank_by <- result$groups[[group]]$layer_specs[[layer]]$rank_by %||% "absolute_score"
    edges <- grn_plot_rank(edges, top_n, "raw_score", rank_by)
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

grn_plot_rank <- function(table, top_n, score, rank_by = "absolute_score") {
  if (!nrow(table) || is.infinite(top_n)) {
    return(table)
  }
  values <- table[[score]]
  if (identical(rank_by, "absolute_score")) {
    values <- abs(values)
  }
  if (!identical(rank_by, "ascending_score")) values <- -values
  table[utils::head(order(values, na.last = TRUE, method = "radix"), top_n), , drop = FALSE]
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
  entity_type <- function(x) {
    if (!layered) x[x == "TF"] <- "gene"
    x
  }
  edges$from <- paste(entity_type(source_types), edges$source, sep = "\r")
  edges$to <- paste(entity_type(target_types), edges$target, sep = "\r")
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
