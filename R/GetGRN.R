#' @title Read a unified gene regulatory network result
#'
#' @description Read complete GRN results without rerunning inference. Standard edge views preserve original scores, direction and sign evidence.
#'
#' @param object A GRNResult or a Seurat object containing one.
#' @param group One fitted group. Required when more than one group is present.
#' @param layer Native evidence layer, such as tf_gene, tf_region, region_gene, or mediated_tf_gene.
#' @param what One of edges, layers, chains, gene_sets, native, parameters, or result.
#' @param tool_name Name of the Seurat \code{tools} entry. If omitted, the only
#'   stored GRN result is used. Specify a name when several results are stored.
#' @param all_groups Return a named list for every fitted group.
#' @param ... Reserved for future extensions.
#' @return The selected table, native result, parameters, gene sets, or complete result.
#'
#' @details Native scientific data are stored once. Edge views map native columns to source,
#'   target, source_type, target_type, layer, raw_score, score_type, directed, sign,
#'   and sign_basis. Additional native scientific columns remain available.
#'   Unknown signs are NA. Unsupported and uncomputed layers give explicit errors;
#'   a computed zero-row table is a successful empty network.
#'   When native gene sets are available, gene_sets returns them. Otherwise direct
#'   TF-gene targets are grouped by source and positive, negative, zero, or unknown
#'   sign. These are predicted target programs and can be passed to CellScoring.
#'   The native selector returns the unchanged method output.
#' @seealso \code{\link{RunGRN}}, \code{\link{GRNPlot}}, \code{\link{CellScoring}}
#' @export
GetGRN <- function(
  object,
  group = NULL,
  layer = "tf_gene",
  what = "edges",
  tool_name = NULL,
  all_groups = FALSE,
  ...
) {
  if (length(list(...))) stop("Unused GetGRN arguments.", call. = FALSE)
  if (inherits(object, "Seurat")) {
    tools <- object@tools
    candidates <- names(tools)[vapply(tools, inherits, logical(1), what = "GRNResult")]
    if (is.null(tool_name)) {
      if (!length(candidates)) stop("No GRNResult is stored in this Seurat object.", call. = FALSE)
      if (length(candidates) > 1L) {
        stop("Multiple GRN results are stored. Specify tool_name: ", paste(candidates, collapse = ", "), ".", call. = FALSE)
      }
      tool_name <- candidates[[1L]]
    }
    if (!is.character(tool_name) || length(tool_name) != 1L || !tool_name %in% candidates) {
      stop("tool_name must identify a stored GRNResult.", call. = FALSE)
    }
    object <- tools[[tool_name]]
  }
  if (!inherits(object, "GRNResult") || !identical(object$schema_version, "1.0")) stop("object must contain a supported GRNResult.", call. = FALSE)
  what <- match.arg(what, c("edges", "layers", "chains", "gene_sets", "native", "parameters", "result"))
  if (identical(what, "result")) {
    return(object)
  }
  available <- names(object$groups)
  if (isTRUE(all_groups)) {
    if (!is.null(group)) stop("Choose group or all_groups, rather than both.", call. = FALSE)
    selected <- available
  } else {
    if (is.null(group) && length(available) != 1L) stop("Select a group or use all_groups = TRUE.", call. = FALSE)
    selected <- group %||% available
    if (length(selected) != 1L || !selected %in% available) stop("Unknown or ambiguous GRN group.", call. = FALSE)
  }
  read_group <- function(name) {
    entry <- object$groups[[name]]
    if (identical(what, "native")) {
      return(entry$native)
    }
    if (identical(what, "parameters")) {
      return(entry$parameters)
    }
    if (identical(what, "edges")) {
      return(grn_layer_edges(entry, layer))
    }
    if (identical(what, "layers")) {
      return(stats::setNames(lapply(names(entry$layer_specs), function(key) {
        spec <- entry$layer_specs[[key]]
        if (identical(spec$status, "computed")) grn_layer_edges(entry, key) else list(status = spec$status)
      }), names(entry$layer_specs)))
    }
    descriptor <- entry[[what]]
    if (identical(what, "gene_sets") && (is.null(descriptor) || !identical(descriptor$status, "computed"))) {
      return(grn_direct_gene_sets(entry))
    }
    if (is.null(descriptor) || !identical(descriptor$status, "computed")) stop(what, " is ", descriptor$status %||% "not_supported", ".", call. = FALSE)
    grn_native_table(entry$native, descriptor$table_path, descriptor$rows)
  }
  if (isTRUE(all_groups)) {
    return(stats::setNames(lapply(selected, read_group), selected))
  }
  read_group(selected)
}

grn_native_table <- function(native, path, rows = NULL) {
  value <- native
  for (key in path) {
    if (!is.list(value) || !key %in% names(value)) stop("The native result does not contain ", paste(path, collapse = "/"), call. = FALSE)
    value <- value[[key]]
  }
  if (is.null(rows)) {
    return(value)
  }
  size <- if (is.data.frame(value)) nrow(value) else length(value)
  if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) ||
    any(rows < 1 | rows > size | rows != as.integer(rows)) || anyDuplicated(rows)) {
    stop("Invalid native GRN row selection.", call. = FALSE)
  }
  rows <- as.integer(rows)
  if (is.data.frame(value)) value[rows, , drop = FALSE] else value[rows]
}

grn_layer_edges <- function(group, layer) {
  spec <- group$layer_specs[[layer]]
  if (is.null(spec) || !identical(spec$status, "computed")) stop("Layer '", layer, "' is ", spec$status %||% "not_supported", ".", call. = FALSE)
  native <- grn_native_table(group$native, spec$table_path, spec$rows)
  if (!is.data.frame(native)) stop("A computed GRN layer must reference a native table.", call. = FALSE)
  required <- c(spec$source, spec$target, spec$score, spec$sign_column)
  if (!all(required %in% names(native))) stop("The native table is missing mapped GRN columns.", call. = FALSE)
  score <- native[[spec$score]]
  if (!is.numeric(score)) stop("Native GRN scores must be numeric.", call. = FALSE)
  direction <- if (is.null(spec$sign_column)) rep(NA_integer_, nrow(native)) else as.integer(sign(native[[spec$sign_column]]))
  edges <- data.frame(
    source = as.character(native[[spec$source]]), target = as.character(native[[spec$target]]),
    source_type = rep(spec$source_type, nrow(native)), target_type = rep(spec$target_type, nrow(native)),
    layer = rep(layer, nrow(native)), raw_score = score, score_type = rep(spec$score_type, nrow(native)),
    directed = rep(isTRUE(spec$directed), nrow(native)), sign = direction,
    sign_basis = rep(spec$sign_basis %||% "unknown", nrow(native)), stringsAsFactors = FALSE
  )
  extra <- native[, setdiff(names(native), c(spec$source, spec$target, spec$score)), drop = FALSE]
  names(extra) <- vapply(names(extra), function(name) {
    while (name %in% c(names(edges), setdiff(names(extra), name))) name <- paste0("native_", name)
    name
  }, character(1))
  rownames(extra) <- NULL
  cbind(edges, extra)
}

grn_direct_gene_sets <- function(group) {
  edges <- grn_layer_edges(group, "tf_gene")
  if (!nrow(edges)) {
    return(list())
  }
  labels <- ifelse(is.na(edges$sign), "unknown", ifelse(edges$sign > 0, "positive", ifelse(edges$sign < 0, "negative", "zero")))
  key <- paste(edges$source, labels, sep = "__")
  lapply(split(edges$target, key), unique)
}
