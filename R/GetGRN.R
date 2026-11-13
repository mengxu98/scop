GetGRN <- function(object, group = NULL, layer = "tf_gene", what = "edges", tool_name = NULL, all_groups = FALSE, ...) {
  if (length(list(...))) stop("Unused GetGRN arguments.", call. = FALSE)
  if (inherits(object, "Seurat")) {
    if (is.null(tool_name) || !tool_name %in% names(object@tools)) stop("Supply a tool_name containing a GRNResult.", call. = FALSE)
    object <- object@tools[[tool_name]]
  }
  if (!inherits(object, "GRNResult") || !identical(object$schema_version, "1.0")) stop("object must contain a supported GRNResult.", call. = FALSE)
  what <- match.arg(what, c("edges", "layers", "chains", "gene_sets", "native", "parameters", "result"))
  if (identical(what, "result")) return(object)
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
    if (identical(what, "native")) return(entry$native)
    if (identical(what, "parameters")) return(entry$parameters)
    if (identical(what, "edges")) return(grn_layer_edges(entry, layer))
    if (identical(what, "layers")) return(stats::setNames(lapply(names(entry$layer_specs), function(key) {
      spec <- entry$layer_specs[[key]]
      if (identical(spec$status, "computed")) grn_layer_edges(entry, key) else list(status = spec$status)
    }), names(entry$layer_specs)))
    descriptor <- entry[[what]]
    if (identical(what, "gene_sets") && (is.null(descriptor) || !identical(descriptor$status, "computed"))) return(grn_direct_gene_sets(entry))
    if (is.null(descriptor) || !identical(descriptor$status, "computed")) stop(what, " is ", descriptor$status %||% "not_supported", ".", call. = FALSE)
    value <- grn_native_table(entry$native, descriptor$table_path)
    if (!is.null(descriptor$rows)) {
      rows <- descriptor$rows
      size <- if (is.data.frame(value)) nrow(value) else length(value)
      if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) || any(rows < 1 | rows > size | rows != as.integer(rows)) || anyDuplicated(rows)) stop("Invalid native GRN row selection.", call. = FALSE)
      value <- if (is.data.frame(value)) value[as.integer(rows), , drop = FALSE] else value[as.integer(rows)]
    }
    value
  }
  if (isTRUE(all_groups)) return(stats::setNames(lapply(selected, read_group), selected))
  read_group(selected)
}
