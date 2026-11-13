grn_native_table <- function(native, path) {
  value <- native
  for (key in path) {
    if (!is.list(value) || !key %in% names(value)) stop("The native result does not contain ", paste(path, collapse = "/"), call. = FALSE)
    value <- value[[key]]
  }
  value
}

grn_layer_edges <- function(group, layer) {
  spec <- group$layer_specs[[layer]]
  if (is.null(spec) || !identical(spec$status, "computed")) stop("Layer '", layer, "' is ", spec$status %||% "not_supported", ".", call. = FALSE)
  native <- grn_native_table(group$native, spec$table_path)
  if (!is.data.frame(native)) stop("A computed GRN layer must reference a native table.", call. = FALSE)
  if (!is.null(spec$rows)) {
    rows <- spec$rows
    if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) || any(rows < 1 | rows > nrow(native) | rows != as.integer(rows)) || anyDuplicated(rows)) stop("Invalid native GRN row selection.", call. = FALSE)
    native <- native[as.integer(rows), , drop = FALSE]
  }
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

grn_validate_result <- function(result) {
  if (!inherits(result, "GRNResult") || !identical(result$schema_version, "1.0") || !length(result$groups) || is.null(names(result$groups)) || anyDuplicated(names(result$groups))) stop("Invalid GRNResult schema or groups.", call. = FALSE)
  for (group in result$groups) {
    if (!identical(group$execution$status, "completed") || is.null(group$cells) || anyDuplicated(group$cells)) stop("A GRN group is incomplete or has invalid cells.", call. = FALSE)
    for (layer in names(group$layer_specs)) {
      spec <- group$layer_specs[[layer]]
      if (!spec$status %in% c("computed", "not_supported", "not_computed")) stop("Invalid GRN layer status.", call. = FALSE)
      if (identical(spec$status, "computed")) {
        edges <- grn_layer_edges(group, layer)
        if (anyNA(edges$source) || anyNA(edges$target) || any(!nzchar(edges$source)) || any(!nzchar(edges$target))) stop("GRN endpoints must have complete identifiers.", call. = FALSE)
      }
    }
  }
  invisible(result)
}

grn_direct_gene_sets <- function(group) {
  edges <- grn_layer_edges(group, "tf_gene")
  if (!nrow(edges)) return(list())
  labels <- ifelse(is.na(edges$sign), "unknown", ifelse(edges$sign > 0, "positive", ifelse(edges$sign < 0, "negative", "zero")))
  key <- paste(edges$source, labels, sep = "__")
  lapply(split(edges$target, key), unique)
}
