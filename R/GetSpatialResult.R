#' Access spatial analysis results
#'
#' @description
#' Read the cluster assignments, parameters, and summary for a spatial method
#' from a Seurat object's stored tool result. If detailed results were not
#' stored, cluster assignments are read from the corresponding metadata column.
#' BANKSY retains a lightweight index for custom result columns and exact sample
#' membership. Sample retrieval requires this index or stored per-sample results;
#' sample identity is never inferred from cluster labels.
#'
#' Results describe cells present in `object`, in its current cell order, rather
#' than the original fit's full cell set. BANKSY domain counts and PRECAST cell,
#' domain, and sample counts are recomputed from the returned assignments.
#' Unfitted PRECAST cells are excluded. Parameters and the native bundles in
#' `object@tools` still describe the original fit; other methods' stored summaries
#' are omitted when their cell set changes. PRECAST can be selected by its tool
#' name or by `"PRECAST"` when that method identifies a unique stored bundle.
#'
#' New BANKSY fits and detailed SmoothClust, CHOIR, RareQ, and PRECAST results
#' retain a private metadata identity column, named by the result index.
#' Preserve this column to retrieve results after `subset()` or `RenameCells()`;
#' names are mapped by saved identities, never by order or sample prefixes.
#' Missing, duplicate, or unknown identities cause an informative error.
#' SmoothClust and CHOIR's explicit `cell` columns also reflect the current IDs.
#' Legacy saved assignments support exact-ID subsets. Contradictory metadata
#' causes an error, but renames cannot always be detected without provenance;
#' rerun the method to establish identity tracking. Legacy PRECAST and
#' metadata-only results use their current metadata assignments.
#'
#' Storage-disabled methods other than BANKSY retain their existing storage
#' behavior. Custom metadata columns require a stored locator; the accessor
#' does not guess a column from its values or prefix.
#'
#' @md
#' @param object A `Seurat` object returned by a spatial method.
#' @param method Name of the method or tool entry, such as `"BANKSY"`.
#' @param sample Optional BANKSY sample label, or a sample label for another
#'   method that stored per-sample results.
#'
#' @return A list with `clusters`, `parameters`, and `summary`. Values that were
#'   not stored by the method are returned as `NULL`; BANKSY and PRECAST summaries are
#'   recomputed from the returned cluster assignments, including empty results.
#'
#' @export
#' @examples
#' data(visium_human_pancreas_sub)
#' \dontrun{
#' spatial <- RunBANKSY(
#'   visium_human_pancreas_sub,
#'   layer = "counts",
#'   features = rownames(visium_human_pancreas_sub)[1:200],
#'   verbose = FALSE
#' )
#' GetSpatialResult(spatial, "BANKSY")
#' }
GetSpatialResult <- function(object, method, sample = NULL) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    log_message("{.arg method} must be one non-empty method name", message_type = "error")
  }
  if (!is.null(sample) && (
    !is.character(sample) || length(sample) != 1L || is.na(sample) || !nzchar(sample)
  )) {
    log_message("{.arg sample} must be one non-empty sample label", message_type = "error")
  }

  bundle <- object@tools[[method]]
  # Integration producers nest each method under a named tool. Resolve only
  # an exact, unique method entry; never guess a metadata column by prefix.
  if (is.null(bundle)) {
    matches <- Filter(function(x) is.list(x) && is.list(x$methods) && !is.null(x$methods[[method]]),
      object@tools)
    if (length(matches) > 1L) {
      log_message("Multiple stored results match {.val {method}}; use the exact tool name",
        message_type = "error")
    }
    if (length(matches) == 1L) bundle <- matches[[1L]]$methods[[method]]
  } else if (!is.null(bundle$active_method) && is.list(bundle$methods)) {
    bundle <- bundle$methods[[bundle$active_method]]
  }
  index <- bundle$result_index
  selected <- bundle
  if (!is.null(sample)) {
    if (is.list(bundle$per_sample) && !is.null(names(bundle$per_sample))) {
      if (!sample %in% names(bundle$per_sample)) {
        log_message(
          "No {.val {method}} results were stored for sample {.val {sample}}",
          message_type = "error"
        )
      }
      selected <- bundle$per_sample[[sample]]
    } else if (is.null(index$samples)) {
      log_message(
        "No exact sample membership was stored for {.val {method}}; rerun the method with {.arg sample.by}",
        message_type = "error"
      )
    } else if (!sample %in% index$samples) {
      log_message(
        "No {.val {method}} results were stored for sample {.val {sample}}",
        message_type = "error"
      )
    }
  }

  clusters <- selected$clusters %||% NULL
  parameters <- selected$parameters %||% NULL
  cluster_colname <- index$cluster_colname %||%
    parameters$cluster_colname %||% paste0(method, "_cluster")
  producer_method <- index$method %||% parameters$method %||% method
  legacy_integration_metadata <- identical(parameters$method, "PRECAST") &&
    is.null(index$cell_id_colname) && is.null(index$object_cells) &&
    cluster_colname %in% colnames(object@meta.data)
  if (is.null(clusters) && !is.null(selected$domains) && !legacy_integration_metadata) {
    domain_cells <- spatial_result_cell_ids(selected$domains)
    clusters <- data.frame(
      stats::setNames(list(as.character(selected$domains)), cluster_colname),
      row.names = domain_cells, stringsAsFactors = FALSE
    )
  }
  cell_map <- if (legacy_integration_metadata) {
    stats::setNames(colnames(object), colnames(object))
  } else {
    spatial_result_cell_map(object, bundle, cluster_colname, method)
  }
  original_cells <- spatial_result_cell_ids(clusters) %||% index$cells %||% bundle$cells

  if (is.null(clusters)) {
    if (!cluster_colname %in% colnames(object@meta.data)) {
      log_message(
        "No cluster results for {.val {method}} were found in stored tools or metadata",
        message_type = "error"
      )
    }
    cells <- names(cell_map)
    if (!is.null(sample)) {
      # Membership comes from the analyzed cells, never a cluster-label prefix.
      fitted_cells <- names(index$samples)[index$samples == sample]
      cells <- cells[cell_map[cells] %in% fitted_cells]
    } else if (!is.null(index$cells)) {
      cells <- cells[cell_map[cells] %in% index$cells]
    }
    values <- as.character(object@meta.data[cells, cluster_colname, drop = TRUE])
    if (legacy_integration_metadata) {
      assigned <- !is.na(values) & nzchar(values)
      cells <- cells[assigned]
      values <- values[assigned]
    }
    if (!is.null(sample) && identical(producer_method, "BANKSY")) {
      values <- substring(values, nchar(sample) + 2L)
    }
    clusters <- data.frame(
      stats::setNames(list(values), cluster_colname),
      row.names = cells,
      stringsAsFactors = FALSE
    )
  } else {
    cells <- names(cell_map)[cell_map %in% original_cells]
    fitted_cells <- unname(cell_map[cells])
    if (is.data.frame(clusters) || is.matrix(clusters)) {
      clusters <- clusters[match(fitted_cells, original_cells), , drop = FALSE]
      rownames(clusters) <- cells
      if (spatial_result_has_cell_column(clusters)) clusters$cell <- cells
    } else {
      clusters <- clusters[fitted_cells]
      names(clusters) <- cells
    }
  }

  summary <- selected$summary %||% NULL
  if (identical(producer_method, "BANKSY")) {
    labels <- spatial_result_cluster_labels(clusters)
    if (!is.null(labels)) {
      summary <- list(
        n_spots = sum(!is.na(labels) & nzchar(labels)),
        domains = spatial_domain_summary(labels)
      )
    } else {
      summary <- NULL
    }
  } else if (identical(producer_method, "PRECAST")) {
    # Fit parameters/native results remain untouched in tools. Only descriptive
    # counts that can be derived from the returned cells belong in this view.
    summary <- list(n_cells = length(cells),
      domains = spatial_domain_summary(spatial_result_cluster_labels(clusters)))
    sample_col <- parameters$sample.by
    if (!is.null(sample_col) && sample_col %in% colnames(object@meta.data)) {
      counts <- table(as.character(object@meta.data[cells, sample_col]))
      samples <- data.frame(sample = as.character(names(counts)),
        count = as.integer(counts), stringsAsFactors = FALSE)
      summary$samples <- samples
    }
  } else if (!is.null(original_cells) &&
    !setequal(original_cells, unname(cell_map[cells]))) {
    # An arbitrary method's fit statistics cannot be recomputed from labels.
    summary <- NULL
  }
  list(clusters = clusters, parameters = parameters, summary = summary)
}

spatial_result_cell_ids <- function(clusters) {
  if (is.null(clusters)) return(NULL)
  ids <- if (spatial_result_has_cell_column(clusters)) {
    cell_ids <- as.character(clusters$cell)
    if (.row_names_info(clusters, type = 1L) > 0L &&
      !identical(rownames(clusters), cell_ids)) {
      log_message("Stored spatial cell columns and row names disagree; rerun the spatial method",
        message_type = "error")
    }
    cell_ids
  } else if (is.data.frame(clusters) || is.matrix(clusters)) {
    rownames(clusters)
  } else if (is.atomic(clusters)) {
    names(clusters)
  } else {
    NULL
  }
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    log_message(
      "Stored spatial clusters lack unique cell identities; rerun the spatial method",
      message_type = "error"
    )
  }
  ids
}

spatial_result_has_cell_column <- function(clusters) {
  is.data.frame(clusters) && all(c("cell", "cluster") %in% colnames(clusters))
}

spatial_result_cell_map <- function(object, bundle, cluster_colname, method) {
  index <- bundle$result_index
  cells <- colnames(object)
  fail <- function() {
    log_message(
      paste0(
        "Cell identities for {.val {method}} results are stale or unverifiable. ",
        "Preserve the result's identity metadata column when subsetting or using RenameCells; ",
        "rerun the method on the current object if identities were removed or changed."
      ),
      message_type = "error"
    )
  }
  if (!is.null(index$cell_id_colname) || !is.null(index$object_cells)) {
    column <- index$cell_id_colname
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !column %in% colnames(object@meta.data)) fail()
    ids <- as.character(object@meta.data[cells, column, drop = TRUE])
    if (is.null(index$object_cells) || anyNA(index$object_cells) ||
      any(!nzchar(index$object_cells)) || anyDuplicated(index$object_cells)) fail()
    for (fitted in list(index$cells, names(index$samples),
      spatial_result_cell_ids(bundle$clusters), names(bundle$domains))) {
      if (!is.null(fitted) && (anyNA(fitted) || any(!nzchar(fitted)) ||
        anyDuplicated(fitted) || any(!fitted %in% index$object_cells))) fail()
    }
    if (length(ids) != length(cells) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids) ||
      any(!ids %in% index$object_cells)) fail()
    return(stats::setNames(ids, cells))
  }

  # Legacy objects can be read by exact cell ID. Metadata-only results have
  # no saved identity provenance and cannot establish a rename history.
  known <- index$cells %||% names(index$samples) %||%
    spatial_result_cell_ids(bundle$clusters) %||% names(bundle$domains)
  if (!is.null(known)) {
    unknown <- setdiff(cells, known)
    if (length(unknown)) {
      if (!cluster_colname %in% colnames(object@meta.data)) fail()
      labels <- as.character(object@meta.data[unknown, cluster_colname, drop = TRUE])
      if (any(!is.na(labels) & nzchar(labels))) fail()
    }
    # Legacy saved labels cannot follow a same-ID-set rename. When metadata
    # contradicts their exact-ID mapping, reject it rather than silently swap.
    saved <- bundle$clusters %||% bundle$domains
    labels <- spatial_result_cluster_labels(saved)
    saved_ids <- spatial_result_cell_ids(saved)
    if (!is.null(labels) && !is.null(saved_ids) &&
      cluster_colname %in% colnames(object@meta.data)) {
      common <- intersect(cells, saved_ids)
      current <- as.character(object@meta.data[common, cluster_colname, drop = TRUE])
      if (!identical(unname(current), unname(labels[match(common, saved_ids)]))) fail()
    }
  }
  stats::setNames(cells, cells)
}

spatial_result_cluster_labels <- function(clusters) {
  if (is.null(clusters)) {
    return(NULL)
  }
  if (is.data.frame(clusters)) {
    if ("cluster" %in% colnames(clusters)) {
      return(as.character(clusters[["cluster"]]))
    }
    if (ncol(clusters) == 1L) {
      return(as.character(clusters[[1L]]))
    }
    return(NULL)
  }
  if (is.atomic(clusters)) {
    return(as.character(clusters))
  }
  NULL
}

# A metadata value follows its cell through Seurat subset/RenameCells. The tool
# index alone does not, so retain both the locator and the original identities.
spatial_record_cell_identity <- function(srt, tool_name, previous_index = NULL,
                                         protected_columns = NULL) {
  index <- srt@tools[[tool_name]]$result_index
  column <- previous_index$cell_id_colname
  other_columns <- spatial_other_identity_columns(srt, tool_name)
  # Reuse only our own tracked field, never a result/user field or another fit's.
  reusable <- length(column) == 1L && !is.na(column) &&
    !is.null(previous_index$object_cells) &&
    column %in% colnames(srt@meta.data) &&
    !column %in% c(index$cluster_colname, protected_columns, other_columns)
  if (reusable) {
    previous_ids <- as.character(srt@meta.data[[column]])
    reusable <- !anyNA(previous_ids) && !anyDuplicated(previous_ids) &&
      all(previous_ids %in% previous_index$object_cells)
  }
  if (!reusable) {
    base <- paste0(".scop_", make.names(tool_name), "_cell_id")
    column <- tail(make.unique(c(colnames(srt@meta.data), base)), 1L)
  }
  srt@meta.data[[column]] <- rownames(srt@meta.data)
  index$cell_id_colname <- column
  index$object_cells <- colnames(srt)
  srt@tools[[tool_name]]$result_index <- index
  srt
}

spatial_other_identity_columns <- function(srt, tool_name) {
  unlist(lapply(srt@tools[setdiff(names(srt@tools), tool_name)], function(tool) {
    if (is.list(tool) && is.list(tool$result_index)) tool$result_index$cell_id_colname
  }), use.names = FALSE)
}

spatial_check_identity_columns <- function(srt, tool_name, columns) {
  if (any(columns %in% spatial_other_identity_columns(srt, tool_name))) {
    log_message("An output column is an identity metadata column for another stored result; choose a different output column",
      message_type = "error")
  }
}
