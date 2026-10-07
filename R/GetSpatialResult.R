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
#' than the original fit's full cell set. BANKSY counts are recomputed from the
#' returned assignments; parameters still describe the original fit. Other
#' methods' stored summaries are omitted when their cell set changes.
#'
#' New BANKSY fits retain a private metadata identity column, named by the
#' result index. Preserve this column to retrieve results after `subset()` or
#' `RenameCells()`; names are mapped by these saved identities, never by order
#' or sample prefixes. Missing, duplicate, or unknown identities cause an
#' informative error. Legacy results support exact-ID subsets, but cannot
#' reliably detect renames that permute an unchanged set of IDs. Rerun BANKSY
#' to establish identity tracking for these older objects.
#'
#' @md
#' @param object A `Seurat` object returned by a spatial method.
#' @param method Name of the method or tool entry, such as `"BANKSY"`.
#' @param sample Optional BANKSY sample label, or a sample label for another
#'   method that stored per-sample results.
#'
#' @return A list with `clusters`, `parameters`, and `summary`. Values that were
#'   not stored by the method are returned as `NULL`; BANKSY summaries are
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
  cell_map <- spatial_result_cell_map(object, bundle, cluster_colname, method)
  original_cells <- spatial_result_cell_ids(clusters)

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
    if (!is.null(sample)) {
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
      clusters <- clusters[fitted_cells, , drop = FALSE]
      rownames(clusters) <- cells
    } else {
      clusters <- clusters[fitted_cells]
      names(clusters) <- cells
    }
  }

  summary <- selected$summary %||% NULL
  if (identical(method, "BANKSY") || identical(index$method, "BANKSY")) {
    labels <- spatial_result_cluster_labels(clusters)
    if (!is.null(labels)) {
      summary <- list(
        n_spots = sum(!is.na(labels) & nzchar(labels)),
        domains = spatial_domain_summary(labels)
      )
    } else {
      summary <- NULL
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
  ids <- if (is.data.frame(clusters) || is.matrix(clusters)) {
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

spatial_result_cell_map <- function(object, bundle, cluster_colname, method) {
  index <- bundle$result_index
  cells <- colnames(object)
  fail <- function() {
    log_message(
      paste0(
        "Cell identities for {.val {method}} results are stale or unverifiable. ",
        "Preserve the BANKSY identity metadata column when subsetting or using RenameCells; ",
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
      spatial_result_cell_ids(bundle$clusters))) {
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
    spatial_result_cell_ids(bundle$clusters)
  if (!is.null(known)) {
    unknown <- setdiff(cells, known)
    if (length(unknown)) {
      if (!cluster_colname %in% colnames(object@meta.data)) fail()
      labels <- as.character(object@meta.data[unknown, cluster_colname, drop = TRUE])
      if (any(!is.na(labels) & nzchar(labels))) fail()
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
