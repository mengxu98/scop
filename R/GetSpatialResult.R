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
#' @md
#' @param object A `Seurat` object returned by a spatial method.
#' @param method Name of the method or tool entry, such as `"BANKSY"`.
#' @param sample Optional BANKSY sample label, or a sample label for another
#'   method that stored per-sample results.
#'
#' @return A list with `clusters`, `parameters`, and `summary`. Values that were
#'   not stored by the method are returned as `NULL`; BANKSY summaries are
#'   derived from the returned cluster assignments when necessary.
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
  if (is.null(clusters)) {
    cluster_colname <- index$cluster_colname %||%
      parameters$cluster_colname %||% paste0(method, "_cluster")
    if (!cluster_colname %in% colnames(object@meta.data)) {
      log_message(
        "No cluster results for {.val {method}} were found in stored tools or metadata",
        message_type = "error"
      )
    }
    values <- object@meta.data[[cluster_colname]]
    names(values) <- rownames(object@meta.data)
    if (!is.null(sample)) {
      # Membership comes from the analyzed cells, never a cluster-label prefix.
      cells <- names(index$samples)[index$samples == sample]
      cells <- cells[cells %in% names(values)]
      values <- values[cells]
      cell_ids <- names(values)
      values <- substring(as.character(values), nchar(sample) + 2L)
      names(values) <- cell_ids
    } else if (!is.null(index$cells)) {
      cells <- index$cells[index$cells %in% names(values)]
      values <- values[cells]
    }
    clusters <- data.frame(
      stats::setNames(list(as.character(values)), cluster_colname),
      row.names = names(values),
      stringsAsFactors = FALSE
    )
  }

  summary <- selected$summary %||% NULL
  if (is.null(summary) && (identical(method, "BANKSY") || identical(index$method, "BANKSY"))) {
    labels <- spatial_result_cluster_labels(clusters)
    if (!is.null(labels)) {
      summary <- list(
        n_spots = sum(!is.na(labels) & nzchar(labels)),
        domains = spatial_domain_summary(labels)
      )
    }
  }
  list(clusters = clusters, parameters = parameters, summary = summary)
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
