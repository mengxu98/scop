#' @title Annotate single cells using SingleR
#'
#' @md
#' @inheritParams RunKNNPredict
#' @inheritParams SingleR::SingleR
#' @inheritParams SingleR::trainSingleR
#' @inheritParams thisutils::log_message
#' @inheritParams thisutils::parallelize_fun
#' @param genes `"genes"` parameter in [SingleR::SingleR] function.
#' @param de.method `"de.method"` parameter in [SingleR::SingleR] function.
#' @param quantile "quantile" parameter in [SingleR::SingleR] function.
#' @param fine.tune `"fine.tune"` parameter in [SingleR::SingleR] function.
#' @param tune.thresh `"tune.thresh"` parameter in [SingleR::SingleR] function.
#' @param prune `"prune"` parameter in [SingleR::SingleR] function.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- CheckDataMerge(
#'   panc8_sub,
#'   batch = "tech"
#' )[["srt_merge"]]
#'
#' # Annotation
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSingleR(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "Standardpca_SNN_res.0.6",
#'   ref_group = "celltype"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = c("singler_annotation", "CellType")
#' )
#'
#' pancreas_sub <- RunSingleR(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = NULL,
#'   ref_group = "celltype"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = c("singler_annotation", "CellType")
#' )
RunSingleR <- function(
    srt_query,
    srt_ref,
    query_group = NULL,
    ref_group = NULL,
    query_assay = "RNA",
    ref_assay = "RNA",
    genes = "de",
    de.method = "wilcox",
    sd.thresh = 1,
    de.n = NULL,
    aggr.ref = FALSE,
    aggr.args = list(),
    quantile = 0.8,
    fine.tune = TRUE,
    tune.thresh = 0.05,
    prune = TRUE,
    cores = 1,
    verbose = TRUE) {
  log_message(
    "Start {.pkg SingleR} annotation",
    verbose = verbose
  )
  check_r(c("SingleR", "scrapper"))
  if (is.null(ref_group)) {
    log_message(
      "'ref_group' must be provided",
      message_type = "error",
      verbose = verbose
    )
  }
  if (length(ref_group) == ncol(srt_ref)) {
    srt_ref[["ref_group"]] <- ref_group
  } else if (length(ref_group) == 1) {
    if (!ref_group %in% colnames(srt_ref@meta.data)) {
      log_message(
        "ref_group must be one of the column names in the meta.data",
        message_type = "error"
      )
    } else {
      srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
    }
  } else {
    log_message(
      "Length of ref_group must be one or length of srt_ref.",
      message_type = "error"
    )
  }
  ref_group <- "ref_group"

  if (!is.null(query_group)) {
    if (length(query_group) == ncol(srt_query)) {
      srt_query[["query_group"]] <- query_group
    } else if (length(query_group) == 1) {
      if (!query_group %in% colnames(srt_query@meta.data)) {
        log_message(
          "query_group must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_query[["query_group"]] <- srt_query[[query_group]]
      }
    } else {
      log_message(
        "Length of query_group must be one or length of srt_query.",
        message_type = "error"
      )
    }
    query_group <- "query_group"
    method <- "SingleRCluster"
  } else {
    method <- "SingleRCell"
  }

  status_query <- CheckDataType(
    object = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected srt_query data type: {.val {status_query}}")
  status_ref <- CheckDataType(
    object = GetAssayData5(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  log_message("Detected srt_ref data type: {.val {status_ref}}")
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    log_message(
      "Data type is unknown or different between query and ref",
      message_type = "warning",
      verbose = verbose
    )
  }

  assays_query <- list(
    counts = GetAssayData5(
      object = srt_query,
      assay = query_assay,
      layer = "counts"
    ),
    logcounts = GetAssayData5(
      object = srt_query,
      assay = query_assay,
      layer = "data"
    )
  )
  sce_query <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_query
  )
  metadata_query <- srt_query[[]]
  SummarizedExperiment::colData(x = sce_query) <- S4Vectors::DataFrame(
    metadata_query
  )

  assays_ref <- list(
    counts = GetAssayData5(
      object = srt_ref,
      assay = ref_assay,
      layer = "counts"
    ),
    logcounts = GetAssayData5(
      object = srt_ref,
      assay = ref_assay,
      layer = "data"
    )
  )
  sce_ref <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_ref
  )
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(
    metadata_ref
  )

  log_message(
    "Perform {.val {method}}",
    verbose = verbose
  )
  if (method == "SingleRCluster") {
    cluster_results <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = factor(
        SummarizedExperiment::colData(sce_ref)[[ref_group]]
      ),
      clusters = factor(
        SummarizedExperiment::colData(sce_query)[[query_group]]
      ),
      de.method = de.method,
      genes = genes,
      sd.thresh = sd.thresh,
      de.n = de.n,
      aggr.ref = aggr.ref,
      aggr.args = aggr.args,
      quantile = quantile,
      fine.tune = fine.tune,
      tune.thresh = tune.thresh,
      prune = prune,
      num.threads = cores
    )
    names(cluster_results$labels) <- levels(
      factor(
        SummarizedExperiment::colData(sce_query)[[query_group]]
      )
    )
    rownames(cluster_results$scores) <- levels(
      factor(
        SummarizedExperiment::colData(sce_query)[[query_group]]
      )
    )
    if (isTRUE(prune)) {
      names(cluster_results$pruned.labels) <- levels(
        factor(
          SummarizedExperiment::colData(sce_query)[[query_group]]
        )
      )
    }

    if (isTRUE(prune)) {
      cluster_labels <- cluster_results$pruned.labels
    } else {
      cluster_labels <- cluster_results$labels
    }

    sce_cluster_ids <- SummarizedExperiment::colData(sce_query)[[query_group]]
    cluster_to_cells <- split(seq_along(sce_cluster_ids), sce_cluster_ids)

    cell_annotations <- character(length(sce_cluster_ids))
    for (cluster_id in names(cluster_to_cells)) {
      if (cluster_id %in% names(cluster_labels)) {
        cell_annotations[cluster_to_cells[[cluster_id]]] <- cluster_labels[cluster_id]
      }
    }

    srt_query$singler_annotation <- cell_annotations
    cell_scores <- numeric(length(sce_cluster_ids))
    for (cluster_id in names(cluster_to_cells)) {
      if (cluster_id %in% names(cluster_labels)) {
        annotation <- cluster_labels[cluster_id]
        if (!is.na(annotation) && annotation %in% colnames(cluster_results$scores)) {
          score <- cluster_results$scores[cluster_id, annotation]
          cell_scores[cluster_to_cells[[cluster_id]]] <- score
        }
      }
    }
    srt_query$singler_score <- cell_scores
  } else if (method == "SingleRCell") {
    cell_results <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = factor(SummarizedExperiment::colData(sce_ref)[[ref_group]]),
      clusters = NULL,
      de.method = de.method,
      genes = genes,
      sd.thresh = sd.thresh,
      de.n = de.n,
      aggr.ref = aggr.ref,
      aggr.args = aggr.args,
      quantile = quantile,
      fine.tune = fine.tune,
      tune.thresh = tune.thresh,
      prune = prune,
      num.threads = cores
    )
    if (isTRUE(prune)) {
      srt_query$singler_annotation <- cell_results$pruned.labels
    } else {
      srt_query$singler_annotation <- cell_results$labels
    }

    cell_scores <- numeric(ncol(srt_query))
    for (i in seq_len(ncol(srt_query))) {
      if (isTRUE(prune)) {
        annotation <- cell_results$pruned.labels[i]
      } else {
        annotation <- cell_results$labels[i]
      }

      if (!is.na(annotation) && annotation %in% colnames(cell_results$scores)) {
        cell_scores[i] <- cell_results$scores[i, annotation]
      } else {
        cell_scores[i] <- NA
      }
    }
    srt_query$singler_score <- cell_scores
  }

  log_message(
    "{.pkg SingleR} annotation completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt_query)
}
