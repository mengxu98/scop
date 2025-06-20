#' Annotate single cells using scmap.
#'
#' @inheritParams RunKNNPredict
#' @param method The method to be used for scmap analysis.
#' Can be any of "scmapCluster" or "scmapCell".
#' The default value is "scmapCluster".
#' @param nfeatures The number of top features to be selected.
#' The default value is 500.
#' @param threshold The threshold value on similarity to determine if a cell is assigned to a cluster.
#' This should be a value between 0 and 1. The default value is 0.5.
#' @param k Number of clusters per group for k-means clustering when method is "scmapCell".
#'
#' @export
#'
#' @examples
#' data("panc8_sub")
#'
#' genenames <- make.unique(
#'   capitalize(rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- check_srt_merge(
#'   panc8_sub,
#'   batch = "tech"
#' )[["srt_merge"]]
#'
#' data("pancreas_sub")
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunScmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   method = "scmapCluster"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "scmap_annotation"
#' )
#'
#' pancreas_sub <- RunScmap(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   method = "scmapCell"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "scmap_annotation"
#' )
RunScmap <- function(
    srt_query,
    srt_ref,
    ref_group = NULL,
    query_assay = "RNA",
    ref_assay = "RNA",
    method = "scmapCluster",
    nfeatures = 500,
    threshold = 0.5,
    k = 10) {
  check_r("scmap")
  if (!is.null(ref_group)) {
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
  } else {
    log_message(
      "'ref_group' must be provided.",
      message_type = "error"
    )
  }

  status_query <- check_data_type(
    data = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected srt_query data type: ", status_query)
  status_ref <- check_data_type(
    data = GetAssayData5(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  log_message("Detected srt_ref data type: ", status_ref)
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    log_message(
      "Data type is unknown or different between query and ref.",
      message_type = "warning"
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
  sce_query <- methods::as(
    SummarizedExperiment::SummarizedExperiment(
      assays = assays_query
    ),
    Class = "SingleCellExperiment"
  )
  SummarizedExperiment::rowData(sce_query)[["feature_symbol"]] <- rownames(
    sce_query
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
  sce_ref <- methods::as(
    SummarizedExperiment::SummarizedExperiment(
      assays = assays_ref
    ),
    Class = "SingleCellExperiment"
  )
  SummarizedExperiment::rowData(sce_ref)[["feature_symbol"]] <- rownames(
    sce_ref
  )
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(
    metadata_ref
  )

  log_message("Perform selectFeatures on the data...")
  sce_ref <- scmap::selectFeatures(
    sce_ref,
    n_features = nfeatures,
    suppress_plot = TRUE
  )
  features <- rownames(sce_ref)[
    SummarizedExperiment::rowData(sce_ref)$scmap_features
  ]

  if (method == "scmapCluster") {
    log_message("Perform indexCluster on the data...")
    sce_ref <- scmap::indexCluster(sce_ref, cluster_col = ref_group)
    log_message("Perform scmapCluster on the data...")
    scmapCluster_results <- scmap::scmapCluster(
      projection = sce_query,
      index_list = list(S4Vectors::metadata(sce_ref)$scmap_cluster_index),
      threshold = threshold
    )
    if (!"scmap_cluster_labs" %in% names(scmapCluster_results)) {
      log_message(
        "scmap failed to run. Please check the warning message.",
        message_type = "error"
      )
    }
    srt_query@tools[["scmapCluster_results"]] <- scmapCluster_results
    srt_query$scmap_annotation <- scmapCluster_results$scmap_cluster_labs[, 1]
    srt_query$scmap_score <- scmapCluster_results$scmap_cluster_siml[, 1]
  } else if (method == "scmapCell") {
    log_message("Perform indexCell on the data...")
    sce_ref <- scmap::indexCell(
      sce_ref,
      M = ifelse(nfeatures <= 1000, nfeatures / 10, 100),
      k = sqrt(ncol(sce_ref))
    )
    log_message("Perform scmapCell on the data...")
    scmapCell_results <- scmap::scmapCell(
      projection = sce_query,
      index_list = list(S4Vectors::metadata(sce_ref)$scmap_cell_index),
      w = k
    )
    if (!"cells" %in% names(scmapCell_results[[1]])) {
      log_message(
        "scmap failed to run. Please check the warning message.",
        message_type = "error"
      )
    }
    srt_query@tools[["scmapCell_results"]] <- scmapCell_results[[1]]
    log_message("Perform scmapCell2Cluster on the data...")
    scmapCell_clusters <- scmap::scmapCell2Cluster(
      scmapCell_results = scmapCell_results,
      cluster_list = list(
        as.character(SummarizedExperiment::colData(sce_ref)[[ref_group]])
      ),
      w = k,
      threshold = threshold
    )
    srt_query@tools[["scmapCell_results"]][[
      "scmapCell2Cluster"
    ]] <- scmapCell_clusters
    srt_query$scmap_annotation <- scmapCell_clusters$scmap_cluster_labs[, 1]
    srt_query$scmap_score <- scmapCell_clusters$scmap_cluster_siml[, 1]
  }
  return(srt_query)
}
