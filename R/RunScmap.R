#' @title Annotate single cells using scmap.
#'
#' @md
#' @inheritParams RunKNNPredict
#' @param method The method to be used for scmap analysis.
#' Can be any of `"scmapCluster"` or `"scmapCell"`.
#' Default is `"scmapCluster"`.
#' @param nfeatures The number of top features to be selected.
#' Default is `500`.
#' @param threshold The threshold value on similarity to determine if a cell is assigned to a cluster.
#' This should be a value between `0` and `1`.
#' Default is `0.5`.
#' @param k Number of clusters per group for k-means clustering when `method` is `"scmapCell"`.
#' Default is `10`.
#'
#' @seealso
#' [RunKNNPredict], [RunKNNMap]
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- standard_scop(panc8_sub)
#'
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
#' data(pancreas_sub)
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
  check_r("scmap", verbose = FALSE)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        log_message(
          "{.arg ref_group} must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      log_message(
        "Length of {.arg ref_group} must be one or length of {.arg srt_ref}",
        message_type = "error"
      )
    }
    ref_group <- "ref_group"
  } else {
    log_message(
      "{.arg ref_group} must be provided",
      message_type = "error"
    )
  }

  status_query <- CheckDataType(
    object = GetAssayData5(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  log_message("Detected {.arg srt_query} data type: {.val {status_query}}")
  status_ref <- CheckDataType(
    object = GetAssayData5(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  log_message("Detected {.arg srt_ref} data type: {.val {status_ref}}")
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
  sce_query <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_query
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
  sce_ref <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_ref
  )
  SummarizedExperiment::rowData(sce_ref)[["feature_symbol"]] <- rownames(
    sce_ref
  )
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(
    metadata_ref
  )

  log_message("Perform selectFeatures")
  sce_ref <- scmap::selectFeatures(
    sce_ref,
    n_features = nfeatures,
    suppress_plot = TRUE
  )
  features <- rownames(sce_ref)[
    SummarizedExperiment::rowData(sce_ref)$scmap_features
  ]

  if (method == "scmapCluster") {
    log_message("Perform indexCluster")
    sce_ref <- scmap::indexCluster(sce_ref, cluster_col = ref_group)
    log_message("Perform scmapCluster")
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
    log_message("Perform indexCell")
    sce_ref <- scmap::indexCell(
      sce_ref,
      M = ifelse(nfeatures <= 1000, nfeatures / 10, 100),
      k = sqrt(ncol(sce_ref))
    )
    log_message("Perform scmapCell")
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
    log_message("Perform scmapCell2Cluster")
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
