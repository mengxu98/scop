
#' Annotate single cells using SingleR
#'
#' @inheritParams RunKNNPredict
#' @inheritParams SingleR::SingleR
#' @inheritParams SingleR::trainSingleR
#' @param genes "genes" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param de.method "de.method" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param quantile "quantile" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param fine.tune "fine.tune" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param tune.thresh "tune.thresh" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param prune "prune" parameter in \code{\link[SingleR]{SingleR}} function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   capitalize(
#'     rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- check_srt_merge(
#'   panc8_sub,
#'   batch = "tech"
#' )[["srt_merge"]]
#'
#' # Annotation
#' data("pancreas_sub")
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSingleR( # bug
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "Standardclusters",
#'   ref_group = "celltype",
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "singler_annotation"
#' )
#'
#' pancreas_sub <- RunSingleR( # bug
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = NULL,
#'   ref_group = "celltype"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "singler_annotation"
#' )
#' }
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
    BPPARAM = BiocParallel::bpparam()) {
  check_r("SingleR")
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
    ref_group <- "ref_group"
  } else {
    stop("'ref_group' must be provided.")
  }
  if (!is.null(query_group)) {
    if (length(query_group) == ncol(srt_query)) {
      srt_query[["query_group"]] <- query_group
    } else if (length(query_group) == 1) {
      if (!query_group %in% colnames(srt_query@meta.data)) {
        stop("query_group must be one of the column names in the meta.data")
      } else {
        srt_query[["query_group"]] <- srt_query[[query_group]]
      }
    } else {
      stop("Length of query_group must be one or length of srt_query.")
    }
    query_group <- "query_group"
    method <- "SingleRCluster"
  } else {
    method <- "SingleRCell"
  }

  status_query <- check_data_type(
    data = SeuratObject::GetAssayData(
      srt_query,
      layer = "data",
      assay = query_assay
    )
  )
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_data_type(
    data = SeuratObject::GetAssayData(
      srt_ref,
      layer = "data",
      assay = ref_assay
    )
  )
  message("Detected srt_ref data type: ", status_ref)
  if (
    status_ref != status_query ||
      any(status_query == "unknown", status_ref == "unknown")
  ) {
    warning(
      "Data type is unknown or different between query and ref.",
      immediate. = TRUE
    )
  }

  assays_query <- list(
    counts = SeuratObject::GetAssayData(
      object = srt_query,
      assay = query_assay,
      layer = "counts"
    ),
    logcounts = SeuratObject::GetAssayData(
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
  metadata_query <- srt_query[[]]
  SummarizedExperiment::colData(x = sce_query) <- S4Vectors::DataFrame(
    metadata_query
  )

  assays_ref <- list(
    counts = SeuratObject::GetAssayData(
      object = srt_ref,
      assay = ref_assay,
      layer = "counts"
    ),
    logcounts = SeuratObject::GetAssayData(
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
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(
    metadata_ref
  )

  if (method == "SingleRCluster") {
    message("Perform ", method, " on the data...")
    SingleRCluster_results <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = factor(SummarizedExperiment::colData(sce_ref)[[ref_group]]),
      clusters = factor(SummarizedExperiment::colData(sce_query)[[
        query_group
      ]]),
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
      BPPARAM = BPPARAM
    )
    names(
      SingleRCluster_results$labels
    ) <- levels(factor(SummarizedExperiment::colData(sce_query)[[query_group]]))
    rownames(
      SingleRCluster_results$scores
    ) <- levels(factor(SummarizedExperiment::colData(sce_query)[[query_group]]))
    if (isTRUE(prune)) {
      names(
        SingleRCluster_results$pruned.labels
      ) <- levels(factor(SummarizedExperiment::colData(sce_query)[[
        query_group
      ]]))
    }
    srt_query$singler_annotation <- if (isTRUE(prune)) {
      SingleRCluster_results$pruned.labels[as.character(srt_query$query_group)]
    } else {
      SingleRCluster_results$labels[as.character(srt_query$query_group)]
    }
    srt_query$singler_score <- sapply(
      as.character(unique(srt_query$query_group)),
      FUN = function(x) {
        if (isTRUE(prune)) {
          y <- SingleRCluster_results$pruned.labels[x]
        } else {
          y <- SingleRCluster_results$labels[x]
        }
        if (is.na(y)) {
          out <- NA
        } else {
          out <- SingleRCluster_results$scores[x, y]
        }
        return(out)
      }
    )[as.character(srt_query$query_group)]
  } else if (method == "SingleRCell") {
    message("Perform ", method, " on the data...")
    SingleRCell_results <- SingleR::SingleR(
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
      BPPARAM = BPPARAM
    )
    srt_query$singler_annotation <- if (isTRUE(prune)) {
      SingleRCell_results$pruned.labels
    } else {
      SingleRCell_results$labels
    }
    srt_query$singler_score <- sapply(
      seq_len(ncol(srt_query)),
      FUN = function(x) {
        if (isTRUE(prune)) {
          y <- SingleRCell_results$pruned.labels[x]
        } else {
          y <- SingleRCell_results$labels[x]
        }
        if (is.na(y)) {
          out <- NA
        } else {
          out <- SingleRCell_results$scores[x, y]
        }
        return(out)
      }
    )
  }
  return(srt_query)
}
