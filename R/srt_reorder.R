#' @title Reorder idents by the gene expression
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @inheritParams CellCorHeatmap
#' @param reorder_by Reorder groups instead of idents.
#' @param log Whether [log1p] transformation needs to be applied.
#' Default is `TRUE`.
#' @param distance_metric Metric to compute distance.
#' Default is `"euclidean"`.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- srt_reorder(
#'   pancreas_sub,
#'   reorder_by = "SubCellType",
#'   layer = "data"
#' )
srt_reorder <- function(
    srt,
    features = NULL,
    reorder_by = NULL,
    layer = "data",
    assay = NULL,
    log = TRUE,
    distance_metric = "euclidean",
    verbose = TRUE) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (is.null(features)) {
    srt <- Seurat::FindVariableFeatures(
      srt,
      assay = SeuratObject::DefaultAssay(srt),
      verbose = FALSE
    )
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  if (is.null(reorder_by)) {
    srt$ident <- SeuratObject::Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  if (length(unique(srt[[reorder_by, drop = TRUE]])) == 1) {
    log_message(
      "Only one cluster found",
      message_type = "warning",
      verbose = verbose
    )
    return(srt)
  }
  simil_methods <- c(
    "cosine",
    "correlation",
    "jaccard",
    "ejaccard",
    "dice",
    "edice",
    "hamman",
    "simple matching",
    "faith"
  )
  dist_methods <- c(
    "euclidean",
    "chisquared",
    "kullback",
    "manhattan",
    "maximum",
    "canberra",
    "minkowski",
    "hamming"
  )
  if (!distance_metric %in% c(simil_methods, dist_methods, "pearson", "spearman")) {
    log_message(
      "{.pkg {distance_metric}} method is invalid",
      message_type = "error"
    )
  }

  data_avg <- Seurat::AggregateExpression(
    object = srt,
    features = features,
    assays = assay,
    group.by = "ident",
    verbose = FALSE
  )[[1]][features, , drop = FALSE]

  if (isTRUE(log)) {
    data_avg <- log1p(data_avg)
  }
  mat <- Matrix::t(data_avg[features, , drop = FALSE])
  if (!inherits(mat, "dgCMatrix")) {
    mat <- SeuratObject::as.sparse(
      mat[seq_len(nrow(mat)), , drop = FALSE]
    )
  }

  if (distance_metric %in% c(simil_methods, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- Matrix::t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 - proxyC::simil(
      SeuratObject::as.sparse(
        mat[seq_len(nrow(mat)), , drop = FALSE]
      ),
      method = distance_metric
    )
  } else if (distance_metric %in% dist_methods) {
    d <- proxyC::dist(
      SeuratObject::as.sparse(
        mat[seq_len(nrow(mat)), , drop = FALSE]
      ),
      method = distance_metric
    )
  }
  data_dist <- stats::as.dist(d)
  hc <- stats::hclust(d = data_dist)
  dd <- stats::as.dendrogram(hc)
  dd_ordered <- stats::reorder(
    dd,
    wts = Matrix::colMeans(data_avg[features, , drop = FALSE]),
    agglo.FUN = mean
  )
  ident_new <- unname(stats::setNames(
    object = seq_along(labels(dd_ordered)),
    nm = labels(dd_ordered)
  )[as.character(srt$ident)])
  ident_new <- factor(ident_new, levels = seq_along(labels(dd_ordered)))
  Idents(srt) <- srt$ident <- ident_new
  return(srt)
}
