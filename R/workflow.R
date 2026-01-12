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

#' @title Append a Seurat object to another
#'
#' @inheritParams thisutils::log_message
#' @param srt_raw A Seurat object to be appended.
#' @param srt_append New Seurat object to append.
#' @param slots slots names.
#' @param pattern A character string containing a regular expression.
#' All data with matching names will be considered for appending.
#' @param overwrite Whether to overwrite.
#'
#' @export
srt_append <- function(
    srt_raw,
    srt_append,
    slots = methods::slotNames(srt_append),
    pattern = NULL,
    overwrite = FALSE,
    verbose = TRUE) {
  if (!inherits(srt_raw, "Seurat") || !inherits(srt_append, "Seurat")) {
    log_message(
      "{.arg srt_raw} or {.arg srt_append} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  pattern <- pattern %||% ""
  for (slot_nm in methods::slotNames(srt_append)) {
    if (!slot_nm %in% slots) {
      log_message(
        "Slot {.val {slot_nm}} is not appended",
        verbose = verbose
      )
      next
    }
    if (identical(slot_nm, "active.ident") && isTRUE(overwrite)) {
      methods::slot(srt_raw, name = "active.ident") <- methods::slot(
        srt_append,
        name = "active.ident"
      )
      next
    }
    slot_data <- methods::slot(srt_append, name = slot_nm)
    for (info in names(slot_data)) {
      if (is.null(info)) {
        if (length(slot_data) > 0 && isTRUE(overwrite)) {
          methods::slot(srt_raw, name = slot_nm) <- slot_data
        }
        next
      }
      if (!grepl(pattern = pattern, x = info)) {
        log_message(
          "{.val {info}} in slot {.val {slot_nm}} is not appended",
          verbose = verbose
        )
        next
      }
      if (!info %in% names(methods::slot(srt_raw, name = slot_nm)) || isTRUE(overwrite)) {
        if (slot_nm %in% c("assays", "graphs", "neighbors", "reductions", "images")) {
          if (identical(slot_nm, "graphs")) {
            srt_raw@graphs[[info]] <- srt_append[[info]]
          } else if (identical(slot_nm, "assays")) {
            if (!info %in% SeuratObject::Assays(srt_raw)) {
              srt_raw[[info]] <- srt_append[[info]]
            } else {
              srt_raw <- SeuratObject::SetAssayData(
                object = srt_raw,
                layer = "counts",
                assay = info,
                new.data = GetAssayData5(
                  srt_append,
                  assay = info,
                  layer = "counts"
                )
              )
              srt_raw <- SeuratObject::SetAssayData(
                object = srt_raw,
                layer = "data",
                assay = info,
                new.data = GetAssayData5(
                  srt_append,
                  assay = info,
                  layer = "data"
                )
              )
              if (inherits(Seurat::GetAssay(srt_raw, assay = info), "Assay5")) {
                srt_raw@assays[[info]]@key <- srt_append@assays[[info]]@key
                srt_raw@assays[[info]]@meta.data$var.features <- srt_append@assays[[info]]@meta.data$var.features
              } else {
                srt_raw[[info]]@key <- srt_append[[info]]@key
                srt_raw[[info]]@var.features <- srt_append[[info]]@var.features
              }
              srt_raw[[info]]@misc <- srt_append[[info]]@misc
              meta_features <- cbind(
                GetFeaturesData(srt_raw, assay = info),
                GetFeaturesData(srt_append, assay = info)[
                  rownames(GetFeaturesData(srt_raw, assay = info)),
                  setdiff(
                    colnames(GetFeaturesData(srt_append, assay = info)),
                    colnames(GetFeaturesData(srt_raw, assay = info))
                  )
                ]
              )
              srt_raw <- AddFeaturesData(
                srt_raw,
                features = meta_features,
                assay = info
              )
            }
          } else {
            srt_raw[[info]] <- srt_append[[info]]
          }
        } else if (identical(slot_nm, "meta.data")) {
          srt_raw@meta.data[, info] <- NULL
          srt_raw@meta.data[[info]] <- srt_append@meta.data[
            colnames(srt_raw),
            info
          ]
        } else {
          methods::slot(srt_raw, name = slot_nm)[[info]] <- methods::slot(
            srt_append,
            name = slot_nm
          )[[info]]
        }
      }
    }
  }
  return(srt_raw)
}
