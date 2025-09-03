#' @title Attempt to recover raw counts from the normalized matrix
#'
#' @param srt A Seurat object.
#' @param assay Name of assay to recover counts.
#' @param trans The transformation function to applied when data is presumed to be log-normalized.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#'  However, due to decimal point preservation during normalization,
#'  the calculated nCount is usually a floating point number close to the integer.
#'  The tolerance is its difference from the integer. Default is 0.1
#' @param sf Set the scaling factor manually.
#' @param verbose Whether to show messages.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' raw_counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#'
#' # Normalized the data
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#'
#' # Now replace counts with the log-normalized data matrix
#' data <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "data"
#' )
#' new_pancreas_sub <- SeuratObject::SetAssayData(
#'   object = pancreas_sub,
#'   layer = "counts",
#'   new.data = data,
#'   assay = "RNA"
#' )
#' # Recover the counts and compare with the raw counts matrix
#' pancreas_sub <- RecoverCounts(new_pancreas_sub)
#' new_counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#' identical(raw_counts, new_counts)
RecoverCounts <- function(
    srt,
    assay = NULL,
    trans = c("expm1", "exp", "none"),
    min_count = c(1, 2, 3),
    tolerance = 0.1,
    sf = NULL,
    verbose = TRUE) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  counts <- GetAssayData5(
    srt,
    layer = "counts",
    assay = assay
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- SeuratObject::as.sparse(
      counts[1:nrow(counts), , drop = FALSE]
    )
  }

  status <- CheckDataType(data = counts)
  if (status == "raw_counts") {
    log_message(
      "The data is already raw counts.",
      verbose = verbose
    )
    return(srt)
  }

  if (status == "log_normalized_counts") {
    log_message(
      "The data is presumed to be log-normalized.",
      verbose = verbose
    )
    trans <- match.arg(trans)
    if (trans %in% c("expm1", "exp")) {
      log_message(
        "Perform ", trans, " on the raw data.",
        verbose = verbose
      )
      counts <- do.call(trans, list(counts))
    }
  }
  if (status == "raw_normalized_counts") {
    log_message(
      "The data is presumed to be normalized without log transformation.",
      verbose = verbose
    )
  }
  if (is.null(sf)) {
    sf <- unique(round(Matrix::colSums(counts)))
    log_message(
      "The presumed scale factor: ",
      paste0(utils::head(sf, 10), collapse = ", "),
      verbose = verbose
    )
  }

  if (length(sf) == 1) {
    counts <- counts / sf
    elements <- split(counts@x, rep(1:ncol(counts), diff(counts@p)))
    min_norm <- sapply(elements, min)
    nCount <- NULL
    for (m in min_count) {
      if (is.null(nCount)) {
        presumed_nCount <- m / min_norm
        diff_value <- abs(presumed_nCount - round(presumed_nCount))
        if (max(diff_value, na.rm = TRUE) < tolerance) {
          nCount <- round(presumed_nCount)
        }
      }
    }

    if (is.null(nCount)) {
      log_message(
        "The presumed nCount of some cells is not valid: ",
        paste0(
          utils::head(colnames(counts)[diff_value < tolerance], 10),
          collapse = ","
        ),
        ", ...",
        message_type = "warning"
      )
      return(srt)
    }
    counts@x <- round(counts@x * rep(nCount, diff(counts@p)))
    srt <- SeuratObject::SetAssayData(
      srt,
      new.data = counts,
      assay = assay,
      layer = "counts"
    )

    nCount <- stats::setNames(
      rep(nCount, length.out = ncol(srt)),
      colnames(srt)
    )
    srt[[paste0("nCount_", assay)]] <- nCount
  } else {
    log_message(
      "Scale factor is not unique. No changes to be made.",
      message_type = "warning"
    )
  }

  return(srt)
}

#' @title Reorder idents by the gene expression
#'
#' @param srt A Seurat object.
#' @param features Features used to reorder idents.
#' @param reorder_by Reorder groups instead of idents.
#' @param layer Specific layer to get data from.
#' @param assay Specific assay to get data from.
#' @param log Whether log1p transformation needs to be applied.
#' Default is \code{TRUE}.
#' @param distance_metric Metric to compute distance.
#' Default is "euclidean".
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- srt_reorder(
#'   srt = pancreas_sub,
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
    distance_metric = "euclidean") {
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
      message_type = "warning"
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

  assay_obj <- Seurat::GetAssay(
    srt,
    assay = assay
  )
  if (inherits(assay_obj, "Assay")) {
    log_message(
      "Using {.fn Seurat::AverageExpression} to calculate pseudo-bulk data for {.cls Assay}"
    )
    data_avg <- Seurat::AverageExpression(
      object = srt,
      features = features,
      layer = layer,
      assays = assay,
      group.by = "ident",
      verbose = FALSE
    )[[1]][features, , drop = FALSE]
  } else if (inherits(assay_obj, "Assay5")) {
    log_message(
      "Using {.fn Seurat::AggregateExpression} to calculate pseudo-bulk data for {.cls Assay5}",
      message_type = "warning"
    )
    data_avg <- Seurat::AggregateExpression(
      object = srt,
      features = features,
      assays = assay,
      group.by = "ident",
      verbose = FALSE
    )[[1]][features, , drop = FALSE]
  } else {
    log_message(
      "Input data in not a Seurat object.",
      message_type = "error"
    )
  }

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
#' @param srt_raw A Seurat object to be appended.
#' @param srt_append New Seurat object to append.
#' @param slots slots names.
#' @param pattern A character string containing a regular expression.
#' All data with matching names will be considered for appending.
#' @param overwrite Whether to overwrite.
#' @param verbose Show messages.
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
      "{.arg srt_raw} or {.arg srt_append} is not a Seurat object",
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
    for (info in names(methods::slot(srt_append, name = slot_nm))) {
      if (is.null(info)) {
        if (length(methods::slot(srt_append, name = slot_nm)) > 0 && isTRUE(overwrite)) {
          methods::slot(srt_raw, name = slot_nm) <- methods::slot(
            srt_append,
            name = slot_nm
          )
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
        if (
          slot_nm %in%
            c("assays", "graphs", "neighbors", "reductions", "images")
        ) {
          if (identical(slot_nm, "graphs")) {
            srt_raw@graphs[[info]] <- srt_append[[info]]
          } else if (identical(slot_nm, "assays")) {
            if (!info %in% SeuratObject::Assays(srt_raw)) {
              srt_raw[[info]] <- srt_append[[info]]
            } else {
              # srt_raw[[info]]@counts <- srt_append[[info]]@counts
              srt_raw <- SeuratObject::SetAssayData(
                object = srt_raw,
                layer = "counts",
                assay = "RNA",
                new.data = GetAssayData5(
                  srt_append,
                  assay = info,
                  layer = "counts"
                )
              )
              # srt_raw[[info]]@data <- srt_append[[info]]@data
              srt_raw <- SeuratObject::SetAssayData(
                object = srt_raw,
                layer = "data",
                assay = "RNA",
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
              meta.features <- cbind(
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
                features = meta.features,
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

#' @title Find the default reduction name in a Seurat object
#'
#' @param srt A Seurat object.
#' @param pattern Character string containing a regular expression to search for.
#' @param min_dim Minimum dimension threshold.
#' @param max_distance Maximum distance allowed for a match.
#'
#' @return Default reduction name.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' names(pancreas_sub@reductions)
#' DefaultReduction(pancreas_sub)
#'
#' # Searches for matches to "pca"
#' DefaultReduction(pancreas_sub, pattern = "pca")
#'
#' # Searches for approximate matches to "pc"
#' DefaultReduction(pancreas_sub, pattern = "pc")
DefaultReduction <- function(
    srt,
    pattern = NULL,
    min_dim = 2,
    max_distance = 0.1) {
  if (length(srt@reductions) == 0) {
    log_message(
      "Unable to find any reductions.",
      message_type = "error"
    )
  }
  pattern_default <- c(
    "umap",
    "tsne",
    "dm",
    "phate",
    "pacmap",
    "trimap",
    "largevis",
    "fr",
    "pca",
    "svd",
    "ica",
    "nmf",
    "mds",
    "glmpca"
  )
  pattern_dim <- c("2D", "3D")
  reduc_all <- names(srt@reductions)
  reduc_all <- reduc_all[unlist(lapply(reduc_all, function(x) {
    dim(srt@reductions[[x]]@cell.embeddings)[2] >= min_dim
  }))]
  if (length(reduc_all) == 0) {
    log_message(
      "No dimensional reduction found in the srt object.",
      message_type = "error"
    )
  }
  if (length(reduc_all) == 1) {
    return(reduc_all)
  }
  if (is.null(pattern)) {
    if (("Default_reduction" %in% names(srt@misc))) {
      pattern <- srt@misc[["Default_reduction"]]
    } else {
      pattern <- pattern_default
    }
  }

  pattern <- c(pattern, paste0(pattern, min_dim, "D"))
  if (any(pattern %in% reduc_all)) {
    return(pattern[pattern %in% reduc_all][1])
  }
  index <- c(unlist(sapply(pattern, function(pat) {
    grep(pattern = pat, x = reduc_all, ignore.case = TRUE)
  })))
  if (length(index) > 0) {
    default_reduc <- reduc_all[index]
  } else {
    index <- c(
      unlist(
        sapply(
          pattern, function(pat) {
            agrep(
              pattern = pat,
              x = reduc_all,
              max.distance = max_distance,
              ignore.case = TRUE
            )
          }
        )
      )
    )
    if (length(index) > 0) {
      default_reduc <- reduc_all[index]
    } else {
      default_reduc <- reduc_all
    }
  }
  if (length(default_reduc) > 1) {
    default_reduc <- default_reduc[unlist(sapply(
      c(pattern_default, pattern_dim),
      function(pat) {
        grep(pattern = pat, x = default_reduc, ignore.case = TRUE)
      }
    ))]
    default_reduc <- default_reduc[which.min(sapply(
      default_reduc,
      function(x) dim(srt@reductions[[x]])[2]
    ))]
  }
  return(default_reduc)
}
