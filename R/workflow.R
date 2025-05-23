#' Attempt to recover raw counts from the normalized matrix.
#'
#' @param srt A Seurat object.
#' @param assay Name of assay to recover counts.
#' @param trans The transformation function to applied when data is presumed to be log-normalized.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#'  However, due to decimal point preservation during normalization, the calculated nCount is usually a floating point number close to the integer.
#'  The tolerance is its difference from the integer. Default is 0.1
#' @param sf Set the scaling factor manually.
#' @param verbose Whether to show messages.
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' raw_counts <- pancreas_sub@assays$RNA@counts
#'
#' # Normalized the data
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#'
#' # Now replace counts with the log-normalized data matrix
#' pancreas_sub@assays$RNA@counts <- pancreas_sub@assays$RNA@data
#'
#' # Recover the counts and compare with the raw counts matrix
#' pancreas_sub <- RecoverCounts(pancreas_sub)
#' identical(raw_counts, pancreas_sub@assays$RNA@counts)
RecoverCounts <- function(
    srt,
    assay = NULL,
    trans = c("expm1", "exp", "none"),
    min_count = c(1, 2, 3),
    tolerance = 0.1,
    sf = NULL,
    verbose = TRUE) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  counts <- SeuratObject::GetAssayData(
    srt,
    layer = "counts",
    assay = assay
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- SeuratObject::as.sparse(
      counts[1:nrow(counts), , drop = FALSE]
    )
  }

  status <- check_data_type(data = counts)
  if (status == "raw_counts") {
    if (isTRUE(verbose)) {
      message("The data is already raw counts.")
    }
    return(srt)
  }

  if (status == "log_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be log-normalized.")
    }
    trans <- match.arg(trans)
    if (trans %in% c("expm1", "exp")) {
      if (isTRUE(verbose)) {
        message("Perform ", trans, " on the raw data.")
      }
      counts <- do.call(trans, list(counts))
    }
  }
  if (status == "raw_normalized_counts") {
    if (isTRUE(verbose)) {
      message(
        "The data is presumed to be normalized without log transformation."
      )
    }
  }
  if (is.null(sf)) {
    sf <- unique(round(colSums(counts)))
    if (isTRUE(verbose)) {
      message(
        "The presumed scale factor: ",
        paste0(head(sf, 10), collapse = ", ")
      )
    }
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
      warning(
        "The presumed nCount of some cells is not valid: ",
        paste0(
          head(colnames(counts)[diff_value < tolerance], 10),
          collapse = ","
        ),
        ", ...",
        immediate. = TRUE
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
    warning(
      "Scale factor is not unique. No changes to be made.",
      immediate. = TRUE
    )
  }

  return(srt)
}

#' Rename features for the Seurat object
#'
#' @param srt A Seurat object.
#' @param newnames A vector with the same length of features in Seurat object,
#' or characters named with old features.
#' @param assays Assays to rename.
#'
#' @export
#'
#' @examples
#' data("panc8_sub")
#' head(rownames(panc8_sub))
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   capitalize(rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' panc8_rename <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' head(rownames(panc8_rename))
RenameFeatures <- function(
    srt,
    newnames = NULL,
    assays = NULL) {
  assays <- assays[assays %in% SeuratObject::Assays(srt)] %||% SeuratObject::Assays(srt)
  if (is.null(names(newnames))) {
    if (!identical(length(newnames), nrow(srt))) {
      stop("'newnames' must be named or the length of features in the srt.")
    }
    if (length(unique(sapply(srt@assays[assays], nrow))) > 1) {
      stop(
        "Assays in the srt object have different number of features. Please use a named vectors."
      )
    }
    names(newnames) <- rownames(srt[[assays[1]]])
  }
  for (assay in assays) {
    message("Rename features for the assay: ", assay)
    assay_obj <- Seurat::GetAssay(srt, assay)
    for (d in c("meta.features", "scale.data", "counts", "data")) {
      index <- which(
        rownames(methods::slot(assay_obj, d)) %in% names(newnames)
      )
      rownames(methods::slot(assay_obj, d))[index] <- newnames[rownames(
        methods::slot(
          assay_obj,
          d
        )
      )[index]]
    }
    if (length(methods::slot(assay_obj, "var.features")) > 0) {
      index <- which(
        methods::slot(
          assay_obj, "var.features"
        ) %in% names(newnames)
      )
      methods::slot(assay_obj, "var.features")[index] <- newnames[methods::slot(
        assay_obj,
        "var.features"
      )[index]]
    }
    srt[[assay]] <- assay_obj
  }
  return(srt)
}

#' Rename clusters for the Seurat object
#'
#' @param srt A Seurat object.
#' @param group.by The old group used to rename cells.
#' @param nameslist A named list of new cluster value.
#' @param name The name of the new cluster stored in the Seurat object.
#' @param keep_levels If the old group is a factor, keep the order of the levels.
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' levels(pancreas_sub@meta.data[["SubCellType"]])
#'
#' # Rename all clusters
#' pancreas_sub <- RenameClusters(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = letters[1:8]
#' )
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Rename specified clusters
#' pancreas_sub <- RenameClusters(pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list("a" = "Alpha", "b" = "Beta")
#' )
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Merge and rename clusters
#' pancreas_sub <- RenameClusters(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list(
#'     "EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")
#'   ),
#'   name = "Merged",
#'   keep_levels = TRUE
#' )
#' CellDimPlot(pancreas_sub, "Merged")
RenameClusters <- function(
    srt,
    group.by,
    nameslist = list(),
    name = "newclusters",
    keep_levels = FALSE) {
  if (missing(group.by)) {
    stop("group.by must be provided")
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    stop(paste0(group.by, " is not in the meta.data of srt object."))
  }
  if (length(nameslist) > 0 && is.null(names(nameslist))) {
    names(nameslist) <- levels(srt@meta.data[[group.by]])
  }
  if (is.list(nameslist) && length(nameslist) > 0) {
    names_assign <- stats::setNames(
      rep(names(nameslist), sapply(nameslist, length)),
      nm = unlist(nameslist)
    )
  } else {
    if (is.null(names(nameslist))) {
      if (!is.factor(srt@meta.data[[group.by]])) {
        stop(
          "'nameslist' must be named when srt@meta.data[[group.by]] is not a factor"
        )
      }
      if (
        !identical(length(nameslist), length(unique(srt@meta.data[[group.by]])))
      ) {
        stop(
          "'nameslist' must be named or the length of ",
          length(unique(srt@meta.data[[group.by]]))
        )
      }
      names(nameslist) <- levels(srt@meta.data[[group.by]])
    }
    names_assign <- nameslist
  }
  if (all(!names(names_assign) %in% srt@meta.data[[group.by]])) {
    stop("No group name mapped.")
  }
  if (is.factor(srt@meta.data[[group.by]])) {
    levels <- levels(srt@meta.data[[group.by]])
  } else {
    levels <- NULL
  }
  index <- which(
    as.character(srt@meta.data[[group.by]]) %in% names(names_assign)
  )
  srt@meta.data[[name]] <- as.character(srt@meta.data[[group.by]])
  srt@meta.data[[name]][index] <- names_assign[srt@meta.data[[name]][index]]
  if (!is.null(levels)) {
    levels[levels %in% names(names_assign)] <- names_assign[levels[
      levels %in% names(names_assign)
    ]]
    if (isFALSE(keep_levels)) {
      levels <- unique(c(names_assign, levels))
    } else {
      levels <- unique(levels)
    }
    srt@meta.data[[name]] <- factor(srt@meta.data[[name]], levels = levels)
  }
  return(srt)
}

#' Reorder idents by the gene expression
#'
#' @param srt A Seurat object.
#' @param features Features used to reorder idents.
#' @param reorder_by Reorder groups instead of idents.
#' @param layer Specific layer to get data from.
#' @param assay Specific assay to get data from.
#' @param log Whether log1p transformation needs to be applied. Default is \code{TRUE}.
#' @param distance_metric Metric to compute distance. Default is "euclidean".
#'
#' @importFrom proxyC simil dist
#' @export
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
    features <- SeuratObject::DefaultAssay(srt, assay = assay)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  if (is.null(reorder_by)) {
    srt$ident <- Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  if (length(unique(srt[[reorder_by, drop = TRUE]])) == 1) {
    warning("Only one cluter found.", immediate. = TRUE)
    return(srt)
  }
  simil_method <- c(
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
  dist_method <- c(
    "euclidean",
    "chisquared",
    "kullback",
    "manhattan",
    "maximum",
    "canberra",
    "minkowski",
    "hamming"
  )
  if (
    !distance_metric %in% c(simil_method, dist_method, "pearson", "spearman")
  ) {
    stop(distance_metric, " method is invalid.")
  }

  data.avg <- Seurat::AverageExpression(
    object = srt,
    features = features,
    layer = layer,
    assays = assay,
    group.by = "ident",
    verbose = FALSE
  )[[1]][features, , drop = FALSE]
  if (isTRUE(log)) {
    data.avg <- log1p(data.avg)
  }
  mat <- Matrix::t(x = data.avg[features, , drop = FALSE])
  if (!inherits(mat, "dgCMatrix")) {
    mat <- SeuratObject::as.sparse(mat[1:nrow(mat), , drop = FALSE])
  }

  if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- Matrix::t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 -
      simil(
        SeuratObject::as.sparse(mat[1:nrow(mat), , drop = FALSE]),
        method = distance_metric
      )
  } else if (distance_metric %in% dist_method) {
    d <- dist(
      SeuratObject::as.sparse(mat[1:nrow(mat), , drop = FALSE]),
      method = distance_metric
    )
  }
  data.dist <- stats::as.dist(d)
  hc <- stats::hclust(d = data.dist)
  dd <- stats::as.dendrogram(hc)
  dd_ordered <- stats::reorder(
    dd,
    wts = Matrix::colMeans(data.avg[features, , drop = FALSE]),
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

#' Append a Seurat object to another
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
SrtAppend <- function(
    srt_raw,
    srt_append,
    slots = slotNames(srt_append),
    pattern = NULL,
    overwrite = FALSE,
    verbose = TRUE) {
  if (!inherits(srt_raw, "Seurat") || !inherits(srt_append, "Seurat")) {
    stop("'srt_raw' or 'srt_append' is not a Seurat object.")
  }

  pattern <- pattern %||% ""
  for (slot_nm in slotNames(srt_append)) {
    if (!slot_nm %in% slots) {
      if (isTRUE(verbose)) {
        message("Slot ", slot_nm, " is not appended.")
      }
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
          methods::slot(srt_raw, name = slot_nm) <- methods::slot(srt_append, name = slot_nm)
        }
        next
      }
      if (!grepl(pattern = pattern, x = info)) {
        if (isTRUE(verbose)) {
          message(info, " in slot ", slot_nm, " is not appended.")
        }
        next
      }
      if (
        !info %in% names(methods::slot(srt_raw, name = slot_nm)) || isTRUE(overwrite)
      ) {
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
              srt_raw[[info]]@counts <- srt_append[[info]]@counts
              srt_raw[[info]]@data <- srt_append[[info]]@data
              srt_raw[[info]]@key <- srt_append[[info]]@key
              srt_raw[[info]]@var.features <- srt_append[[info]]@var.features
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

#' Run dimensionality reduction
#'
#' @param srt A Seurat object.
#' @param prefix The prefix used to name the result.
#' @param features Use features expression data to run linear or nonlinear dimensionality reduction.
#' @param assay Specific assay to get data from.
#' @param layer Specific layer to get data from.
#' @param linear_reduction Method of linear dimensionality reduction. Options are "pca", "ica", "nmf", "mds", "glmpca".
#' @param linear_reduction_dims Total number of dimensions to compute and store for \code{linear_reduction}.
#' @param linear_reduction_params Other parameters passed to the \code{linear_reduction} method.
#' @param force_linear_reduction Whether force to do linear dimensionality reduction.
#' @param nonlinear_reduction Method of nonlinear dimensionality reduction. Options are "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"
#' @param nonlinear_reduction_dims Total number of dimensions to compute and store for \code{nonlinear_reduction}.
#' @param reduction_use Which dimensional reduction to use as input for \code{nonlinear_reduction}.
#' @param reduction_dims Which dimensions to use as input for \code{nonlinear_reduction}, used only if \code{features} is \code{NULL}.
#' @param neighbor_use Name of neighbor to use for the \code{nonlinear_reduction}.
#' @param graph_use Name of graph to use for the \code{nonlinear_reduction}.
#' @param nonlinear_reduction_params  Other parameters passed to the \code{nonlinear_reduction} method.
#' @param force_nonlinear_reduction Whether force to do nonlinear dimensionality reduction.
#' @param verbose Show messages.
#' @param seed Set a seed.
#'
#' @importFrom Signac RunSVD
#' @export
RunDimReduction <- function(
    srt,
    prefix = "",
    features = NULL,
    assay = NULL,
    layer = "data",
    linear_reduction = NULL,
    linear_reduction_dims = 50,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = NULL,
    nonlinear_reduction_dims = 2,
    reduction_use = NULL,
    reduction_dims = NULL,
    graph_use = NULL,
    neighbor_use = NULL,
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    verbose = TRUE,
    seed = 11) {
  set.seed(seed)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    type <- "Chromatin"
  } else {
    type <- "RNA"
  }
  linear_reduction_dims <- min(
    linear_reduction_dims,
    nrow(srt[[assay]]) - 1,
    ncol(srt[[assay]]) - 1,
    na.rm = TRUE
  )
  nonlinear_reduction_dims <- min(
    nonlinear_reduction_dims,
    nrow(srt[[assay]]) - 1,
    ncol(srt[[assay]]) - 1,
    na.rm = TRUE
  )
  if (!is.null(linear_reduction)) {
    if (
      any(
        !linear_reduction %in%
          c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt))
      ) ||
        length(linear_reduction) > 1
    ) {
      stop(
        "'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'."
      )
    }
  }
  if (!is.null(nonlinear_reduction)) {
    if (
      any(
        !nonlinear_reduction %in%
          c(
            "umap",
            "umap-naive",
            "tsne",
            "dm",
            "phate",
            "pacmap",
            "trimap",
            "largevis",
            "fr",
            Reductions(srt)
          )
      ) ||
        length(nonlinear_reduction) > 1
    ) {
      stop(
        "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
      )
    }
    if (
      is.null(features) &&
        is.null(reduction_use) &&
        is.null(neighbor_use) &&
        is.null(graph_use)
    ) {
      stop(
        "'features', 'reduction_use', 'neighbor_use', or 'graph_use' must be provided when running non-linear dimensionality reduction."
      )
    }
    if (nonlinear_reduction %in% c("fr")) {
      if (!is.null(graph_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Graphs(",
          graph_use,
          ") as input"
        )
      } else if (!is.null(neighbor_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Neighbors(",
          neighbor_use,
          ") as input"
        )
      } else if (!is.null(features)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Features(length:",
          length(features),
          ") as input"
        )
      } else if (!is.null(reduction_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Reduction(",
          reduction_use,
          ", dims:",
          min(reduction_dims),
          "-",
          max(reduction_dims),
          ") as input"
        )
      }
    } else {
      if (!is.null(features)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Features(length:",
          length(features),
          ") as input"
        )
      } else if (!is.null(reduction_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Reduction(",
          reduction_use,
          ", dims:",
          min(reduction_dims),
          "-",
          max(reduction_dims),
          ") as input"
        )
      } else if (!is.null(neighbor_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Neighbors(",
          neighbor_use,
          ") as input"
        )
      } else if (!is.null(graph_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Graphs(",
          graph_use,
          ") as input"
        )
      }
    }
  }
  if (!is.null(linear_reduction)) {
    if (!isTRUE(force_linear_reduction)) {
      if (linear_reduction %in% Reductions(srt)) {
        if (srt[[linear_reduction]]@assay.used == assay) {
          message(
            "linear_reduction(",
            linear_reduction,
            ") is already existed. Skip calculation."
          )
          reduc <- srt[[linear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, linear_reduction, "_")
          srt[[paste0(prefix, linear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
          return(srt)
        } else {
          message(
            "assay.used is ",
            srt[[linear_reduction]]@assay.used,
            ", which is not the same as the ",
            assay,
            " specified. Recalculate the linear reduction(pca)"
          )
          linear_reduction <- "pca"
        }
      }
    }
    if (is.null(features) || length(features) == 0) {
      message("No features provided. Use variable features.")
      if (length(SeuratObject::DefaultAssay(srt, assay = assay)) == 0) {
        srt <- FindVariableFeatures(srt, assay = assay, verbose = FALSE)
      }
      features <- SeuratObject::DefaultAssay(srt, assay = assay)
    }
    fun_use <- switch(linear_reduction,
      "pca" = "RunPCA",
      "svd" = "RunSVD",
      "ica" = "RunICA",
      "nmf" = "RunNMF",
      "mds" = "RunMDS",
      "glmpca" = "RunGLMPCA"
    )
    key_use <- switch(linear_reduction,
      "pca" = "PC_",
      "svd" = "LSI_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    components_nm <- switch(linear_reduction,
      "pca" = "npcs",
      "svd" = "n",
      "ica" = "nics",
      "nmf" = "nbes",
      "mds" = "nmds",
      "glmpca" = "L"
    )
    params <- list(
      object = srt,
      assay = assay,
      layer = layer,
      features = features,
      components_nm = linear_reduction_dims,
      reduction.name = paste0(prefix, linear_reduction),
      reduction.key = paste0(prefix, key_use),
      verbose = verbose,
      seed.use = seed
    )
    if (fun_use %in% c("RunSVD", "RunICA")) {
      params <- params[!names(params) %in% "layer"]
    }
    if (fun_use == "RunGLMPCA") {
      params[["layer"]] <- "counts"
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(linear_reduction_params)) {
      params[[nm]] <- linear_reduction_params[[nm]]
    }
    srt <- invoke_fun(.fn = fun_use, .args = params)

    if (is.null(rownames(srt[[paste0(prefix, linear_reduction)]]))) {
      rownames(
        srt[[paste0(prefix, linear_reduction)]]@cell.embeddings
      ) <- colnames(srt)
    }
    if (linear_reduction == "pca") {
      pca.out <- srt[[paste0(prefix, linear_reduction)]]
      center <- rowMeans(
        Seurat::GetAssayData(
          object = srt,
          layer = "scale.data",
          assay = assay
        )[features, , drop = FALSE]
      )
      model <- list(
        sdev = pca.out@stdev,
        rotation = pca.out@feature.loadings,
        center = center,
        scale = FALSE,
        x = pca.out@cell.embeddings
      )
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[[
        "model"
      ]] <- model
    }
    if (linear_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:linear_reduction_dims
    } else {
      dim_est <- tryCatch(
        expr = {
          min(
            intrinsicDimension::maxLikGlobalDimEst(
              data = Embeddings(
                srt,
                reduction = paste0(prefix, linear_reduction)
              ),
              k = 20
            )[["dim.est"]],
            ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction)))
          )
        },
        error = function(e) {
          message(
            "Can not estimate intrinsic dimensions with maxLikGlobalDimEst."
          )
          return(NA)
        }
      )
      if (!is.na(dim_est)) {
        dims_estimate <- seq_len(max(
          min(
            ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))),
            10
          ),
          ceiling(dim_est)
        ))
      } else {
        dims_estimate <- seq_len(min(
          ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))),
          30
        ))
      }
    }
    srt@reductions[[paste0(prefix, linear_reduction)]]@misc[[
      "dims_estimate"
    ]] <- dims_estimate
    srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
  } else if (!is.null(nonlinear_reduction)) {
    if (!isTRUE(force_nonlinear_reduction)) {
      if (nonlinear_reduction %in% Reductions(srt)) {
        if (srt[[nonlinear_reduction]]@assay.used == assay) {
          message(
            "nonlinear_reduction(",
            nonlinear_reduction,
            ") is already existed. Skip calculation."
          )
          reduc <- srt[[nonlinear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, nonlinear_reduction, "_")
          srt[[paste0(prefix, nonlinear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction)
          return(srt)
        } else {
          message(
            "assay.used is ",
            srt[[nonlinear_reduction]]@assay.used,
            ", which is not the same as the ",
            assay,
            " specified. Recalculate the nonlinear reduction(umap)"
          )
          nonlinear_reduction <- "umap"
        }
      }
    }
    # if (!is.null(neighbor_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'neighbor_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    # if (!is.null(graph_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'graph_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    fun_use <- switch(nonlinear_reduction,
      "umap" = "RunUMAP2",
      "umap-naive" = "RunUMAP2",
      "tsne" = "RunTSNE",
      "dm" = "RunDM",
      "phate" = "RunPHATE",
      "pacmap" = "RunPaCMAP",
      "trimap" = "RunTriMap",
      "largevis" = "RunLargeVis",
      "fr" = "RunFR"
    )
    components_nm <- switch(nonlinear_reduction,
      "umap" = "n.components",
      "umap-naive" = "n.components",
      "tsne" = "dim.embed",
      "dm" = "ndcs",
      "phate" = "n_components",
      "pacmap" = "n_components",
      "trimap" = "n_components",
      "largevis" = "n_components",
      "fr" = "ndim"
    )
    other_params <- switch(nonlinear_reduction,
      "umap" = list(umap.method = "uwot", return.model = TRUE),
      "umap-naive" = list(umap.method = "naive", return.model = TRUE),
      "tsne" = list(
        tsne.method = "Rtsne",
        num_threads = 0,
        check_duplicates = FALSE
      ),
      "dm" = list(),
      "phate" = list(),
      "pacmap" = list(),
      "trimap" = list(),
      "largevis" = list(),
      "fr" = list()
    )
    nonlinear_reduction_sim <- toupper(gsub(
      pattern = "-.*",
      replacement = "",
      x = nonlinear_reduction
    ))
    params <- list(
      object = srt,
      assay = assay,
      layer = layer,
      components_nm = nonlinear_reduction_dims,
      features = features,
      reduction = reduction_use,
      dims = reduction_dims,
      reduction.name = paste0(
        prefix,
        nonlinear_reduction_sim,
        nonlinear_reduction_dims,
        "D"
      ),
      reduction.key = paste0(
        prefix,
        nonlinear_reduction_sim,
        nonlinear_reduction_dims,
        "D_"
      ),
      verbose = verbose,
      seed.use = seed
    )
    if (!is.null(neighbor_use)) {
      params[["neighbor"]] <- neighbor_use
    }
    if (!is.null(graph_use)) {
      params[["graph"]] <- graph_use
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(other_params)) {
      params[[nm]] <- other_params[[nm]]
    }
    for (nm in names(nonlinear_reduction_params)) {
      params[[nm]] <- nonlinear_reduction_params[[nm]]
    }
    srt <- invoke_fun(.fn = fun_use, .args = params)

    srt@reductions[[paste0(
      prefix,
      nonlinear_reduction_sim,
      nonlinear_reduction_dims,
      "D"
    )]]@misc[["reduction_dims"]] <- reduction_dims
    srt@reductions[[paste0(
      prefix,
      nonlinear_reduction_sim,
      nonlinear_reduction_dims,
      "D"
    )]]@misc[["reduction_use"]] <- reduction_use
    srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction_sim)
  }
  return(srt)
}

#' Find the default reduction name in a Seurat object.
#'
#' @param srt A Seurat object.
#' @param pattern Character string containing a regular expression to search for.
#' @param min_dim Minimum dimension threshold.
#' @param max_distance Maximum distance allowed for a match.
#'
#' @examples
#' data("pancreas_sub")
#' names(pancreas_sub@reductions)
#' DefaultReduction(pancreas_sub)
#'
#' # Searches for matches to "pca"
#' DefaultReduction(pancreas_sub, pattern = "pca")
#'
#' # Searches for approximate matches to "pc"
#' DefaultReduction(pancreas_sub, pattern = "pc")
#'
#' @return Default reduction name.
#'
#' @export
DefaultReduction <- function(
    srt,
    pattern = NULL,
    min_dim = 2,
    max_distance = 0.1) {
  if (length(srt@reductions) == 0) {
    stop("Unable to find any reductions.")
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
    stop("No dimensional reduction found in the srt object.")
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
    index <- c(unlist(sapply(pattern, function(pat) {
      agrep(
        pattern = pat,
        x = reduc_all,
        max.distance = max_distance,
        ignore.case = TRUE
      )
    })))
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

#' Uncorrected_integrate
#'
#' @inheritParams Integration_scop
#'
#' @export
Uncorrected_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(Uncorrected) on the data...\n"
  ))

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(Seurat::GetAssayData(
              srt_merge,
              layer = "scale.data",
              assay = SeuratObject::DefaultAssay(srt_merge)
            ))
        ))
  ) {
    if (normalization_method != "SCT") {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srt_merge <- Seurat::ScaleData(
        object = srt_merge,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srt_merge),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "Uncorrected",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "Uncorrected",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srt_merge <- tryCatch(
    {
      srt_merge <- FindNeighbors(
        object = srt_merge,
        reduction = paste0("Uncorrected", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("Uncorrected_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srt_merge <- FindClusters(
        object = srt_merge,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "Uncorrected_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srt_merge <- srt_reorder(
        srt_merge,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt_merge[["seurat_clusters"]] <- NULL
      srt_merge[[paste0("Uncorrected", linear_reduction, "clusters")]] <- Idents(
        srt_merge
      )
      srt_merge
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srt_merge)
    }
  )

  srt_merge <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srt_merge <- RunDimReduction(
            srt_merge,
            prefix = "Uncorrected",
            reduction_use = paste0("Uncorrected", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "Uncorrected_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srt_merge
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srt_merge)
    }
  )

  SeuratObject::DefaultAssay(srt_merge) <- assay
  SeuratObject::VariableFeatures(srt_merge) <- srt_merge@misc[["Uncorrected_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(assay, "|Uncorrected|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_merge)
  }
}

#' Seurat_integrate
#'
#' @inheritParams Integration_scop
#' @param FindIntegrationAnchors_params A list of parameters for the Seurat::FindIntegrationAnchors function, default is an empty list.
#' @param IntegrateData_params A list of parameters for the Seurat::IntegrateData function, default is an empty list.
#' @param IntegrateEmbeddings_params A list of parameters for the Seurat::IntegrateEmbeddings function, default is an empty list.
#'
#' @importFrom Signac RunTFIDF
#' @export
Seurat_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    FindIntegrationAnchors_params = list(),
    IntegrateData_params = list(),
    IntegrateEmbeddings_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "svd", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    if (normalization_method == "TFIDF") {
      srt_merge <- Reduce(merge, srt_list)
      SeuratObject::VariableFeatures(srt_merge) <- HVF
    }
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srt_list, ncol)) < 50) {
    warning(
      "The cell count in some batches is lower than 50, which may not be suitable for the current integration method.",
      immediate. = TRUE
    )
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srt_merge)
    }
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'rlsi' integration workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      FindIntegrationAnchors_params[["dims"]] <- 2:min(
        linear_reduction_dims,
        30
      )
    }
    srt_merge <- RunTFIDF(
      object = srt_merge,
      assay = SeuratObject::DefaultAssay(srt_merge),
      verbose = FALSE
    )
    srt_merge <- RunDimReduction(
      srt_merge,
      prefix = "",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt_merge),
      linear_reduction = "svd",
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    srt_merge[["lsi"]] <- srt_merge[["svd"]]
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform linear dimension reduction (svd) on the data ",
        i,
        " ...\n"
      ))
      srt <- RunDimReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "svd",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = FALSE,
        seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srt_list[[i]] <- srt
    }
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    cat(paste0("[", Sys.time(), "]", " Use 'rpca' integration workflow...\n"))
    for (i in seq_along(srt_list)) {
      srt <- srt_list[[i]]
      if (
        isTRUE(do_scaling) ||
          (is.null(do_scaling) &&
            any(
              !HVF %in%
                rownames(
                  Seurat::GetAssayData(
                    srt,
                    layer = "scale.data",
                    assay = SeuratObject::DefaultAssay(srt)
                  )
                )
            ))
      ) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform ScaleData on the data ",
          i,
          " ...\n"
        ))
        srt <- Seurat::ScaleData(
          object = srt,
          assay = SeuratObject::DefaultAssay(srt),
          features = HVF,
          vars.to.regress = vars_to_regress,
          model.use = regression_model,
          verbose = FALSE
        )
      }
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform linear dimension reduction (pca) on the data ",
        i,
        " ...\n"
      ))
      srt <- RunDimReduction(
        srt,
        prefix = "",
        features = HVF,
        assay = SeuratObject::DefaultAssay(srt),
        linear_reduction = "pca",
        linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params,
        force_linear_reduction = force_linear_reduction,
        verbose = FALSE,
        seed = seed
      )
      srt_list[[i]] <- srt
    }
  }

  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (normalization_method %in% c("LogNormalize", "SCT")) {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform FindIntegrationAnchors on the data...\n"
    ))
    params1 <- list(
      object.list = srt_list,
      normalization.method = normalization_method,
      anchor.features = HVF,
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(.fn = Seurat::FindIntegrationAnchors, .args = params1)

    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform integration(Seurat) on the data...\n"
    ))
    params2 <- list(
      anchorset = srt_anchors,
      new.assay.name = "Seuratcorrected",
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srtIntegrated <- invoke_fun(.fn = Seurat::IntegrateData, .args = params2)

    SeuratObject::DefaultAssay(srtIntegrated) <- "Seuratcorrected"
    SeuratObject::VariableFeatures(srtIntegrated[["Seuratcorrected"]]) <- HVF

    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) &&
          any(
            !HVF %in%
              rownames(
                Seurat::GetAssayData(
                  srtIntegrated,
                  layer = "scale.data",
                  assay = SeuratObject::DefaultAssay(srtIntegrated)
                )
              )
          ))
    ) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srtIntegrated <- Seurat::ScaleData(
        object = srtIntegrated,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srtIntegrated),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }

    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform linear dimension reduction (",
      linear_reduction,
      ") on the data...\n"
    ))
    srtIntegrated <- RunDimReduction(
      srtIntegrated,
      prefix = "Seurat",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srtIntegrated),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srtIntegrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]] %||%
        1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform FindIntegrationAnchors on the data...\n"
    ))
    params1 <- list(
      object.list = srt_list,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke_fun(.fn = Seurat::FindIntegrationAnchors, .args = params1)

    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform integration(Seurat) on the data...\n"
    ))
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srt_merge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srtIntegrated <- invoke_fun(.fn = IntegrateEmbeddings, .args = params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:max(srtIntegrated@reductions[[paste0(
        "Seurat",
        linear_reduction
      )]]@misc[["dims_estimate"]]) %||%
        2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = paste0("Seurat", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("Seurat_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "Seurat_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Seurat",
            reduction_use = paste0("Seurat", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "Seurat_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|Seurat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' scVI_integrate
#'
#' @inheritParams Integration_scop
#' @param scVI_dims_use A vector specifying the dimensions returned by scVI that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param model A string indicating the scVI model to be used.
#' Options are "SCVI" and "PEAKVI".
#' Default is "SCVI".
#' @param SCVI_params A list of parameters for the SCVI model.
#' Default is an empty list.
#' @param PEAKVI_params A list of parameters for the PEAKVI model.
#' Default is an empty list.
#' @param num_threads An integer setting the number of threads for scVI.
#' Default is 8.
#'
#' @export
scVI_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    scVI_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    model = "SCVI",
    SCVI_params = list(),
    PEAKVI_params = list(),
    num_threads = 8,
    seed = 11) {
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (
    .Platform$OS.type == "windows" &&
      !exist_Python_pkgs(packages = "scvi-tools")
  ) {
    suppressWarnings(system2(
      command = conda_python(),
      args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver",
      stdout = TRUE
    ))
  }

  check_python("scvi-tools")
  scvi <- reticulate::import("scvi")
  scipy <- reticulate::import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(num_threads)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  adata <- srt_to_adata(
    srt_merge,
    features = HVF,
    assay_x = SeuratObject::DefaultAssay(srt_merge),
    assay_y = NULL,
    verbose = FALSE
  )
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  if (model == "SCVI") {
    scvi$model$SCVI$setup_anndata(adata, batch_key = batch)
    params <- list(
      adata = adata
    )
    for (nm in names(SCVI_params)) {
      params[[nm]] <- SCVI_params[[nm]]
    }
    model <- invoke_fun(.fn = scvi$model$SCVI, .args = params)
    model$train()
    srtIntegrated <- srt_merge
    srt_merge <- NULL
    corrected <- Matrix::t(
      Matrix::as.matrix(
        model$get_normalized_expression()
      )
    )
    srtIntegrated[["scVIcorrected"]] <- SeuratObject::CreateAssayObject(
      counts = corrected
    )
    SeuratObject::DefaultAssay(srtIntegrated) <- "scVIcorrected"
    SeuratObject::VariableFeatures(srtIntegrated[["scVIcorrected"]]) <- HVF
  } else if (model == "PEAKVI") {
    message("Assay is ChromatinAssay. Using PeakVI workflow.")
    scvi$model$PEAKVI$setup_anndata(adata, batch_key = batch)
    params <- list(
      adata = adata
    )
    for (nm in names(PEAKVI_params)) {
      params[[nm]] <- PEAKVI_params[[nm]]
    }
    model <- invoke_fun(.fn = scvi$model$PEAKVI, .args = params)
    model$train()
    srtIntegrated <- srt_merge
    srt_merge <- NULL
  }

  latent <- Matrix::as.matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
  srtIntegrated[["scVI"]] <- CreateDimReducObject(
    embeddings = latent,
    key = "scVI_",
    assay = SeuratObject::DefaultAssay(srtIntegrated)
  )
  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srtIntegrated[["scVI"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "scVI",
        dims = scVI_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("scVI_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "scVI_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["scVIclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "scVI",
            reduction_use = "scVI",
            reduction_dims = scVI_dims_use,
            graph_use = "scVI_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["scVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|scVI|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' MNN_integrate
#'
#' @inheritParams Integration_scop
#' @param mnnCorrect_params A list of parameters for the batchelor::mnnCorrect function, default is an empty list.
#'
#' @export
MNN_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    mnnCorrect_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- SplitObject(object = srt_merge, split.by = batch)
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  sceList <- lapply(srt_list, function(srt) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(
      counts = Seurat::GetAssayData(
        srt,
        layer = "data",
        assay = SeuratObject::DefaultAssay(srt)
      )[
        HVF, ,
        drop = FALSE
      ]
    ))
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- Matrix::as.matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(MNN) on the data...\n"
  ))
  params <- list(
    sceList,
    cos.norm.out = FALSE
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- invoke_fun(.fn = batchelor::mnnCorrect, .args = params)

  srtIntegrated <- srt_merge
  srt_merge <- NULL
  srtIntegrated[["MNNcorrected"]] <- CreateAssayObject(
    counts = out@assays@data$corrected
  )
  SeuratObject::VariableFeatures(srtIntegrated[["MNNcorrected"]]) <- HVF
  SeuratObject::DefaultAssay(srtIntegrated) <- "MNNcorrected"

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(
              Seurat::GetAssayData(
                srtIntegrated,
                layer = "scale.data",
                assay = SeuratObject::DefaultAssay(srtIntegrated)
              )
            )
        ))
  ) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- Seurat::ScaleData(
      object = srtIntegrated,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srtIntegrated),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srtIntegrated <- RunDimReduction(
    srtIntegrated,
    prefix = "MNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0(
      "MNN",
      linear_reduction
    )]]@misc[["dims_estimate"]] %||%
      1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = paste0("MNN", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("MNN_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "MNN_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["MNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "MNN",
            reduction_use = paste0("MNN", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "MNN_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["MNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|MNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' fastMNN_integrate
#'
#' @inheritParams Integration_scop
#' @param fastMNN_dims_use A vector specifying the dimensions returned by fastMNN that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param fastMNN_params A list of parameters for the batchelor::fastMNN function, default is an empty list.
#'
#' @export
fastMNN_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    fastMNN_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    fastMNN_params = list(),
    seed = 11) {
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("batchelor")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- SplitObject(object = srt_merge, split.by = batch)
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  sceList <- lapply(srt_list, function(srt) {
    sce <- Seurat::as.SingleCellExperiment(
      Seurat::CreateSeuratObject(
        counts = Seurat::GetAssayData(
          srt,
          layer = "data",
          assay = SeuratObject::DefaultAssay(srt)
        )[HVF, , drop = FALSE]
      )
    )
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- Matrix::as.matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(fastMNN) on the data...\n"
  ))
  params <- list(
    sceList
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke_fun(.fn = batchelor::fastMNN, .args = params)

  srtIntegrated <- srt_merge
  srt_merge <- NULL
  srtIntegrated[["fastMNNcorrected"]] <- CreateAssayObject(
    counts = Matrix::as.matrix(out@assays@data$reconstructed)
  )
  SeuratObject::DefaultAssay(srtIntegrated) <- "fastMNNcorrected"
  SeuratObject::VariableFeatures(srtIntegrated[["fastMNNcorrected"]]) <- HVF
  reduction <- out@int_colData$reducedDims$corrected
  colnames(reduction) <- paste0("fastMNN_", seq_len(ncol(reduction)))
  srtIntegrated[["fastMNN"]] <- CreateDimReducObject(
    embeddings = reduction,
    key = "fastMNN_",
    assay = "fastMNNcorrected"
  )

  if (is.null(fastMNN_dims_use)) {
    fastMNN_dims_use <- 1:ncol(srtIntegrated[["fastMNN"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "fastMNN",
        dims = fastMNN_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("fastMNN", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "fastMNN_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["fastMNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "fastMNN",
            reduction_use = "fastMNN",
            reduction_dims = fastMNN_dims_use,
            graph_use = "fastMNN_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["fastMNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|fastMNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Harmony_integrate
#'
#' @inheritParams Integration_scop
#' @param Harmony_dims_use A vector specifying the dimensions returned by RunHarmony that will be utilized for downstream cell cluster finding and non-linear reduction.
#' If set to NULL, all the returned dimensions will be used by default.
#' @param RunHarmony_params A list of parameters for the harmony::RunHarmony function, default is an empty list.
#'
#' @export
Harmony_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    Harmony_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    RunHarmony_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("harmony@1.1.0")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(
              Seurat::GetAssayData(
                srt_merge,
                layer = "scale.data",
                assay = SeuratObject::DefaultAssay(srt_merge)
              )
            )
        ))
  ) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "Harmony",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "Harmony",
      linear_reduction
    )]]@misc[["dims_estimate"]] %||%
      1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(Harmony) on the data...\n"
  ))
  message(
    "Harmony integration using Reduction(",
    paste0("Harmony", linear_reduction),
    ", dims:",
    min(linear_reduction_dims_use),
    "-",
    max(linear_reduction_dims_use),
    ") as input"
  )
  params <- list(
    object = srt_merge,
    group.by.vars = batch,
    assay.use = SeuratObject::DefaultAssay(srt_merge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  if (
    nrow(
      Seurat::GetAssayData(
        srt_merge,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt_merge)
      )
    ) ==
      0
  ) {
    params[["project.dim"]] <- FALSE
  }
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srtIntegrated <- invoke_fun(.fn = RunHarmony2, .args = params)

  if (is.null(Harmony_dims_use)) {
    Harmony_dims_use <- 1:ncol(srtIntegrated[["Harmony"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "Harmony",
        dims = Harmony_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("Harmony", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "Harmony_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Harmonyclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Harmony",
            reduction_use = "Harmony",
            reduction_dims = Harmony_dims_use,
            graph_use = "Harmony_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            force_nonlinear_reduction = force_nonlinear_reduction,
            nonlinear_reduction_params = nonlinear_reduction_params,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Harmony_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|Harmony|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Scanorama_integrate
#'
#' @inheritParams Integration_scop
#' @param Scanorama_dims_use  A vector specifying the dimensions returned by Scanorama that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param return_corrected Logical indicating whether to return the corrected data. Default is FALSE.
#' @param Scanorama_params A list of parameters for the scanorama.correct function, default is an empty list.
#'
#' @importFrom stats sd
#' @export
Scanorama_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    Scanorama_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    return_corrected = FALSE,
    Scanorama_params = list(),
    seed = 11) {
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_python("scanorama")
  scanorama <- reticulate::import("scanorama")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    srt_list <- SplitObject(object = srt_merge, split.by = batch)
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  srtIntegrated <- Reduce(merge, srt_list)

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(Scanorama) on the data...\n"
  ))
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srt_list)) {
    assaylist[[i]] <- Matrix::t(
      Matrix::as.matrix(
        Seurat::GetAssayData(
          object = srt_list[[i]],
          layer = "data",
          assay = SeuratObject::DefaultAssay(srt_list[[i]])
        )[HVF, , drop = FALSE]
      )
    )
    genelist[[i]] <- HVF
  }
  if (isTRUE(return_corrected)) {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      return_dimred = TRUE,
      return_dense = TRUE,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    corrected <- invoke_fun(.fn = scanorama$correct, .args = params)

    cor_value <- Matrix::t(do.call(rbind, corrected[[2]]))
    rownames(cor_value) <- corrected[[3]]
    colnames(cor_value) <- unlist(sapply(assaylist, rownames))
    srtIntegrated[["Scanoramacorrected"]] <- CreateAssayObject(data = cor_value)
    SeuratObject::VariableFeatures(srtIntegrated[["Scanoramacorrected"]]) <- HVF

    dim_reduction <- do.call(rbind, corrected[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0(
      "Scanorama_",
      seq_len(ncol(dim_reduction))
    )
  } else {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      verbose = FALSE
    )
    for (nm in names(Scanorama_params)) {
      params[[nm]] <- Scanorama_params[[nm]]
    }
    integrated <- invoke_fun(.fn = scanorama$integrate, .args = params)

    dim_reduction <- do.call(rbind, integrated[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0(
      "Scanorama_",
      seq_len(ncol(dim_reduction))
    )
  }
  srtIntegrated[["Scanorama"]] <- CreateDimReducObject(
    embeddings = dim_reduction,
    key = "Scanorama_",
    assay = SeuratObject::DefaultAssay(srtIntegrated)
  )

  if (is.null(Scanorama_dims_use)) {
    Scanorama_dims_use <- 1:ncol(srtIntegrated[["Scanorama"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "Scanorama",
        dims = Scanorama_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("Scanorama_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "Scanorama_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Scanoramaclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Scanorama",
            reduction_use = "Scanorama",
            reduction_dims = Scanorama_dims_use,
            graph_use = "Scanorama_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[[
    "Scanorama_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|Scanorama|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' BBKNN_integrate
#'
#' @inheritParams Integration_scop
#' @param bbknn_params A list of parameters for the bbknn.matrix.bbknn function, default is an empty list.
#'
#' @export
BBKNN_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    bbknn_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_python("bbknn")
  bbknn <- reticulate::import("bbknn")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(
              Seurat::GetAssayData(
                srt_merge,
                layer = "scale.data",
                assay = SeuratObject::DefaultAssay(srt_merge)
              )
            )
        ))
  ) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "BBKNN",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "BBKNN",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(BBKNN) on the data...\n"
  ))
  message(
    "BBKNN integration using Reduction(",
    paste0("BBKNN", linear_reduction),
    ", dims:",
    min(linear_reduction_dims_use),
    "-",
    max(linear_reduction_dims_use),
    ") as input"
  )
  emb <- Embeddings(srt_merge, reduction = paste0("BBKNN", linear_reduction))[,
    linear_reduction_dims_use,
    drop = FALSE
  ]
  params <- list(
    pca = emb,
    batch_list = srt_merge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke_fun(.fn = bbknn$matrix$bbknn, .args = params)
  n.neighbors <- bem[[3]]$n_neighbors
  srtIntegrated <- srt_merge

  bbknn_graph <- SeuratObject::as.sparse(bem[[2]][1:nrow(bem[[2]]), , drop = FALSE])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- SeuratObject::as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- SeuratObject::DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- Matrix::t(
    SeuratObject::as.sparse(
      bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]
    )
  )
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  bbknn_dist <- SeuratObject::as.Graph(bbknn_dist)
  bbknn_dist@assay.used <- SeuratObject::DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["BBKNN_dist"]] <- bbknn_dist

  val <- split(bbknn_dist@x, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  pos <- split(bbknn_dist@i + 1, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  idx <- Matrix::t(mapply(
    function(x, y) {
      out <- y[head(order(x, decreasing = F), n.neighbors - 1)]
      length(out) <- n.neighbors - 1
      return(out)
    },
    x = val,
    y = pos
  ))
  idx[is.na(idx)] <- sample(
    seq_len(nrow(idx)),
    size = sum(is.na(idx)),
    replace = TRUE
  )
  idx <- cbind(seq_len(nrow(idx)), idx)
  dist <- Matrix::t(mapply(
    function(x, y) {
      out <- y[head(order(x, decreasing = F), n.neighbors - 1)]
      length(out) <- n.neighbors - 1
      out[is.na(out)] <- 0
      return(out)
    },
    x = val,
    y = val
  ))
  dist <- cbind(0, dist)
  srtIntegrated[["BBKNN_neighbors"]] <- new(
    Class = "Neighbor",
    nn.idx = idx,
    nn.dist = dist,
    alg.info = list(),
    cell.names = rownames(emb)
  )
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

  srtIntegrated <- tryCatch(
    {
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        graph.name = "BBKNN",
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["BBKNNclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(
          "Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n",
          sep = ""
        )
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors
        }
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "BBKNN",
            neighbor_use = "BBKNN_neighbors",
            graph_use = "BBKNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["BBKNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|BBKNN|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' CSS_integrate
#'
#' @inheritParams Integration_scop
#' @param CSS_dims_use A vector specifying the dimensions returned by CSS that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param CSS_params A list of parameters for the simspec::cluster_sim_spectrum function, default is an empty list.
#' @export
CSS_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    CSS_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    CSS_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r(c("quadbiolab/simspec", "qlcMatrix"))
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(
              Seurat::GetAssayData(
                srt_merge,
                layer = "scale.data",
                assay = SeuratObject::DefaultAssay(srt_merge)
              )
            )
        ))
  ) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srt_merge <- Seurat::ScaleData(
      object = srt_merge,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srt_merge),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srt_merge <- RunDimReduction(
    srt_merge,
    prefix = "CSS",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srt_merge@reductions[[paste0(
      "CSS",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(CSS) on the data...\n"
  ))
  message(
    "CSS integration using Reduction(",
    paste0("CSS", linear_reduction),
    ", dims:",
    min(linear_reduction_dims_use),
    "-",
    max(linear_reduction_dims_use),
    ") as input"
  )
  params <- list(
    object = srt_merge,
    use_dr = paste0("CSS", linear_reduction),
    dims_use = linear_reduction_dims_use,
    var_genes = HVF,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  for (nm in names(CSS_params)) {
    params[[nm]] <- CSS_params[[nm]]
  }
  srtIntegrated <- invoke_fun(
    .fn = get("cluster_sim_spectrum",
      envir = getNamespace("simspec")
    ),
    .args = params
  )

  if (any(is.na(srtIntegrated@reductions[["CSS"]]@cell.embeddings))) {
    stop(
      "NA detected in the CSS embeddings. You can try to use a lower resolution value in the CSS_param."
    )
  }
  if (is.null(CSS_dims_use)) {
    CSS_dims_use <- 1:ncol(srtIntegrated[["CSS"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "CSS",
        dims = CSS_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("CSS", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "CSS_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["CSSclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(
          "Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n",
          sep = ""
        )
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "CSS",
            reduction_use = "CSS",
            reduction_dims = CSS_dims_use,
            graph_use = "CSS_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["CSS_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|CSS|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' LIGER_integrate
#'
#' @inheritParams Integration_scop
#' @param LIGER_dims_use A vector specifying the dimensions returned by LIGER that will be utilized for downstream cell cluster finding and non-linear reduction. If set to NULL, all the returned dimensions will be used by default.
#' @param optimizeALS_params A list of parameters for the rliger::optimizeALS function, default is an empty list.
#' @param quantilenorm_params A list of parameters for the rliger::quantile_norm function, default is an empty list.
#'
#' @export
LIGER_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    LIGER_dims_use = NULL,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    optimizeALS_params = list(),
    quantilenorm_params = list(),
    seed = 11) {
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("rliger")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srt_list, ncol)) < 30) {
    warning(
      "The cell count in some batches is lower than 30, which may not be suitable for the current integration method.",
      immediate. = TRUE
    )
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srt_merge)
    }
  }

  scale.data <- list()
  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) &&
          any(
            !HVF %in%
              rownames(
                Seurat::GetAssayData(
                  srt,
                  layer = "scale.data",
                  assay = SeuratObject::DefaultAssay(srt)
                )
              )
          ))
    ) {
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform ScaleData on the data ",
        i,
        " ...\n"
      ))
      srt <- Seurat::ScaleData(
        object = srt,
        assay = SeuratObject::DefaultAssay(srt),
        features = HVF,
        do.center = FALSE,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
    scale.data[[i]] <- Matrix::t(
      x = Seurat::GetAssayData(
        object = srt,
        layer = "scale.data",
        assay = SeuratObject::DefaultAssay(srt)
      )
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(LIGER) on the data...\n"
  ))
  params1 <- list(
    object = scale.data,
    k = 20,
    verbose = FALSE
  )
  for (nm in names(optimizeALS_params)) {
    params1[[nm]] <- optimizeALS_params[[nm]]
  }
  out1 <- invoke_fun(.fn = rliger::optimizeALS, .args = params1)
  cat("\n")
  colnames(x = out1$W) <- colnames(scale.data[[1]])
  reduction1 <- do.call(what = "rbind", args = out1$H)
  colnames(reduction1) <- paste0("riNMF_", seq_len(ncol(reduction1)))
  loadings1 <- Matrix::t(x = out1$W)
  rownames(loadings1) <- colnames(scale.data[[1]])
  colnames(loadings1) <- paste0("riNMF_", seq_len(ncol(loadings1)))
  srt_merge[["iNMF_raw"]] <- CreateDimReducObject(
    embeddings = reduction1,
    loadings = loadings1,
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "riNMF_"
  )

  embeddings <- sapply(
    X = SplitObject(object = srt_merge, split.by = batch),
    FUN = function(x) {
      return(Embeddings(object = x[["iNMF_raw"]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  num.samples <- vapply(
    X = embeddings,
    FUN = nrow,
    FUN.VALUE = integer(length = 1L)
  )
  ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  params2 <- list(
    object = embeddings,
    ref_dataset = ref_dataset
  )
  for (nm in names(quantilenorm_params)) {
    params2[[nm]] <- quantilenorm_params[[nm]]
  }
  out2 <- invoke_fun(.fn = rliger::quantile_norm, .args = params2)
  srt_merge[["LIGER"]] <- CreateDimReducObject(
    embeddings = out2$H.norm,
    assay = SeuratObject::DefaultAssay(srt_merge),
    key = "LIGER_"
  )
  srtIntegrated <- srt_merge
  srt_merge <- NULL
  if (is.null(LIGER_dims_use)) {
    LIGER_dims_use <- 1:ncol(srtIntegrated[["LIGER"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = "LIGER",
        dims = LIGER_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("LIGER", "_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "LIGER_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["LIGERclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "LIGER",
            reduction_use = "LIGER",
            reduction_dims = LIGER_dims_use,
            graph_use = "LIGER_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["LIGER_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|LIGER|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Conos_integrate
#'
#' @inheritParams Integration_scop
#' @param buildGraph_params A list of parameters for the buildGraph function, default is an empty list.
#' @param num_threads  An integer setting the number of threads for Conos, default is 2.
#'
#' @importFrom igraph as_adjacency_matrix
#' @export
Conos_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    buildGraph_params = list(),
    num_threads = 2,
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "fr"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive', 'fr'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("conos")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (min(sapply(srt_list, ncol)) < 30) {
    warning(
      "The cell count in some batches is lower than 30, which may not be suitable for the current integration method.",
      immediate. = TRUE
    )
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(srt_merge)
    }
  }

  srtIntegrated <- srt_merge
  srt_merge <- NULL

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    if (
      isTRUE(do_scaling) ||
        (is.null(do_scaling) &&
          any(
            !HVF %in%
              rownames(
                Seurat::GetAssayData(
                  srt,
                  layer = "scale.data",
                  assay = SeuratObject::DefaultAssay(srt)
                )
              )
          ))
    ) {
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform ScaleData on the data ",
        i,
        " ...\n"
      ))
      srt <- Seurat::ScaleData(
        object = srt,
        assay = SeuratObject::DefaultAssay(srt),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " Perform linear dimension reduction (",
      linear_reduction,
      ") on the data ",
      i,
      " ...\n"
    ))
    srt <- RunDimReduction(
      srt,
      prefix = "Conos",
      features = HVF,
      assay = SeuratObject::DefaultAssay(srt),
      linear_reduction = linear_reduction,
      linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params,
      force_linear_reduction = force_linear_reduction,
      verbose = FALSE,
      seed = seed
    )
    srt[["pca"]] <- srt[[paste0("Conos", linear_reduction)]]
    srt_list[[i]] <- srt
  }
  if (is.null(names(srt_list))) {
    names(srt_list) <- paste0("srt_", seq_along(srt_list))
  }

  if (is.null(linear_reduction_dims_use)) {
    maxdims <- max(unlist(sapply(
      srt_list,
      function(srt) {
        max(srt@reductions[[paste0("Conos", linear_reduction)]]@misc[[
          "dims_estimate"
        ]])
      }
    )))
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(Conos) on the data...\n"
  ))
  message(
    "Conos integration using Reduction(",
    linear_reduction,
    ", dims_max:",
    maxdims,
    ") as input"
  )
  srt_list_con <- conos::Conos$new(srt_list, n.cores = num_threads)
  params <- list(
    ncomps = maxdims,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invoke_fun(.fn = srt_list_con[["buildGraph"]], .args = params)
  conos_graph <- as_adjacency_matrix(
    srt_list_con$graph,
    type = "both",
    attr = "weight",
    names = TRUE,
    sparse = TRUE
  )
  conos_graph <- as.Graph(conos_graph)
  conos_graph@assay.used <- SeuratObject::DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

  srtIntegrated <- tryCatch(
    {
      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        graph.name = "Conos",
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Conosclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(
          "Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n",
          sep = ""
        )
        if (nr %in% c("fr")) {
          nonlinear_reduction_params[["n.neighbors"]] <- NULL
        } else {
          nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]
        }
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Conos",
            graph_use = "Conos",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Conos_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|Conos|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Combat_integrate
#'
#' @inheritParams Integration_scop
#' @param ComBat_params A list of parameters for the sva::ComBat function, default is an empty list.
#'
#' @export
ComBat_integrate <- function(
    srt_merge = NULL,
    batch = NULL,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    ComBat_params = list(),
    seed = 11) {
  if (length(linear_reduction) > 1) {
    warning(
      "Only the first method in the 'linear_reduction' will be used.",
      immediate. = TRUE
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop(
      "'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'."
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (
    any(
      !nonlinear_reduction %in%
        c(
          "umap",
          "umap-naive",
          "tsne",
          "dm",
          "phate",
          "pacmap",
          "trimap",
          "largevis",
          "fr"
        )
    )
  ) {
    stop(
      "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
    )
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_r("sva")
  set.seed(seed)

  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("srt_list and srt_merge were all empty.")
  }
  if (!is.null(srt_list) && !is.null(srt_merge)) {
    cell1 <- sort(unique(unlist(lapply(srt_list, colnames))))
    cell2 <- sort(unique(colnames(srt_merge)))
    if (!identical(cell1, cell2)) {
      stop("srt_list and srt_merge have different cells.")
    }
  }
  if (!is.null(srt_merge)) {
    srt_merge_raw <- srt_merge
  } else {
    srt_merge_raw <- NULL
  }
  if (!is.null(srt_list)) {
    checked <- check_srt_list(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- check_srt_merge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0(
      "[",
      Sys.time(),
      "]",
      " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"
    ))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform integration(Combat) on the data...\n"
  ))
  dat <- Seurat::GetAssayData(
    srt_merge,
    layer = "data",
    assay = SeuratObject::DefaultAssay(srt_merge)
  )[
    HVF, ,
    drop = FALSE
  ]
  batch <- srt_merge[[batch, drop = TRUE]]
  params <- list(
    dat = dat,
    batch = batch
  )
  for (nm in names(ComBat_params)) {
    params[[nm]] <- ComBat_params[[nm]]
  }
  corrected <- suppressMessages(
    invoke_fun(
      .fn = sva::ComBat,
      .args = params
    )
  )

  srtIntegrated <- srt_merge
  srt_merge <- NULL
  srtIntegrated[["ComBatcorrected"]] <- CreateAssayObject(data = corrected)
  SeuratObject::DefaultAssay(srtIntegrated) <- "ComBatcorrected"
  SeuratObject::VariableFeatures(srtIntegrated[["ComBatcorrected"]]) <- HVF

  if (
    isTRUE(do_scaling) ||
      (is.null(do_scaling) &&
        any(
          !HVF %in%
            rownames(
              Seurat::GetAssayData(
                srtIntegrated,
                layer = "scale.data",
                assay = SeuratObject::DefaultAssay(srtIntegrated)
              )
            )
        ))
  ) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- Seurat::ScaleData(
      srtIntegrated,
      split.by = if (isTRUE(scale_within_batch)) batch else NULL,
      assay = SeuratObject::DefaultAssay(srtIntegrated),
      features = HVF,
      vars.to.regress = vars_to_regress,
      model.use = regression_model,
      verbose = FALSE
    )
  }

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Perform linear dimension reduction (",
    linear_reduction,
    ") on the data...\n"
  ))
  srtIntegrated <- RunDimReduction(
    srtIntegrated,
    prefix = "ComBat",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = FALSE,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0(
      "ComBat",
      linear_reduction
    )]]@misc[["dims_estimate"]]
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated,
        reduction = paste0("ComBat", linear_reduction),
        dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric,
        k.param = neighbor_k,
        # force.recalc = TRUE,
        graph.name = paste0("ComBat_", c("KNN", "SNN")),
        verbose = FALSE
      )

      cat(paste0(
        "[",
        Sys.time(),
        "]",
        " Perform FindClusters (",
        cluster_algorithm,
        ") on the data...\n"
      ))
      srtIntegrated <- FindClusters(
        object = srtIntegrated,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        method = "igraph",
        graph.name = "ComBat_SNN",
        verbose = FALSE
      )
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- srt_reorder(
        srtIntegrated,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["ComBatclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0(
          "[",
          Sys.time(),
          "]",
          " Perform nonlinear dimension reduction (",
          nr,
          ") on the data...\n"
        ))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "ComBat",
            reduction_use = paste0("ComBat", linear_reduction),
            reduction_dims = linear_reduction_dims_use,
            graph_use = "ComBat_SNN",
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE,
            seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message(
        "Error when performing nonlinear dimension reduction. Skip this step..."
      )
      return(srtIntegrated)
    }
  )

  SeuratObject::DefaultAssay(srtIntegrated) <- assay
  SeuratObject::VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["ComBat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- SrtAppend(
      srt_raw = srt_merge_raw,
      srt_append = srtIntegrated,
      pattern = paste0(assay, "|ComBat|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Integration_scop
#'
#' Integrate single-cell RNA-seq data using various integration methods.
#'
#' @inheritParams check_srt_list
#' @inheritParams check_srt_merge
#' @inheritParams standard_scop
#' @param scale_within_batch  Whether to scale data within each batch. Only valid when the \code{integration_method} is one of \code{"Uncorrected"}, \code{"Seurat"}, \code{"MNN"}, \code{"Harmony"}, \code{"BBKNN"}, \code{"CSS"}, \code{"ComBat"}.
#' @param integration_method  A character string specifying the integration method to use.
#'   Supported methods are: \code{"Uncorrected"}, \code{"Seurat"}, \code{"scVI"}, \code{"MNN"}, \code{"fastMNN"}, \code{"Harmony"},
#'   \code{"Scanorama"}, \code{"BBKNN"}, \code{"CSS"}, \code{"LIGER"}, \code{"Conos"}, \code{"ComBat"}. Default is \code{"Uncorrected"}.
#' @param append Logical, if TRUE, the integrated data will be appended to the original Seurat object (srt_merge).
#' @param ... Additional arguments to be passed to the integration method function.
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link{Seurat_integrate}} \code{\link{scVI_integrate}} \code{\link{MNN_integrate}} \code{\link{fastMNN_integrate}} \code{\link{Harmony_integrate}} \code{\link{Scanorama_integrate}} \code{\link{BBKNN_integrate}} \code{\link{CSS_integrate}} \code{\link{LIGER_integrate}} \code{\link{Conos_integrate}} \code{\link{ComBat_integrate}} \code{\link{standard_scop}}
#'
#' @examples
#' data("panc8_sub")
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Uncorrected",
#'   HVF_min_intersection = 5, scale_within_batch = TRUE
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat",
#'   FindIntegrationAnchors_params = list(reduction = "rpca")
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' \dontrun{
#' integration_methods <- c(
#'   "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony",
#'   "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat"
#' )
#' for (method in integration_methods) {
#'   panc8_sub <- Integration_scop(
#'     srt_merge = panc8_sub, batch = "tech",
#'     integration_method = method,
#'     linear_reduction_dims_use = 1:50,
#'     nonlinear_reduction = "umap"
#'   )
#'   print(CellDimPlot(panc8_sub,
#'     group.by = c("tech", "celltype"),
#'     reduction = paste0(method, "UMAP2D"),
#'     xlab = "", ylab = "", title = method,
#'     legend.position = "none", theme_use = "theme_blank"
#'   ))
#' }
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"
#' )
#' panc8_sub <- Integration_scop(
#'   srt_merge = panc8_sub, batch = "tech",
#'   integration_method = "Seurat",
#'   linear_reduction_dims_use = 1:50,
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' for (nr in nonlinear_reductions) {
#'   print(CellDimPlot(panc8_sub,
#'     group.by = c("tech", "celltype"),
#'     reduction = paste0("Seurat", nr, "2D"),
#'     xlab = "", ylab = "", title = nr,
#'     legend.position = "none", theme_use = "theme_blank"
#'   ))
#' }
#' }
#'
#' @export
Integration_scop <- function(
    srt_merge = NULL,
    batch,
    append = TRUE,
    srt_list = NULL,
    assay = NULL,
    integration_method = "Uncorrected",
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    do_scaling = TRUE,
    vars_to_regress = NULL,
    regression_model = "linear",
    scale_within_batch = FALSE,
    linear_reduction = "pca",
    linear_reduction_dims = 50,
    linear_reduction_dims_use = NULL,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = c(2, 3),
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    neighbor_metric = "euclidean",
    neighbor_k = 20L,
    cluster_algorithm = "louvain",
    cluster_resolution = 0.6,
    seed = 11,
    ...) {
  if (is.null(srt_list) && is.null(srt_merge)) {
    stop("Neither 'srt_list' nor 'srt_merge' was found.")
  }
  if (
    length(integration_method) == 1 &&
      integration_method %in%
        c(
          "Uncorrected",
          "Seurat",
          "scVI",
          "MNN",
          "fastMNN",
          "Harmony",
          "Scanorama",
          "BBKNN",
          "CSS",
          "LIGER",
          "Conos",
          "ComBat"
        )
  ) {
    # Convert the arguments of the function call to a list and remove the function itself
    args <- as.list(match.call())[-1]

    # Create a new environment, the parent of which is the environment that called the foo function
    new_env <- new.env(parent = parent.frame())

    # Evaluate the arguments in the new environment to get the correct values
    args <- lapply(args, function(x) eval(x, envir = new_env))

    # Keep srt_merge and srt_list as type of 'symbol' when use `do.call` function
    # args[!names(args) %in% c("srt_merge", "srt_list")] <- lapply(args[!names(args) %in% c("srt_merge", "srt_list")], function(x) eval(x, envir = new_env))

    # print("================ args ================ ")
    # print(args)

    # Get the function's formal arguments and their default values
    formals <- mget(names(formals()))
    formals <- formals[names(formals) != "..."]

    # print("================ formals ================ ")
    # print(formals)

    # Merge the formal arguments with the actual arguments, so that all arguments are included
    args <- modifyList(formals, args)

    time_start <- Sys.time()
    cat(paste0(
      "[",
      time_start,
      "] ",
      paste0("Start ", integration_method, "_integrate"),
      "\n"
    ))
    srtIntegrated <- invoke_fun(
      .fn = paste0(integration_method, "_integrate"),
      .args = args[
        names(args) %in% formalArgs(paste0(integration_method, "_integrate"))
      ]
    )
    time_end <- Sys.time()
    cat(paste0(
      "[",
      time_end,
      "] ",
      paste0(integration_method, "_integrate done\n")
    ))
    cat(
      "Elapsed time:",
      format(
        round(difftime(time_end, time_start), 2),
        format = "%Y-%m-%d %H:%M:%S"
      ),
      "\n"
    )

    return(srtIntegrated)
  } else {
    stop(paste(integration_method, "is not a supported integration method!"))
  }
}
