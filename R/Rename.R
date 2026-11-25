#' @title Rename features for the Seurat object
#'
#' @inheritParams RunStandardWorkflow
#' @param newnames One name per feature in the first selected assay, or a character
#'   vector named by old feature IDs. Unnamed vectors require the selected assays
#'   to have the same feature IDs.
#' @param assays Assays to rename (all assays by default). Standard Seurat v3 and
#'   v5 assays are supported. Rename features before fitting reductions that rely
#'   on those features.
#' @details Previously computed results stored in the object are not rewritten.
#'
#' @export
#'
#' @examples
#' data(panc8_sub)
#' head(rownames(panc8_sub))
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   thisutils::capitalize(rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_rename <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' head(rownames(panc8_rename))
RenameFeatures <- function(
  object,
  newnames = NULL,
  assays = NULL,
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  available_assays <- SeuratObject::Assays(srt)
  assays <- intersect(assays %||% available_assays, available_assays)
  if (!length(assays)) {
    return(srt)
  }
  if (!is.character(newnames) || !length(newnames) ||
    anyNA(newnames) || any(!nzchar(newnames))) {
    log_message("{.arg newnames} must contain non-empty feature names.", message_type = "error")
  }

  if (is.null(names(newnames))) {
    features <- rownames(srt[[assays[[1L]]]])
    if (length(newnames) != length(features)) {
      log_message(
        "Supply one name per feature or a named old-to-new mapping.",
        message_type = "error"
      )
    }
    if (!all(vapply(assays, function(assay) {
      setequal(features, rownames(srt[[assay]]))
    }, logical(1)))) {
      log_message(
        "Selected assays have different feature IDs. Supply a named old-to-new mapping.",
        message_type = "error"
      )
    }
    names(newnames) <- features
  } else if (anyNA(names(newnames)) || any(!nzchar(names(newnames))) ||
    anyDuplicated(names(newnames))) {
    log_message(
      "Old feature names must be non-empty and unique.",
      message_type = "error"
    )
  }

  rename <- function(ids) {
    mapping <- match(ids, names(newnames))
    found <- !is.na(mapping)
    ids[found] <- unname(newnames[mapping[found]])
    ids
  }
  for (assay in assays) {
    assay_obj <- srt[[assay]]
    features <- rownames(assay_obj)
    renamed <- rename(features)
    changed <- features[features != renamed]
    if (!length(changed)) {
      next
    }
    if (anyDuplicated(renamed)) {
      log_message(
        "Renamed features must be unique within assay {.val {assay}}.",
        message_type = "error"
      )
    }

    assay_class <- class(assay_obj)[[1L]]
    if (!assay_class %in% c("Assay", "Assay5")) {
      log_message(
        "Renaming features in assay class {.val {assay_class}} is not supported.",
        message_type = "error"
      )
    }
    for (reduction in SeuratObject::Reductions(srt)) {
      reduction_obj <- srt[[reduction]]
      if (!identical(SeuratObject::DefaultAssay(reduction_obj), assay)) {
        next
      }
      loading_features <- c(
        rownames(SeuratObject::Loadings(reduction_obj)),
        rownames(SeuratObject::Loadings(reduction_obj, projected = TRUE))
      )
      if (any(changed %in% loading_features)) {
        log_message(
          "Reduction {.val {reduction}} contains feature loadings for renamed features. Rename before computing the reduction.",
          message_type = "error"
        )
      }
    }
    log_message("Rename features for the assay: ", assay, verbose = verbose)
    if (identical(assay_class, "Assay")) {
      for (layer in c("meta.features", "scale.data", "counts", "data")) {
        data <- methods::slot(assay_obj, layer)
        if (length(rownames(data))) {
          rownames(data) <- rename(rownames(data))
          methods::slot(assay_obj, layer) <- data
        }
      }
      methods::slot(assay_obj, "var.features") <- rename(
        methods::slot(assay_obj, "var.features")
      )
    } else {
      # The assay-wide feature map supplies names to each layer independently.
      dimnames(assay_obj) <- list(renamed, colnames(assay_obj))
      metadata <- methods::slot(assay_obj, "meta.data")
      if ("var.features" %in% names(metadata)) {
        metadata$var.features <- rename(metadata$var.features)
        methods::slot(assay_obj, "meta.data") <- metadata
      }
    }
    srt[[assay]] <- assay_obj
  }
  return(srt)
}

#' @title Rename clusters for the Seurat object
#'
#' @inheritParams RunStandardWorkflow
#' @param group.by The old group used to rename cells.
#' @param nameslist A named list of new cluster value.
#' @param name New cluster stored in the Seurat object.
#' @param keep_levels If the old group is a factor, keep the order of the levels.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#'
#' # Rename all clusters
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' levels(pancreas_sub@meta.data[["SubCellType"]]) <- unique(
#'   pancreas_sub@meta.data[["SubCellType"]]
#' )
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
  object,
  group.by,
  nameslist = list(),
  name = "newclusters",
  keep_levels = FALSE,
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (missing(group.by)) {
    log_message(
      "group.by must be provided",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      paste0(group.by, " is not in the meta.data of srt object."),
      message_type = "error"
    )
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
        log_message(
          "'nameslist' must be named when srt@meta.data[[group.by]] is not a factor",
          message_type = "error"
        )
      }
      if (
        !identical(length(nameslist), length(unique(srt@meta.data[[group.by]])))
      ) {
        log_message(
          "'nameslist' must be named or the length of ",
          length(unique(srt@meta.data[[group.by]])),
          message_type = "error"
        )
      }
      names(nameslist) <- levels(srt@meta.data[[group.by]])
    }
    names_assign <- nameslist
  }
  if (all(!names(names_assign) %in% srt@meta.data[[group.by]])) {
    log_message(
      "No group name mapped.",
      message_type = "error"
    )
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
