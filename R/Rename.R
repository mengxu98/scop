#' @title Rename features for the Seurat object
#'
#' @inheritParams standard_scop
#' @param newnames A vector with the same length of features in Seurat object,
#' or characters named with old features.
#' @param assays Assays to rename.
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
    srt,
    newnames = NULL,
    assays = NULL) {
  assays <- assays[assays %in% SeuratObject::Assays(srt)] %||% SeuratObject::Assays(srt)
  if (is.null(names(newnames))) {
    if (length(newnames) == nrow(srt)) {
      log_message(
        "'newnames' must be named or the length of features in the srt.",
        message_type = "error"
      )
    }
    if (length(unique(sapply(srt@assays[assays], nrow))) > 1) {
      log_message(
        "Assays in the srt object have different number of features. Please use a named vectors.",
        message_type = "error"
      )
    }
    names(newnames) <- rownames(srt[[assays[1]]])
  }
  for (assay in assays) {
    log_message("Rename features for the assay: ", assay)
    assay_obj <- Seurat::GetAssay(srt, assay = assay)
    if (inherits(assay_obj, "Assay")) {
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
    }

    if (inherits(assay_obj, "Assay5")) {
      meta.data <- methods::slot(assay_obj, "meta.data")
      if ("var.features" %in% names(meta.data)) {
        index <- which(
          meta.data$var.features %in% names(newnames)
        )
        meta.data$var.features[index] <- newnames[meta.data$var.features[index]]
        methods::slot(assay_obj, "meta.data") <- meta.data
      }

      index_all <- which(
        rownames(assay_obj) %in% names(newnames)
      )
      rownames(assay_obj) <- newnames[rownames(assay_obj)[index_all]]
    }

    srt[[assay]] <- assay_obj
  }
  return(srt)
}

#' Rename clusters for the Seurat object
#'
#' @inheritParams standard_scop
#' @param group.by The old group used to rename cells.
#' @param nameslist A named list of new cluster value.
#' @param name The name of the new cluster stored in the Seurat object.
#' @param keep_levels If the old group is a factor, keep the order of the levels.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#'
#' # Rename all clusters
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
    srt,
    group.by,
    nameslist = list(),
    name = "newclusters",
    keep_levels = FALSE) {
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
