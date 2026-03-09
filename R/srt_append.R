#' @title Append a Seurat object to another
#'
#' @md
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
