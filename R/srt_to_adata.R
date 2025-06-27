#' Convert a seurat object to an anndata object using reticulate
#'
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @param srt A Seurat object.
#' @param assay_x Assay to convert as the main data matrix (X) in the anndata object.
#' @param layer_x Layer name for assay_x in the Seurat object.
#' @param assay_y Assays to convert as layers in the anndata object.
#' @param layer_y Layer names for the assay_y in the Seurat object.
#' @param convert_tools Logical indicating whether to convert the tool-specific data.
#' @param convert_misc Logical indicating whether to convert the miscellaneous data.
#' @param features Optional vector of features to include in the anndata object.
#' Defaults to all features in assay_x.
#' @param verbose Logical indicating whether to print verbose messages during the conversion process.
#'
#' @return A \code{anndata} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ## Or save as an h5ad file or a loom file
#' # adata$write_h5ad(
#' #   "pancreas_sub.h5ad"
#' # )
#' # adata$write_loom(
#' #   "pancreas_sub.loom",
#' #   write_obsm_varm = TRUE
#' # )
#' }
srt_to_adata <- function(
    srt,
    features = NULL,
    assay_x = "RNA",
    layer_x = "counts",
    assay_y = c("spliced", "unspliced"),
    layer_y = "counts",
    convert_tools = FALSE,
    convert_misc = FALSE,
    verbose = TRUE) {
  check_python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    log_message(
      "'srt' is not a Seurat object.",
      message_type = "error"
    )
  }

  if (is.null(features)) {
    features <- rownames(srt[[assay_x]])
  }
  if (length(layer_y) == 1) {
    layer_y <- rep(layer_y, length(assay_y))
    names(layer_y) <- assay_y
  } else if (length(layer_y) != length(assay_y)) {
    log_message(
      "layer_y must be one character or the same length of the assay_y",
      message_type = "error"
    )
  }

  sc <- reticulate::import("scanpy", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(
          as.character(obs[, i]),
          levels = c("TRUE", "FALSE")
        )
      }
    }
  }

  var <- GetFeaturesData(srt, assay = assay_x)[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (
        is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")
      ) {
        var[, i] <- factor(
          as.character(var[, i]),
          levels = c("TRUE", "FALSE")
        )
      }
    }
  }
  if (length(SeuratObject::VariableFeatures(srt, assay = assay_x) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in%
      SeuratObject::VariableFeatures(srt, assay = assay_x)
  }

  X <- Matrix::t(
    GetAssayData5(
      srt,
      assay = assay_x,
      layer = layer_x
    )[features, , drop = FALSE]
  )
  adata <- sc$AnnData(
    X = reticulate::np_array(X, dtype = np$float32),
    obs = obs,
    var = cbind(
      data.frame(features = features),
      var
    )
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_x]) {
    if (assay %in% assay_y) {
      layer <- Matrix::t(
        GetAssayData5(
          srt,
          assay = assay,
          layer = layer_y[assay]
        )
      )
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          log_message(
            "The following features in the '",
            assay_x,
            "' assay can not be found in the '",
            assay,
            "' assay:\n  ",
            paste0(
              utils::head(colnames(X)[!colnames(X) %in% colnames(layer)], 10),
              collapse = ","
            ),
            "...",
            message_type = "warning"
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      log_message(
        paste0("Assay '", assay, "' is in the srt object but not converted."),
        message_type = "warning",
        verbose = verbose
      )
    }
  }
  if (length(layer_list) > 0) {
    adata$layers <- layer_list
  }

  reduction_list <- list()
  for (reduction in names(srt@reductions)) {
    reduction_list[[paste0(reduction)]] <- srt[[reduction]]@cell.embeddings
    reduction_list[[paste0("X_", reduction)]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  obsp_list <- list()
  for (graph in names(srt@graphs)) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in names(srt@neighbors)) {
    obsp_list[[neighbor]] <- srt[[neighbor]]
  }
  if (length(obsp_list) > 0) {
    adata$obsp <- obsp_list
  }

  uns_list <- list()
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@misc[[nm]]
      }
    }
  } else {
    log_message(
      "'misc' slot is not converted.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    log_message(
      "'tools' slot is not converted.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}
