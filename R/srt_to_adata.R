#' @title Convert a Seurat object to an AnnData object
#'
#' @description
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @md
#' @param srt A Seurat object.
#' @param assay_x Assay to convert as the main data matrix (X) in the anndata object.
#' @param layer_x Layer name for assay_x in the Seurat object.
#' @param assay_y Assays to convert as layers in the anndata object.
#' @param layer_y Layer names for the assay_y in the Seurat object.
#' @param convert_tools Whether to convert the tool-specific data. Default is `FALSE`.
#' @param convert_misc Whether to convert the miscellaneous data. Default is `FALSE`.
#' @param features Optional vector of features to include in the anndata object.
#' Defaults to all features in assay_x.
#' @param verbose Whether to print verbose messages during the conversion process. Default is `TRUE`.
#'
#' @return A `anndata` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ## Or save as a h5ad/loom file
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
      "{.arg srt} is not a Seurat object",
      message_type = "error"
    )
  }
  log_message(
    "Converting {.cls Seurat} object to {.cls AnnData} object...",
    verbose = verbose
  )

  if (is.null(features)) {
    features <- rownames(srt[[assay_x]])
  }
  if (length(layer_y) == 1) {
    layer_y <- rep(layer_y, length(assay_y))
    names(layer_y) <- assay_y
  } else if (length(layer_y) != length(assay_y)) {
    log_message(
      "{.arg layer_y} must be one character or the same length of the {.arg assay_y}",
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
      layer = layer_x,
      verbose = FALSE
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
          layer = layer_y[assay],
          verbose = FALSE
        )
      )
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          features_null <- colnames(X)[!colnames(X) %in% colnames(layer)]
          log_message(
            "The following features in the {.val {assay_x}} are not found in the {.val {assay}}: {.val {features_null}}",
            message_type = "warning",
            verbose = verbose
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      log_message(
        "{.val {assay}} is in the srt object but not converted",
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
      "{.val misc} slot is not converted",
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
      "{.val tools} slot is not converted",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }
  log_message(
    "Convert {.cls Seurat} object to {.cls AnnData} object completed",
    message_type = "success",
    verbose = verbose
  )

  return(adata)
}
