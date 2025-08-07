#' @title Convert a Seurat object to an AnnData object
#'
#' @description
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @md
#' @param srt A Seurat object.
#' @param assay_x Assay to convert as the main data matrix in the anndata object.
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
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(
          as.character(var[, i]),
          levels = c("TRUE", "FALSE")
        )
      }
    }
  }
  var_features <- SeuratObject::VariableFeatures(
    srt,
    assay = assay_x
  )
  if (length(var_features) > 0) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in% var_features
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

#' @title Convert an anndata object to a seurat object using reticulate
#'
#' @param adata A connected python anndata object.
#' @param verbose Whether to print messages.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(
#'   adata = adata,
#'   group_by = "SubCellType",
#'   linear_reduction = "PCA",
#'   nonlinear_reduction = "UMAP"
#' )
#' srt <- adata_to_srt(adata)
#' srt
#'
#' ## Or convert a h5ad file to Seurat object
#' # library(reticulate)
#' # check_python("scanpy")
#' # sc <- import("scanpy")
#' # adata <- sc$read_h5ad("pancreas.h5ad")
#' # srt <- adata_to_srt(adata)
#' # srt
#' }
adata_to_srt <- function(
    adata,
    verbose = TRUE) {
  if (!inherits(adata, "python.builtin.object")) {
    log_message(
      "{.val adata} is not a python.builtin.object.",
      message_type = "error"
    )
  }
  log_message(
    "Converting {.cls AnnData} object to {.cls Seurat} object...",
    verbose = verbose
  )

  x <- Matrix::t(py_to_r2(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- SeuratObject::as.sparse(x)
  }
  rownames(x) <- py_to_r2(adata$var_names$values)
  colnames(x) <- py_to_r2(adata$obs_names$values)
  rownames(x) <- as.character(rownames(x))
  colnames(x) <- as.character(colnames(x))

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(py_to_r2(adata$obs))
    colnames(metadata) <- make.names(colnames(metadata))
  }

  srt <- Seurat::CreateSeuratObject(
    counts = x,
    meta.data = metadata
  )

  if (inherits(adata$layers, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$layers$keys())
  } else {
    keys <- names(adata$layers)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      layer <- py_to_r2(adata$layers[[k]])
      if (!inherits(layer, c("Matrix", "matrix"))) {
        log_message(
          paste0(
            "The object in {.val {k}} layers is not a matrix: ",
            paste0(class(adata$layers[[k]]), collapse = ",")
          ),
          message_type = "error"
        )
      }
      layer <- Matrix::t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- SeuratObject::as.sparse(layer)
      }
      rownames(layer) <- py_to_r2(adata$var_names$values)
      colnames(layer) <- py_to_r2(adata$obs_names$values)
      srt[[py_to_r2(k)]] <- Seurat::CreateAssayObject(counts = layer)
    }
  }

  if (inherits(adata$obsm, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$obsm$keys())
  } else {
    keys <- names(adata$obsm)
  }
  if (length(keys) > 0) {
    processed_reductions <- character(0)
    for (k in keys) {
      k_clean <- gsub(pattern = "^X_", replacement = "", x = py_to_r2(k))

      if (k_clean %in% processed_reductions) {
        next
      }

      processed_reductions <- c(processed_reductions, k_clean)
      obsm <- tryCatch(py_to_r2(adata$obsm[[k]]), error = identity)
      if (inherits(obsm, "error")) {
        log_message(
          "{.val obsm}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- Matrix::as.matrix(obsm)
      }
      colnames(obsm) <- paste0(k_clean, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- py_to_r2(adata$obs_names$values)
      srt[[py_to_r2(k)]] <- Seurat::CreateDimReducObject(
        embeddings = obsm,
        assay = "RNA",
        key = paste0(gsub(pattern = "_", replacement = "", x = k_clean), "_")
      )
    }
  }

  if (inherits(adata$obsp, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$obsp$keys())
  } else {
    keys <- names(adata$obsp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsp <- tryCatch(py_to_r2(adata$obsp[[k]]), error = identity)
      if (inherits(obsp, "error")) {
        log_message(
          "{.val obsp}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- SeuratObject::as.sparse(obsp)
      }
      colnames(obsp) <- py_to_r2(adata$obs_names$values)
      rownames(obsp) <- py_to_r2(adata$obs_names$values)
      obsp <- SeuratObject::as.Graph(obsp)
      SeuratObject::DefaultAssay(obsp) <- "RNA"
      srt[[py_to_r2(k)]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- Seurat::AddMetaData(
      srt[["RNA"]],
      metadata = as.data.frame(py_to_r2(adata$var))
    )
  }

  if (inherits(adata$varm, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$varm$keys())
  } else {
    keys <- names(adata$varm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varm <- tryCatch(py_to_r2(adata$varm[[k]]), error = identity)
      if (inherits(varm, "error")) {
        log_message(
          "{.val varm}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- Matrix::as.matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r2(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- py_to_r2(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.loadings"]][[py_to_r2(k)]] <- varm
    }
  }

  if (inherits(adata$varp, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$varp$keys())
  } else {
    keys <- names(adata$varp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varp <- tryCatch(py_to_r2(adata$varp[[k]]), error = identity)
      if (inherits(varp, "error")) {
        log_message(
          "{.val varp}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- Matrix::as.matrix(varp)
      }
      colnames(varp) <- py_to_r2(adata$var_names$values)
      rownames(varp) <- py_to_r2(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.graphs"]][[py_to_r2(k)]] <- varp
    }
  }

  if (inherits(adata$uns, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$uns$keys())
  } else {
    keys <- names(adata$uns)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      uns <- tryCatch(py_to_r2(adata$uns[[k]]), error = identity)
      if (inherits(uns, "error")) {
        log_message(
          "{.val uns}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        log_message(
          "{.val uns}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r2(k)]] <- uns
      } else {
        log_message(
          "{.val uns}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
    }
  }
  log_message(
    "Convert {.cls AnnData} object to {.cls Seurat} object completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}

py_to_r2 <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    reticulate::py_to_r(x)
  } else {
    x
  }
}
