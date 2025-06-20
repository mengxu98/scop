#' Convert an anndata object to a seurat object using reticulate
#'
#' @param adata a connected python anndata object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
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
#' ### Or convert a h5ad file to Seurat object
#' # library(reticulate)
#' # check_python("scanpy")
#' # sc <- import("scanpy")
#' # adata <- sc$read_h5ad("pancreas.h5ad")
#' # srt <- adata_to_srt(adata)
#' # srt
#' }
adata_to_srt <- function(adata) {
  if (!inherits(adata, "python.builtin.object")) {
    log_message(
      "'adata' is not a python.builtin.object.",
      message_type = "error"
    )
  }
  sc <- reticulate::import("scanpy", convert = TRUE)
  np <- reticulate::import("numpy", convert = TRUE)

  x <- Matrix::t(py_to_r_auto(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- SeuratObject::as.sparse(x[1:nrow(x), , drop = FALSE])
  }
  rownames(x) <- py_to_r_auto(adata$var_names$values)
  colnames(x) <- py_to_r_auto(adata$obs_names$values)
  rownames(x) <- as.character(rownames(x))
  colnames(x) <- as.character(colnames(x))

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(py_to_r_auto(adata$obs))
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
      layer <- py_to_r_auto(adata$layers[[k]])
      if (!inherits(layer, c("Matrix", "matrix"))) {
        log_message(
          paste0(
            "The object in '",
            k,
            "' layers is not a matrix: ",
            paste0(class(adata$layers[[k]]), collapse = ",")
          ),
          message_type = "error"
        )
      }
      layer <- Matrix::t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- SeuratObject::as.sparse(layer[1:nrow(layer), , drop = FALSE])
      }
      rownames(layer) <- py_to_r_auto(adata$var_names$values)
      colnames(layer) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- Seurat::CreateAssayObject(counts = layer)
    }
  }

  if (inherits(adata$obsm, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$obsm$keys())
  } else {
    keys <- names(adata$obsm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsm <- tryCatch(py_to_r_auto(adata$obsm[[k]]), error = identity)
      if (inherits(obsm, "error")) {
        log_message(
          "'obsm: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- Matrix::as.matrix(obsm)
      }
      k <- gsub(pattern = "^X_", replacement = "", x = py_to_r_auto(k))
      colnames(obsm) <- paste0(k, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- Seurat::CreateDimReducObject(
        embeddings = obsm,
        assay = "RNA",
        key = paste0(gsub(pattern = "_", replacement = "", x = k), "_")
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
      obsp <- tryCatch(py_to_r_auto(adata$obsp[[k]]), error = identity)
      if (inherits(obsp, "error")) {
        log_message(
          "'obsp: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- SeuratObject::as.sparse(obsp[1:nrow(obsp), , drop = FALSE])
      }
      colnames(obsp) <- py_to_r_auto(adata$obs_names$values)
      rownames(obsp) <- py_to_r_auto(adata$obs_names$values)
      obsp <- SeuratObject::as.Graph(obsp[seq_len(nrow(obsp)), , drop = FALSE])
      DefaultAssay(object = obsp) <- "RNA"
      srt[[py_to_r_auto(k)]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- Seurat::AddMetaData(
      srt[["RNA"]],
      metadata = as.data.frame(py_to_r_auto(adata$var))
    )
  }

  if (inherits(adata$varm, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$varm$keys())
  } else {
    keys <- names(adata$varm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varm <- tryCatch(py_to_r_auto(adata$varm[[k]]), error = identity)
      if (inherits(varm, "error")) {
        log_message(
          "'varm: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- Matrix::as.matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r_auto(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.loadings"]][[py_to_r_auto(k)]] <- varm
    }
  }

  if (inherits(adata$varp, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$varp$keys())
  } else {
    keys <- names(adata$varp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varp <- tryCatch(py_to_r_auto(adata$varp[[k]]), error = identity)
      if (inherits(varp, "error")) {
        log_message(
          "'varp: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- Matrix::as.matrix(varp)
      }
      colnames(varp) <- py_to_r_auto(adata$var_names$values)
      rownames(varp) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.graphs"]][[py_to_r_auto(k)]] <- varp
    }
  }

  if (inherits(adata$uns, "python.builtin.object")) {
    keys <- reticulate::iterate(adata$uns$keys())
  } else {
    keys <- names(adata$uns)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      uns <- tryCatch(py_to_r_auto(adata$uns[[k]]), error = identity)
      if (inherits(uns, "error")) {
        log_message(
          "'uns: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        log_message(
          "'uns: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r_auto(k)]] <- uns
      } else {
        log_message(
          "'uns: ",
          k,
          "' will not be converted. You may need to convert it manually.",
          message_type = "warning"
        )
        next
      }
    }
  }
  return(srt)
}

py_to_r_auto <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    x <- reticulate::py_to_r(x)
  }
  return(x)
}
