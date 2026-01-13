#' @title Convert a Seurat object to an AnnData object
#'
#' @description
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @param assay_x Assay to convert as the main data matrix in the anndata object.
#' Default is `"RNA"`.
#' @param layer_x Layer name for assay_x in the Seurat object.
#' Default is `"counts"`.
#' @param assay_y Assays to convert as layers in the anndata object.
#' Default is `c("spliced", "unspliced")`.
#' @param layer_y Layer names for the assay_y in the Seurat object.
#' Default is `"counts"`.
#' @param convert_tools Whether to convert the tool-specific data.
#' Default is `FALSE`.
#' @param convert_misc Whether to convert the miscellaneous data.
#' Default is `FALSE`.
#' @param features Optional vector of features to include in the anndata object.
#' Default is all features in `assay_x`.
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
#' # Or save as a h5ad/loom file
#' adata$write_h5ad(
#'   "pancreas_sub.h5ad"
#' )
#' adata$write_loom(
#'   "pancreas_sub.loom",
#'   write_obsm_varm = TRUE
#' )
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
  PrepareEnv()
  check_python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  log_message(
    "Converting {.cls Seurat} to {.cls AnnData} ...",
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
    reduction_list[[reduction]] <- srt[[reduction]]@cell.embeddings
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
    "Convert {.cls Seurat} to {.cls AnnData} object completed",
    message_type = "success",
    verbose = verbose
  )

  return(adata)
}

#' @title Convert an anndata object to a seurat object using reticulate
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param adata An AnnData object. Can be a Python AnnData object (from scanpy/reticulate),
#'   an R6 AnnData object from the `anndata` package (AnnDataR6), or an R6 AnnData object
#'   from the `anndataR` package (InMemoryAnnData).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(
#'   adata = adata,
#'   group_by = "SubCellType",
#'   linear_reduction = "X_pca",
#'   nonlinear_reduction = "X_umap"
#' )
#' srt <- adata_to_srt(adata)
#' srt
#'
#' # Or convert a h5ad file to Seurat object
#' sc <- reticulate::import("scanpy")
#' adata <- sc$read_h5ad("pancreas.h5ad")
#' srt <- adata_to_srt(adata)
#' srt
#' }
adata_to_srt <- function(
    adata,
    verbose = TRUE) {
  PrepareEnv()
  data_types <- c(
    "python.builtin.object", "AnnDataR6", "InMemoryAnnData", "AbstractAnnData"
  )
  if (!inherits(adata, data_types)) {
    log_message(
      "{.val adata} must be one of the following classes: {.val {data_types}}",
      message_type = "error"
    )
  }
  log_message(
    "Converting {.cls {class(adata)}} object to {.cls Seurat}...",
    verbose = verbose
  )

  x <- Matrix::t(py_to_r2(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- SeuratObject::as.sparse(x)
  }
  rownames(x) <- get_adata_names(adata, "var")
  colnames(x) <- get_adata_names(adata, "obs")
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
      layer <- py_to_r2(get_adata_element(adata$layers, k))
      if (!inherits(layer, c("Matrix", "matrix"))) {
        log_message(
          "The object in {.val {k}} layers is not a matrix: {.val {class(get_adata_element(adata$layers, k))}}",
          message_type = "error"
        )
      }
      layer <- Matrix::t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- SeuratObject::as.sparse(layer)
      }
      rownames(layer) <- get_adata_names(adata, "var")
      colnames(layer) <- get_adata_names(adata, "obs")
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
      k_clean <- py_to_r2(k)

      if (k_clean %in% processed_reductions) {
        next
      }

      processed_reductions <- c(processed_reductions, k_clean)
      obsm <- tryCatch(
        py_to_r2(get_adata_element(adata$obsm, k)),
        error = identity
      )
      if (inherits(obsm, "error")) {
        log_message(
          "{.val obsm}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- as_matrix(obsm)
      }
      colnames(obsm) <- paste0(k_clean, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- get_adata_names(adata, "obs")
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
      obsp <- tryCatch(
        py_to_r2(get_adata_element(adata$obsp, k)),
        error = identity
      )
      if (inherits(obsp, "error")) {
        log_message(
          "{.val obsp}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- SeuratObject::as.sparse(obsp)
      }
      colnames(obsp) <- get_adata_names(adata, "obs")
      rownames(obsp) <- get_adata_names(adata, "obs")
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
      varm <- tryCatch(
        py_to_r2(get_adata_element(adata$varm, k)),
        error = identity
      )
      if (inherits(varm, "error")) {
        log_message(
          "{.val varm}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- as_matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r2(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- get_adata_names(adata, "var")
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
      varp <- tryCatch(
        py_to_r2(get_adata_element(adata$varp, k)),
        error = identity
      )
      if (inherits(varp, "error")) {
        log_message(
          "{.val varp}: {.val {k}} will not be converted.",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- as_matrix(varp)
      }
      colnames(varp) <- get_adata_names(adata, "var")
      rownames(varp) <- get_adata_names(adata, "var")
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
      uns <- tryCatch(
        py_to_r2(get_adata_element(adata$uns, k)),
        error = identity
      )
      if (inherits(uns, "error")) {
        log_message(
          "{.val uns}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        log_message(
          "{.val uns}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r2(k)]] <- uns
      } else {
        log_message(
          "{.val uns}: {.val {k}} will not be converted",
          message_type = "warning",
          verbose = verbose
        )
        next
      }
    }
  }
  log_message(
    "Convert {.cls AnnData} object to {.cls Seurat} completed",
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

get_adata_element <- function(container, key) {
  if (is.function(container$get)) {
    result <- tryCatch(
      {
        container$get(key)
      },
      error = function(e) NULL
    )
    if (!is.null(result)) {
      return(result)
    }
  }
  tryCatch(
    {
      return(container[[key]])
    },
    error = function(e) {
      log_message(
        "Cannot access element: {.val {key}}",
        message_type = "error"
      )
    }
  )
}

get_adata_names <- function(adata, name_type = c("var", "obs")) {
  name_type <- match.arg(name_type)
  attr_name <- paste0(name_type, "_names")

  names_obj <- adata[[attr_name]]
  if (is.null(names_obj)) {
    log_message(
      "Cannot access element: {.val {attr_name}}",
      message_type = "error"
    )
  }

  names_r <- py_to_r2(names_obj)

  if (is.vector(names_r) || is.character(names_r)) {
    return(names_r)
  }

  if (is.list(names_r) && "values" %in% names(names_r)) {
    return(py_to_r2(names_r$values))
  }

  if (inherits(names_obj, "python.builtin.object")) {
    tryCatch(
      {
        return(py_to_r2(names_obj$values))
      },
      error = function(e) {
        return(names_r)
      }
    )
  }

  return(names_r)
}

check_python_element <- function(
    x,
    depth = max_depth(x)) {
  if (depth == 0 || !is.list(x) || !inherits(x, "python.builtin.object")) {
    if (inherits(x, "python.builtin.object")) {
      x_r <- tryCatch(
        py_to_r2(x),
        error = identity
      )
      if (inherits(x_r, "error")) {
        return(x)
      } else {
        return(x_r)
      }
    } else {
      return(x)
    }
  } else {
    raw_depth <- max_depth(x)
    x <- lapply(
      x, function(element) {
        if (inherits(element, "python.builtin.object")) {
          element_r <- tryCatch(
            py_to_r2(element),
            error = identity
          )
          if (inherits(element_r, "error")) {
            return(element)
          } else {
            return(element_r)
          }
        } else {
          return(element)
        }
      }
    )
    cur_depth <- max_depth(x)
    if (cur_depth > raw_depth) {
      depth <- depth + 1
    }
    x_checked <- lapply(x, check_python_element, depth - 1)
    return(x_checked)
  }
}
