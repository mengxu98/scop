#' @title Convert a Seurat object to an AnnData object
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
#' @param reductions Character vector specifying which Seurat reductions to
#' convert into `obsm`. Default is `NULL`, which converts all available
#' reductions.
#' @param graphs Character vector specifying which Seurat graphs to convert into
#' `obsp`. Default is `NULL`, which converts all available graphs.
#' @param neighbors Character vector specifying which Seurat neighbor objects to
#' convert into `obsp`. Default is `NULL`, which converts all available neighbor
#' objects.
#' @param convert_tools Whether to convert the tool-specific data.
#' Default is `FALSE`.
#' @param convert_misc Whether to convert the miscellaneous data.
#' Default is `FALSE`.
#' @param features Optional vector of features to include in the anndata object.
#' Default is all features in `assay_x`.
#'
#' @return A `anndata` object.
#'
#' @seealso [adata_to_srt]
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
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  convert_tools = FALSE,
  convert_misc = FALSE,
  verbose = TRUE
) {
  if (!isTRUE(getOption("scop_skip_python_prepare", FALSE))) {
    old_log_verbose <- getOption("log_message.verbose", TRUE)
    if (!isTRUE(verbose)) {
      options(log_message.verbose = FALSE)
      on.exit(
        options(log_message.verbose = old_log_verbose),
        add = TRUE
      )
    }
    PrepareEnv(modules = "scanpy")
    check_python(
      c("anndata", "numpy"),
      verbose = FALSE
    )
  }

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

  ad <- reticulate::import("anndata", convert = FALSE)
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
  adata <- ad$AnnData(
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

  reduction_names <- reductions %||% names(srt@reductions)
  reduction_names <- intersect(reduction_names, names(srt@reductions))
  reduction_list <- list()
  for (reduction in reduction_names) {
    reduction_list[[reduction]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  graph_names <- graphs %||% names(srt@graphs)
  graph_names <- intersect(graph_names, names(srt@graphs))
  neighbor_names <- neighbors %||% names(srt@neighbors)
  neighbor_names <- intersect(neighbor_names, names(srt@neighbors))
  obsp_list <- list()
  for (graph in graph_names) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in neighbor_names) {
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

#' @title Convert an anndata object to a seurat object
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param adata An AnnData object.
#' Can be a Python AnnData object (from `scanpy`/`reticulate``),
#' an `AnnDataR6` object from the `anndata` package,
#' or an `InMemoryAnnData` object from the `anndataR` package.
#'
#' @export
#'
#' @seealso [srt_to_adata]
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(
#'   adata = adata,
#'   group.by = "SubCellType",
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
  verbose = TRUE
) {
  if (!isTRUE(getOption("scop_skip_python_prepare", FALSE))) {
    old_log_verbose <- getOption("log_message.verbose", TRUE)
    if (!isTRUE(verbose)) {
      options(log_message.verbose = FALSE)
      on.exit(
        options(log_message.verbose = old_log_verbose),
        add = TRUE
      )
    }
    PrepareEnv(modules = "scanpy")
  }
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
  skipped_layers <- list()
  if (length(keys) > 0) {
    for (k in keys) {
      k_clean <- py_to_r2(k)
      err <- tryCatch(
        {
          raw <- get_adata_element(adata$layers, k, missing = "null")
          if (is.null(raw)) {
            stop("cannot access layer", call. = FALSE)
          }
          layer <- py_to_r2(raw)
          if (!inherits(layer, c("Matrix", "matrix"))) {
            stop(
              "not a matrix: ",
              paste(class(layer), collapse = ", "),
              call. = FALSE
            )
          }
          layer <- Matrix::t(layer)
          if (!inherits(layer, "dgCMatrix")) {
            layer <- SeuratObject::as.sparse(layer)
          }
          rownames(layer) <- get_adata_names(adata, "var")
          colnames(layer) <- get_adata_names(adata, "obs")
          srt[[k_clean]] <- Seurat::CreateAssayObject(counts = layer)
          NULL
        },
        error = function(e) e
      )
      if (inherits(err, "error")) {
        skipped_layers[[k_clean]] <- conditionMessage(err)
      }
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
  if (length(skipped_layers) > 0) {
    log_message(
      "Some layers were skipped: {.val {names(skipped_layers)}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  return(srt)
}

#' @title Read an `.h5ad` file and convert to a `Seurat`
#'
#' @md
#' @inheritParams adata_to_srt
#' @param path Path to an `.h5ad` file (passed to `anndata.read_h5ad()`).
#' @param prepare_for_reticulate If `TRUE` (default), coerces `X` and each layer
#'   matrix to CSR `float64` in Python (avoids invalid `dgRMatrix` conversion
#'   via reticulate). Layers that still fail in [adata_to_srt()] are skipped and
#'   reported. Set to `FALSE` for a plain `read_h5ad` then convert.
#'
#' @return A `Seurat` object.
#'
#' @export
#'
#' @seealso [adata_to_srt], [srt_to_adata]
#'
#' @examples
#' \dontrun{
#' srt <- h5ad_to_srt("path/to/data.h5ad")
#' srt
#' }
h5ad_to_srt <- function(
  path,
  verbose = TRUE,
  prepare_for_reticulate = TRUE
) {
  old_log_verbose <- getOption("log_message.verbose", TRUE)
  if (!isTRUE(verbose)) {
    options(log_message.verbose = FALSE)
    on.exit(
      options(log_message.verbose = old_log_verbose),
      add = TRUE
    )
  }
  PrepareEnv(modules = "scanpy")
  if (isTRUE(prepare_for_reticulate)) {
    check_python(
      c("anndata", "numpy", "scipy"),
      verbose = FALSE
    )
  } else {
    check_python(
      c("anndata", "numpy"),
      verbose = FALSE
    )
  }

  path <- normalizePath(path.expand(path), mustWork = TRUE, winslash = "/")
  if (!length(path) || !nzchar(path[1L])) {
    log_message(
      "{.arg path} must be a non-empty path",
      message_type = "error"
    )
  }

  if (isTRUE(prepare_for_reticulate)) {
    path_py <- gsub("\"", "\\\\\"", path, fixed = TRUE)
    reticulate::py_run_string(paste0(
      "
import numpy as np
import anndata as ad
import __main__
adata = ad.read_h5ad(r\"",
      path_py,
      "\")

def as_csr_f64(x):
    import scipy.sparse as sp
    if sp.issparse(x):
        return x.tocsr().astype(np.float64)
    return np.asarray(x, dtype=np.float64)

adata.X = as_csr_f64(adata.X)
for k in list(adata.layers.keys()):
    adata.layers[k] = as_csr_f64(adata.layers[k])
__main__.adata = adata
"
    ))
    main <- reticulate::import("__main__", convert = FALSE)
    adata <- main$adata
  } else {
    ad <- reticulate::import("anndata", convert = FALSE)
    adata <- ad$read_h5ad(path)
  }

  old_skip_prepare <- getOption("scop_skip_python_prepare", FALSE)
  options(scop_skip_python_prepare = TRUE)
  on.exit(
    options(scop_skip_python_prepare = old_skip_prepare),
    add = TRUE
  )
  adata_to_srt(adata, verbose = verbose)
}

#' @title Read a `.loom` file as an AnnData object
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param path Path to a `.loom` file (passed to `scanpy.read_loom()`).
#' @param ... Additional arguments passed to `scanpy.read_loom()`.
#'
#' @return A Python `anndata.AnnData` object.
#'
#' @details
#' This is a Python-backed wrapper and requires `reticulate` plus a Python
#' environment with `scanpy` and `loompy` available. It is independent from
#' [loom_to_srt()], which reads loom files directly in R without initializing
#' Python.
#'
#' @export
#'
#' @seealso [loom_to_srt], [adata_to_srt], [srt_to_adata]
#'
#' @examples
#' \dontrun{
#' adata <- loom_to_adata("path/to/data.loom")
#' adata
#' }
loom_to_adata <- function(
  path,
  verbose = TRUE,
  ...
) {
  old_log_verbose <- getOption("log_message.verbose", TRUE)
  if (!isTRUE(verbose)) {
    options(log_message.verbose = FALSE)
    on.exit(
      options(log_message.verbose = old_log_verbose),
      add = TRUE
    )
  }
  PrepareEnv(modules = "scanpy")
  check_python(
    c("scanpy", "loompy"),
    verbose = FALSE
  )

  path <- normalizePath(path.expand(path), mustWork = TRUE, winslash = "/")
  if (!length(path) || !nzchar(path[1L])) {
    log_message(
      "{.arg path} must be a non-empty path",
      message_type = "error"
    )
  }

  log_message(
    "Reading {.file {path}} as {.cls AnnData} via {.pkg scanpy}",
    verbose = verbose
  )
  sc <- reticulate::import("scanpy", convert = FALSE)
  adata <- sc$read_loom(path, ...)
  log_message(
    "Read {.file {path}} as {.cls AnnData} completed",
    message_type = "success",
    verbose = verbose
  )
  adata
}

#' @title Read a `.loom` file and convert to a `Seurat`
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param path Path to a `.loom` file.
#' @param layers Character vector of loom layers to import as additional Seurat
#' assays. Missing layers are skipped with a warning. Default is
#' `c("spliced", "unspliced")`.
#' @param chunk_rows Number of feature rows to read from each matrix dataset per
#' chunk. Larger values can be faster but use more memory.
#'
#' @return A `Seurat` object.
#'
#' @details
#' This function reads loom/HDF5 files directly in R using `rhdf5`; it does not
#' call `reticulate`, `scanpy`, or [loom_to_adata()]. The loom `/matrix` dataset
#' is imported as the `RNA` assay. Requested `/layers/*` datasets are imported
#' as additional assays, which is useful for velocity-style loom files with
#' `spliced` and `unspliced` layers.
#'
#' @export
#'
#' @seealso [loom_to_adata], [adata_to_srt], [srt_to_adata]
#'
#' @examples
#' \dontrun{
#' srt <- loom_to_srt("path/to/data.loom")
#' srt
#' }
loom_to_srt <- function(
  path,
  layers = c("spliced", "unspliced"),
  verbose = TRUE,
  chunk_rows = 1000
) {
  old_log_verbose <- getOption("log_message.verbose", TRUE)
  if (!isTRUE(verbose)) {
    options(log_message.verbose = FALSE)
    on.exit(
      options(log_message.verbose = old_log_verbose),
      add = TRUE
    )
  }
  check_r("rhdf5", verbose = FALSE)

  path <- normalizePath(path.expand(path), mustWork = TRUE, winslash = "/")
  if (!length(path) || !nzchar(path[1L])) {
    log_message(
      "{.arg path} must be a non-empty path",
      message_type = "error"
    )
  }
  if (!is.numeric(chunk_rows) || length(chunk_rows) != 1 || is.na(chunk_rows) || chunk_rows < 1) {
    log_message(
      "{.arg chunk_rows} must be a positive number",
      message_type = "error"
    )
  }
  chunk_rows <- as.integer(chunk_rows)
  layers <- unique(as.character(layers %||% character(0)))
  layers <- layers[nzchar(layers)]

  log_message(
    "Reading {.file {path}} as {.cls Seurat} directly from loom/HDF5",
    verbose = verbose
  )

  paths <- loom_h5_paths(path)
  if (!loom_h5_exists(paths, "/matrix")) {
    log_message(
      "{.file {path}} does not contain the required {.val /matrix} dataset",
      message_type = "error"
    )
  }

  matrix_dim <- loom_h5_dim(paths, "/matrix")
  if (length(matrix_dim) != 2) {
    log_message(
      "{.val /matrix} must be a two-dimensional dataset",
      message_type = "error"
    )
  }
  n_features <- loom_h5_length(
    paths,
    c("/row_attrs/Gene", "/row_attrs/gene", "/row_attrs/Accession")
  )
  n_cells <- loom_h5_length(
    paths,
    c("/col_attrs/CellID", "/col_attrs/cell_id", "/col_attrs/CellID_1")
  )
  if (is.na(n_features) || is.na(n_cells)) {
    n_features <- matrix_dim[[1]]
    n_cells <- matrix_dim[[2]]
  }
  matrix_orientation <- loom_matrix_orientation(
    matrix_dim = matrix_dim,
    n_features = n_features,
    n_cells = n_cells
  )

  features <- loom_read_names(
    path = path,
    paths = paths,
    candidates = c("/row_attrs/Gene", "/row_attrs/gene", "/row_attrs/Accession"),
    n = n_features,
    fallback_prefix = "feature"
  )
  cells <- loom_read_names(
    path = path,
    paths = paths,
    candidates = c("/col_attrs/CellID", "/col_attrs/cell_id", "/col_attrs/CellID_1"),
    n = n_cells,
    fallback_prefix = "cell"
  )

  counts <- loom_read_matrix_sparse(
    path = path,
    paths = paths,
    dataset = "/matrix",
    orientation = matrix_orientation,
    features = features,
    cells = cells,
    chunk_rows = chunk_rows
  )
  col_metadata <- loom_read_attrs_dataframe(
    path = path,
    paths = paths,
    group = "/col_attrs",
    n = n_cells,
    rownames = cells
  )
  feature_metadata <- loom_read_attrs_dataframe(
    path = path,
    paths = paths,
    group = "/row_attrs",
    n = n_features,
    rownames = features
  )

  srt <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = col_metadata
  )
  if (ncol(feature_metadata) > 0) {
    srt[["RNA"]] <- Seurat::AddMetaData(
      object = srt[["RNA"]],
      metadata = feature_metadata
    )
  }

  for (layer in layers) {
    layer_path <- paste0("/layers/", layer)
    if (!loom_h5_exists(paths, layer_path)) {
      log_message(
        "Loom layer {.val {layer}} was not found and will be skipped",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    layer_dim <- loom_h5_dim(paths, layer_path)
    if (!identical(as.integer(layer_dim), as.integer(matrix_dim))) {
      log_message(
        "Loom layer {.val {layer}} has dimensions {.val {paste(layer_dim, collapse = ' x ')}} and will be skipped",
        message_type = "warning",
        verbose = verbose
      )
      next
    }
    layer_counts <- loom_read_matrix_sparse(
      path = path,
      paths = paths,
      dataset = layer_path,
      orientation = matrix_orientation,
      features = features,
      cells = cells,
      chunk_rows = chunk_rows
    )
    srt[[layer]] <- Seurat::CreateAssayObject(counts = layer_counts)
    if (ncol(feature_metadata) > 0) {
      srt[[layer]] <- Seurat::AddMetaData(
        object = srt[[layer]],
        metadata = feature_metadata
      )
    }
  }

  log_message(
    "Read {.file {path}} as {.cls Seurat} completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

py_to_r2 <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    reticulate::py_to_r(x)
  } else {
    x
  }
}

get_adata_element <- function(
  container,
  key,
  missing = c("error", "null")
) {
  missing <- match.arg(missing)
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
      if (missing == "null") {
        return(NULL)
      }
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
  depth = max_depth(x)
) {
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

loom_h5_paths <- function(path) {
  h5 <- rhdf5::h5ls(path, recursive = TRUE)
  h5[["full_path"]] <- ifelse(
    h5[["group"]] == "/",
    paste0("/", h5[["name"]]),
    paste0(h5[["group"]], "/", h5[["name"]])
  )
  h5[["full_path"]] <- gsub("//+", "/", h5[["full_path"]])
  h5
}

loom_h5_exists <- function(paths, name) {
  name <- gsub("//+", "/", name)
  name %in% paths[["full_path"]]
}

loom_h5_dim <- function(paths, name) {
  name <- gsub("//+", "/", name)
  i <- match(name, paths[["full_path"]])
  if (is.na(i)) {
    return(integer(0))
  }
  dim_raw <- paths[["dim"]][[i]]
  if (is.null(dim_raw) || length(dim_raw) == 0 || anyNA(dim_raw)) {
    return(integer(0))
  }
  if (is.numeric(dim_raw)) {
    return(as.integer(dim_raw))
  }
  dim_raw <- gsub(",", "", as.character(dim_raw))
  dim_raw <- trimws(unlist(strsplit(dim_raw, " x ", fixed = TRUE)))
  as.integer(dim_raw[nzchar(dim_raw)])
}

loom_h5_length <- function(paths, candidates) {
  for (candidate in candidates) {
    dims <- loom_h5_dim(paths, candidate)
    if (length(dims) > 0 && all(!is.na(dims))) {
      return(prod(dims))
    }
  }
  NA_integer_
}

loom_matrix_orientation <- function(
  matrix_dim,
  n_features,
  n_cells
) {
  if (identical(as.integer(matrix_dim), as.integer(c(n_features, n_cells)))) {
    return("features_by_cells")
  }
  if (identical(as.integer(matrix_dim), as.integer(c(n_cells, n_features)))) {
    return("cells_by_features")
  }
  log_message(
    "{.val /matrix} dimensions do not match loom row/column attributes",
    message_type = "error"
  )
}

loom_read_names <- function(
  path,
  paths,
  candidates,
  n,
  fallback_prefix
) {
  for (candidate in candidates) {
    if (!loom_h5_exists(paths, candidate)) {
      next
    }
    values <- tryCatch(
      as.vector(rhdf5::h5read(path, name = candidate)),
      error = function(e) NULL
    )
    if (!is.null(values) && length(values) == n) {
      values <- as.character(values)
      values[is.na(values) | !nzchar(values)] <- paste0(
        fallback_prefix,
        "_",
        which(is.na(values) | !nzchar(values))
      )
      return(make.unique(values))
    }
  }
  make.unique(paste0(fallback_prefix, "_", seq_len(n)))
}

loom_read_attrs_dataframe <- function(
  path,
  paths,
  group,
  n,
  rownames
) {
  group <- sub("/+$", "", group)
  rows <- paths[paths[["group"]] == group & paths[["otype"]] == "H5I_DATASET", , drop = FALSE]
  if (nrow(rows) == 0) {
    out <- data.frame(row.names = rownames)
    return(out)
  }

  attrs <- list()
  for (i in seq_len(nrow(rows))) {
    attr_path <- rows[["full_path"]][[i]]
    attr_name <- make.names(rows[["name"]][[i]])
    values <- tryCatch(
      as.vector(rhdf5::h5read(path, name = attr_path)),
      error = function(e) NULL
    )
    if (is.null(values) || length(values) != n) {
      next
    }
    if (is.raw(values) || is.list(values)) {
      next
    }
    attrs[[attr_name]] <- values
  }

  if (length(attrs) == 0) {
    out <- data.frame(row.names = rownames)
    return(out)
  }
  out <- as.data.frame(attrs, stringsAsFactors = FALSE, check.names = TRUE)
  rownames(out) <- rownames
  out
}

loom_read_matrix_sparse <- function(
  path,
  paths,
  dataset,
  orientation,
  features,
  cells,
  chunk_rows
) {
  dims <- loom_h5_dim(paths, dataset)
  if (length(dims) != 2) {
    log_message(
      "{.val {dataset}} must be a two-dimensional dataset",
      message_type = "error"
    )
  }
  n_features <- dims[[1]]
  n_cells <- dims[[2]]
  if (identical(orientation, "cells_by_features")) {
    n_features <- dims[[2]]
    n_cells <- dims[[1]]
  }
  if (length(features) != n_features || length(cells) != n_cells) {
    log_message(
      "The supplied dimnames do not match {.val {dataset}}",
      message_type = "error"
    )
  }

  triplet_i <- integer(0)
  triplet_j <- integer(0)
  triplet_x <- numeric(0)
  feature_starts <- seq.int(1L, n_features, by = chunk_rows)
  for (start in feature_starts) {
    end <- min(start + chunk_rows - 1L, n_features)
    if (identical(orientation, "features_by_cells")) {
      block <- rhdf5::h5read(
        file = path,
        name = dataset,
        index = list(start:end, seq_len(n_cells))
      )
      if (!is.matrix(block)) {
        block <- matrix(block, nrow = end - start + 1L, ncol = n_cells)
      }
      sparse_block <- Matrix::Matrix(block, sparse = TRUE)
    } else if (identical(orientation, "cells_by_features")) {
      block <- rhdf5::h5read(
        file = path,
        name = dataset,
        index = list(seq_len(n_cells), start:end)
      )
      if (!is.matrix(block)) {
        block <- matrix(block, nrow = n_cells, ncol = end - start + 1L)
      }
      sparse_block <- Matrix::t(Matrix::Matrix(block, sparse = TRUE))
    } else {
      log_message(
        "{.arg orientation} is invalid",
        message_type = "error"
      )
    }
    sparse_block <- Matrix::drop0(sparse_block)
    if (length(sparse_block@x) == 0) {
      next
    }
    summary_block <- Matrix::summary(sparse_block)
    triplet_i <- c(triplet_i, as.integer(summary_block[["i"]]) + start - 1L)
    triplet_j <- c(triplet_j, as.integer(summary_block[["j"]]))
    triplet_x <- c(triplet_x, as.numeric(summary_block[["x"]]))
  }

  mat <- Matrix::sparseMatrix(
    i = triplet_i,
    j = triplet_j,
    x = triplet_x,
    dims = c(n_features, n_cells),
    dimnames = list(features, cells)
  )
  methods::as(mat, "dgCMatrix")
}
