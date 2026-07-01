#' @title Run Giotto nearest-network clustering
#'
#' @description
#' Run Giotto as a temporary backend for nearest-network clustering and return
#' the complete Giotto object together with extracted cluster results. The input
#' `Seurat` object is not modified.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used as the expression matrix.
#' @param features Features used for PCA and clustering. If `NULL`, current
#' variable features are used, falling back to all assay features.
#' @param method Giotto clustering method.
#' @param dims Dimensions used to build the Giotto nearest-neighbor network.
#' @param k Number of nearest neighbors used by Giotto.
#' @param resolution Resolution passed to Giotto clustering.
#' @param cluster_colname Result column name recorded in returned parameters.
#' This function does not write to `srt@meta.data`.
#' @param tool_name Result name recorded in returned parameters. This function
#' does not write to `srt@tools`.
#' @param store_giotto Deprecated compatibility argument. The complete Giotto
#' object is always returned in the `giotto` element.
#' @param conversion_params Additional parameters passed to
#' `Giotto::createGiottoObject()`.
#' @param preprocess_params Additional parameters passed to `Giotto::runPCA()`.
#' @param network_params Additional parameters passed to
#' `Giotto::createNearestNetwork()`.
#' @param cluster_params Additional parameters passed to Giotto clustering.
#'
#' @return A `giotto2_result` list containing the full Giotto object,
#' cluster assignments, Giotto metadata, parameters, features, and cells.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' giotto_clusters <- list(
#'   clusters = data.frame(
#'     cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
#'     row.names = colnames(spatial)
#'   ),
#'   parameters = list(
#'     cluster_colname = "Giotto_cluster",
#'     coord.cols = c("x", "y"),
#'     k = 8,
#'     resolution = 0.4
#'   )
#' )
#' class(giotto_clusters) <- c("giotto2_cluster", "giotto2_result", "list")
#' head(giotto_clusters$clusters)
#' GiottoPlot(
#'   giotto_clusters,
#'   srt = spatial,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   requireNamespace("Giotto", quietly = TRUE) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' spatial <- Seurat::FindVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 300,
#'   verbose = FALSE
#' )
#' giotto_clusters <- RunGiottoCluster(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "data",
#'   dims = 1:10,
#'   k = 8,
#'   resolution = 0.4,
#'   coord.cols = c("x", "y")
#' )
#'
#' head(giotto_clusters$clusters)
#' }
#'
#' @export
RunGiottoCluster <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = NULL,
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  cluster_colname = "Giotto_cluster",
  tool_name = "GiottoCluster",
  store_giotto = TRUE,
  conversion_params = list(),
  preprocess_params = list(),
  network_params = list(),
  cluster_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  giotto_validate_scalar_string(cluster_colname, "cluster_colname")
  giotto_validate_scalar_string(tool_name, "tool_name")
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_validate_named_list(preprocess_params, "preprocess_params")
  giotto_validate_named_list(network_params, "network_params")
  giotto_validate_named_list(cluster_params, "cluster_params")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    log_message(
      "{.arg k} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.numeric(resolution) || length(resolution) != 1L || is.na(resolution) || resolution <= 0) {
    log_message(
      "{.arg resolution} must be a positive number",
      message_type = "error"
    )
  }
  dims <- unique(as.integer(dims))
  dims <- dims[is.finite(dims) & dims > 0L]
  if (length(dims) == 0L) {
    log_message(
      "{.arg dims} must contain at least one positive dimension",
      message_type = "error"
    )
  }

  giotto_require(verbose = verbose)
  giotto_old_options <- options(
    giotto.has_conda = FALSE,
    giotto.use_conda = FALSE,
    giotto.update_param = FALSE,
    giotto.no_python_warn = TRUE
  )
  on.exit(options(giotto_old_options), add = TRUE)

  log_message(
    "Create temporary {.pkg Giotto} object from {.cls Seurat}",
    verbose = verbose
  )
  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols
  )
  max_dims <- min(nrow(input$expr), ncol(input$expr)) - 1L
  dims <- dims[dims <= max_dims]
  if (length(dims) == 0L) {
    log_message(
      "No valid {.arg dims} remain for Giotto PCA with the current expression matrix",
      message_type = "error"
    )
  }
  k <- min(as.integer(k), ncol(input$expr) - 1L)
  gobject <- giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )

  log_message(
    "Run {.pkg Giotto} PCA and nearest-network clustering",
    verbose = verbose
  )
  expression_values <- giotto_expression_values(layer)
  if (identical(expression_values, "raw")) {
    normalizeGiotto <- giotto_get_fun("normalizeGiotto")
    gobject <- giotto_call(
      normalizeGiotto,
      giotto_merge_args(
        list(gobject = gobject, verbose = FALSE),
        preprocess_params[["normalize_params"]] %||% list()
      )
    )
    expression_values <- "normalized"
  }

  runPCA <- giotto_get_fun("runPCA")
  pca_args <- giotto_merge_args(
    list(
      gobject = gobject,
      expression_values = expression_values,
      feats_to_use = input$features,
      name = "pca",
      ncp = max(dims),
      scale_unit = TRUE,
      center = TRUE,
      set_seed = TRUE,
      seed_number = seed
    ),
    preprocess_params[setdiff(names(preprocess_params), "normalize_params")],
    arg_name = "preprocess_params"
  )
  gobject <- giotto_call(runPCA, pca_args)
  pca_name <- pca_args[["name"]] %||% "pca"
  giotto_validate_scalar_string(pca_name, "preprocess_params$name")

  network_name <- network_params[["name"]] %||% "scop_NN"
  createNearestNetwork <- giotto_get_fun("createNearestNetwork")
  network_args <- giotto_merge_args(
    list(
      gobject = gobject,
      dim_reduction_to_use = "pca",
      dim_reduction_name = pca_name,
      dimensions_to_use = dims,
      k = as.integer(k),
      name = network_name
    ),
    network_params,
    arg_name = "network_params",
    reserved = c("gobject", "return_gobject")
  )
  network_name <- network_args[["name"]] %||% network_name
  network_type <- network_args[["type"]] %||% "sNN"
  gobject <- giotto_call(createNearestNetwork, network_args)

  giotto_cluster_name <- cluster_params[["name"]] %||% paste0(method, "_clus")
  cluster_fun <- giotto_get_fun(
    if (identical(method, "leiden")) "doLeidenCluster" else "doLouvainCluster"
  )
  cluster_args <- giotto_merge_args(
    list(
      gobject = gobject,
      network_name = network_name,
      nn_network_to_use = network_type,
      resolution = resolution,
      name = giotto_cluster_name,
      return_gobject = TRUE,
      set_seed = TRUE,
      seed_number = seed
    ),
    cluster_params,
    arg_name = "cluster_params",
    reserved = c("gobject", "return_gobject")
  )
  gobject <- giotto_call(cluster_fun, cluster_args)

  giotto_metadata <- giotto_get_cell_metadata(gobject)
  clusters <- giotto_extract_clusters(
    giotto_metadata = giotto_metadata,
    cells = input$cells,
    cluster_name = giotto_cluster_name,
    method = method
  )

  result <- giotto_result(
    result_type = "cluster",
    giotto = gobject,
    clusters = data.frame(
      cluster = as.character(clusters),
      row.names = names(clusters),
      stringsAsFactors = FALSE
    ),
    cluster_vector = clusters,
    giotto_metadata = giotto_metadata,
    parameters = list(
      assay = assay,
      layer = layer,
      method = method,
      image = image,
      coord.cols = coord.cols,
      dims = dims,
      k = as.integer(k),
      resolution = resolution,
      cluster_colname = cluster_colname,
      giotto_cluster_name = giotto_cluster_name,
      conversion_params = conversion_params,
      preprocess_params = preprocess_params,
      network_params = network_params,
      cluster_params = cluster_params,
      tool_name = tool_name,
      store_giotto = store_giotto,
      seed = seed
    ),
    features = input$features,
    cells = input$cells
  )

  log_message(
    "{.pkg Giotto} clustering completed; returning standalone result without modifying {.arg srt}",
    message_type = "success",
    verbose = verbose
  )
  result
}

giotto_require <- function(verbose = TRUE) {
  if (!giotto_namespace_available()) {
    status <- check_r("giotto-suite/Giotto", verbose = FALSE)
    if (!isTRUE(unname(unlist(status))[1]) || !giotto_namespace_available()) {
      log_message(
        "Please install {.pkg Giotto} before running {.fn RunGiottoCluster}: {.code check_r('giotto-suite/Giotto')}",
        message_type = "error",
        verbose = verbose
      )
    }
  }
  invisible(TRUE)
}

giotto_prepare_input <- function(
  srt,
  assay,
  layer,
  features = NULL,
  image = NULL,
  coord.cols = NULL
) {
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  coords <- giotto_spatial_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols
  )
  cells <- intersect(colnames(srt), colnames(expr))
  cells <- intersect(cells, rownames(coords))
  coords <- coords[cells, c("x", "y"), drop = FALSE]
  finite_coords <- is.finite(coords$x) & is.finite(coords$y)
  cells <- cells[finite_coords]
  if (length(cells) < 3L) {
    log_message(
      "At least three Seurat cells/spots with finite spatial coordinates in assay layer {.val {layer}} are required",
      message_type = "error"
    )
  }
  coords <- coords[cells, , drop = FALSE]
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
    if (length(features) == 0L) {
      features <- rownames(expr)
    }
  }
  features <- intersect(unique(features), rownames(expr))
  if (length(features) == 0L) {
    log_message(
      "No requested {.arg features} are present in assay {.val {assay}}",
      message_type = "error"
    )
  }
  expr <- expr[features, cells, drop = FALSE]
  keep_features <- Matrix::rowSums(expr != 0) > 0
  if (!any(keep_features)) {
    log_message(
      "No non-zero features remain for {.fn RunGiottoCluster}",
      message_type = "error"
    )
  }
  expr <- expr[keep_features, , drop = FALSE]
  features <- rownames(expr)

  counts <- tryCatch(
    GetAssayData5(srt, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  raw_expr <- if (
    !is.null(counts) &&
      all(features %in% rownames(counts)) &&
      all(cells %in% colnames(counts))
  ) {
    counts[features, cells, drop = FALSE]
  } else {
    expr
  }
  spatial_locs <- data.frame(
    cell_ID = cells,
    sdimx = coords$x,
    sdimy = coords$y,
    stringsAsFactors = FALSE
  )
  metadata <- srt@meta.data[cells, , drop = FALSE]
  metadata <- data.frame(
    cell_ID = cells,
    metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    expr = expr,
    raw_expr = raw_expr,
    spatial_locs = spatial_locs,
    metadata = metadata,
    features = features,
    cells = cells
  )
}

giotto_spatial_coords <- function(srt, image = NULL, coord.cols = NULL) {
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  if (length(images) > 0L) {
    image <- image %||% images[1L]
    if (!image %in% images) {
      log_message(
        "{.arg image} {.val {image}} is not present in {.cls Seurat}",
        message_type = "error"
      )
    }
    coords <- as.data.frame(SeuratObject::GetTissueCoordinates(srt[[image]]))
    cell_col <- if ("cell" %in% colnames(coords)) "cell" else NULL
    cells <- if (is.null(cell_col)) rownames(coords) else coords[[cell_col]]
    rownames(coords) <- cells
    x_col <- spatial_dim_pick_col(
      coords,
      c("x", "pxl_col_in_fullres", "imagecol")
    )
    y_col <- spatial_dim_pick_col(
      coords,
      c("y", "pxl_row_in_fullres", "imagerow")
    )
    return(data.frame(
      x = coords[[x_col]],
      y = coords[[y_col]],
      row.names = cells,
      stringsAsFactors = FALSE
    ))
  }

  coord.cols <- scop_spatial_resolve_coord_cols(srt, coord.cols = coord.cols)
  data.frame(
    x = srt@meta.data[[coord.cols[1L]]],
    y = srt@meta.data[[coord.cols[2L]]],
    row.names = rownames(srt@meta.data),
    stringsAsFactors = FALSE
  )
}

giotto_create_object <- function(
  input,
  layer,
  conversion_params = list(),
  verbose = TRUE
) {
  createGiottoObject <- giotto_get_fun("createGiottoObject")
  args <- list(
    expression = input$raw_expr,
    spatial_locs = input$spatial_locs,
    cell_metadata = input$metadata
  )
  args <- giotto_merge_args(args, conversion_params, arg_name = "conversion_params")
  gobject <- giotto_suppress_known_warnings(
    giotto_call(createGiottoObject, args)
  )
  if (is.null(gobject)) {
    log_message(
      "{.pkg Giotto} did not return a valid object",
      message_type = "error",
      verbose = verbose
    )
  }
  if (!identical(layer, "counts")) {
    createExprObj <- giotto_get_fun("createExprObj")
    expr_obj <- giotto_suppress_known_warnings(
      giotto_call(
        createExprObj,
        list(
          expression_data = input$expr,
          name = giotto_expression_values(layer),
          spat_unit = "cell",
          feat_type = "rna",
          provenance = "cell"
        )
      )
    )
    setExpression <- giotto_get_fun("setExpression")
    gobject <- giotto_call(
      setExpression,
      list(
        gobject = gobject,
        x = expr_obj,
        verbose = FALSE
      )
    )
  }
  gobject
}

giotto_get_cell_metadata <- function(gobject) {
  pDataDT <- giotto_get_fun("pDataDT")
  meta <- giotto_call(pDataDT, list(gobject = gobject))
  as.data.frame(meta, stringsAsFactors = FALSE)
}

giotto_extract_clusters <- function(
  giotto_metadata,
  cells,
  cluster_name,
  method
) {
  candidates <- unique(c(cluster_name, paste0(method, "_clus")))
  cluster_col <- candidates[candidates %in% colnames(giotto_metadata)][1]
  if (is.na(cluster_col)) {
    log_message(
      "{.pkg Giotto} did not return a cluster column among {.val {candidates}}",
      message_type = "error"
    )
  }
  ids <- giotto_metadata_cell_ids(giotto_metadata, cells = cells)
  clusters <- giotto_metadata[[cluster_col]]
  names(clusters) <- ids
  missing_cells <- setdiff(cells, names(clusters))
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg Giotto} metadata could not be matched to Seurat cells/spots: {.val {missing_cells}}",
      message_type = "error"
    )
  }
  clusters[cells]
}

giotto_metadata_cell_ids <- function(giotto_metadata, cells) {
  id_candidates <- c("cell_ID", "cell", "spat_ID", "spot", "barcode")
  id_col <- id_candidates[id_candidates %in% colnames(giotto_metadata)][1]
  if (!is.na(id_col)) {
    return(as.character(giotto_metadata[[id_col]]))
  }
  rn <- rownames(giotto_metadata)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(giotto_metadata))))) {
    return(as.character(rn))
  }
  if (nrow(giotto_metadata) == length(cells)) {
    return(cells)
  }
  log_message(
    "Unable to identify Giotto cell/spot IDs in returned metadata",
    message_type = "error"
  )
}

giotto_expression_slot <- function(layer) {
  if (identical(layer, "scale.data")) {
    return("norm_scaled_expr")
  }
  if (identical(layer, "counts")) {
    return("raw_exprs")
  }
  "norm_expr"
}

giotto_expression_values <- function(layer) {
  if (identical(layer, "scale.data")) {
    return("scaled")
  }
  if (identical(layer, "counts")) {
    return("raw")
  }
  "normalized"
}

giotto_get_fun <- function(name) {
  pkg <- "Giotto"
  if (!giotto_namespace_available()) {
    log_message(
      "Please install {.pkg Giotto} before running {.fn RunGiottoCluster}",
      message_type = "error"
    )
  }
  exported <- tryCatch(getExportedValue(pkg, name), error = function(e) NULL)
  if (is.function(exported)) {
    return(exported)
  }
  ns <- asNamespace(pkg)
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  pkg_env <- tryCatch(as.environment(paste0("package:", pkg)), error = function(e) NULL)
  if (!is.null(pkg_env) && exists(name, envir = pkg_env, inherits = FALSE)) {
    return(get(name, envir = pkg_env, inherits = FALSE))
  }
  log_message(
    "Function {.val {name}} not found in {.pkg Giotto} namespace",
    message_type = "error"
  )
}

giotto_call <- function(fun, args) {
  fmls <- giotto_formal_names(fun)
  if (!is.null(fmls) && !"..." %in% fmls) {
    dropped <- setdiff(names(args), fmls)
    if (length(dropped) > 0L) {
      log_message(
        "Unsupported Giotto argument(s): {.val {dropped}}",
        message_type = "error"
      )
    }
  }
  call <- as.call(c(list(fun), args))
  eval(call, envir = parent.frame())
}

giotto_suppress_known_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("\\[createExprObj\\].*deprecated", msg)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

giotto_formal_names <- function(fun) {
  tryCatch(names(formals(fun)), error = function(e) NULL)
}

giotto_merge_args <- function(defaults, extra, arg_name = "params", reserved = character()) {
  if (length(extra) == 0L) {
    return(defaults)
  }
  blocked <- intersect(names(extra), reserved)
  if (length(blocked) > 0L) {
    log_message(
      "{.arg {arg_name}} cannot override required Giotto argument(s): {.val {blocked}}",
      message_type = "error"
    )
  }
  utils::modifyList(defaults, extra)
}

giotto_validate_scalar_string <- function(x, arg_name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg_name}} must be a non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

giotto_validate_named_list <- function(x, arg_name) {
  if (!is.list(x) || (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x)))))) {
    log_message(
      "{.arg {arg_name}} must be a named list",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

giotto_namespace_available <- function() {
  do.call("requireNamespace", list(package = "Giotto", quietly = TRUE))
}

giotto_result <- function(result_type, giotto, ...) {
  structure(
    list(
      result_type = result_type,
      giotto = giotto,
      ...
    ),
    class = c(paste0("giotto2_", result_type), "giotto2_result", "list")
  )
}
