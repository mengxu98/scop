#' @title Run Giotto nearest-network clustering
#'
#' @description
#' Run Giotto as a temporary backend for nearest-network clustering and write
#' the resulting spot or cell clusters back to a `Seurat` object.
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
#' @param cluster_colname Metadata column written back to `srt`.
#' @param tool_name Name of the `srt@tools` entry used to store results.
#' @param store_giotto Whether to store the full Giotto object in `srt@tools`.
#' @param conversion_params Additional parameters passed to
#' `Giotto::createGiottoObject()`.
#' @param preprocess_params Additional parameters passed to `Giotto::runPCA()`.
#' @param network_params Additional parameters passed to
#' `Giotto::createNearestNetwork()`.
#' @param cluster_params Additional parameters passed to Giotto clustering.
#'
#' @return A `Seurat` object with Giotto clusters in metadata and raw results in
#' `srt@tools[[tool_name]]`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(visium_human_pancreas_sub)
#' spatial <- Seurat::NormalizeData(
#'   visium_human_pancreas_sub,
#'   assay = "Spatial",
#'   verbose = FALSE
#' )
#' spatial <- RunGiottoCluster(
#'   spatial,
#'   assay = "Spatial",
#'   method = "leiden",
#'   k = 15,
#'   resolution = 0.4
#' )
#' SpatialSpotPlot(spatial, group.by = "Giotto_cluster")
#' }
RunGiottoCluster <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  cluster_colname = "Giotto_cluster",
  tool_name = "GiottoCluster",
  store_giotto = FALSE,
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
      center = TRUE
    ),
    preprocess_params[setdiff(names(preprocess_params), "normalize_params")]
  )
  gobject <- giotto_call(runPCA, pca_args)

  network_name <- network_params[["name"]] %||% "scop_NN"
  createNearestNetwork <- giotto_get_fun("createNearestNetwork")
  network_args <- giotto_merge_args(
    list(
      gobject = gobject,
      dim_reduction_to_use = "pca",
      dim_reduction_name = "pca",
      dimensions_to_use = dims,
      k = as.integer(k),
      name = network_name
    ),
    network_params
  )
  gobject <- giotto_call(createNearestNetwork, network_args)

  giotto_cluster_name <- cluster_params[["name"]] %||% paste0(method, "_clus")
  cluster_fun <- giotto_get_fun(
    if (identical(method, "leiden")) "doLeidenCluster" else "doLouvainCluster"
  )
  cluster_args <- giotto_merge_args(
    list(
      gobject = gobject,
      network_name = network_name,
      resolution = resolution,
      name = giotto_cluster_name,
      return_gobject = TRUE,
      set_seed = TRUE,
      seed_number = seed
    ),
    cluster_params
  )
  gobject <- giotto_call(cluster_fun, cluster_args)

  giotto_metadata <- giotto_get_cell_metadata(gobject)
  clusters <- giotto_extract_clusters(
    giotto_metadata = giotto_metadata,
    cells = input$cells,
    cluster_name = giotto_cluster_name,
    method = method
  )
  srt <- giotto_add_cluster_metadata(
    srt = srt,
    clusters = clusters,
    cluster_colname = cluster_colname
  )

  tool <- list(
    clusters = data.frame(
      cluster = as.character(clusters),
      row.names = names(clusters),
      stringsAsFactors = FALSE
    ),
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
      seed = seed
    ),
    features = input$features,
    cells = input$cells
  )
  if (isTRUE(store_giotto)) {
    tool$giotto <- gobject
  }
  srt@tools[[tool_name]] <- tool

  log_message(
    "{.pkg Giotto} clusters stored in metadata column {.val {cluster_colname}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

giotto_require <- function(verbose = TRUE) {
  if (!requireNamespace("Giotto", quietly = TRUE)) {
    status <- check_r("giotto-suite/Giotto", verbose = FALSE)
    if (!isTRUE(unname(unlist(status))[1]) || !requireNamespace("Giotto", quietly = TRUE)) {
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
  coord.cols = c("x", "y")
) {
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
  cells <- intersect(colnames(srt), rownames(coords))
  if (length(cells) < 3L) {
    log_message(
      "At least three Seurat cells/spots with spatial coordinates are required",
      message_type = "error"
    )
  }
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
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

  raw_expr <- tryCatch(
    GetAssayData5(srt, assay = assay, layer = "counts")[features, cells, drop = FALSE],
    error = function(e) expr
  )
  coords <- coords[cells, c("x", "y"), drop = FALSE]
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
  args <- giotto_merge_args(args, conversion_params)
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

giotto_add_cluster_metadata <- function(
  srt,
  clusters,
  cluster_colname
) {
  cluster_df <- data.frame(
    cluster = as.character(clusters),
    row.names = names(clusters),
    stringsAsFactors = FALSE
  )
  colnames(cluster_df) <- cluster_colname
  Seurat::AddMetaData(srt, metadata = cluster_df)
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
  if (!requireNamespace("Giotto", quietly = TRUE)) {
    log_message(
      "Please install {.pkg Giotto} before running {.fn RunGiottoCluster}",
      message_type = "error"
    )
  }
  exported <- tryCatch(getExportedValue("Giotto", name), error = function(e) NULL)
  if (is.function(exported)) {
    return(exported)
  }
  ns <- asNamespace("Giotto")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  pkg_env <- tryCatch(as.environment("package:Giotto"), error = function(e) NULL)
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
    args <- args[names(args) %in% fmls]
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

giotto_merge_args <- function(defaults, extra) {
  if (length(extra) == 0L) {
    return(defaults)
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
