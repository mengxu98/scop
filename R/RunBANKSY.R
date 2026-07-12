#' @title Run BANKSY spatial clustering
#'
#' @description
#' Build neighborhood-augmented BANKSY features from a spatial `Seurat` object
#' and store spatial domain or microenvironment clusters in metadata.
#'
#' @md
#' @inheritParams RunBayesSpace
#' @param layer Assay layer used as BANKSY input.
#' @param features Optional features to use. If `NULL`, all assay features are
#' used after zero-count filtering.
#' @param coord.cols Metadata coordinate columns used when no image coordinate
#' source is available.
#' @param lambda BANKSY spatial weighting parameter.
#' @param k_geom Number of spatial neighbors used by BANKSY.
#' @param M Highest azimuthal Fourier harmonic passed to BANKSY.
#' @param npcs Number of principal components to compute.
#' @param use_agf Whether to use azimuthal Gabor filters.
#' @param algo Clustering algorithm passed to `Banksy::clusterBanksy()`.
#' @param k_neighbors Number of neighbors for graph clustering.
#' @param resolution Graph clustering resolution.
#' @param group Optional metadata column used by BANKSY for multi-sample
#' scaling. It is copied into the `SpatialExperiment` colData.
#' @param seed Optional seed for PCA and clustering.
#' @param compute_banksy_params Additional parameters passed to
#' `Banksy::computeBanksy()`.
#' @param run_pca_params Additional parameters passed to
#' `Banksy::runBanksyPCA()`.
#' @param cluster_banksy_params Additional parameters passed to
#' `Banksy::clusterBanksy()`.
#' @param cluster_source Optional BANKSY `colData` column to copy. If `NULL`,
#' the first cluster name reported by `Banksy::clusterNames()` is used when
#' available.
#' @param cluster_colname Metadata column used for BANKSY clusters.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param store_results Whether to store detailed BANKSY results in
#' `srt@tools`.
#' @param coordinate_space Coordinate space used for BANKSY spatial input.
#'   The default preserves the historical coordinate behavior.
#'
#' @return A `Seurat` object with BANKSY clusters in metadata. When
#' `store_results = TRUE`, detailed results are stored in
#' `srt@tools[[tool_name]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial$BANKSY_cluster <- factor(
#'   paste0("BANKSY", (seq_len(ncol(spatial)) - 1) %% 3 + 1)
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "BANKSY_cluster",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' spatial <- RunBANKSY(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "counts",
#'   coord.cols = c("x", "y"),
#'   features = rownames(spatial)[1:300],
#'   lambda = 0.2,
#'   k_geom = 8,
#'   resolution = 0.6,
#'   verbose = FALSE
#' )
RunBANKSY <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  lambda = 0.2,
  k_geom = 15,
  M = 1,
  npcs = 20,
  use_agf = FALSE,
  algo = "leiden",
  k_neighbors = 50,
  resolution = 0.6,
  group = NULL,
  seed = 1,
  compute_banksy_params = list(),
  run_pca_params = list(),
  cluster_banksy_params = list(),
  cluster_source = NULL,
  cluster_colname = "BANKSY_cluster",
  tool_name = "BANKSY",
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("legacy_display", "raw")
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  banksy_assert_scalar_string(cluster_colname, "cluster_colname")
  banksy_assert_scalar_string(tool_name, "tool_name")
  banksy_validate_param_list(compute_banksy_params, "compute_banksy_params")
  banksy_validate_param_list(run_pca_params, "run_pca_params")
  banksy_validate_param_list(cluster_banksy_params, "cluster_banksy_params")
  coordinate_space <- match.arg(coordinate_space)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  features_use <- features %||% rownames(srt[[assay]])
  features_use <- intersect(features_use, rownames(srt[[assay]]))
  if (length(features_use) == 0L) {
    log_message(
      "No features are available for {.fn RunBANKSY}",
      message_type = "error"
    )
  }
  expr <- banksy_get_matrix(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features_use
  )
  keep_features <- Matrix::rowSums(expr) > 0
  keep_spots <- Matrix::colSums(expr) > 0
  expr <- expr[keep_features, keep_spots, drop = FALSE]
  if (nrow(expr) == 0L || ncol(expr) == 0L) {
    log_message(
      "No non-zero features or spots remain for {.fn RunBANKSY}",
      message_type = "error"
    )
  }
  coords <- rctd_get_spatial_coords(
    srt = srt,
    spot_ids = colnames(expr),
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )
  coldata <- banksy_coldata(
    srt = srt,
    spot_ids = colnames(expr),
    coords = coords,
    group = group
  )

  log_message(
    "Run {.pkg Banksy} with {.val {nrow(expr)}} features and {.val {ncol(expr)}} spatial spots",
    verbose = verbose
  )
  backend <- banksy_run_backend(
    expr = expr,
    coords = coords,
    coldata = coldata,
    assay_name = "scop_input",
    lambda = lambda,
    k_geom = k_geom,
    M = M,
    npcs = npcs,
    use_agf = use_agf,
    algo = algo,
    k_neighbors = k_neighbors,
    resolution = resolution,
    group = group,
    seed = seed,
    compute_banksy_params = compute_banksy_params,
    run_pca_params = run_pca_params,
    cluster_banksy_params = cluster_banksy_params
  )

  cluster_source <- banksy_resolve_cluster_source(
    se = backend$se,
    before_cols = backend$before_cols,
    cluster_source = cluster_source
  )
  cdata <- as.data.frame(SummarizedExperiment::colData(backend$se))
  clusters <- as.character(cdata[colnames(expr), cluster_source, drop = TRUE])
  cluster_df <- data.frame(
    BANKSY_cluster = clusters,
    row.names = colnames(expr),
    stringsAsFactors = FALSE
  )
  colnames(cluster_df) <- cluster_colname
  srt <- Seurat::AddMetaData(srt, metadata = cluster_df)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      clusters = cluster_df,
      cluster_source = cluster_source,
      colData = cdata,
      coords = coords,
      features = rownames(expr),
      se = backend$se,
      summary = list(
        n_spots = nrow(cluster_df),
        domains = scop_spatial_domain_summary(cluster_df[[cluster_colname]])
      ),
      parameters = list(
        assay = assay,
        layer = layer,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        lambda = lambda,
        k_geom = k_geom,
        M = M,
        npcs = npcs,
        use_agf = use_agf,
        algo = algo,
        k_neighbors = k_neighbors,
        resolution = resolution,
        group = group,
        seed = seed,
        compute_banksy_params = compute_banksy_params,
        run_pca_params = run_pca_params,
        cluster_banksy_params = cluster_banksy_params,
        cluster_source = cluster_source,
        cluster_colname = cluster_colname,
        tool_name = tool_name
      )
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "BANKSY",
      result_type = "domain",
      source = c(
        attr(coords, "spatial_source") %||% list(),
        list(transform = attr(coords, "spatial_transform"))
      ),
      provenance = list(producer = "RunBANKSY", backend_id = "banksy")
    )
  }

  log_message(
    "{.pkg BANKSY} clusters stored in metadata column {.val {cluster_colname}}",
    verbose = verbose
  )
  srt
}

banksy_run_backend <- function(
  expr,
  coords,
  coldata,
  assay_name,
  lambda,
  k_geom,
  M,
  npcs,
  use_agf,
  algo,
  k_neighbors,
  resolution,
  group,
  seed,
  compute_banksy_params,
  run_pca_params,
  cluster_banksy_params
) {
  check_r(
    c("Banksy", "SpatialExperiment", "SummarizedExperiment", "S4Vectors"),
    verbose = FALSE
  )
  spatial_experiment <- get_namespace_fun("SpatialExperiment", "SpatialExperiment")
  se <- spatial_experiment(
    assays = stats::setNames(list(expr), assay_name),
    colData = S4Vectors::DataFrame(coldata),
    spatialCoords = as.matrix(coords)
  )
  before_cols <- colnames(as.data.frame(SummarizedExperiment::colData(se)))
  compute_fun <- get_namespace_fun("Banksy", "computeBanksy")
  pca_fun <- get_namespace_fun("Banksy", "runBanksyPCA")
  cluster_fun <- get_namespace_fun("Banksy", "clusterBanksy")

  compute_args <- c(
    list(
      assay_name = assay_name,
      coord_names = c("x", "y"),
      compute_agf = use_agf,
      M = M,
      k_geom = k_geom
    ),
    compute_banksy_params
  )
  se <- banksy_do_call(compute_fun, se, compute_args)

  pca_args <- c(
    list(
      assay_name = assay_name,
      M = M,
      lambda = lambda,
      npcs = npcs,
      use_agf = use_agf,
      group = group,
      seed = seed
    ),
    run_pca_params
  )
  se <- banksy_do_call(pca_fun, se, pca_args)

  cluster_args <- c(
    list(
      assay_name = assay_name,
      M = M,
      lambda = lambda,
      use_agf = use_agf,
      npcs = npcs,
      algo = algo,
      k_neighbors = k_neighbors,
      resolution = resolution,
      group = group,
      seed = seed
    ),
    cluster_banksy_params
  )
  se <- banksy_do_call(cluster_fun, se, cluster_args)
  list(se = se, before_cols = before_cols)
}

banksy_get_matrix <- function(srt, assay, layer, features) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  mat <- mat[features, , drop = FALSE]
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(
      if (is.data.frame(mat)) as.matrix(mat) else mat,
      sparse = TRUE
    )
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  mat@x[!is.finite(mat@x)] <- 0
  Matrix::drop0(mat)
}

banksy_coldata <- function(srt, spot_ids, coords, group = NULL) {
  coldata <- data.frame(
    x = coords$x,
    y = coords$y,
    row.names = spot_ids,
    stringsAsFactors = FALSE
  )
  if (!is.null(group)) {
    if (length(group) != 1L || !group %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group} must be a single metadata column in {.arg srt}",
        message_type = "error"
      )
    }
    coldata[[group]] <- srt@meta.data[spot_ids, group, drop = TRUE]
  }
  coldata
}

banksy_resolve_cluster_source <- function(se, before_cols, cluster_source = NULL) {
  cdata <- as.data.frame(SummarizedExperiment::colData(se))
  if (!is.null(cluster_source)) {
    if (length(cluster_source) != 1L || !cluster_source %in% colnames(cdata)) {
      log_message(
        "{.arg cluster_source} must be a single BANKSY colData column",
        message_type = "error"
      )
    }
    return(cluster_source)
  }
  cluster_names <- tryCatch(
    get_namespace_fun("Banksy", "clusterNames")(se),
    error = function(e) character()
  )
  cluster_names <- cluster_names[cluster_names %in% colnames(cdata)]
  if (length(cluster_names) > 0L) {
    return(cluster_names[1L])
  }
  new_cols <- setdiff(colnames(cdata), before_cols)
  candidate <- new_cols[vapply(cdata[, new_cols, drop = FALSE], function(x) {
    is.factor(x) || is.character(x) || is.numeric(x)
  }, logical(1))]
  if (length(candidate) == 0L) {
    log_message(
      "{.pkg Banksy} did not add a detectable cluster column",
      message_type = "error"
    )
  }
  candidate[1L]
}

banksy_do_call <- function(fun, se, args) {
  args <- args[!vapply(args, is.null, logical(1))]
  fmls <- names(formals(fun))
  if (!"..." %in% fmls) {
    args <- args[names(args) %in% fmls]
  }
  do.call(fun, c(list(se), args))
}

banksy_validate_param_list <- function(x, arg_name) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg_name}} must be a list",
      message_type = "error"
    )
  }
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      "{.arg {arg_name}} must contain named arguments only",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

banksy_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
}
