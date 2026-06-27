#' @title Run smoothclust spatial domain clustering
#'
#' @description
#' Smooth expression across spatial neighborhoods with the optional
#' `smoothclust` package, then cluster the smoothed profiles into spatial
#' domains with PCA and k-means.
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay used for expression. If `NULL`, the default assay is used.
#' @param layer Assay layer used for expression values.
#' @param image Name of the Seurat spatial image. If `NULL`, the first image is
#' used when present.
#' @param coord.cols Metadata coordinate columns used when no Seurat image is
#' available.
#' @param features Features to use. If `NULL`, current variable features are
#' used; if no variable features are present, the top `nfeatures` by variance
#' are used.
#' @param nfeatures Number of variance-ranked features to use when
#' `features = NULL` and no variable features are present.
#' @param min_spots Minimum number of spots with non-zero expression required
#' for a feature to be used.
#' @param smooth_method Smoothing method passed to `smoothclust::smoothclust()`.
#' @param bandwidth,k,truncate,n_threads Smoothing parameters passed to
#' `smoothclust::smoothclust()`.
#' @param n_clusters Number of spatial domains for k-means clustering. This
#' must be supplied explicitly.
#' @param n_pcs Number of principal components used for k-means.
#' @param center,scale Whether to center and scale features before PCA.
#' @param nstart,iter.max,algorithm Parameters passed to `stats::kmeans()`.
#' @param cluster_colname Metadata column used for smoothclust clusters.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param store_results Whether to store detailed results in `srt@tools`.
#' @param store_smoothed Whether to store the smoothed expression matrix in
#' `srt@tools[[tool_name]]`. This can be large.
#' @param seed Random seed used for k-means.
#' @param verbose Whether to print progress messages.
#' @param ... Additional arguments passed to `smoothclust::smoothclust()`.
#'
#' @return A `Seurat` object with smoothclust clusters in metadata. When
#' `store_results = TRUE`, detailed outputs are stored in
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:120],
#'   features = rownames(visium_human_pancreas_sub)[1:400]
#' )
#' spatial$SmoothClust_cluster <- factor(
#'   paste0("SmoothClust", (seq_len(ncol(spatial)) - 1) %% 3 + 1)
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "SmoothClust_cluster",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   requireNamespace("smoothclust", quietly = TRUE) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' spatial <- Seurat::FindVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 200,
#'   verbose = FALSE
#' )
#'
#' spatial <- RunSmoothClust(
#'   spatial,
#'   assay = "Spatial",
#'   n_clusters = 3,
#'   smooth_method = "knn",
#'   coord.cols = c("x", "y"),
#'   k = 6,
#'   verbose = FALSE
#' )
#'
#' table(spatial$SmoothClust_cluster)
#' }
RunSmoothClust <- function(
  srt,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("col", "row"),
  features = NULL,
  nfeatures = 2000,
  min_spots = 5,
  smooth_method = c("uniform", "kernel", "knn"),
  bandwidth = 0.05,
  k = 18,
  truncate = 0.05,
  n_threads = 1,
  n_clusters,
  n_pcs = 15,
  center = TRUE,
  scale = TRUE,
  nstart = 10,
  iter.max = 100,
  algorithm = "Hartigan-Wong",
  cluster_colname = "SmoothClust_cluster",
  tool_name = "SmoothClust",
  store_results = TRUE,
  store_smoothed = FALSE,
  seed = 11,
  verbose = TRUE,
  ...
) {
  log_message(
    "Running smoothclust spatial domain clustering",
    message_type = "running",
    verbose = verbose
  )
  smoothclust_validate_srt(srt)
  smooth_method <- match.arg(smooth_method)
  smoothclust_assert_string(cluster_colname, "cluster_colname")
  smoothclust_assert_string(tool_name, "tool_name")
  nfeatures <- smoothclust_assert_positive_integer(nfeatures, "nfeatures")
  min_spots <- smoothclust_assert_positive_integer(min_spots, "min_spots")
  k <- smoothclust_assert_positive_integer(k, "k")
  n_threads <- smoothclust_assert_positive_integer(n_threads, "n_threads")
  if (missing(n_clusters)) {
    log_message(
      "{.arg n_clusters} must be supplied explicitly",
      message_type = "error"
    )
  }
  n_clusters <- smoothclust_assert_positive_integer(
    n_clusters,
    "n_clusters",
    allow_null = FALSE
  )
  n_pcs <- smoothclust_assert_positive_integer(n_pcs, "n_pcs")
  nstart <- smoothclust_assert_positive_integer(nstart, "nstart")
  iter.max <- smoothclust_assert_positive_integer(iter.max, "iter.max")
  smoothclust_assert_flag(center, "center")
  smoothclust_assert_flag(scale, "scale")
  smoothclust_assert_flag(store_results, "store_results")
  smoothclust_assert_flag(store_smoothed, "store_smoothed")
  smoothclust_assert_positive_number(bandwidth, "bandwidth")
  smoothclust_assert_positive_number(truncate, "truncate")
  if (truncate >= 1) {
    log_message(
      "{.arg truncate} must be less than 1",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }

  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
  spots <- intersect(colnames(srt), rownames(coords))
  if (length(spots) == 0L) {
    log_message(
      "No spatial coordinates match cells or spots in {.arg srt}",
      message_type = "error"
    )
  }
  coords <- coords[spots, , drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  if (!all(keep_coords)) {
    log_message(
      "Drop {.val {sum(!keep_coords)}} spots with missing or non-finite spatial coordinates",
      verbose = verbose
    )
  }
  coords <- coords[keep_coords, , drop = FALSE]
  spots <- rownames(coords)
  if (length(spots) < 3L) {
    log_message(
      "At least three spots with finite coordinates are required",
      message_type = "error"
    )
  }
  if (identical(smooth_method, "knn") && k >= length(spots)) {
    log_message(
      "{.arg k} must be smaller than the number of spots for {.arg smooth_method = 'knn'}",
      message_type = "error"
    )
  }
  if (n_clusters > length(spots)) {
    log_message(
      "{.arg n_clusters} must be no larger than the number of clustered spots",
      message_type = "error"
    )
  }

  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  selected <- smoothclust_select_features(
    srt = srt,
    expr = expr,
    assay = assay,
    spots = spots,
    features = features,
    nfeatures = nfeatures,
    min_spots = min_spots
  )
  expr_use <- as.matrix(expr[selected$features, spots, drop = FALSE])
  expr_use[!is.finite(expr_use)] <- 0

  check_r("smoothclust", verbose = FALSE)
  smooth_fun <- smoothclust_get_fun("smoothclust")
  smooth_args <- c(
    list(
      input = expr_use,
      spatial_coords = as.matrix(coords[, c("x", "y"), drop = FALSE]),
      method = smooth_method,
      bandwidth = bandwidth,
      k = k,
      truncate = truncate,
      n_threads = n_threads
    ),
    list(...)
  )
  smooth_args <- smoothclust_filter_backend_args(
    fun = smooth_fun,
    args = smooth_args
  )
  log_message(
    "Run {.pkg smoothclust} with {.val {nrow(expr_use)}} features and {.val {ncol(expr_use)}} spots",
    verbose = verbose
  )
  smoothed <- do.call(smooth_fun, smooth_args)
  smoothed <- smoothclust_normalize_smoothed(
    smoothed = smoothed,
    features = rownames(expr_use),
    spots = colnames(expr_use)
  )

  pca <- smoothclust_run_pca(
    smoothed = smoothed,
    n_pcs = n_pcs,
    center = center,
    scale = scale,
    verbose = verbose
  )
  set.seed(seed)
  km <- stats::kmeans(
    x = pca$embedding,
    centers = n_clusters,
    nstart = nstart,
    iter.max = iter.max,
    algorithm = algorithm
  )
  clusters <- factor(
    paste0("SmoothClust", km$cluster),
    levels = paste0("SmoothClust", sort(unique(km$cluster)))
  )
  names(clusters) <- rownames(pca$embedding)

  cluster_full <- rep(NA_character_, ncol(srt))
  names(cluster_full) <- colnames(srt)
  cluster_full[names(clusters)] <- as.character(clusters)
  cluster_df <- data.frame(
    SmoothClust_cluster = factor(cluster_full, levels = levels(clusters)),
    row.names = colnames(srt)
  )
  colnames(cluster_df) <- cluster_colname
  srt <- Seurat::AddMetaData(srt, metadata = cluster_df)

  smoothness <- smoothclust_smoothness_metric(
    coords = coords[names(clusters), , drop = FALSE],
    labels = as.integer(clusters),
    k = min(6L, length(clusters) - 1L)
  )

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      clusters = data.frame(
        cell = names(clusters),
        cluster = as.character(clusters),
        stringsAsFactors = FALSE
      ),
      coords = coords,
      features = selected$features,
      feature_selection = selected$summary,
      pca = pca,
      kmeans = km,
      smoothness = smoothness,
      parameters = list(
        assay = assay,
        layer = layer,
        image = image,
        coord.cols = coord.cols,
        nfeatures = nfeatures,
        min_spots = min_spots,
        smooth_method = smooth_method,
        bandwidth = bandwidth,
        k = k,
        truncate = truncate,
        n_threads = n_threads,
        n_clusters = n_clusters,
        n_pcs = pca$n_pcs,
        center = center,
        scale = scale,
        nstart = nstart,
        iter.max = iter.max,
        algorithm = algorithm,
        cluster_colname = cluster_colname,
        tool_name = tool_name,
        store_smoothed = store_smoothed,
        seed = seed
      )
    )
    if (isTRUE(store_smoothed)) {
      srt@tools[[tool_name]][["smoothed"]] <- smoothed
    }
  }

  log_message(
    "{.pkg smoothclust} clusters stored in metadata column {.val {cluster_colname}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

smoothclust_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

smoothclust_select_features <- function(
  srt,
  expr,
  assay,
  spots,
  features = NULL,
  nfeatures = 2000,
  min_spots = 5
) {
  source <- "features"
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
    source <- "variable_features"
    if (length(features) == 0L) {
      expr_spots <- expr[, spots, drop = FALSE]
      feature_vars <- fast_row_vars(expr_spots)
      feature_vars[!is.finite(feature_vars)] <- -Inf
      ordered <- order(feature_vars, decreasing = TRUE)
      features <- rownames(expr_spots)[utils::head(ordered, nfeatures)]
      source <- "variance"
    }
  }
  features <- unique(as.character(features))
  features <- features[!is.na(features) & nzchar(features)]
  features <- intersect(features, rownames(expr))
  if (length(features) == 0L) {
    log_message(
      "No requested {.arg features} are present in assay {.val {assay}}",
      message_type = "error"
    )
  }
  if (length(features) > nfeatures) {
    features <- utils::head(features, nfeatures)
  }

  expr_use <- expr[features, spots, drop = FALSE]
  expressed_spots <- Matrix::rowSums(expr_use > 0)
  keep <- expressed_spots >= min_spots
  if (!any(keep)) {
    log_message(
      "No features are expressed in at least {.val {min_spots}} spots",
      message_type = "error"
    )
  }
  features <- features[keep]
  list(
    features = features,
    summary = data.frame(
      feature = features,
      source = source,
      n_spots = as.integer(expressed_spots[keep]),
      stringsAsFactors = FALSE
    )
  )
}

smoothclust_normalize_smoothed <- function(smoothed, features, spots) {
  if (inherits(smoothed, "SpatialExperiment")) {
    log_message(
      "{.pkg smoothclust} returned a SpatialExperiment for matrix input",
      message_type = "error"
    )
  }
  smoothed <- as.matrix(smoothed)
  if (nrow(smoothed) != length(features) || ncol(smoothed) != length(spots)) {
    log_message(
      "{.pkg smoothclust} returned a smoothed matrix with unexpected dimensions",
      message_type = "error"
    )
  }
  rownames(smoothed) <- features
  colnames(smoothed) <- spots
  smoothed[!is.finite(smoothed)] <- 0
  smoothed
}

smoothclust_run_pca <- function(
  smoothed,
  n_pcs = 15,
  center = TRUE,
  scale = TRUE,
  verbose = TRUE
) {
  vars <- fast_row_vars(smoothed)
  keep <- is.finite(vars) & vars > 0
  if (!any(keep)) {
    log_message(
      "Smoothed expression has no variable features for PCA",
      message_type = "error"
    )
  }
  if (!all(keep)) {
    log_message(
      "Drop {.val {sum(!keep)}} features with zero variance after smoothing",
      verbose = verbose
    )
  }
  mat <- t(smoothed[keep, , drop = FALSE])
  max_pcs <- min(nrow(mat), ncol(mat))
  if (max_pcs < 1L) {
    log_message(
      "At least one feature and one spot are required for PCA",
      message_type = "error"
    )
  }
  n_pcs_use <- min(n_pcs, max_pcs)
  pca <- stats::prcomp(
    mat,
    center = center,
    scale. = scale
  )
  embedding <- pca$x[, seq_len(n_pcs_use), drop = FALSE]
  list(
    embedding = embedding,
    rotation = pca$rotation[, seq_len(n_pcs_use), drop = FALSE],
    sdev = pca$sdev[seq_len(n_pcs_use)],
    center = pca$center,
    scale = pca$scale,
    features = colnames(mat),
    n_pcs = n_pcs_use
  )
}

smoothclust_smoothness_metric <- function(coords, labels, k = 6L) {
  if (length(labels) < 2L || k < 1L) {
    return(list(n_discordant = rep(NA_real_, length(labels)), mean_discordant = NA_real_))
  }
  smooth_metric <- tryCatch(
    smoothclust_get_fun("smoothness_metric"),
    error = function(e) NULL
  )
  if (!is.null(smooth_metric)) {
    return(tryCatch(
      smooth_metric(
        spatial_coords = as.matrix(coords[, c("x", "y"), drop = FALSE]),
        labels = labels,
        k = k
      ),
      error = function(e) smoothclust_native_smoothness(coords, labels, k)
    ))
  }
  smoothclust_native_smoothness(coords, labels, k)
}

smoothclust_native_smoothness <- function(coords, labels, k = 6L) {
  coord_mat <- as.matrix(coords[, c("x", "y"), drop = FALSE])
  dmat <- as.matrix(stats::dist(coord_mat))
  diag(dmat) <- Inf
  nn <- t(apply(dmat, 1L, function(x) utils::head(order(x), k)))
  if (k == 1L) {
    nn <- matrix(nn, ncol = 1L)
  }
  n_discordant <- rowSums(matrix(labels[nn], nrow = nrow(nn)) != labels)
  list(
    n_discordant = as.integer(n_discordant),
    mean_discordant = mean(n_discordant)
  )
}

smoothclust_get_fun <- function(fun) {
  tryCatch(
    getExportedValue("smoothclust", fun),
    error = function(e) {
      log_message(
        "{.pkg smoothclust} does not export required function {.fn {fun}}",
        message_type = "error"
      )
    }
  )
}

smoothclust_filter_backend_args <- function(fun, args) {
  formals_use <- tryCatch(names(formals(fun)), error = function(e) NULL)
  if (is.null(formals_use) || "..." %in% formals_use) {
    return(args)
  }
  args[names(args) %in% formals_use]
}

smoothclust_assert_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a non-empty character string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

smoothclust_assert_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

smoothclust_assert_positive_integer <- function(
  x,
  arg,
  allow_null = FALSE
) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }
  if (
    is.null(x) ||
      length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      x < 1 ||
      x != floor(x)
  ) {
    log_message(
      "{.arg {arg}} must be a single positive integer",
      message_type = "error"
    )
  }
  as.integer(x)
}

smoothclust_assert_positive_number <- function(x, arg) {
  if (
    length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0
  ) {
    log_message(
      "{.arg {arg}} must be a single positive number",
      message_type = "error"
    )
  }
  invisible(TRUE)
}
