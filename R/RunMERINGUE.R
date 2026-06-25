#' @title Run MERINGUE spatial autocorrelation analysis
#'
#' @description
#' Run `MERINGUE` spatial autocorrelation, spatial cross-correlation, and
#' spatial module analysis for a spatial `Seurat` object.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @param mode MERINGUE analysis modes to run. `"autocorrelation"` computes
#' spatial autocorrelation, `"cross_correlation"` computes pairwise spatial
#' cross-correlation, and `"modules"` detects spatial gene modules.
#' @param filterDist Euclidean distance cutoff passed to
#' `MERINGUE::getSpatialNeighbors()`.
#' @param binary Whether to binarize the MERINGUE spatial neighbor matrix.
#' @param alternative Alternative hypothesis passed to MERINGUE Moran tests.
#' @param ncores Number of cores passed to MERINGUE permutation tests.
#' @param pairwise_features Features used for spatial cross-correlation. If
#' `NULL`, top spatially autocorrelated features are used.
#' @param neighbor_params,moran_params,cross_cor_params,module_params Named
#' lists of additional arguments passed to the corresponding MERINGUE steps.
#'
#' @return A `Seurat` object with MERINGUE results stored in
#' `srt@tools[["MERINGUE"]]` and top autocorrelated features stored in
#' `srt@misc[["MERINGUEFeatures"]]`.
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
#'
#' spatial <- RunMERINGUE(
#'   spatial,
#'   assay = "Spatial",
#'   mode = c("autocorrelation", "cross_correlation"),
#'   nfeatures = 50
#' )
#'
#' head(spatial@tools[["MERINGUE"]]$autocorrelation)
#' SpatialSpotPlot(
#'   spatial,
#'   features = spatial@misc[["MERINGUEFeatures"]][1:2]
#' )
#' }
RunMERINGUE <- function(
  srt,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("col", "row"),
  features = NULL,
  mode = c("autocorrelation", "cross_correlation", "modules"),
  nfeatures = 2000,
  min_spots = 5,
  filterDist = NA_real_,
  binary = TRUE,
  alternative = "greater",
  nperm = 0,
  ncores = 1,
  pairwise_features = NULL,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  seed = 11,
  neighbor_params = list(),
  moran_params = list(),
  cross_cor_params = list(),
  module_params = list()
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
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
  mode <- match.arg(
    mode,
    choices = c("autocorrelation", "cross_correlation", "modules"),
    several.ok = TRUE
  )
  nfeatures <- meringue_check_positive_integer(nfeatures, "nfeatures")
  min_spots <- meringue_check_positive_integer(min_spots, "min_spots")
  nperm <- meringue_check_nonnegative_integer(nperm, "nperm")
  ncores <- meringue_check_positive_integer(ncores, "ncores")
  meringue_check_scalar_logical(binary, "binary")
  meringue_check_scalar_logical(set_variable_features, "set_variable_features")
  meringue_check_scalar_logical(store_results, "store_results")
  if (length(filterDist) != 1L || (!is.numeric(filterDist) && !is.na(filterDist))) {
    log_message("{.arg filterDist} must be a single numeric value", message_type = "error")
  }
  filterDist <- as.numeric(filterDist)
  meringue_validate_named_param_list(neighbor_params, "neighbor_params")
  meringue_validate_named_param_list(moran_params, "moran_params")
  meringue_validate_named_param_list(cross_cor_params, "cross_cor_params")
  meringue_validate_named_param_list(module_params, "module_params")

  log_message(
    "Running MERINGUE spatial analysis",
    message_type = "running",
    verbose = verbose
  )
  meringue_require_package("MERINGUE")
  set.seed(seed)

  inputs <- meringue_prepare_inputs(
    srt = srt,
    assay = assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols,
    features = features,
    min_spots = min_spots
  )
  expr <- inputs$expr
  coords <- inputs$coords
  expressed_spots <- inputs$expressed_spots

  weight <- meringue_run_neighbors(
    coords = coords,
    filterDist = filterDist,
    binary = binary,
    verbose = verbose,
    neighbor_params = neighbor_params
  )

  need_autocorrelation <- any(mode %in% c("autocorrelation", "modules")) ||
    (is.null(pairwise_features) && "cross_correlation" %in% mode)
  autocorrelation <- meringue_empty_autocorrelation()
  if (isTRUE(need_autocorrelation)) {
    autocorrelation <- meringue_run_autocorrelation(
      expr = expr,
      weight = weight,
      expressed_spots = expressed_spots,
      alternative = alternative,
      nperm = nperm,
      ncores = ncores,
      seed = seed,
      moran_params = moran_params
    )
  }

  pairwise_features <- meringue_resolve_pairwise_features(
    pairwise_features = pairwise_features,
    autocorrelation = autocorrelation,
    expr = expr,
    nfeatures = nfeatures
  )
  cross_correlation <- meringue_empty_cross_correlation()
  cross_correlation_matrix <- NULL
  if ("cross_correlation" %in% mode) {
    cross_out <- meringue_run_cross_correlation(
      expr = expr,
      weight = weight,
      features = pairwise_features,
      ncores = ncores,
      cross_cor_params = cross_cor_params
    )
    cross_correlation <- cross_out$result
    cross_correlation_matrix <- cross_out$matrix
  }

  modules <- meringue_empty_modules()
  if ("modules" %in% mode) {
    module_features <- meringue_resolve_module_features(
      autocorrelation = autocorrelation,
      expr = expr,
      nfeatures = nfeatures
    )
    module_cross <- cross_correlation_matrix
    if (is.null(module_cross) || !all(module_features %in% rownames(module_cross))) {
      module_cross <- meringue_run_cross_correlation_matrix(
        expr = expr,
        weight = weight,
        features = module_features
      )
    } else {
      module_cross <- module_cross[module_features, module_features, drop = FALSE]
    }
    modules <- meringue_run_modules(
      coords = coords,
      expr = expr,
      features = module_features,
      scc = module_cross,
      verbose = verbose,
      module_params = module_params
    )
  }

  top_features <- utils::head(
    autocorrelation$feature[is.finite(autocorrelation$score)],
    nfeatures
  )
  if (length(top_features) > 0L) {
    srt@misc[["MERINGUEFeatures"]] <- top_features
    if (isTRUE(set_variable_features)) {
      SeuratObject::VariableFeatures(srt, assay = assay) <- top_features
    }
  }

  parameters <- meringue_parameters_df(list(
    assay = assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols,
    mode = mode,
    nfeatures = nfeatures,
    min_spots = min_spots,
    filterDist = filterDist,
    binary = binary,
    alternative = alternative,
    nperm = nperm,
    ncores = ncores,
    pairwise_features = pairwise_features,
    set_variable_features = set_variable_features,
    store_results = store_results,
    seed = seed
  ))
  if (isTRUE(store_results)) {
    srt@tools[["MERINGUE"]] <- list(
      autocorrelation = autocorrelation,
      cross_correlation = cross_correlation,
      modules = modules,
      coords = coords,
      weight = weight,
      features = rownames(expr),
      pairwise_features = pairwise_features,
      parameters = parameters
    )
  }
  log_message(
    "Stored {.val {length(top_features)}} MERINGUE spatial features",
    message_type = "success",
    verbose = verbose
  )
  srt
}

meringue_prepare_inputs <- function(
  srt,
  assay,
  layer,
  image,
  coord.cols,
  features = NULL,
  min_spots = 5
) {
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
  spots <- intersect(colnames(srt), rownames(coords))
  if (length(spots) == 0L) {
    log_message(
      "No spatial coordinates match spots in {.arg srt}",
      message_type = "error"
    )
  }
  coords <- coords[spots, , drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  coords <- coords[keep_coords, , drop = FALSE]
  spots <- rownames(coords)
  if (length(spots) < 3L) {
    log_message(
      "At least three spots with finite coordinates are required",
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
  features <- unique(as.character(features))
  features <- intersect(features, rownames(expr))
  if (length(features) == 0L) {
    log_message(
      "No requested {.arg features} are present in assay {.val {assay}}",
      message_type = "error"
    )
  }
  expr <- expr[features, spots, drop = FALSE]
  expressed_spots <- Matrix::rowSums(expr > 0)
  keep_features <- expressed_spots >= min_spots
  if (!any(keep_features)) {
    log_message(
      "No features are expressed in at least {.val {min_spots}} spots",
      message_type = "error"
    )
  }
  expr <- expr[keep_features, , drop = FALSE]
  expressed_spots <- expressed_spots[keep_features]
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  expr[!is.finite(expr)] <- 0
  list(
    expr = expr,
    coords = coords,
    expressed_spots = expressed_spots
  )
}

meringue_run_neighbors <- function(
  coords,
  filterDist = NA,
  binary = TRUE,
  verbose = TRUE,
  neighbor_params = list()
) {
  args <- utils::modifyList(
    list(
      pos = as.matrix(coords[, c("x", "y"), drop = FALSE]),
      filterDist = filterDist,
      binary = binary,
      verbose = verbose
    ),
    neighbor_params
  )
  do.call(meringue_get_fun("getSpatialNeighbors"), args)
}

meringue_run_autocorrelation <- function(
  expr,
  weight,
  expressed_spots,
  alternative = "greater",
  nperm = 0,
  ncores = 1,
  seed = 11,
  moran_params = list()
) {
  rows <- lapply(seq_len(nrow(expr)), function(i) {
    feature <- rownames(expr)[[i]]
    out <- if (nperm > 0L) {
      args <- utils::modifyList(
        list(
          z = expr[i, ],
          w = weight,
          alternative = alternative,
          N = nperm,
          seed = seed,
          ncores = ncores,
          plot = FALSE
        ),
        moran_params
      )
      do.call(meringue_get_fun("moranPermutationTest"), args)
    } else {
      args <- utils::modifyList(
        list(
          x = expr[i, ],
          weight = weight,
          alternative = alternative
        ),
        moran_params
      )
      do.call(meringue_get_fun("moranTest"), args)
    }
    meringue_normalize_moran_result(out, feature = feature)
  })
  result <- do.call(rbind, rows)
  result$q_value <- if (all(is.na(result$p_value))) {
    NA_real_
  } else {
    stats::p.adjust(result$p_value, method = "BH")
  }
  result$score <- result$statistic
  if (!any(is.finite(result$score)) && any(is.finite(result$p_value))) {
    result$score <- -log10(pmax(result$p_value, .Machine$double.xmin))
  }
  result$mean <- rowMeans(expr[result$feature, , drop = FALSE])
  result$variance <- fast_row_vars(expr[result$feature, , drop = FALSE])
  result$n_spots <- as.integer(expressed_spots[result$feature])
  score_order <- ifelse(is.finite(result$score), result$score, -Inf)
  p_order <- ifelse(is.finite(result$p_value), result$p_value, Inf)
  q_order <- ifelse(is.finite(result$q_value), result$q_value, Inf)
  result <- result[order(-score_order, p_order, q_order), , drop = FALSE]
  result$rank <- seq_len(nrow(result))
  rownames(result) <- NULL
  meringue_reorder_cols(
    result,
    c(
      "feature", "rank", "statistic", "expected", "sd", "p_value",
      "q_value", "score", "mean", "variance", "n_spots"
    )
  )
}

meringue_normalize_moran_result <- function(out, feature) {
  values <- meringue_flatten_numeric(out)
  names_lower <- tolower(names(values) %||% character(length(values)))
  statistic <- meringue_pick_named_numeric(
    values,
    names_lower,
    c("observed", "statistic", "moran.i", "morani", "i"),
    fallback_index = 1L
  )
  expected <- meringue_pick_named_numeric(
    values,
    names_lower,
    c("expected", "expectation", "null", "ei"),
    fallback_index = 2L
  )
  sd <- meringue_pick_named_numeric(
    values,
    names_lower,
    c("sd", "standard.deviation", "standard_deviation", "se", "std"),
    fallback_index = 3L
  )
  p_value <- meringue_pick_named_numeric(
    values,
    names_lower,
    c("p.value", "p_value", "pval", "pv", "p"),
    fallback_index = 4L
  )
  data.frame(
    feature = feature,
    statistic = statistic,
    expected = expected,
    sd = sd,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

meringue_flatten_numeric <- function(x) {
  if (is.data.frame(x)) {
    x <- as.list(x[1L, , drop = TRUE])
  }
  if (is.list(x)) {
    values <- unlist(x, recursive = TRUE, use.names = TRUE)
  } else {
    values <- x
  }
  suppressWarnings(as.numeric(values))
}

meringue_pick_named_numeric <- function(values, names_lower, candidates, fallback_index) {
  hit <- which(names_lower %in% tolower(candidates))
  if (length(hit) > 0L) {
    return(values[[hit[[1L]]]])
  }
  if (length(values) >= fallback_index) {
    return(values[[fallback_index]])
  }
  NA_real_
}

meringue_resolve_pairwise_features <- function(
  pairwise_features,
  autocorrelation,
  expr,
  nfeatures
) {
  if (is.null(pairwise_features)) {
    if (is.data.frame(autocorrelation) && nrow(autocorrelation) > 0L) {
      pairwise_features <- utils::head(autocorrelation$feature, nfeatures)
    } else {
      pairwise_features <- utils::head(rownames(expr), nfeatures)
    }
  }
  pairwise_features <- unique(as.character(pairwise_features))
  pairwise_features <- intersect(pairwise_features, rownames(expr))
  if (length(pairwise_features) < 2L) {
    log_message(
      "{.arg pairwise_features} must contain at least two assay features",
      message_type = "error"
    )
  }
  pairwise_features
}

meringue_run_cross_correlation <- function(
  expr,
  weight,
  features,
  ncores = 1,
  cross_cor_params = list()
) {
  scc <- meringue_run_cross_correlation_matrix(
    expr = expr,
    weight = weight,
    features = features
  )
  result <- meringue_cross_matrix_to_df(scc)
  if (isTRUE(cross_cor_params$test)) {
    test_params <- cross_cor_params[setdiff(names(cross_cor_params), "test")]
    result$p_value <- meringue_run_cross_correlation_tests(
      expr = expr,
      weight = weight,
      result = result,
      ncores = ncores,
      test_params = test_params
    )
  } else {
    result$p_value <- NA_real_
  }
  list(result = result, matrix = scc)
}

meringue_run_cross_correlation_matrix <- function(expr, weight, features) {
  features <- intersect(features, rownames(expr))
  if (length(features) < 2L) {
    log_message(
      "At least two features are required for MERINGUE spatial cross-correlation",
      message_type = "error"
    )
  }
  spatial_cross_cor <- meringue_get_fun("spatialCrossCorMatrix")
  scc <- spatial_cross_cor(
    mat = expr[features, , drop = FALSE],
    weight = weight
  )
  scc <- as.matrix(scc)
  rownames(scc) <- rownames(scc) %||% features
  colnames(scc) <- colnames(scc) %||% features
  scc
}

meringue_cross_matrix_to_df <- function(scc) {
  if (nrow(scc) < 2L || ncol(scc) < 2L) {
    return(meringue_empty_cross_correlation())
  }
  idx <- which(upper.tri(scc), arr.ind = TRUE)
  result <- data.frame(
    feature1 = rownames(scc)[idx[, "row"]],
    feature2 = colnames(scc)[idx[, "col"]],
    correlation = as.numeric(scc[idx]),
    stringsAsFactors = FALSE
  )
  score_order <- ifelse(is.finite(result$correlation), abs(result$correlation), -Inf)
  result <- result[order(-score_order), , drop = FALSE]
  result$rank <- seq_len(nrow(result))
  rownames(result) <- NULL
  meringue_reorder_cols(result, c("feature1", "feature2", "rank", "correlation"))
}

meringue_run_cross_correlation_tests <- function(
  expr,
  weight,
  result,
  ncores = 1,
  test_params = list()
) {
  cross_test <- meringue_get_fun("spatialCrossCorTest")
  vapply(seq_len(nrow(result)), function(i) {
    args <- utils::modifyList(
      list(
        x = expr[result$feature1[[i]], ],
        y = expr[result$feature2[[i]], ],
        w = weight,
        ncores = ncores,
        plot = FALSE
      ),
      test_params
    )
    out <- do.call(cross_test, args)
    val <- meringue_flatten_numeric(out)
    if (length(val) == 0L) {
      NA_real_
    } else {
      val[[1L]]
    }
  }, numeric(1))
}

meringue_resolve_module_features <- function(autocorrelation, expr, nfeatures) {
  features <- if (is.data.frame(autocorrelation) && nrow(autocorrelation) > 0L) {
    autocorrelation$feature[is.finite(autocorrelation$score)]
  } else {
    rownames(expr)
  }
  features <- utils::head(unique(features), nfeatures)
  features <- intersect(features, rownames(expr))
  if (length(features) < 2L) {
    log_message(
      "At least two features are required for MERINGUE module analysis",
      message_type = "error"
    )
  }
  features
}

meringue_run_modules <- function(
  coords,
  expr,
  features,
  scc,
  verbose = TRUE,
  module_params = list()
) {
  args <- utils::modifyList(
    list(
      pos = as.matrix(coords[, c("x", "y"), drop = FALSE]),
      mat = expr[features, , drop = FALSE],
      scc = scc,
      plot = FALSE,
      verbose = verbose
    ),
    module_params
  )
  out <- tryCatch(
    do.call(meringue_get_fun("groupSigSpatialPatterns"), args),
    error = function(e) {
      log_message(
        "{.pkg MERINGUE} module analysis failed: {.val {conditionMessage(e)}}",
        message_type = "error"
      )
    }
  )
  meringue_normalize_modules(out, features = features)
}

meringue_normalize_modules <- function(out, features) {
  if (is.data.frame(out)) {
    df <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    df <- meringue_rename_first(df, "feature", c("feature", "features", "gene", "genes"))
    df <- meringue_rename_first(df, "module", c("module", "cluster", "clusters", "pattern", "group"))
  } else if (is.list(out) && any(c("groups", "clusters", "cluster", "modules", "module") %in% names(out))) {
    nm <- intersect(c("groups", "clusters", "cluster", "modules", "module"), names(out))[[1L]]
    df <- meringue_module_vector_to_df(out[[nm]], features = features)
  } else {
    df <- meringue_module_vector_to_df(out, features = features)
  }
  if (!"feature" %in% colnames(df)) {
    df$feature <- features[seq_len(nrow(df))]
  }
  if (!"module" %in% colnames(df)) {
    df$module <- NA
  }
  df$feature <- as.character(df$feature)
  df <- df[df$feature %in% features, , drop = FALSE]
  df$module <- as.character(df$module)
  module_size <- table(df$module)
  df$module_size <- as.integer(module_size[df$module])
  df$rank <- seq_len(nrow(df))
  rownames(df) <- NULL
  meringue_reorder_cols(df, c("feature", "module", "module_size", "rank"))
}

meringue_module_vector_to_df <- function(x, features) {
  x <- unlist(x, recursive = TRUE, use.names = TRUE)
  if (length(x) == 0L) {
    return(meringue_empty_modules())
  }
  feature <- names(x)
  if (is.null(feature) || all(!nzchar(feature))) {
    feature <- utils::head(features, length(x))
  }
  data.frame(
    feature = feature,
    module = as.character(x),
    stringsAsFactors = FALSE
  )
}

meringue_empty_autocorrelation <- function() {
  data.frame(
    feature = character(),
    rank = integer(),
    statistic = numeric(),
    expected = numeric(),
    sd = numeric(),
    p_value = numeric(),
    q_value = numeric(),
    score = numeric(),
    mean = numeric(),
    variance = numeric(),
    n_spots = integer(),
    stringsAsFactors = FALSE
  )
}

meringue_empty_cross_correlation <- function() {
  data.frame(
    feature1 = character(),
    feature2 = character(),
    rank = integer(),
    correlation = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
}

meringue_empty_modules <- function() {
  data.frame(
    feature = character(),
    module = character(),
    module_size = integer(),
    rank = integer(),
    stringsAsFactors = FALSE
  )
}

meringue_require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_message(
      "Please install required package before running {.fn RunMERINGUE}: {.val {pkg}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

meringue_get_fun <- function(fun) {
  meringue_require_package("MERINGUE")
  tryCatch(
    getExportedValue("MERINGUE", fun),
    error = function(e) {
      log_message(
        "{.pkg MERINGUE} does not export required function {.fn {fun}}",
        message_type = "error"
      )
    }
  )
}

meringue_parameters_df <- function(parameters) {
  parameters$scop_version <- meringue_package_version("scop")
  parameters$MERINGUE_version <- meringue_package_version("MERINGUE")
  data.frame(
    key = names(parameters),
    value = vapply(parameters, meringue_collapse_value, character(1)),
    stringsAsFactors = FALSE
  )
}

meringue_package_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

meringue_collapse_value <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (length(x) == 0L) {
    return("")
  }
  paste(as.character(x), collapse = ",")
}

meringue_rename_first <- function(df, target, candidates) {
  hit <- intersect(candidates, colnames(df))
  hit <- setdiff(hit, target)
  if (length(hit) > 0L && !target %in% colnames(df)) {
    colnames(df)[match(hit[[1L]], colnames(df))] <- target
  }
  df
}

meringue_reorder_cols <- function(df, first_cols) {
  first_cols <- intersect(first_cols, colnames(df))
  df[, c(first_cols, setdiff(colnames(df), first_cols)), drop = FALSE]
}

meringue_check_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1) {
    log_message(
      "{.arg {arg_name}} must be a positive number",
      message_type = "error"
    )
  }
  as.integer(x)
}

meringue_check_nonnegative_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
    log_message(
      "{.arg {arg_name}} must be a non-negative number",
      message_type = "error"
    )
  }
  as.integer(x)
}

meringue_check_scalar_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg_name}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

meringue_validate_named_param_list <- function(x, arg_name) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg_name}} must be a named list",
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
