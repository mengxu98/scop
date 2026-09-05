#' @title Run spatial gradient feature screening
#'
#' @description
#' Run native spatial trajectory or annotation gradient screening for Seurat
#' objects. The compiled C++ backend computes distance-based screening and
#' stores validated result tables in `srt@tools` when requested.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @inheritParams SpatialSpotPlot
#' @inheritParams scop-params
#' @param reference Spatial reference type: `"trajectory"` for STS or
#' `"annotation"` for SAS.
#' @param backend Computation backend. The supported backend is `"cpp"`, which
#' uses the compiled native spatial implementation.
#' @param result_name Name used to store this result. If `NULL`, a name is
#' generated from `reference`.
#' @param variables Numeric variables or genes tested by the native backend. If `NULL`,
#' `srt@tools[["SpatialVariableFeatures"]]` is used first, then variable
#' features, then all assay features.
#' @param coord.cols Metadata coordinate columns used by the `"cpp"`
#' backend when no image coordinates are available.
#' @param coordinate_space Coordinate system used by the C++ distance
#' calculations. The default is raw acquisition coordinates, so `start`,
#' `end`, trajectory positions, and C++ distances share raw
#' coordinate units. Use `"legacy_display"` explicitly for pre-0.9.0 display
#' coordinates. Results retain the raw coordinate units used by the native
#' backend.
#' @param start,end,traj_df Trajectory geometry used by the native backend when
#' `reference = "trajectory"`.
#' @param annotation_ids Reserved compatibility input. It is not supported by
#' the native backend; use `annotation.by` with `annotation.groups`, or
#' `annotation.variable` with `annotation.threshold`.
#' @param annotation.by,annotation.groups Metadata grouping used to create the
#' native annotation mask.
#' @param annotation.variable,annotation.threshold Numeric variable and
#' threshold used to create a native numeric annotation mask. Numeric
#' thresholds are interpreted as `">{threshold}"`.
#' @param sample_name,platform,img_scale_fct,assay_modality,trajectory_id,width,annotation_id,core,distance,angle_span,resolution,model_add,model_subset,model_remove,control
#' Legacy compatibility inputs. The native backend ignores these values and
#' does not include them in the stored effective-parameter summary.
#' @param unit Optional coordinate-unit label stored in the result source. It
#' does not transform coordinates.
#' @param sign_var,sign_threshold Column and cutoff used to select top gradient
#' variables from the native significance table.
#' @param n_random,seed Permutation count and random seed used by the native
#' screening backend.
#' @param n_bins Number of distance bins used for the `"cpp"` backend
#' screening curve.
#' @param min_spots Minimum number of non-zero spots required for a variable in
#' the `"cpp"` backend.
#' @param nfeatures Number of top gradient variables retained in
#' `top_variables` and optionally set as Seurat variable features.
#' @param set_variable_features Whether to set top gradient variables as Seurat
#' variable features independently of `store_results`. If explicitly `TRUE`,
#' an empty selection clears the target assay's variable features.
#' @param store_results Whether to store the normalized result in `srt@tools`.
#' @param ... Additional named arguments accepted for compatibility. The native
#' backend ignores them and does not store them as effective parameters.
#'
#' @return A `Seurat` object. When `store_results = TRUE`, spatial gradient
#' screening results are stored in `srt@tools[["SpatialGradientFeatures"]]`.
#' Otherwise existing stored results are unchanged. Variable features are
#' independently updated only when `set_variable_features = TRUE`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial <- RunSpatialGradientFeatures(
#'   spatial,
#'   reference = "trajectory",
#'   backend = "cpp",
#'   result_name = "ductal_axis",
#'   variables = rownames(spatial)[1:8],
#'   start = c(min(spatial$x), min(spatial$y)),
#'   end = c(max(spatial$x), max(spatial$y)),
#'   layer = "counts",
#'   coord.cols = c("x", "y"),
#'   n_random = 0,
#'   n_bins = 5,
#'   min_spots = 3,
#'   sign_threshold = 1,
#'   nfeatures = 4,
#'   verbose = FALSE
#' )
#'
#' SpatialGradientPlot(
#'   spatial,
#'   plot_type = "surface",
#'   features = rownames(spatial)[1:2],
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y"),
#'   pt.size = 1.2
#' )
RunSpatialGradientFeatures <- function(
  srt,
  reference = c("trajectory", "annotation"),
  backend = "cpp",
  result_name = NULL,
  assay = NULL,
  layer = "data",
  variables = NULL,
  sample_name = NULL,
  platform = "Undefined",
  image = NULL,
  coord.cols = c("x", "y"),
  img_scale_fct = "lowres",
  assay_modality = "gene",
  trajectory_id = "scop_gradient",
  start = NULL,
  end = NULL,
  traj_df = NULL,
  width = NULL,
  annotation_ids = NULL,
  annotation.by = NULL,
  annotation.groups = NULL,
  annotation.variable = NULL,
  annotation.threshold = NULL,
  annotation_id = "scop_gradient",
  core = FALSE,
  distance = "dte",
  angle_span = c(0, 360),
  resolution = NULL,
  unit = NULL,
  sign_var = "fdr",
  sign_threshold = 0.05,
  model_add = NULL,
  model_subset = NULL,
  model_remove = NULL,
  n_random = 10000,
  seed = 123,
  control = NULL,
  n_bins = 50,
  min_spots = 3,
  nfeatures = 2000,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
) {
  for (flag in c("store_results", "set_variable_features")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(paste0(flag, " must be a single non-missing logical value"), call. = FALSE)
    }
  }
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  reference <- match.arg(reference)
  backend <- match.arg(backend)
  coordinate_space <- match.arg(coordinate_space)
  requested_image <- image
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.numeric(nfeatures) || length(nfeatures) != 1L || is.na(nfeatures) || nfeatures < 1) {
    log_message("{.arg nfeatures} must be a positive number", message_type = "error")
  }
  if (!is.numeric(n_random) || length(n_random) != 1L || is.na(n_random) || n_random < 0) {
    log_message("{.arg n_random} must be a non-negative number", message_type = "error")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    log_message("{.arg seed} must be a single numeric value", message_type = "error")
  }
  nfeatures <- as.integer(nfeatures)
  n_random <- as.integer(n_random)
  seed <- as.integer(seed)
  if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) || n_bins < 1) {
    log_message("{.arg n_bins} must be a positive number", message_type = "error")
  }
  if (!is.numeric(min_spots) || length(min_spots) != 1L || is.na(min_spots) || min_spots < 1) {
    log_message("{.arg min_spots} must be a positive number", message_type = "error")
  }
  n_bins <- as.integer(n_bins)
  min_spots <- as.integer(min_spots)

  log_message(
    "Running spatial gradient screening with {.val {backend}} backend",
    message_type = "running",
    verbose = verbose
  )

  variables <- sgf_resolve_variables(
    srt = srt,
    assay = assay,
    layer = layer,
    variables = variables
  )
  result_name <- result_name %||% paste0(reference, "_", format(Sys.time(), "%Y%m%d%H%M%S"))

  if (identical(backend, "cpp")) {
    result <- sgf_run_cpp_gradient(
      srt = srt,
      reference = reference,
      assay = assay,
      layer = layer,
      variables = variables,
      image = image,
      coord.cols = coord.cols,
      coordinate_space = coordinate_space,
      start = start,
      end = end,
      traj_df = traj_df,
      annotation_ids = annotation_ids,
      annotation.by = annotation.by,
      annotation.groups = annotation.groups,
      annotation.variable = annotation.variable,
      annotation.threshold = annotation.threshold,
      n_random = n_random,
      seed = seed,
      n_bins = n_bins,
      min_spots = min_spots,
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      nfeatures = nfeatures,
      parameters = list(
        result_name = result_name,
        reference = reference,
        backend = backend,
        assay = assay,
        layer = layer,
        image = image,
        coord.cols = paste(coord.cols, collapse = ","),
        coordinate_space = coordinate_space,
        annotation.by = annotation.by,
        annotation.groups = paste(annotation.groups %||% character(0), collapse = ","),
        annotation.variable = annotation.variable,
        annotation.threshold = annotation.threshold,
        n_random = n_random,
        seed = seed,
        n_bins = n_bins,
        min_spots = min_spots
      )
    )
  }

  vars <- sgf_validate_result(result)
  if (isTRUE(set_variable_features) && any(!vars %in% rownames(srt[[assay]]))) {
    stop("Spatial gradient result contains variables absent from the target assay", call. = FALSE)
  }
  if (isTRUE(store_results)) {
    source <- result$source %||% list(
      image = image,
      coord.cols = coord.cols[seq_len(min(2L, length(coord.cols)))],
      image_policy = "strict",
      coordinate_space = coordinate_space
    )
    source$selection_strategy <- source$selection_strategy %||% if (is.null(source$image)) {
      "metadata_coord_cols"
    } else if (!is.null(requested_image)) {
      "explicit_image"
    } else {
      "single_available_image"
    }
    source$units <- unit %||% NA_character_
    source$layer <- layer
    srt <- sgf_store_result(
      srt = srt,
      result_name = result_name,
      result = result,
      assay = assay,
      backend = backend,
      source = source
    )
  }
  if (isTRUE(set_variable_features)) {
    srt <- spatial_set_active_variable_features(srt, assay, vars)
  }
  saved_message <- if (store_results) "results stored" else "this run's results were not stored"
  feature_message <- if (!set_variable_features) "VariableFeatures unchanged" else if (length(vars)) "VariableFeatures updated" else "VariableFeatures cleared"
  log_message(
    "Completed spatial gradient screening: {.val {length(vars)}} features; {saved_message}; {feature_message}",
    message_type = "success",
    verbose = verbose
  )
  srt
}


#' @title Plot spatial gradient screening results
#'
#' @description
#' Visualize normalized results produced by `RunSpatialGradientFeatures()`
#'
#' @md
#' @inheritParams SpatialSpotPlot
#' @param result_name Stored spatial gradient result name. If `NULL`, the latest
#' stored result is used.
#' @param plot_type Plot type: `"summary"`, `"surface"`, `"line"`, `"model"`, or
#' `"combined"`.
#' @param features Variables to plot. If `NULL`, top variables from the stored
#' result are used.
#' @param nfeatures Number of top variables used when `features = NULL`.
#' @param palette,palcolor Color palette passed to SCOP plotting helpers.
#' @param legend.position Legend position for surface, line, and model plots.
#' @param nrow,ncol,byrow Layout controls for multi-feature plots.
#' @param line_size Size of fitted gradient lines.
#' @param line_alpha Alpha for raw value points.
#' @param line_fit Gradient line source. `"stored"` uses the saved
#' `screening$estimate` values produced by the selected backend. `"lm"` draws a
#' fresh linear fit from `screening$value`, which is useful for showing a simple
#' monotonic trend even when the backend stores a smoothed curve.
#'
#' @return A `ggplot` or `patchwork` object.
#'
#' @examples
#' data(visium_human_pancreas_results_sub)
#' SpatialGradientPlot(
#'   visium_human_pancreas_results_sub,
#'   result_name = "scop_gradient_fixture",
#'   plot_type = "surface",
#'   features = rownames(visium_human_pancreas_results_sub)[1:2],
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y"),
#'   pt.size = 1.2
#' )
#' @export
SpatialGradientPlot <- function(
  srt,
  result_name = NULL,
  plot_type = c("summary", "surface", "line", "model", "combined"),
  features = NULL,
  nfeatures = 4,
  assay = NULL,
  layer = "data",
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  pt.size = NULL,
  pt.alpha = 0.9,
  stroke = 0.1,
  palette = "Spectral",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  line_size = 1,
  line_alpha = 0.35,
  line_fit = c("stored", "lm"),
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  image.scale = c("lowres", "hires")
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  image.scale <- match.arg(image.scale)
  line_fit <- match.arg(line_fit)
  result <- sgf_get_result(srt, result_name = result_name)
  spatial_require_coordinate_contract(result, "RunSpatialGradientFeatures()")
  features <- sgf_plot_features(result, features = features, nfeatures = nfeatures)
  layer <- sgf_plot_layer(srt = srt, result = result, assay = assay, layer = layer, features = features)

  if (identical(plot_type, "surface")) {
    return(sgf_surface_plot(
      srt = srt,
      result = result,
      features = features,
      assay = assay,
      layer = layer,
      image = image,
      image.scale = image.scale,
      overlay_image = overlay_image,
      image.alpha = image.alpha,
      coord.cols = coord.cols,
      flip.y = flip.y,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      stroke = stroke,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    ))
  }
  if (identical(plot_type, "line")) {
    return(sgf_line_plot(
      result = result,
      features = features,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args,
      line_size = line_size,
      line_alpha = line_alpha,
      line_fit = line_fit,
      nrow = nrow,
      ncol = ncol
    ))
  }
  if (identical(plot_type, "model")) {
    return(sgf_model_plot(
      result = result,
      features = features,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
  if (identical(plot_type, "summary")) {
    return(sgf_summary_plot(
      result = result,
      features = features,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  check_r("patchwork", verbose = FALSE)
  surface <- sgf_surface_plot(
    srt = srt,
    result = result,
    features = features,
    assay = assay,
    layer = layer,
    image = image,
    image.scale = image.scale,
    overlay_image = overlay_image,
    image.alpha = image.alpha,
    coord.cols = coord.cols,
    flip.y = flip.y,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    stroke = stroke,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow
  )
  line <- sgf_line_plot(
    result = result,
    features = features,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args,
    line_size = line_size,
    line_alpha = line_alpha,
    line_fit = line_fit,
    nrow = nrow,
    ncol = ncol
  )
  check_r("patchwork", verbose = FALSE)
  wrap_plots <- get_namespace_fun("patchwork", "wrap_plots")
  wrap_plots(surface, line, ncol = 1)
}


sgf_resolve_variables <- function(srt, assay, layer, variables = NULL) {
  if (is.null(variables)) {
    variables <- srt@tools[["SpatialVariableFeatures"]]$summary$top_features
  }
  if (is.null(variables) || length(variables) == 0L) {
    variables <- SeuratObject::VariableFeatures(srt, assay = assay)
  }
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(variables) || length(variables) == 0L) {
    variables <- rownames(expr)
  }
  variables <- unique(as.character(variables))
  variables <- intersect(variables, rownames(expr))
  if (length(variables) == 0L) {
    log_message(
      "No requested {.arg variables} are present in assay {.val {assay}}",
      message_type = "error"
    )
  }
  variables
}


sgf_run_cpp_gradient <- function(
  srt,
  reference,
  assay,
  layer,
  variables,
  image,
  coord.cols,
  coordinate_space,
  start,
  end,
  traj_df,
  annotation_ids,
  annotation.by,
  annotation.groups,
  annotation.variable,
  annotation.threshold,
  n_random,
  seed,
  n_bins,
  min_spots,
  sign_var,
  sign_threshold,
  nfeatures,
  parameters
) {
  if (!is.null(annotation_ids) && length(annotation_ids) > 0L) {
    log_message(
      "{.arg annotation_ids} is not supported by the native {.arg backend = 'cpp'} path; provide annotation metadata through {.arg annotation.by}/{.arg annotation.groups} or {.arg annotation.variable}/{.arg annotation.threshold}",
      message_type = "error"
    )
  }

  coords <- sgf_cpp_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )
  coordinate_source <- attr(coords, "spatial_source") %||% list(
    image = image,
    coord.cols = coord.cols[seq_len(min(2L, length(coord.cols)))],
    image_policy = "strict",
    coordinate_space = coordinate_space
  )
  spots <- intersect(colnames(srt), rownames(coords))
  if (length(spots) == 0L) {
    log_message("No spatial coordinates match spots in {.arg srt}", message_type = "error")
  }
  coords <- coords[spots, c("x", "y"), drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  coords <- coords[keep_coords, , drop = FALSE]
  spots <- rownames(coords)
  if (length(spots) < 3L) {
    log_message("At least three spots with finite coordinates are required", message_type = "error")
  }

  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  expr <- expr[variables, spots, drop = FALSE]
  if (!methods::is(expr, "dgCMatrix")) {
    expr <- methods::as(Matrix::Matrix(expr, sparse = TRUE), "dgCMatrix")
  }

  reference_spots <- sgf_cpp_reference_spots(
    srt = srt,
    spots = spots,
    reference = reference,
    annotation.by = annotation.by,
    annotation.groups = annotation.groups,
    annotation.variable = annotation.variable,
    annotation.threshold = annotation.threshold
  )
  trajectory <- sgf_cpp_trajectory(
    reference = reference,
    start = start,
    end = end,
    traj_df = traj_df
  )

  out <- spatial_gradient_screening_cpp(
    expr = expr,
    coords = as.matrix(coords),
    reference_spots = reference_spots,
    trajectory = trajectory,
    variables = rownames(expr),
    mode = reference,
    n_bins = n_bins,
    n_random = as.integer(n_random),
    seed = as.integer(seed),
    min_spots = min_spots
  )
  significance <- sgf_standardize_significance(out$significance)
  if ("p_value" %in% colnames(significance) && nrow(significance) > 0L) {
    significance$fdr <- if (all(is.na(significance$p_value))) {
      NA_real_
    } else {
      stats::p.adjust(significance$p_value, method = "BH")
    }
  }
  model_fits <- sgf_standardize_model_fits(out$model_fits)
  top_variables <- sgf_top_variables(
    screening_out = list(),
    significance = significance,
    model_fits = model_fits,
    nfeatures = nfeatures,
    sign_var = sign_var,
    sign_threshold = sign_threshold
  )
  list(
    screening = sgf_standardize_screening_df(out$screening, reference = reference),
    significance = significance,
    model_fits = model_fits,
    top_variables = top_variables,
    parameters = parameters_summary_df(parameters),
    source = coordinate_source
  )
}

sgf_cpp_coords <- function(srt, image, coord.cols, coordinate_space = "raw") {
  if (length(coord.cols) < 2L) {
    log_message("{.arg coord.cols} must contain at least two coordinate columns", message_type = "error")
  }
  coord.cols <- coord.cols[seq_len(2L)]
  resolved <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space,
    image_policy = "strict"
  )
  out <- resolved$data
  attr(out, "spatial_source") <- resolved$source
  out
}

sgf_cpp_reference_spots <- function(
  srt,
  spots,
  reference,
  annotation.by,
  annotation.groups,
  annotation.variable,
  annotation.threshold
) {
  if (!identical(reference, "annotation")) {
    return(rep(FALSE, length(spots)))
  }
  meta <- srt@meta.data[spots, , drop = FALSE]
  if (!is.null(annotation.by) && !is.null(annotation.groups)) {
    if (!annotation.by %in% colnames(meta)) {
      log_message("{.arg annotation.by} {.val {annotation.by}} is not present in metadata", message_type = "error")
    }
    out <- as.character(meta[[annotation.by]]) %in% as.character(annotation.groups)
  } else if (!is.null(annotation.variable) && !is.null(annotation.threshold)) {
    if (!annotation.variable %in% colnames(meta)) {
      log_message(
        "{.arg annotation.variable} {.val {annotation.variable}} is not present in metadata",
        message_type = "error"
      )
    }
    values <- meta[[annotation.variable]]
    if (!is.numeric(values)) {
      values <- suppressWarnings(as.numeric(as.character(values)))
    }
    out <- sgf_eval_annotation_threshold(values, annotation.threshold)
  } else {
    log_message(
      paste(
        "For {.arg reference = 'annotation'} with {.arg backend = 'cpp'}, provide",
        "{.arg annotation.by} with {.arg annotation.groups}, or",
        "{.arg annotation.variable} with {.arg annotation.threshold}"
      ),
      message_type = "error"
    )
  }
  out[is.na(out)] <- FALSE
  if (!any(out)) {
    log_message("The annotation reference contains no matching spots", message_type = "error")
  }
  as.logical(out)
}

sgf_eval_annotation_threshold <- function(values, threshold) {
  threshold <- sgf_format_annotation_threshold(threshold)
  finite <- is.finite(values)
  out <- rep(FALSE, length(values))
  if (grepl("^kmeans_(high|low)$", threshold)) {
    if (sum(finite) < 2L || length(unique(values[finite])) < 2L) {
      log_message("k-means annotation threshold requires at least two finite values", message_type = "error")
    }
    set.seed(123)
    km <- stats::kmeans(values[finite], centers = 2)
    means <- tapply(values[finite], km$cluster, mean)
    keep <- if (identical(threshold, "kmeans_high")) {
      names(which.max(means))
    } else {
      names(which.min(means))
    }
    out[finite] <- as.character(km$cluster) == keep
    return(out)
  }
  op <- sub("^\\s*(>=|<=|>|<).*$", "\\1", threshold)
  cutoff <- suppressWarnings(as.numeric(sub("^\\s*(>=|<=|>|<)\\s*", "", threshold)))
  if (!is.finite(cutoff)) {
    log_message("{.arg annotation.threshold} contains a non-finite cutoff", message_type = "error")
  }
  out[finite] <- switch(op,
    ">" = values[finite] > cutoff,
    ">=" = values[finite] >= cutoff,
    "<" = values[finite] < cutoff,
    "<=" = values[finite] <= cutoff,
    log_message("{.arg annotation.threshold} must start with >, >=, <, or <=", message_type = "error")
  )
  out
}

sgf_cpp_trajectory <- function(reference, start, end, traj_df) {
  if (!identical(reference, "trajectory")) {
    return(matrix(numeric(0), ncol = 2L))
  }
  if (!is.null(traj_df)) {
    traj_df <- as.data.frame(traj_df, check.names = FALSE)
    x_col <- intersect(c("x", "X", "col", "imagecol", "pxl_col_in_fullres"), colnames(traj_df))[1L]
    y_col <- intersect(c("y", "Y", "row", "imagerow", "pxl_row_in_fullres"), colnames(traj_df))[1L]
    if (is.na(x_col) || is.na(y_col)) {
      numeric_cols <- names(traj_df)[vapply(traj_df, is.numeric, logical(1))]
      if (length(numeric_cols) < 2L) {
        log_message("{.arg traj_df} must contain x/y columns or at least two numeric columns", message_type = "error")
      }
      x_col <- numeric_cols[1L]
      y_col <- numeric_cols[2L]
    }
    out <- as.matrix(traj_df[, c(x_col, y_col), drop = FALSE])
  } else if (!is.null(start) && !is.null(end)) {
    out <- rbind(as.numeric(start), as.numeric(end))
  } else {
    log_message(
      "For {.arg reference = 'trajectory'} with {.arg backend = 'cpp'}, provide {.arg start}/{.arg end} or {.arg traj_df}",
      message_type = "error"
    )
  }
  if (!is.numeric(out) || ncol(out) != 2L || nrow(out) < 2L || any(!is.finite(out))) {
    log_message("{.arg start}/{.arg end} or {.arg traj_df} must define at least two finite x/y coordinates", message_type = "error")
  }
  out
}


sgf_format_annotation_threshold <- function(threshold) {
  if (length(threshold) != 1L || is.na(threshold)) {
    log_message(
      "{.arg annotation.threshold} must be a single non-missing value",
      message_type = "error"
    )
  }
  if (is.numeric(threshold)) {
    return(paste0(">", format(threshold, scientific = FALSE, trim = TRUE)))
  }
  threshold <- trimws(as.character(threshold))
  if (!nzchar(threshold)) {
    log_message("{.arg annotation.threshold} must not be empty", message_type = "error")
  }
  if (grepl("^kmeans_(high|low)$", threshold)) {
    return(threshold)
  }
  if (grepl("^(>=|<=|>|<)\\s*[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?$", threshold)) {
    return(gsub("\\s+", "", threshold))
  }
  if (grepl("^[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?$", threshold)) {
    return(paste0(">", threshold))
  }
  log_message(
    paste(
      "{.arg annotation.threshold} must be numeric, a comparison string such as {.val >0},",
      "{.val <=1}, or one of {.val kmeans_high}/{.val kmeans_low}"
    ),
    message_type = "error"
  )
}


sgf_standardize_significance <- function(df) {
  df <- sgf_as_df(df)
  df <- rename_first_column(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- rename_first_column(df, "p_value", c("p_value", "p.val", "pval", "p.value"))
  df <- rename_first_column(df, "fdr", c("fdr", "q_value", "q.value", "padj", "adjustedPval"))
  if (!"variable" %in% colnames(df)) {
    df$variable <- character(nrow(df))
  }
  reorder_first_columns(df, c("variable", "tot_var", "p_value", "fdr"))
}

sgf_standardize_model_fits <- function(df) {
  df <- sgf_as_df(df)
  df <- rename_first_column(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- rename_first_column(df, "model", c("model", "models"))
  if (!"variable" %in% colnames(df)) {
    df$variable <- character(nrow(df))
  }
  if (!"model" %in% colnames(df)) {
    df$model <- character(nrow(df))
  }
  reorder_first_columns(df, c("variable", "model", "mae", "rmse", "r2", "R2"))
}

sgf_standardize_screening_df <- function(df, reference) {
  df <- sgf_as_df(df)
  df <- rename_first_column(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- rename_first_column(df, "distance", c("distance", "dist", "x", "x_orig", "x_original", "bin"))
  df <- rename_first_column(df, "value", c("value", "values", "expr", "expression", "y"))
  df <- rename_first_column(df, "estimate", c("estimate", "estimates", "fit", "fitted", "smooth", "smoothed"))
  for (col in c("variable", "distance", "value", "estimate")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  df$reference <- if ("reference" %in% colnames(df)) df$reference else reference
  df$mode <- if ("mode" %in% colnames(df)) df$mode else reference
  reorder_first_columns(df, c("variable", "distance", "value", "estimate", "reference", "mode"))
}

sgf_top_variables <- function(
  screening_out,
  significance,
  model_fits,
  nfeatures,
  sign_var,
  sign_threshold
) {
  vars <- NULL
  vars <- if (is.null(vars)) character(0) else unique(as.character(vars))
  if (length(vars) == 0L && nrow(significance) > 0L) {
    sign_use <- significance
    if (sign_var %in% colnames(sign_use) && is.finite(sign_threshold)) {
      sign_use <- sign_use[is.na(sign_use[[sign_var]]) | sign_use[[sign_var]] <= sign_threshold, , drop = FALSE]
    }
    order_cols <- intersect(c("fdr", "p_value"), colnames(sign_use))
    if (length(order_cols) > 0L) {
      ord <- do.call(order, c(sign_use[order_cols], list(na.last = TRUE)))
      sign_use <- sign_use[ord, , drop = FALSE]
    }
    vars <- unique(sign_use$variable)
  }
  vars <- utils::head(vars[nzchar(vars)], nfeatures)
  top <- data.frame(variable = vars, rank = seq_along(vars), stringsAsFactors = FALSE)
  top <- merge(top, significance, by = "variable", all.x = TRUE, sort = FALSE)
  best <- sgf_best_model_fits(model_fits)
  top <- merge(top, best, by = "variable", all.x = TRUE, sort = FALSE)
  top <- top[order(top$rank), , drop = FALSE]
  rownames(top) <- NULL
  reorder_first_columns(top, c("variable", "rank", "fdr", "p_value", "best_model", "mae", "rmse"))
}


sgf_best_model_fits <- function(model_fits) {
  if (!is.data.frame(model_fits) || nrow(model_fits) == 0L) {
    return(data.frame(variable = character(), best_model = character(), mae = numeric(), rmse = numeric()))
  }
  df <- model_fits
  if (!"mae" %in% colnames(df)) {
    df$mae <- NA_real_
  }
  if (!"rmse" %in% colnames(df)) {
    df$rmse <- NA_real_
  }
  split_df <- split(df, df$variable)
  best <- lapply(split_df, function(x) {
    ord <- order(x$mae, x$rmse, na.last = TRUE)
    x <- x[ord[1L], , drop = FALSE]
    data.frame(
      variable = x$variable[1L],
      best_model = x$model[1L],
      mae = x$mae[1L],
      rmse = x$rmse[1L],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, best)
}

sgf_validate_result <- function(result) {
  expected <- c("screening", "significance", "model_fits", "top_variables", "parameters")
  missing <- setdiff(expected, names(result))
  if (length(missing) > 0L) {
    log_message(
      "Spatial gradient result is missing required table(s): {.val {paste(missing, collapse = ', ')}}",
      message_type = "error"
    )
  }
  for (nm in expected) {
    if (!is.data.frame(result[[nm]])) {
      log_message(
        "Spatial gradient result item {.val {nm}} must be a data.frame",
        message_type = "error"
      )
    }
  }
  vars <- result$top_variables$variable
  if (!is.character(vars) || anyNA(vars) || any(!nzchar(vars))) {
    stop("Spatial gradient top_variables must contain a character variable column without missing or empty names", call. = FALSE)
  }
  if (!is.null(result$source) && !is.list(result$source)) {
    stop("Spatial gradient result source must be a list", call. = FALSE)
  }
  unique(vars)
}

sgf_store_result <- function(srt, result_name, result, assay, backend = "cpp", source = NULL) {
  vars <- sgf_validate_result(result)
  expected <- c("screening", "significance", "model_fits", "top_variables", "parameters")
  if (is.null(srt@tools[["SpatialGradientFeatures"]])) {
    srt@tools[["SpatialGradientFeatures"]] <- list()
  }
  source <- source %||% result$source %||% list()
  if (!is.list(source)) {
    log_message(
      "Spatial gradient result source must be a list",
      message_type = "error"
    )
  }
  stored_result <- c(result[expected], list(source = source))
  attr(stored_result, "coordinate_contract_version") <-
    .spatial_coordinate_contract_version
  srt@tools[["SpatialGradientFeatures"]][[result_name]] <- stored_result
  srt@tools[["SpatialGradientFeatures"]]$parameters <- list(
    result_name = result_name,
    assay = assay,
    backend = backend,
    coordinate_contract_version = .spatial_coordinate_contract_version
  )
  srt@tools[["SpatialGradientFeatures"]]$summary <- list(
    active_result = result_name,
    top_features = vars,
    n_results = length(setdiff(
      names(srt@tools[["SpatialGradientFeatures"]]),
      c("parameters", "summary")
    ))
  )
  srt
}

sgf_get_result <- function(srt, result_name = NULL) {
  all_results <- srt@tools[["SpatialGradientFeatures"]]
  if (is.null(all_results) || length(all_results) == 0L) {
    log_message(
      "No spatial gradient results are stored in {.code srt@tools[['SpatialGradientFeatures']]}",
      message_type = "error"
    )
  }
  metadata_names <- c("parameters", "summary")
  result_names <- setdiff(names(all_results), metadata_names)
  if (is.null(result_name)) {
    result_name <- all_results$summary$active_result %||% result_names[length(result_names)]
  }
  if (!result_name %in% result_names) {
    log_message(
      "{.arg result_name} {.val {result_name}} is not present in stored spatial gradient results",
      message_type = "error"
    )
  }
  all_results[[result_name]]
}

sgf_plot_features <- function(result, features = NULL, nfeatures = 4) {
  if (is.null(features)) {
    if (!is.data.frame(result$top_variables) || nrow(result$top_variables) == 0L) {
      log_message("No top spatial gradient variables are available", message_type = "error")
    }
    features <- utils::head(result$top_variables$variable, nfeatures)
  }
  features <- unique(as.character(features))
  features <- features[!is.na(features) & nzchar(features)]
  if (length(features) == 0L) {
    log_message("No {.arg features} are available for plotting", message_type = "error")
  }
  features
}

sgf_plot_layer <- function(srt, result, assay, layer, features) {
  stored_layer <- sgf_parameter_value(result, "layer")
  if (is.null(layer) || identical(layer, "") || identical(layer, NA_character_)) {
    return(stored_layer %||% "data")
  }
  if (!identical(layer, "data") || is.null(stored_layer) || identical(stored_layer, layer)) {
    return(layer)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  expr <- tryCatch(
    suppressWarnings(GetAssayData5(srt, assay = assay, layer = layer)),
    error = function(e) NULL
  )
  has_features <- !is.null(expr) && all(features %in% rownames(expr))
  if (isTRUE(has_features)) {
    return(layer)
  }
  stored_expr <- tryCatch(
    suppressWarnings(GetAssayData5(srt, assay = assay, layer = stored_layer)),
    error = function(e) NULL
  )
  if (!is.null(stored_expr) && all(features %in% rownames(stored_expr))) {
    return(stored_layer)
  }
  layer
}

sgf_parameter_value <- function(result, key) {
  params <- result$parameters
  if (!is.data.frame(params) || !"key" %in% colnames(params) || !"value" %in% colnames(params)) {
    return(NULL)
  }
  idx <- which(params$key %in% key)
  if (length(idx) == 0L) {
    return(NULL)
  }
  val <- params$value[[idx[[1L]]]]
  if (is.na(val) || !nzchar(val)) {
    return(NULL)
  }
  as.character(val)
}

sgf_surface_plot <- function(
  srt,
  result,
  features,
  assay,
  layer,
  image,
  image.scale,
  overlay_image,
  image.alpha,
  coord.cols,
  flip.y,
  pt.size,
  pt.alpha,
  stroke,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args,
  nrow,
  ncol,
  byrow
) {
  if (is.null(theme_use)) {
    theme_use <- ggplot2::theme_minimal
  }
  SpatialSpotPlot(
    srt = srt,
    features = features,
    assay = assay,
    layer = layer,
    image = image,
    image.scale = image.scale,
    overlay_image = overlay_image,
    image.alpha = image.alpha,
    coord.cols = coord.cols,
    flip.y = flip.y,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    stroke = stroke,
    palette = palette,
    palcolor = palcolor,
    legend.position = legend.position,
    theme_use = theme_use,
    theme_args = theme_args,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow
  )
}

sgf_line_plot <- function(
  result,
  features,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args,
  line_size,
  line_alpha,
  line_fit,
  nrow,
  ncol
) {
  df <- result$screening
  if (!is.data.frame(df) || nrow(df) == 0L) {
    log_message(
      "Stored spatial gradient result does not contain screening data for line plotting",
      message_type = "error"
    )
  }
  df <- df[df$variable %in% features, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message("No screening records are available for requested {.arg features}", message_type = "error")
  }
  df$distance <- suppressWarnings(as.numeric(df$distance))
  df$value <- suppressWarnings(as.numeric(df$value))
  df$estimate <- suppressWarnings(as.numeric(df$estimate))
  y_col <- if (any(is.finite(df$estimate))) "estimate" else "value"
  cols <- sgf_feature_colors(features, palette = palette, palcolor = palcolor)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["distance"]], color = .data[["variable"]])) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data[["value"]]),
      alpha = line_alpha,
      size = 0.8,
      na.rm = TRUE
    )
  if (identical(line_fit, "lm")) {
    p <- p + ggplot2::geom_smooth(
      ggplot2::aes(y = .data[["value"]]),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = line_size,
      na.rm = TRUE
    )
  } else {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[[y_col]]),
        linewidth = line_size,
        na.rm = TRUE
      )
  }
  p <- p +
    ggplot2::facet_wrap(~variable, scales = "free_y", nrow = nrow, ncol = ncol) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::labs(x = "Gradient distance", y = "Expression") +
    apply_plot_theme(theme_use = theme_use, theme_args = theme_args, fallback = ggplot2::theme_minimal) +
    ggplot2::theme(legend.position = legend.position)
  p
}

sgf_summary_plot <- function(
  result,
  features,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args
) {
  df <- result$top_variables
  if (!is.data.frame(df) || nrow(df) == 0L) {
    log_message("No top spatial gradient variables are available", message_type = "error")
  }
  df <- df[df$variable %in% features, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message("No top variable records are available for requested {.arg features}", message_type = "error")
  }
  if (!"rmse" %in% colnames(df)) {
    df$rmse <- NA_real_
  }
  df <- df[order(df$rank, na.last = TRUE), , drop = FALSE]
  df$variable <- factor(df$variable, levels = rev(as.character(df$variable)))
  use_rmse <- any(is.finite(df$rmse))
  df$.metric <- if (isTRUE(use_rmse)) df$rmse else df$rank
  cols <- sgf_feature_colors(as.character(df$variable), palette = palette, palcolor = palcolor)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[".metric"]], y = .data[["variable"]], color = .data[["variable"]])) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[[".metric"]], yend = .data[["variable"]]),
      linewidth = 0.35,
      alpha = 0.45,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = if (isTRUE(use_rmse)) "RMSE" else "Rank", y = NULL) +
    apply_plot_theme(theme_use = theme_use, theme_args = theme_args, fallback = ggplot2::theme_minimal) +
    ggplot2::theme(legend.position = legend.position)
}

sgf_model_plot <- function(
  result,
  features,
  palette,
  palcolor,
  legend.position,
  theme_use,
  theme_args
) {
  df <- result$model_fits
  if (!is.data.frame(df) || nrow(df) == 0L) {
    log_message("No spatial gradient model fit table is available", message_type = "error")
  }
  df <- df[df$variable %in% features, , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message("No model fit records are available for requested {.arg features}", message_type = "error")
  }
  df$variable <- factor(df$variable, levels = rev(features))
  if (!"rmse" %in% colnames(df)) {
    df$rmse <- NA_real_
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["model"]], y = .data[["variable"]], fill = .data[["rmse"]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.finite(.data[["rmse"]]), formatC(.data[["rmse"]], format = "f", digits = 2), "")),
      size = 3,
      color = "grey15"
    ) +
    ggplot2::scale_fill_gradientn(colors = rev({
      .inline0 <- palette
      .inline1 <- palcolor
      palette_colors(palette = .inline0, palcolor = .inline1, n = 9)
    }), na.value = "grey85") +
    ggplot2::labs(x = "Model", y = NULL, fill = "RMSE") +
    apply_plot_theme(theme_use = theme_use, theme_args = theme_args, fallback = ggplot2::theme_minimal) +
    ggplot2::theme(
      legend.position = legend.position,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

sgf_feature_colors <- function(features, palette = "Spectral", palcolor = NULL) {
  features <- unique(as.character(features))
  cols <- palette_colors(features, palette = palette, palcolor = palcolor)
  stats::setNames(cols[seq_along(features)], features)
}


sgf_as_df <- function(df) {
  if (!is.data.frame(df)) {
    df <- data.frame()
  }
  as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
}
