#' @title Run spatial gradient feature screening
#'
#' @description
#' Run spatial trajectory or annotation gradient screening for Seurat objects.
#' The native `"cpp"` backend avoids SPATA2 object construction for fast
#' distance-based screening, while the `"r"` backend keeps full upstream
#' SPATA2 SAS/STS behavior. Results are normalized into plain data.frames and
#' stored in `srt@tools[["SpatialGradientFeatures"]]`; the SPATA2 object itself
#' is never stored.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @inheritParams SpatialSpotPlot
#' @param reference Spatial reference type: `"trajectory"` for STS or
#' `"annotation"` for SAS.
#' @param backend Computation backend. `"cpp"` uses SCOP's native fast spatial
#' gradient implementation and avoids SPATA2 object construction. `"r"`
#' uses SPATA2 directly for full upstream SAS/STS behavior.
#' @param result_name Name used to store this result. If `NULL`, a name is
#' generated from `reference`.
#' @param spata_object Optional pre-built SPATA2 object. If `NULL`, `srt` is
#' converted with `SPATA2::asSPATA2()`.
#' @param variables Numeric variables or genes passed to SPATA2. If `NULL`,
#' `srt@misc[["SpatialVariableFeatures"]]` is used first, then variable
#' features, then all assay features.
#' @param sample_name,platform,img_scale_fct,assay_modality Arguments forwarded
#' to `SPATA2::asSPATA2()` when `spata_object` is `NULL`.
#' @param coord.cols Metadata coordinate columns used by the native `"cpp"`
#' backend when no image coordinates are available.
#' @param trajectory_id,start,end,traj_df,width Trajectory setup passed to
#' `SPATA2::addSpatialTrajectory()` and `SPATA2::spatialTrajectoryScreening()`.
#' @param annotation_ids Existing SPATA2 spatial annotation ids. If `NULL`,
#' annotations are created from `annotation.by` and `annotation.groups`, or from
#' `annotation.variable` and `annotation.threshold`.
#' @param annotation.by,annotation.groups Metadata grouping used to create
#' SPATA2 group annotations.
#' @param annotation.variable,annotation.threshold Numeric variable and
#' threshold used to create SPATA2 numeric annotations. Numeric thresholds are
#' interpreted as `">{threshold}"`.
#' @param annotation_id Base id used when creating annotations.
#' @param core,distance,angle_span SAS parameters forwarded to SPATA2.
#' @param resolution,unit,sign_var,sign_threshold,model_add,model_subset,model_remove,n_random,seed,control
#' SPATA2 screening parameters.
#' @param n_bins Number of distance bins used for the native `"cpp"` backend
#' screening curve.
#' @param min_spots Minimum number of non-zero spots required for a variable in
#' the native `"cpp"` backend.
#' @param nfeatures Number of top gradient variables retained in
#' `top_variables` and optionally set as Seurat variable features.
#' @param set_variable_features Whether to set top gradient variables as Seurat
#' variable features.
#' @param store_results Whether to store the normalized result in `srt@tools`.
#' @param ... Additional arguments forwarded to the SPATA2 screening function.
#'
#' @return A `Seurat` object with spatial gradient screening results stored in
#' `srt@tools[["SpatialGradientFeatures"]]`.
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
#' SpatialGradientPlot(spatial, plot_type = "summary", nfeatures = 4)
#' SpatialGradientPlot(spatial, plot_type = "line", nfeatures = 2)
#' SpatialGradientPlot(spatial, plot_type = "model", nfeatures = 2)
#' SpatialGradientPlot(
#'   spatial,
#'   plot_type = "surface",
#'   nfeatures = 2,
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
RunSpatialGradientFeatures <- function(
  srt,
  reference = c("trajectory", "annotation"),
  backend = c("cpp", "r"),
  result_name = NULL,
  spata_object = NULL,
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
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  reference <- match.arg(reference)
  backend <- match.arg(backend)
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
        trajectory_id = trajectory_id,
        annotation.by = annotation.by,
        annotation.groups = paste(annotation.groups %||% character(0), collapse = ","),
        annotation.variable = annotation.variable,
        annotation.threshold = annotation.threshold,
        distance = distance,
        n_random = n_random,
        seed = seed,
        n_bins = n_bins,
        min_spots = min_spots
      )
    )
  } else {
    sgf_require_spata2()
    spata_object <- spata_object %||% sgf_as_spata2(
      srt = srt,
      sample_name = sample_name %||% "scop_sample",
      platform = platform,
      assay = assay,
      image = image,
      img_scale_fct = img_scale_fct,
      assay_modality = assay_modality,
      verbose = verbose
    )

    if (identical(reference, "trajectory")) {
      spata_object <- sgf_prepare_trajectory(
        object = spata_object,
        trajectory_id = trajectory_id,
        start = start,
        end = end,
        traj_df = traj_df,
        width = width,
        verbose = verbose
      )
      screening_out <- sgf_run_trajectory_screening(
        object = spata_object,
        trajectory_id = trajectory_id,
        variables = variables,
        resolution = resolution,
        width = width,
        unit = unit,
        sign_var = sign_var,
        sign_threshold = sign_threshold,
        model_add = model_add,
        model_subset = model_subset,
        model_remove = model_remove,
        n_random = n_random,
        seed = seed,
        control = control,
        verbose = verbose,
        ...
      )
      annotation_ids_use <- character(0)
    } else {
      prep <- sgf_prepare_annotations(
        object = spata_object,
        annotation_ids = annotation_ids,
        annotation.by = annotation.by,
        annotation.groups = annotation.groups,
        annotation.variable = annotation.variable,
        annotation.threshold = annotation.threshold,
        annotation_id = annotation_id,
        verbose = verbose
      )
      spata_object <- prep$object
      annotation_ids_use <- prep$annotation_ids
      screening_out <- sgf_run_annotation_screening(
        object = spata_object,
        annotation_ids = annotation_ids_use,
        variables = variables,
        core = core,
        distance = distance,
        resolution = resolution,
        angle_span = angle_span,
        unit = unit,
        sign_var = sign_var,
        sign_threshold = sign_threshold,
        model_add = model_add,
        model_subset = model_subset,
        model_remove = model_remove,
        n_random = n_random,
        seed = seed,
        control = control,
        verbose = verbose,
        ...
      )
    }

    result <- sgf_normalize_screening_result(
      screening_out = screening_out,
      spata_object = spata_object,
      reference = reference,
      variables = variables,
      nfeatures = nfeatures,
      trajectory_id = trajectory_id,
      annotation_ids = annotation_ids_use,
      distance = distance,
      width = width,
      unit = unit,
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      parameters = list(
        result_name = result_name,
        reference = reference,
        backend = backend,
        assay = assay,
        layer = layer,
        sample_name = sample_name %||% "scop_sample",
        platform = platform,
        image = image,
        img_scale_fct = img_scale_fct,
        assay_modality = assay_modality,
        trajectory_id = trajectory_id,
        annotation_ids = paste(annotation_ids_use, collapse = ","),
        annotation.by = annotation.by,
        annotation.groups = paste(annotation.groups %||% character(0), collapse = ","),
        annotation.variable = annotation.variable,
        annotation.threshold = annotation.threshold,
        core = core,
        distance = distance,
        resolution = resolution,
        unit = unit,
        sign_var = sign_var,
        sign_threshold = sign_threshold,
        n_random = n_random,
        seed = seed
      )
    )
  }

  if (isTRUE(store_results)) {
    srt <- sgf_store_result(
      srt = srt,
      result_name = result_name,
      result = result,
      assay = assay,
      set_variable_features = set_variable_features
    )
  }
  log_message(
    "Stored {.val {nrow(result$top_variables)}} spatial gradient features",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Plot spatial gradient screening results
#'
#' @description
#' Visualize normalized results produced by `RunSpatialGradientFeatures()`
#' without requiring the original SPATA2 object.
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
#' counts <- matrix(
#'   c(4, 1, 0, 2, 1, 3, 2, 0),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("REG1A", "COL1A1")
#' colnames(counts) <- paste0("spot", 1:4)
#' srt <- Seurat::CreateSeuratObject(counts)
#' srt <- Seurat::NormalizeData(srt, verbose = FALSE)
#' srt$col <- c(0, 1, 0, 1)
#' srt$row <- c(0, 0, 1, 1)
#'
#' gradient_result <- list(
#'   screening = data.frame(
#'     variable = rep(c("REG1A", "COL1A1"), each = 4),
#'     distance = rep(seq(0, 1, length.out = 4), 2),
#'     value = c(0.1, 0.4, 0.8, 1.1, 1.0, 0.7, 0.3, 0.1),
#'     estimate = c(0.15, 0.45, 0.75, 1.05, 0.95, 0.65, 0.35, 0.05)
#'   ),
#'   significance = data.frame(
#'     variable = c("REG1A", "COL1A1"),
#'     p_value = c(0.004, 0.018),
#'     q_value = c(0.008, 0.024)
#'   ),
#'   model_fits = data.frame(
#'     variable = rep(c("REG1A", "COL1A1"), each = 2),
#'     model = rep(c("linear", "spline"), 2),
#'     rmse = c(0.12, 0.08, 0.18, 0.11)
#'   ),
#'   top_variables = data.frame(
#'     variable = c("REG1A", "COL1A1"),
#'     rank = 1:2,
#'     rmse = c(0.08, 0.11)
#'   ),
#'   parameters = data.frame(
#'     key = c("assay", "layer", "reference"),
#'     value = c("RNA", "data", "ductal_axis")
#'   )
#' )
#' srt@tools[["SpatialGradientFeatures"]] <- list(ductal_axis = gradient_result)
#' srt@misc[["SpatialGradientFeaturesResult"]] <- "ductal_axis"
#'
#' SpatialGradientPlot(srt, plot_type = "summary", nfeatures = 2)
#' SpatialGradientPlot(srt, plot_type = "line", nfeatures = 2)
#' SpatialGradientPlot(srt, plot_type = "model", nfeatures = 2)
#' SpatialGradientPlot(
#'   srt,
#'   plot_type = "surface",
#'   nfeatures = 2,
#'   overlay_image = FALSE,
#'   coord.cols = c("col", "row"),
#'   pt.size = 4
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
  byrow = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  line_fit <- match.arg(line_fit)
  result <- sgf_get_result(srt, result_name = result_name)
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

  sgf_require_package("patchwork")
  surface <- sgf_surface_plot(
    srt = srt,
    result = result,
    features = features,
    assay = assay,
    layer = layer,
    image = image,
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
  sgf_require_package("patchwork")
  wrap_plots <- get_namespace_fun("patchwork", "wrap_plots")
  wrap_plots(surface, line, ncol = 1)
}

sgf_require_package <- function(pkg) {
  repo <- if (identical(pkg, "SPATA2")) "theMILOlab/SPATA2" else pkg
  status <- tryCatch(check_r(repo, verbose = FALSE), error = function(e) FALSE)
  if (!isTRUE(unname(unlist(status))[1])) {
    log_message(
      "Please install required package before running this function: {.val {pkg}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

sgf_spata2_pkg <- function() {
  "SPATA2"
}

sgf_require_spata2 <- function() {
  status <- tryCatch(check_r("theMILOlab/SPATA2", verbose = FALSE), error = function(e) FALSE)
  if (!isTRUE(unname(unlist(status))[1])) {
    log_message(
      paste(
        "Please install SPATA2 before running spatial gradient screening.",
        "Official installation uses devtools::install_github('theMILOlab/SPATA2')."
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

sgf_spata_fun <- function(fun, required = TRUE) {
  pkg <- sgf_spata2_pkg()
  sgf_require_spata2()
  out <- tryCatch(get_namespace_fun(pkg, fun), error = function(e) NULL)
  if (is.null(out) && isTRUE(required)) {
    log_message(
      "Installed SPATA2 does not provide required function {.fn {fun}}",
      message_type = "error"
    )
  }
  out
}

sgf_try_spata_call <- function(fun, args, required = FALSE) {
  f <- sgf_spata_fun(fun, required = required)
  if (is.null(f)) {
    return(NULL)
  }
  tryCatch(do.call(f, args), error = function(e) NULL)
}

sgf_resolve_variables <- function(srt, assay, layer, variables = NULL) {
  if (is.null(variables)) {
    variables <- srt@misc[["SpatialVariableFeatures"]]
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

sgf_as_spata2 <- function(
  srt,
  sample_name,
  platform,
  assay,
  image,
  img_scale_fct,
  assay_modality,
  verbose
) {
  args <- sgf_drop_nulls(list(
    object = srt,
    sample_name = sample_name,
    platform = platform,
    assay_name = assay,
    assay_modality = assay_modality,
    img_name = image,
    img_scale_fct = img_scale_fct,
    transfer_meta_data = TRUE,
    transfer_dim_red = FALSE,
    verbose = verbose
  ))
  do.call(sgf_spata_fun("asSPATA2"), args)
}

sgf_run_cpp_gradient <- function(
  srt,
  reference,
  assay,
  layer,
  variables,
  image,
  coord.cols,
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
      "{.arg annotation_ids} requires {.arg backend = 'r'} because SPATA2 annotation ids are not stored in Seurat metadata",
      message_type = "error"
    )
  }

  coords <- sgf_cpp_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols
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
    parameters = sgf_parameters_df(parameters)
  )
}

sgf_cpp_coords <- function(srt, image, coord.cols) {
  if (length(coord.cols) < 2L) {
    log_message("{.arg coord.cols} must contain at least two coordinate columns", message_type = "error")
  }
  coord.cols <- coord.cols[seq_len(2L)]
  if (is.null(image) && all(coord.cols %in% colnames(srt@meta.data))) {
    return(data.frame(
      x = srt@meta.data[[coord.cols[1L]]],
      y = srt@meta.data[[coord.cols[2L]]],
      row.names = rownames(srt@meta.data),
      stringsAsFactors = FALSE
    ))
  }
  spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
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
  out[finite] <- switch(
    op,
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

sgf_prepare_trajectory <- function(object, trajectory_id, start, end, traj_df, width, verbose) {
  if (is.null(start) && is.null(end) && is.null(traj_df)) {
    return(object)
  }
  invisible(verbose)
  args <- sgf_drop_nulls(list(
    object = object,
    id = trajectory_id,
    width = width,
    traj_df = traj_df,
    start = start,
    end = end,
    overwrite = TRUE
  ))
  do.call(sgf_spata_fun("addSpatialTrajectory"), args)
}

sgf_prepare_annotations <- function(
  object,
  annotation_ids,
  annotation.by,
  annotation.groups,
  annotation.variable,
  annotation.threshold,
  annotation_id,
  verbose
) {
  if (!is.null(annotation_ids) && length(annotation_ids) > 0L) {
    return(list(object = object, annotation_ids = as.character(annotation_ids)))
  }
  if (!is.null(annotation.by) && !is.null(annotation.groups)) {
    object <- do.call(sgf_spata_fun("createGroupAnnotations"), sgf_drop_nulls(list(
      object = object,
      grouping = annotation.by,
      group = annotation.groups,
      id = annotation_id,
      tags = annotation_id,
      overwrite = TRUE,
      verbose = verbose
    )))
    return(list(
      object = object,
      annotation_ids = sgf_get_created_annotation_ids(object, annotation_id)
    ))
  }
  if (!is.null(annotation.variable) && !is.null(annotation.threshold)) {
    annotation.threshold <- sgf_format_annotation_threshold(annotation.threshold)
    object <- do.call(sgf_spata_fun("createNumericAnnotations"), sgf_drop_nulls(list(
      object = object,
      variable = annotation.variable,
      threshold = annotation.threshold,
      id = annotation_id,
      tags = annotation_id,
      overwrite = TRUE,
      verbose = verbose
    )))
    return(list(
      object = object,
      annotation_ids = sgf_get_created_annotation_ids(object, annotation_id)
    ))
  }
  log_message(
    paste(
      "For {.arg reference = 'annotation'}, provide {.arg annotation_ids},",
      "or {.arg annotation.by} with {.arg annotation.groups},",
      "or {.arg annotation.variable} with {.arg annotation.threshold}"
    ),
    message_type = "error"
  )
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

sgf_get_created_annotation_ids <- function(object, annotation_id) {
  ids <- sgf_try_spata_call(
    "getSpatAnnIds",
    list(object = object, tags = annotation_id, test = "any"),
    required = FALSE
  )
  if (is.null(ids) || length(ids) == 0L) {
    ids <- annotation_id
  }
  as.character(ids)
}

sgf_run_trajectory_screening <- function(
  object,
  trajectory_id,
  variables,
  resolution,
  width,
  unit,
  sign_var,
  sign_threshold,
  model_add,
  model_subset,
  model_remove,
  n_random,
  seed,
  control,
  verbose,
  ...
) {
  args <- sgf_drop_nulls(c(
    list(
      object = object,
      id = trajectory_id,
      variables = variables,
      resolution = resolution,
      width = width,
      unit = unit,
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      model_add = model_add,
      model_subset = model_subset,
      model_remove = model_remove,
      n_random = n_random,
      seed = seed,
      control = control,
      verbose = verbose
    ),
    list(...)
  ))
  do.call(sgf_spata_fun("spatialTrajectoryScreening"), args)
}

sgf_run_annotation_screening <- function(
  object,
  annotation_ids,
  variables,
  core,
  distance,
  resolution,
  angle_span,
  unit,
  sign_var,
  sign_threshold,
  model_add,
  model_subset,
  model_remove,
  n_random,
  seed,
  control,
  verbose,
  ...
) {
  args <- sgf_drop_nulls(c(
    list(
      object = object,
      ids = annotation_ids,
      variables = variables,
      core = core,
      distance = distance,
      resolution = resolution,
      angle_span = angle_span,
      unit = unit,
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      model_add = model_add,
      model_subset = model_subset,
      model_remove = model_remove,
      n_random = n_random,
      seed = seed,
      control = control,
      verbose = verbose
    ),
    list(...)
  ))
  do.call(sgf_spata_fun("spatialAnnotationScreening"), args)
}

sgf_normalize_screening_result <- function(
  screening_out,
  spata_object,
  reference,
  variables,
  nfeatures,
  trajectory_id,
  annotation_ids,
  distance,
  width,
  unit,
  sign_var,
  sign_threshold,
  parameters
) {
  significance <- sgf_extract_significance(screening_out)
  model_fits <- sgf_extract_model_fits(screening_out)
  top_variables <- sgf_top_variables(
    screening_out = screening_out,
    significance = significance,
    model_fits = model_fits,
    nfeatures = nfeatures,
    sign_var = sign_var,
    sign_threshold = sign_threshold
  )
  screening_variables <- if (nrow(top_variables) > 0L) top_variables$variable else variables
  screening <- sgf_extract_screening_df(
    spata_object = spata_object,
    screening_out = screening_out,
    reference = reference,
    variables = screening_variables,
    trajectory_id = trajectory_id,
    annotation_ids = annotation_ids,
    distance = distance,
    width = width,
    unit = unit
  )
  list(
    screening = screening,
    significance = significance,
    model_fits = model_fits,
    top_variables = top_variables,
    parameters = sgf_parameters_df(parameters)
  )
}

sgf_extract_significance <- function(screening_out) {
  df <- sgf_slot_result(screening_out, "significance")
  if (!is.data.frame(df) || nrow(df) == 0L) {
    df <- sgf_try_spata_call("getSgsResultsDf", list(screening_out), required = FALSE)
  }
  sgf_standardize_significance(df)
}

sgf_extract_model_fits <- function(screening_out) {
  df <- sgf_slot_result(screening_out, "model_fits")
  sgf_standardize_model_fits(df)
}

sgf_slot_result <- function(screening_out, name) {
  out <- tryCatch({
    res <- methods::slot(screening_out, "results")
    res[[name]]
  }, error = function(e) NULL)
  if (!is.data.frame(out)) {
    out <- data.frame()
  }
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

sgf_extract_screening_df <- function(
  spata_object,
  screening_out,
  reference,
  variables,
  trajectory_id,
  annotation_ids,
  distance,
  width,
  unit
) {
  variables <- unique(as.character(variables))
  if (identical(reference, "trajectory")) {
    raw <- sgf_try_spata_candidates(
      "getStsDf",
      list(
        sgf_drop_nulls(list(object = spata_object, variables = variables, id = trajectory_id, width = width, unit = unit)),
        sgf_drop_nulls(list(spata_object, variables = variables, id = trajectory_id)),
        list(screening_out)
      )
    )
  } else {
    raw <- sgf_try_spata_candidates(
      "getSasDf",
      list(
        sgf_drop_nulls(list(object = spata_object, variables = variables, ids = annotation_ids, distance = distance, unit = unit)),
        sgf_drop_nulls(list(spata_object, variables = variables, ids = annotation_ids)),
        list(screening_out)
      )
    )
  }
  sgf_standardize_screening_df(raw, reference = reference)
}

sgf_try_spata_candidates <- function(fun, candidates) {
  f <- sgf_spata_fun(fun, required = FALSE)
  if (is.null(f)) {
    return(data.frame())
  }
  for (args in candidates) {
    out <- tryCatch(do.call(f, args), error = function(e) NULL)
    if (is.data.frame(out)) {
      return(as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE))
    }
  }
  data.frame()
}

sgf_standardize_significance <- function(df) {
  df <- sgf_as_df(df)
  df <- sgf_rename_first(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- sgf_rename_first(df, "p_value", c("p_value", "p.val", "pval", "p.value"))
  df <- sgf_rename_first(df, "fdr", c("fdr", "q_value", "q.value", "padj", "adjustedPval"))
  if (!"variable" %in% colnames(df)) {
    df$variable <- character(nrow(df))
  }
  sgf_reorder_cols(df, c("variable", "tot_var", "p_value", "fdr"))
}

sgf_standardize_model_fits <- function(df) {
  df <- sgf_as_df(df)
  df <- sgf_rename_first(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- sgf_rename_first(df, "model", c("model", "models"))
  if (!"variable" %in% colnames(df)) {
    df$variable <- character(nrow(df))
  }
  if (!"model" %in% colnames(df)) {
    df$model <- character(nrow(df))
  }
  sgf_reorder_cols(df, c("variable", "model", "mae", "rmse", "r2", "R2"))
}

sgf_standardize_screening_df <- function(df, reference) {
  df <- sgf_as_df(df)
  df <- sgf_rename_first(df, "variable", c("variable", "variables", "gene", "genes", "feature", "features"))
  df <- sgf_rename_first(df, "distance", c("distance", "dist", "x", "x_orig", "x_original", "bin"))
  df <- sgf_rename_first(df, "value", c("value", "values", "expr", "expression", "y"))
  df <- sgf_rename_first(df, "estimate", c("estimate", "estimates", "fit", "fitted", "smooth", "smoothed"))
  for (col in c("variable", "distance", "value", "estimate")) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  df$reference <- if ("reference" %in% colnames(df)) df$reference else reference
  df$mode <- if ("mode" %in% colnames(df)) df$mode else reference
  sgf_reorder_cols(df, c("variable", "distance", "value", "estimate", "reference", "mode"))
}

sgf_top_variables <- function(
  screening_out,
  significance,
  model_fits,
  nfeatures,
  sign_var,
  sign_threshold
) {
  vars <- sgf_try_spata_call("getSgsResultsVec", list(screening_out), required = FALSE)
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
  sgf_reorder_cols(top, c("variable", "rank", "fdr", "p_value", "best_model", "mae", "rmse"))
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

sgf_store_result <- function(srt, result_name, result, assay, set_variable_features) {
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
  if (is.null(srt@tools[["SpatialGradientFeatures"]])) {
    srt@tools[["SpatialGradientFeatures"]] <- list()
  }
  srt@tools[["SpatialGradientFeatures"]][[result_name]] <- result[expected]
  vars <- result$top_variables$variable
  vars <- vars[!is.na(vars) & nzchar(vars)]
  srt@misc[["SpatialGradientFeatures"]] <- vars
  srt@misc[["SpatialGradientFeaturesResult"]] <- result_name
  if (isTRUE(set_variable_features) && length(vars) > 0L) {
    SeuratObject::VariableFeatures(srt, assay = assay) <- vars
  }
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
  if (is.null(result_name)) {
    result_name <- srt@misc[["SpatialGradientFeaturesResult"]] %||% names(all_results)[length(all_results)]
  }
  if (!result_name %in% names(all_results)) {
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
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = "Gradient distance", y = "Expression", color = "Variable") +
    sgf_plot_theme(theme_use = theme_use, theme_args = theme_args) +
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
  df$variable <- factor(df$variable, levels = rev(df$variable))
  cols <- sgf_feature_colors(as.character(df$variable), palette = palette, palcolor = palcolor)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["rank"]], y = .data[["variable"]], color = .data[["variable"]])) +
    ggplot2::geom_point(ggplot2::aes(size = .data[["rmse"]]), na.rm = TRUE) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = "Rank", y = NULL, color = "Variable", size = "RMSE") +
    sgf_plot_theme(theme_use = theme_use, theme_args = theme_args) +
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
    ggplot2::scale_fill_gradientn(colors = rev(sgf_gradient_colors(palette, palcolor)), na.value = "grey85") +
    ggplot2::labs(x = "Model", y = NULL, fill = "RMSE") +
    sgf_plot_theme(theme_use = theme_use, theme_args = theme_args) +
    ggplot2::theme(
      legend.position = legend.position,
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

sgf_plot_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(ggplot2::theme_minimal())
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_fun <- if (is.character(theme_use)) {
    get(theme_use, mode = "function", inherits = TRUE)
  } else {
    theme_use
  }
  do.call(theme_fun, theme_args)
}

sgf_feature_colors <- function(features, palette = "Spectral", palcolor = NULL) {
  features <- unique(as.character(features))
  cols <- palette_colors(features, palette = palette, palcolor = palcolor)
  stats::setNames(cols[seq_along(features)], features)
}

sgf_gradient_colors <- function(palette = "Spectral", palcolor = NULL) {
  palette_colors(palette = palette, palcolor = palcolor, n = 9)
}

sgf_parameters_df <- function(parameters) {
  parameters$scop_version <- sgf_package_version("scop")
  parameters$spata2_version <- sgf_package_version(sgf_spata2_pkg())
  data.frame(
    key = names(parameters),
    value = vapply(parameters, sgf_collapse_value, character(1)),
    stringsAsFactors = FALSE
  )
}

sgf_package_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

sgf_collapse_value <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (length(x) == 0L) {
    return("")
  }
  paste(as.character(x), collapse = ",")
}

sgf_drop_nulls <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

sgf_as_df <- function(df) {
  if (!is.data.frame(df)) {
    df <- data.frame()
  }
  as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
}

sgf_rename_first <- function(df, target, candidates) {
  hit <- intersect(candidates, colnames(df))
  hit <- setdiff(hit, target)
  if (length(hit) > 0L && !target %in% colnames(df)) {
    colnames(df)[match(hit[1L], colnames(df))] <- target
  }
  df
}

sgf_reorder_cols <- function(df, first_cols) {
  first_cols <- intersect(first_cols, colnames(df))
  df[, c(first_cols, setdiff(colnames(df), first_cols)), drop = FALSE]
}
