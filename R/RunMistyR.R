#' @title Run mistyR multiview spatial modeling
#'
#' @description
#' Build a small `mistyR` view composition from a spatial `Seurat` object,
#' train MISTy models, collect results, and store a SCOP-style result bundle in
#' `srt@tools`. The intraview is always created from the selected assay layer;
#' optional juxtaview and paraview components describe local and broader spatial
#' context. `mistyR` is an optional Bioconductor dependency installable with
#' `BiocManager::install("mistyR")`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param assay Assay used for expression. If `NULL`, the default assay is used.
#' @param layer Assay layer used for expression values.
#' @param features Features used by MISTy. If `NULL`, variable features are used
#' when available; otherwise all assay features are used.
#' @param image Name of the Seurat spatial image. Required when multiple images
#' are present; a single image is selected automatically when `NULL`.
#' @param coord.cols Metadata coordinate columns used when no Seurat image
#' coordinates are available.
#' @param coordinate_space Coordinate system used to build MISTy views. The
#' default is raw acquisition coordinates; `"legacy_display"` remains an
#' explicit compatibility option.
#' @param views Spatial views to add besides the required intraview. One or both
#' of `"para"` and `"juxta"`.
#' @param para_l,para_zoi,para_family,para_approx,para_nn Parameters passed to
#' `mistyR::add_paraview()`. `para_l` and `para_zoi` use the selected coordinate
#' units; `para_nn` is a unitless neighbor count.
#' @param juxta_neighbor_thr Neighbor threshold passed to
#' `mistyR::add_juxtaview()`, expressed in the selected coordinate units.
#' @param view_cached Whether generated mistyR views should use cache.
#' @param results_folder Folder passed to `mistyR::run_misty()`. If `NULL`, a
#' temporary folder is used.
#' @param seed,target_subset,bypass_intra,cv_folds,model_cached,append
#' Parameters passed to `mistyR::run_misty()`.
#' @param tool_name Name used to store results in `srt@tools`.
#' @param store_results Whether to store results in `srt@tools`.
#' @param store_views Whether to store the mistyR view composition in
#' `srt@tools`. This can be large.
#' @param ... Additional named arguments passed to `mistyR::run_misty()`.
#'
#' @return A `Seurat` object with results stored in `srt@tools[[tool_name]]`
#' when `store_results = TRUE`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#'
#' spatial <- RunMistyR(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "data",
#'   features = rownames(spatial)[1:10],
#'   coord.cols = c("x", "y"),
#'   views = "para",
#'   para_l = 5,
#'   cv_folds = 3,
#'   verbose = FALSE
#' )
#' spatial@tools$MistyR$summary
RunMistyR <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  views = "para",
  para_l = 10,
  para_zoi = 0,
  para_family = c("gaussian", "exponential", "linear", "constant"),
  para_approx = 1,
  para_nn = NULL,
  juxta_neighbor_thr = 15,
  view_cached = FALSE,
  results_folder = NULL,
  seed = 42,
  target_subset = NULL,
  bypass_intra = FALSE,
  cv_folds = 10,
  model_cached = FALSE,
  append = FALSE,
  tool_name = "MistyR",
  store_results = TRUE,
  store_views = FALSE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
) {
  coordinate_space <- match.arg(coordinate_space)
  log_message(
    "Running mistyR multiview spatial modeling",
    message_type = "running",
    verbose = verbose
  )
  mistyr_validate_srt(srt)
  mistyr_assert_string(tool_name, "tool_name")
  mistyr_assert_flag(view_cached, "view_cached")
  mistyr_assert_flag(bypass_intra, "bypass_intra")
  mistyr_assert_flag(model_cached, "model_cached")
  mistyr_assert_flag(append, "append")
  mistyr_assert_flag(store_results, "store_results")
  mistyr_assert_flag(store_views, "store_views")
  para_family <- match.arg(para_family)
  views <- mistyr_validate_views(views)
  para_l <- mistyr_assert_positive_number(para_l, "para_l")
  para_zoi <- mistyr_assert_nonnegative_number(para_zoi, "para_zoi")
  para_approx <- mistyr_assert_positive_number(para_approx, "para_approx")
  if (!is.null(para_nn)) {
    para_nn <- mistyr_validate_positive_integer(para_nn, "para_nn")
  }
  juxta_neighbor_thr <- mistyr_assert_positive_number(juxta_neighbor_thr, "juxta_neighbor_thr")
  cv_folds <- mistyr_validate_positive_integer(cv_folds, "cv_folds")
  extra_args <- list(...)
  mistyr_validate_named_list(extra_args, "...")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  input <- mistyr_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )

  check_r("mistyR", verbose = FALSE)
  create_initial_view <- get_namespace_fun("mistyR", "create_initial_view")
  add_paraview <- get_namespace_fun("mistyR", "add_paraview")
  add_juxtaview <- get_namespace_fun("mistyR", "add_juxtaview")
  run_misty <- get_namespace_fun("mistyR", "run_misty")
  collect_results <- get_namespace_fun("mistyR", "collect_results")

  misty_views <- create_initial_view(input$expression)
  if ("juxta" %in% views) {
    misty_views <- add_juxtaview(
      misty_views,
      input$positions,
      neighbor.thr = juxta_neighbor_thr,
      cached = view_cached,
      verbose = verbose
    )
  }
  if ("para" %in% views) {
    misty_views <- add_paraview(
      misty_views,
      input$positions,
      l = para_l,
      zoi = para_zoi,
      family = para_family,
      approx = para_approx,
      nn = para_nn,
      cached = view_cached,
      verbose = verbose
    )
  }

  results_folder <- results_folder %||% tempfile("scop-mistyr-")
  run_args <- list(
    views = misty_views,
    results.folder = results_folder,
    seed = seed,
    target.subset = target_subset,
    bypass.intra = bypass_intra,
    cv.folds = cv_folds,
    cached = model_cached,
    append = append
  )
  run_args <- utils::modifyList(run_args, extra_args)
  run_path <- do.call(run_misty, run_args)
  misty_results <- collect_results(run_path)
  summary <- mistyr_summary(misty_results, views = names(misty_views))

  if (isTRUE(store_results)) {
    bundle <- list(
      results = misty_results,
      summary = summary,
      results_folder = run_path,
      features = input$features,
      cells = input$cells,
      parameters = list(
        assay = assay,
        layer = layer,
        image = image,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space,
        views = views,
        para_l = para_l,
        para_zoi = para_zoi,
        para_family = para_family,
        para_approx = para_approx,
        para_nn = para_nn,
        juxta_neighbor_thr = juxta_neighbor_thr,
        view_cached = view_cached,
        results_folder = results_folder,
        seed = seed,
        target_subset = target_subset,
        bypass_intra = bypass_intra,
        cv_folds = cv_folds,
        model_cached = model_cached,
        append = append,
        tool_name = tool_name,
        store_results = store_results,
        store_views = store_views,
        backend_args = extra_args
      )
    )
    if (isTRUE(store_views)) {
      bundle$views <- misty_views
    }
    bundle <- spatial_result_build(
      bundle = bundle,
      method = "MistyR",
      result_type = "neighborhood",
      provenance = list(producer = "RunMistyR", backend_id = "mistyr")
    )
    srt@tools[[tool_name]] <- bundle
  }

  log_message(
    "{.pkg mistyR} completed for {.val {length(input$features)}} feature{?s}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

mistyr_prepare_input <- function(
  srt,
  assay,
  layer,
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("raw", "legacy_display")
) {
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  cells <- colnames(srt)[colnames(srt) %in% rownames(coords)]
  if (length(cells) == 0L) {
    log_message("No cells or spots match spatial coordinates", message_type = "error")
  }
  coords <- coords[cells, , drop = FALSE]
  keep <- is.finite(coords$x) & is.finite(coords$y)
  if (!all(keep)) {
    cells <- cells[keep]
    coords <- coords[cells, , drop = FALSE]
  }
  if (length(cells) < 3L) {
    log_message("At least three cells or spots with finite coordinates are required", message_type = "error")
  }
  features <- mistyr_resolve_features(srt, assay = assay, features = features)
  expr <- expr[features, cells, drop = FALSE]
  expr <- as.matrix(expr)
  expr[!is.finite(expr)] <- 0
  expression <- as.data.frame(t(expr), check.names = FALSE)
  positions <- data.frame(
    x = as.numeric(coords$x),
    y = as.numeric(coords$y),
    row.names = cells,
    stringsAsFactors = FALSE
  )
  list(
    expression = expression,
    positions = positions,
    features = features,
    cells = cells
  )
}

mistyr_resolve_features <- function(srt, assay, features = NULL) {
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(srt, assay = assay)
    if (length(features) == 0L) {
      features <- rownames(srt[[assay]])
    }
  }
  features <- unique(as.character(features))
  missing <- setdiff(features, rownames(srt[[assay]]))
  if (length(missing) > 0L) {
    log_message(
      "{.arg features} not found in assay {.val {assay}}: {.val {missing}}",
      message_type = "error"
    )
  }
  if (length(features) < 2L) {
    log_message(
      "At least two features are required for {.fn RunMistyR}",
      message_type = "error"
    )
  }
  features
}

mistyr_summary <- function(results, views = character()) {
  improvements <- results$improvements %||% data.frame()
  contributions <- results$contributions %||% data.frame()
  importances <- results$importances.aggregated %||% list()
  list(
    views = views,
    n_targets = if ("target" %in% colnames(improvements)) length(unique(improvements$target)) else NA_integer_,
    n_improvement_records = if (is.data.frame(improvements)) nrow(improvements) else NA_integer_,
    n_contribution_records = if (is.data.frame(contributions)) nrow(contributions) else NA_integer_,
    importance_views = names(importances)
  )
}

#' @title Plot stored MISTy results
#'
#' @description Plot model improvements or view contributions from a result
#' produced by [RunMistyR()] without rerunning the backend.
#'
#' @param object Optional `Seurat` object containing `MistyR` results.
#' @param res Optional result list, usually `object@tools$MistyR`.
#' @param type Result table to plot.
#' @param top_n Maximum number of records shown after ranking by absolute value.
#' @param target Optional target feature filter.
#' @return A `ggplot` object.
#' @export
MistyRPlot <- function(
  object = NULL,
  res = NULL,
  type = c("improvements", "contributions"),
  top_n = 20,
  target = NULL
) {
  type <- match.arg(type)
  if (is.null(res)) {
    if (is.null(object) || !inherits(object, "Seurat")) {
      log_message("Provide a {.cls Seurat} {.arg object} or a MistyR {.arg res}", message_type = "error")
    }
    res <- object@tools[["MistyR"]]
  }
  tab <- res$results[[type]] %||% NULL
  if (!is.data.frame(tab) || nrow(tab) == 0L) {
    log_message("MistyR result {.val {type}} is empty", message_type = "error")
  }
  if (!is.null(target) && "target" %in% colnames(tab)) {
    tab <- tab[as.character(tab$target) %in% as.character(target), , drop = FALSE]
  }
  value_col <- intersect(c("value", "gain.R2", "contribution", "importance"), colnames(tab))[1L]
  if (is.na(value_col)) {
    numeric_cols <- colnames(tab)[vapply(tab, is.numeric, logical(1))]
    value_col <- if (length(numeric_cols)) numeric_cols[[1L]] else NA_character_
  }
  if (is.na(value_col) || nrow(tab) == 0L) {
    log_message("MistyR {.val {type}} has no plottable numeric records", message_type = "error")
  }
  tab$.value <- as.numeric(tab[[value_col]])
  tab <- tab[is.finite(tab$.value), , drop = FALSE]
  label_cols <- intersect(c("target", "view", "measure"), colnames(tab))
  tab$.label <- if (length(label_cols)) do.call(paste, c(tab[label_cols], sep = " / ")) else rownames(tab)
  thisplot::StatPlot(
    meta.data = tab,
    stat.by = ".label",
    value.by = ".value",
    stat_type = "value",
    plot_type = "bar",
    top_n = top_n,
    rank.by = "absolute",
    bar_fill = "#3C8DBC",
    bar_width = 0.75,
    title = paste("MISTy", type),
    ylab = value_col,
    flip = TRUE,
    theme_use = theme_scop
  )
}

mistyr_validate_views <- function(views) {
  if (is.null(views)) {
    return(character())
  }
  views <- unique(as.character(views))
  allowed <- c("para", "juxta")
  missing <- setdiff(views, allowed)
  if (length(missing) > 0L) {
    log_message(
      "{.arg views} must contain only {.val para} and/or {.val juxta}. Invalid: {.val {missing}}",
      message_type = "error"
    )
  }
  views
}

mistyr_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  invisible(TRUE)
}

mistyr_assert_string <- function(x, arg) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    log_message("{.arg {arg}} must be a single non-empty string", message_type = "error")
  }
  invisible(TRUE)
}

mistyr_assert_flag <- function(x, arg) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    log_message("{.arg {arg}} must be TRUE or FALSE", message_type = "error")
  }
  invisible(TRUE)
}

mistyr_assert_positive_number <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    log_message("{.arg {arg}} must be a positive finite number", message_type = "error")
  }
  as.numeric(x)
}

mistyr_assert_nonnegative_number <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    log_message("{.arg {arg}} must be a non-negative finite number", message_type = "error")
  }
  as.numeric(x)
}

mistyr_validate_positive_integer <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1) {
    log_message("{.arg {arg}} must be a positive integer", message_type = "error")
  }
  as.integer(x)
}

mistyr_validate_named_list <- function(x, arg) {
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message("{.arg {arg}} must contain named arguments only", message_type = "error")
  }
  invisible(TRUE)
}
