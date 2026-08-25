#' @title Plot stored spatial deconvolution proportions
#'
#' @description
#' Plot spot-by-cell-type proportions stored by [RunRCTD()], [RunCARD()],
#' [RunSPOTlight()], or [RunSpatialDWLS()]. The plot reads the stored result
#' directly from `srt@tools[[tool_name]]` and never reruns a deconvolution
#' backend. [RunCSIDE()] is intentionally excluded because its output
#' represents differential or context effects rather than cell-type
#' proportions.
#'
#' @md
#' @param srt A spatial `Seurat` object containing a stored deconvolution result.
#' @param tool_name Explicit non-empty key in `srt@tools`. Results are never
#' discovered implicitly.
#' @param cell_types Optional cell types to display. The default uses all stored
#' cell types.
#' @param plot_type Plot proportions as separate point maps, one dominant-type
#' map derived from the stored proportions, or one spot-level pie map.
#' @param combine Whether to combine point maps. If `FALSE`, return a named list.
#' @param nrow,ncol,byrow Point-map layout controls. When both dimensions are
#' `NULL`, a near-square layout with at most three columns is used.
#' @param image.scale Image scale factor matching the selected raster.
#' @param ... Additional arguments passed to [SpatialSpotPlot()].
#'
#' @return A `ggplot`, `patchwork`, or named list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' data(panc8_sub)
#' keep_spots <- unique(round(seq(1, ncol(visium_human_pancreas_sub), length.out = 200)))
#' spatial <- visium_human_pancreas_sub[, keep_spots]
#' reference <- panc8_sub[, panc8_sub$celltype %in% c("ductal", "alpha", "beta")]
#' reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
#' shared <- intersect(
#'   SeuratObject::VariableFeatures(reference),
#'   rownames(spatial)
#' )
#' spatial <- RunSpatialDWLS(
#'   spatial[shared, ],
#'   reference = reference,
#'   reference_label = "celltype",
#'   features = shared,
#'   coord.cols = c("x", "y"),
#'   normalize = FALSE,
#'   verbose = FALSE
#' )
#' SpatialDeconvolutionPlot(
#'   spatial,
#'   tool_name = "SpatialDWLS",
#'   cell_types = colnames(spatial@tools$SpatialDWLS$proportions)[1],
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
SpatialDeconvolutionPlot <- function(
  srt,
  tool_name = NULL,
  cell_types = NULL,
  plot_type = c("point", "dominant", "pie"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...,
  image.scale = c("lowres", "hires")
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  image.scale <- match.arg(image.scale)
  if (
    is.null(tool_name) || !is.character(tool_name) || length(tool_name) != 1L ||
      is.na(tool_name) || !nzchar(tool_name)
  ) {
    log_message(
      "{.arg tool_name} must be an explicit non-empty key in {.code srt@tools}",
      message_type = "error"
    )
  }
  stored <- srt@tools[[tool_name]]
  if (!is.list(stored)) {
    log_message(
      "Stored result {.val {tool_name}} is not a spatial deconvolution result",
      message_type = "error"
    )
  }
  if (
    !is.null(stored$parameters$coordinate_space) ||
      !is.null(stored$coords) || !is.null(stored$coordinates)
  ) {
    spatial_require_coordinate_contract(stored, paste0(tool_name, " producer"))
  }
  proportions <- spatial_deconvolution_proportions(
    stored$proportions %||% stored$weights,
    spot_ids = colnames(srt),
    tool_name = tool_name
  )
  if (!is.null(cell_types)) {
    if (
      !is.character(cell_types) || length(cell_types) == 0L ||
        anyNA(cell_types) || any(!nzchar(cell_types)) || anyDuplicated(cell_types)
    ) {
      log_message("{.arg cell_types} must contain unique non-empty names", message_type = "error")
    }
    missing_types <- setdiff(cell_types, colnames(proportions))
    if (length(missing_types) > 0L) {
      log_message("Unknown {.arg cell_types}: {.val {missing_types}}", message_type = "error")
    }
    proportions <- proportions[, cell_types, drop = FALSE]
  }
  if (identical(plot_type, "pie")) {
    return(SpatialSpotPlot(
      srt,
      values = proportions,
      plot_type = "pie",
      image.scale = image.scale,
      ...
    ))
  }
  if (identical(plot_type, "dominant")) {
    dominant <- spatial_deconvolution_dominant(proportions)
    return(SpatialSpotPlot(
      srt,
      values = stats::setNames(dominant, rownames(proportions)),
      plot_type = "point",
      image.scale = image.scale,
      ...
    ))
  }
  plots <- SpatialSpotPlot(
    srt,
    values = proportions,
    plot_type = "point",
    image.scale = image.scale,
    combine = FALSE,
    ...
  )
  legend_title <- list(...)$legend.title %||% "Proportion"
  plots <- Map(
    function(plot, cell_type) {
      set_continuous_color_scale(
        plot = plot,
        limits = c(0, 1),
        title = legend_title,
        context = "proportion"
      ) + ggplot2::labs(title = cell_type)
    },
    plots,
    colnames(proportions)
  )
  if (isFALSE(combine)) {
    return(plots)
  }
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- min(3L, ceiling(sqrt(length(plots))))
  }
  patchwork::wrap_plots(
    plots,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    guides = "collect"
  ) + patchwork::plot_annotation(title = paste0(tool_name, " proportions"))
}

spatial_deconvolution_proportions <- function(x, spot_ids, tool_name) {
  if (is.null(x) || (!is.matrix(x) && !inherits(x, "Matrix") && !is.data.frame(x))) {
    log_message(
      "Stored result {.val {tool_name}} does not contain matrix-like proportions",
      message_type = "error"
    )
  }
  x <- as.matrix(x)
  if (nrow(x) == 0L || ncol(x) == 0L) {
    log_message("Stored proportions are empty", message_type = "error")
  }
  if (!is.numeric(x)) {
    log_message("Stored proportions must be numeric", message_type = "error")
  }
  valid_names <- function(nm) {
    !is.null(nm) && !anyNA(nm) && all(nzchar(nm)) && !anyDuplicated(nm)
  }
  if (!valid_names(rownames(x))) {
    log_message("Stored proportions must have unique, non-missing spot names", message_type = "error")
  }
  if (!valid_names(colnames(x))) {
    log_message("Stored proportions must have unique, non-missing cell-type names", message_type = "error")
  }
  missing_spots <- setdiff(spot_ids, rownames(x))
  extra_spots <- setdiff(rownames(x), spot_ids)
  if (length(missing_spots) > 0L || length(extra_spots) > 0L) {
    log_message(
      "Stored proportions are stale or incomplete: {.val {length(missing_spots)}} missing and {.val {length(extra_spots)}} unknown spots",
      message_type = "error"
    )
  }
  finite <- x[is.finite(x)]
  if (length(finite) > 0L && (any(finite < 0) || any(finite > 1 + sqrt(.Machine$double.eps)))) {
    log_message("Stored proportions must lie between zero and one", message_type = "error")
  }
  x[spot_ids, , drop = FALSE]
}

spatial_deconvolution_dominant <- function(proportions) {
  out <- rep(NA_character_, nrow(proportions))
  for (i in seq_len(nrow(proportions))) {
    values <- proportions[i, ]
    values[!is.finite(values)] <- NA_real_
    if (all(is.na(values)) || max(values, na.rm = TRUE) <= 0) {
      next
    }
    out[[i]] <- colnames(proportions)[which.max(replace(values, is.na(values), -Inf))]
  }
  factor(out, levels = colnames(proportions))
}
