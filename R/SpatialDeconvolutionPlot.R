#' @title Plot stored spatial deconvolution proportions
#'
#' @description
#' Plot spot-by-cell-type proportions stored by [RunRCTD()], [RunCARD()],
#' [RunSPOTlight()], or [RunSpatialDWLS()]. The plot reads a schema-v1 result
#' through [GetSpatialResult()] and never reruns a deconvolution backend.
#' [RunCSIDE()] is intentionally excluded because its output represents
#' differential or context effects rather than cell-type proportions.
#'
#' @md
#' @param srt A spatial `Seurat` object containing a stored deconvolution result.
#' @param tool_name Exact key in `srt@tools`. If `NULL`, exactly one compatible
#' stored result must be discoverable.
#' @param cell_types Optional cell types to display. The default uses all stored
#' cell types.
#' @param plot_type Plot proportions as separate point maps, one dominant-type
#' map derived from the stored proportions, or one spot-level pie map.
#' @param combine Whether to combine point maps. If `FALSE`, return a named list.
#' @param nrow,ncol,byrow Point-map layout controls. When both dimensions are
#' `NULL`, a near-square layout with at most three columns is used.
#' @param ... Additional arguments passed to [SpatialSpotPlot()].
#'
#' @return A `ggplot`, `patchwork`, or named list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' data(pancreas_sub)
#' shared <- head(intersect(
#'   rownames(visium_human_pancreas_sub),
#'   rownames(pancreas_sub)
#' ), 40)
#' spatial <- RunSpatialDWLS(
#'   visium_human_pancreas_sub[shared, 1:20],
#'   reference = pancreas_sub,
#'   reference_label = "CellType",
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
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  plot_type <- match.arg(plot_type)
  tool_name <- spatial_deconvolution_resolve_tool(srt, tool_name)
  stored <- GetSpatialResult(srt, tool_name = tool_name)
  allowed <- c("RunRCTD", "RunCARD", "RunSPOTlight", "RunSpatialDWLS")
  producer <- stored$provenance$producer %||% NA_character_
  if (length(producer) != 1L || is.na(producer) || !producer %in% allowed) {
    log_message(
      "Stored result {.val {tool_name}} is not a supported spatial proportion result",
      message_type = "error"
    )
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
    return(SpatialSpotPlot(srt, values = proportions, plot_type = "pie", ...))
  }
  if (identical(plot_type, "dominant")) {
    dominant <- spatial_deconvolution_dominant(proportions)
    return(SpatialSpotPlot(
      srt,
      values = stats::setNames(dominant, rownames(proportions)),
      plot_type = "point",
      ...
    ))
  }
  plots <- SpatialSpotPlot(
    srt,
    values = proportions,
    plot_type = "point",
    combine = FALSE,
    ...
  )
  legend_title <- list(...)$legend.title %||% "Proportion"
  plots <- Map(
    function(plot, cell_type) {
      spatial_deconvolution_set_point_scale(
        plot = plot,
        limits = c(0, 1),
        title = legend_title
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

spatial_deconvolution_resolve_tool <- function(srt, tool_name = NULL) {
  if (!is.null(tool_name)) {
    if (length(tool_name) != 1L || !is.character(tool_name) || is.na(tool_name) || !nzchar(tool_name)) {
      log_message("{.arg tool_name} must be a single non-empty string", message_type = "error")
    }
    return(tool_name)
  }
  index <- spatial_result_index(srt)
  allowed <- c("RunRCTD", "RunCARD", "RunSPOTlight", "RunSpatialDWLS")
  matches <- index$tool_name[index$registry_method %in% allowed]
  if (length(matches) != 1L) {
    log_message(
      "Select exactly one compatible stored result with {.arg tool_name}; found {.val {length(matches)}}",
      message_type = "error"
    )
  }
  matches[[1L]]
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

spatial_deconvolution_set_point_scale <- function(plot, limits, title) {
  scale_index <- which(vapply(
    plot$scales$scales,
    function(scale) any(scale$aesthetics %in% c("colour", "color")),
    logical(1)
  ))
  if (length(scale_index) != 1L) {
    log_message(
      "Unable to identify the continuous proportion color scale",
      message_type = "error"
    )
  }
  plot$scales$scales[[scale_index]]$limits <- limits
  plot$scales$scales[[scale_index]]$name <- title
  plot$labels$colour <- title
  plot
}
