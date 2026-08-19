# Inherited via @inheritParams.

#' @title Shared parameters
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay to use. `NULL` uses the default assay.
#' @param layer Assay layer to use.
#' @param verbose Print progress messages.
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param group.by Metadata column(s) to group cells by.
#' @param split.by Metadata column to split the analysis or plot by.
#' @param reduction Dimensionality reduction to use. `NULL` uses [DefaultReduction].
#' @param features Features to use (assay names or numeric metadata columns).
#' @param cells Cell names to include.
#' @param palette Color palette name. See [thisplot::show_palettes].
#' @param palcolor Custom colors used to build the palette.
#' @param species Species name, e.g. `"Homo_sapiens"` or `"Mus_musculus"`.
#' @param prefix Prefix for stored result names.
#' @param image Spatial image name. Required when multiple images are present;
#' a single image is selected automatically when `NULL`.
#' @param coord.cols Metadata coordinate columns used when no image is available.
#' @param coordinate_space Coordinate space for spatial distances. `"raw"` is
#' the default; `"legacy_display"` remains available for compatibility.
#' @param backend Implementation backend.
#' @param tool_name Name of the `srt@tools` slot used to store results.
#' @param overwrite Overwrite existing results.
#' @param ... Additional arguments passed to the underlying method.
#'
#' @keywords internal
#' @name scop-params
NULL
