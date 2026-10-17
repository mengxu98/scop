#' @md
#' @param show_points Draw cell points or hexagonal bins and highlighted cell
#' points. Set to `FALSE` for contour-only or surface-only plots. Labels and
#' graph or lineage overlays are controlled separately. Defaults to `TRUE`.
#' @param add_density,density_color,density_filled,density_filled_palette,density_filled_palcolor
#' Add a smooth feature-expression layer. For unblended features, density
#' contours are colored by local Gaussian-weighted mean feature values,
#' including zero values. Contour shape represents cell density; color
#' represents local mean expression and may vary along one contour.
#' `density_filled = TRUE` draws the local mean expression surface with fine density contours instead.
#' With points visible, the original point-expression scale is preserved and
#' local means use that scale. Without points, the scale is based on the local
#' means. Cutoffs and `keep_scale` apply to the selected scale. `bg_cutoff` only
#' affects cell marks. When neither fill-color option is supplied, expression
#' surfaces inherit `palette` and `palcolor`; explicit fill options override
#' this default.
#' For compatibility, explicitly supplying `density_color` retains the old
#' fixed-color cell-density layer. Blended multi-feature plots using
#' `compare_features = TRUE` also retain the original cell-density layer,
#' including filled density bands and their original default colors.
#' Use [CellDimPlot()] for cell-distribution density layers.
#' @param force Draw even when more than 100 features are requested.
#' @details
#' Expression-colored contours use a 240 by 240 grid, Gaussian bandwidths from
#' [stats::bw.nrd()], and nine density levels from 8% to 92% of peak density.
#' Major contours use a wider stroke every third level. Grid limits extend three bandwidths beyond the finite coordinates to avoid
#' clipping outer contours. Local means are interpolated onto contour vertices;
#' filled surfaces omit grid points below 2% of peak kernel support, fade in
#' between 2% and 3.5% support, and interpolate raster pixels for smooth edges.
#' Cells with non-finite coordinates or feature values are omitted from this
#' calculation. Panels with fewer than three finite cells or no variation in
#' either coordinate return no contours with a warning. No additional public
#' plotting function or dependency is required.
