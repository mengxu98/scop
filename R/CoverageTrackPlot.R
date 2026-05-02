#' @title Coverage track plot for ATAC data
#'
#' @description
#' Thin wrapper around `Signac::CoveragePlot` with scop defaults.
#'
#' @md
#' @inheritParams CellDimPlot
#' @param region Genomic region passed to `Signac::CoveragePlot`.
#' @param assay ATAC assay used for plotting.
#' @param extend.upstream,extend.downstream Distance to extend around the region.
#' @param annotation,peaks,links Whether to show gene annotation, peaks and links.
#' @param tile Whether to show fragment tiles in the coverage plot.
#' @param ranges Optional genomic ranges added as external tracks.
#' @param ranges.group.by Optional grouping variable used for `ranges`.
#' @param region.highlight Optional genomic ranges highlighted in the locus panel.
#' @param verbose Whether to print progress messages.
#' @param ... Additional parameters passed to `Signac::CoveragePlot`.
#'
#' @return A coverage plot object.
#' @export
#' @examples
#' \dontrun{
#' data("pbmcmultiome_sub", package = "scop")
#' # Coverage plotting requires an ATAC object with valid fragment information.
#' CoverageTrackPlot(
#'   pbmcmultiome_sub,
#'   region = rownames(pbmcmultiome_sub[["peaks"]])[1],
#'   assay = "peaks",
#'   group.by = "CellType"
#' )
#' }
CoverageTrackPlot <- function(
  srt,
  region,
  assay = NULL,
  group.by = NULL,
  palette = "Chinese",
  palcolor = NULL,
  extend.upstream = 1000,
  extend.downstream = 1000,
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  tile = FALSE,
  ranges = NULL,
  ranges.group.by = NULL,
  region.highlight = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    log_message(
      "{.arg assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }

  cols <- NULL
  if (!is.null(group.by)) {
    groups <- unique(as.character(stats::na.omit(srt[[group.by, drop = TRUE]])))
    if (length(groups) > 0) {
      cols <- palette_colors(
        groups,
        palette = palette,
        palcolor = palcolor
      )
    }
  }

  log_message(
    "Drawing ATAC coverage plot for {.val {region}}",
    verbose = verbose
  )

  Signac::CoveragePlot(
    object = srt,
    region = region,
    assay = assay,
    group.by = group.by,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream,
    annotation = annotation,
    peaks = peaks,
    links = links,
    tile = tile,
    ranges = ranges,
    ranges.group.by = ranges.group.by,
    region.highlight = region.highlight,
    cols = cols,
    ...
  )
}
