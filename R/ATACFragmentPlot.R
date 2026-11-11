#' @title Plot ATAC fragment lengths
#' @md
#' @description Plot fragment lengths in an explicitly selected genomic region
#' using `Signac::FragmentHistogram()`. This is a regional diagnostic, not a
#' whole-library fragment-length estimate.
#' @param object A Seurat object containing a ChromatinAssay.
#' @param region Genomic region(s), as a GRanges or strings accepted by Signac.
#'   Required explicitly to avoid silently assuming chromosome naming or genome.
#' @param assay Name of the ChromatinAssay; defaults to the active assay.
#' @param group.by One metadata column used for grouping, or NULL for identities.
#' @param cells Cell names in the selected assay; NULL uses all assay cells.
#' @param log.scale Whether to use a logarithmic count axis.
#' @param ... Arguments passed to other functions.
#' @return A ggplot object. The input object is not modified.
#' @details Requires local bgzip fragment files with a `.tbi` or `.csi` index
#'   and explicit cell/barcode mappings. The histogram counts fragment
#'   records, not column-five read multiplicities, and displays lengths between
#'   0 and 800 bp by default.
#'   As in Signac, a fragment returned for several query intervals can contribute
#'   more than once; use a single region for an unambiguous regional histogram.
#' @seealso [RunATACQC()]
#' @export
#' @examples
#' \dontrun{
#' ATACFragmentPlot(object, region = "chr1-1-2000000", assay = "peaks")
#' }
ATACFragmentPlot <- function(object, region, assay = NULL, group.by = NULL,
                             cells = NULL, log.scale = FALSE, ...) {
  checked <- atac_plot_input(object, assay, cells, group.by)
  object <- checked$object
  assay <- checked$assay
  if (missing(region) || is.null(region) || length(region) == 0L) {
    stop("Supply an explicit non-empty genomic region for the fragment histogram.", call. = FALSE)
  }
  if (!inherits(region, "GRanges") &&
    (!is.character(region) || anyNA(region) || any(!nzchar(region)))) {
    stop("region must be genomic region strings or a GRanges object.", call. = FALSE)
  }
  if (!is.logical(log.scale) || length(log.scale) != 1L || is.na(log.scale)) {
    stop("log.scale must be TRUE or FALSE.", call. = FALSE)
  }
  fragments <- Signac::Fragments(object[[assay]])
  get_data <- Signac::GetFragmentData
  active <- vapply(fragments, function(fragment) {
    any(names(get_data(fragment, "cells")) %in% checked$cells)
  }, logical(1))
  fragments <- fragments[active]
  if (!length(fragments)) {
    stop("Fragment plotting requires mapped fragment files for the selected cells.", call. = FALSE)
  }
  maps <- lapply(fragments, function(fragment) {
    mapping <- get_data(fragment, "cells")
    mapping[names(mapping) %in% checked$cells]
  })
  mapped <- unlist(lapply(maps, names), use.names = FALSE)
  if (anyDuplicated(mapped) || !setequal(mapped, checked$cells) ||
    !all(vapply(maps, function(x) atac_unique_ids(x) && atac_unique_ids(names(x)), logical(1)))) {
    stop("Fragment mappings must cover each selected cell exactly once.", call. = FALSE)
  }
  paths <- vapply(fragments, get_data, character(1), slot = "path")
  indexed <- file.exists(paste0(paths, ".tbi")) | file.exists(paste0(paths, ".csi"))
  if (any(!file.exists(paths)) || any(!indexed)) {
    stop("Fragment plotting requires accessible local fragment files and .tbi or .csi indexes.", call. = FALSE)
  }
  if (anyDuplicated(normalizePath(paths, mustWork = TRUE))) {
    stop("Fragment objects must refer to different files, not aliases of the same path.", call. = FALSE)
  }
  # fastmatch may cache a hash directly on the barcode vector by reference.
  # Use the selected mapping copies so plotting cannot alter input attributes.
  for (i in seq_along(fragments)) fragments[[i]]@cells <- maps[[i]]
  object[[assay]]@fragments <- fragments
  plot <- Signac::FragmentHistogram(
    object = object, assay = assay, region = region, group.by = group.by,
    cells = checked$cells, log.scale = log.scale, ...
  )
  if (is.null(plot$data) || !nrow(plot$data)) {
    stop("No mapped fragments overlap the requested region; check chromosome names and coordinates.", call. = FALSE)
  }
  plot
}
