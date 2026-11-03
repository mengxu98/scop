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
#' @param ... Additional arguments to `Signac::FragmentHistogram()`.
#' @return A ggplot object. The input object is not modified.
#' @details Requires local bgzip fragment files with a `.tbi` or `.csi` index
#'   and explicit cell/barcode mappings. The native histogram counts fragment
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
  checked <- atac_fragment_plot_input(object, assay, cells, group.by)
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
  fragments <- atac_backend_fun("Signac", "Fragments")(object[[assay]])
  get_data <- atac_backend_fun("Signac", "GetFragmentData")
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
  plot <- atac_backend_fun("Signac", "FragmentHistogram")(
    object = object, assay = assay, region = region, group.by = group.by,
    cells = checked$cells, log.scale = log.scale, ...
  )
  if (is.null(plot$data) || !nrow(plot$data)) {
    stop("No mapped fragments overlap the requested region; check chromosome names and coordinates.", call. = FALSE)
  }
  plot
}

atac_fragment_plot_input <- function(object, assay = NULL, cells = NULL, group.by = NULL) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.", call. = FALSE)
  assay <- if (is.null(assay)) SeuratObject::DefaultAssay(object) else assay
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) ||
      !assay %in% SeuratObject::Assays(object) || !inherits(object[[assay]], "ChromatinAssay")) {
    stop("assay must name a ChromatinAssay.", call. = FALSE)
  }
  assay.cells <- colnames(object[[assay]])
  cells <- if (is.null(cells)) assay.cells else cells
  if (!atac_unique_ids(cells) || !all(cells %in% assay.cells)) {
    stop("cells must be unique cell names belonging to the selected assay.", call. = FALSE)
  }
  if (!is.null(group.by) && (!is.character(group.by) || length(group.by) != 1L ||
      is.na(group.by) || !group.by %in% colnames(object[[]]))) {
    stop("group.by must name one existing metadata column.", call. = FALSE)
  }
  groups <- if (is.null(group.by)) SeuratObject::Idents(object)[cells] else object[[]][cells, group.by]
  if (anyNA(groups)) stop("Selected cells must have non-missing group labels.", call. = FALSE)
  list(object = object, assay = assay, cells = cells)
}
