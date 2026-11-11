#' @title Plot ATAC TSS enrichment profiles
#' @md
#' @description Plot previously computed per-position TSS enrichment with
#' `Signac::TSSPlot()`. A scalar TSS enrichment score alone is insufficient.
#' @param object A Seurat object containing a ChromatinAssay.
#' @param assay Name of the ChromatinAssay; defaults to the active assay.
#' @param cells Cell names in the selected assay; NULL uses all assay cells.
#' @param group.by One metadata column used for grouping, or NULL for identities.
#' @param idents Group labels to retain, or NULL for all groups.
#' @return A ggplot object. The input object is not modified.
#' @details First calculate positional enrichment with
#'   `Signac::TSSEnrichment(..., fast = FALSE)` or
#'   `RunATACQC(..., tss.positions = tss, fast = FALSE)`. Supplying explicit TSS
#'   positions to RunATACQC forces recalculation even if imported scores exist.
#'   Use TSS coordinates from the matching reference and recompute positional
#'   enrichment after changing reference coordinates or preprocessing.
#'   This function never silently recomputes enrichment or invents a curve from
#'   scalar scores. Only cells in the selected assay are included. Signac 1.17
#'   marks TSSPlot deprecated; any upstream lifecycle warning is retained.
#' @seealso [RunATACQC()]
#' @export
#' @examples
#' \dontrun{
#' object <- RunATACQC(object, tss.positions = tss, fast = FALSE)
#' ATACTSSPlot(object, assay = "peaks", group.by = "celltype")
#' }
ATACTSSPlot <- function(object, assay = NULL, group.by = NULL, idents = NULL,
                        cells = NULL) {
  checked <- atac_plot_input(object, assay, cells, group.by)
  object <- checked$object
  assay <- checked$assay
  groups <- if (is.null(group.by)) SeuratObject::Idents(object)[checked$cells] else object[[]][checked$cells, group.by]
  if (!is.null(idents) &&
    (!is.character(idents) || !length(idents) || anyNA(idents) ||
      any(!idents %in% as.character(groups)))) {
    stop("idents must contain existing non-missing group labels.", call. = FALSE)
  }
  matrix <- object[[assay]]@positionEnrichment[["TSS"]]
  if (is.null(matrix) || length(dim(matrix)) != 2L || ncol(matrix) < 2L) {
    stop("TSS position enrichment is unavailable; run TSSEnrichment(fast = FALSE) or RunATACQC(tss.positions = tss, fast = FALSE) first.", call. = FALSE)
  }
  if (is.null(rownames(matrix)) || anyDuplicated(rownames(matrix)) ||
    !all(checked$cells %in% rownames(matrix))) {
    stop("TSS position enrichment must have unique rows for every selected assay cell.", call. = FALSE)
  }
  positions <- suppressWarnings(as.numeric(colnames(matrix)))
  if (length(positions) != ncol(matrix) || any(!is.finite(positions)) || anyDuplicated(positions)) {
    stop("TSS position enrichment must have unique numeric position column names.", call. = FALSE)
  }
  matrix <- matrix[checked$cells, , drop = FALSE]
  if (any(!is.finite(matrix))) {
    stop("TSS position enrichment contains non-finite values for selected cells.", call. = FALSE)
  }
  if (!setequal(checked$cells, colnames(object))) {
    # Avoid native one-cell matrix dimension dropping during assay subsetting.
    # These caches are not needed to subset counts and are restored below.
    object[[assay]]@positionEnrichment <- list()
    object <- subset(object, cells = checked$cells)
  }
  object[[assay]]@positionEnrichment[["TSS"]] <- matrix
  Signac::TSSPlot(
    object = object, assay = assay, group.by = group.by, idents = idents
  )
}
