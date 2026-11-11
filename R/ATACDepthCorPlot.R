#' @title Plot ATAC reduction-depth correlations
#' @md
#' @description Compute and plot correlations between peak-count totals (or
#' explicit depth metadata) and reduced components, aligning cells by name.
#' @param object A Seurat object containing a ChromatinAssay.
#' @param assay Name of the ChromatinAssay; defaults to the active assay.
#' @param cells Cell names in the selected assay; NULL uses all assay cells.
#' @param reduction Name of a dimensional reduction, normally `"lsi"`.
#' @param n Maximum number of leading components; NULL uses all components.
#' @param depth.column Optional numeric non-negative metadata column. NULL
#'   recomputes totals from the selected assay's counts layer, avoiding stale
#'   nCount metadata. Peak-matrix totals are not necessarily total fragments or
#'   sequencing reads; select a verified fragment-count column explicitly when
#'   that is the desired depth measure.
#' @param method Correlation method: `"pearson"`, `"spearman"`, or `"kendall"`.
#' @param return.data Return the correlation table instead of a plot.
#' @return A ggplot, or a data.frame with component, correlation and n.cells when
#'   return.data is TRUE. Constant components have NA correlations. Nothing is
#'   written into the object and no component is automatically excluded.
#' @details All selected assay cells must occur in the reduction and have finite
#'   depth and embedding values; missing cells are not silently dropped. Depth
#'   must vary across at least three cells. Component numbers are signed
#'   correlations; SVD signs are arbitrary.
#'   This diagnostic supports, but does not by itself determine, LSI selection.
#' @seealso [RunATACQC()]
#' @export
#' @examples
#' \dontrun{
#' ATACDepthCorPlot(object, assay = "peaks", reduction = "lsi")
#' ATACDepthCorPlot(object, assay = "peaks", return.data = TRUE)
#' }
ATACDepthCorPlot <- function(object, assay = NULL, reduction = "lsi", n = 10L,
                             depth.column = NULL, method = "pearson",
                             cells = NULL, return.data = FALSE) {
  checked <- atac_plot_input(object, assay, cells)
  assay <- checked$assay
  cells <- checked$cells
  if (!is.character(reduction) || length(reduction) != 1L || is.na(reduction) ||
    !reduction %in% SeuratObject::Reductions(object)) {
    stop("reduction must name an existing dimensional reduction.", call. = FALSE)
  }
  embedding <- SeuratObject::Embeddings(object[[reduction]])
  if (is.null(rownames(embedding)) || anyDuplicated(rownames(embedding)) ||
    !all(cells %in% rownames(embedding))) {
    stop("The reduction must contain every selected assay cell exactly once.", call. = FALSE)
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n))) {
    stop("n must be a positive integer or NULL.", call. = FALSE)
  }
  n <- if (is.null(n)) ncol(embedding) else min(n, ncol(embedding))
  if (n < 1L) stop("The reduction has no components.", call. = FALSE)
  embedding <- embedding[cells, seq_len(n), drop = FALSE]
  if (is.null(depth.column)) {
    counts <- SeuratObject::GetAssayData(object, assay = assay, layer = "counts")
    if (!all(cells %in% colnames(counts))) {
      stop("The counts layer must contain every selected assay cell.", call. = FALSE)
    }
    depth <- Matrix::colSums(counts[, cells, drop = FALSE])
    depth.label <- "peak counts (counts layer)"
  } else {
    if (!is.character(depth.column) || length(depth.column) != 1L || is.na(depth.column) ||
      !depth.column %in% colnames(object[[]])) {
      stop("depth.column must name an existing numeric metadata column.", call. = FALSE)
    }
    depth <- object[[]][cells, depth.column]
    depth.label <- depth.column
  }
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  result <- atac_depth_cor_data(embedding, depth, method)
  if (!is.logical(return.data) || length(return.data) != 1L || is.na(return.data)) {
    stop("return.data must be TRUE or FALSE.", call. = FALSE)
  }
  if (return.data) {
    return(result)
  }
  ggplot2::ggplot(result, ggplot2::aes(x = .data$component, y = .data$correlation)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = result$component) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      x = "Component", y = paste(method, "correlation"),
      title = "Correlation between depth and reduced components",
      subtitle = paste("Assay:", assay, "| Reduction:", reduction, "| Depth:", depth.label)
    ) +
    ggplot2::theme_classic()
}

atac_depth_cor_data <- function(embedding, depth, method) {
  if (!is.numeric(depth) || length(depth) != nrow(embedding) || length(depth) < 3L ||
    any(!is.finite(depth)) || any(depth < 0) || any(!is.finite(embedding))) {
    stop("Depth and embeddings must be finite and aligned for at least three cells; depth must be non-negative.", call. = FALSE)
  }
  if (stats::sd(depth) == 0) stop("Depth must vary across selected cells.", call. = FALSE)
  correlations <- vapply(seq_len(ncol(embedding)), function(i) {
    if (stats::sd(embedding[, i]) == 0) {
      return(NA_real_)
    }
    stats::cor(embedding[, i], depth, method = method)
  }, numeric(1))
  data.frame(
    component = seq_len(ncol(embedding)), correlation = correlations,
    n.cells = nrow(embedding), row.names = NULL
  )
}
