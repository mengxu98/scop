#' Profile cell neighborhoods at several distances
#'
#' Count labelled neighbors of selected cells without removing the surrounding
#' tissue. A native spatial index directly accumulates counts at all distances;
#' no cell-cell edge table is stored. This is an observed composition summary,
#' not a colocalization test, an optimal-radius selector, or a communication model.
#'
#' @param object A Seurat object. The object is not modified.
#' @param group.by Metadata column containing nonmissing cell or spot labels.
#' @param radii Positive, strictly increasing distances in raw coordinate units.
#'   These are not necessarily micrometers. For spot-based assays, counts refer
#'   to labelled spots, not inferred numbers of cells.
#' @param cells Unique target cell IDs, in the desired output order. NULL uses
#'   all cells in the resolved image/context. Neighbors always come from the
#'   complete resolved context, not just these targets.
#' @param sample.by Metadata column identifying independent samples. NULL treats
#'   the resolved context as one sample; supply this for pooled metadata-only
#'   coordinates to prevent cross-sample neighbors.
#' @param image Spatial image name. Required when multiple images are present.
#' @param coord.cols Metadata coordinate columns when no image is available.
#' @param cumulative FALSE returns shells [0, r1], (r1, r2], ...; TRUE returns
#'   cumulative neighborhoods [0, r] at each radius.
#' @inheritParams thisutils::log_message
#'
#' @return A data.frame with cell_id, sample, group, lower, radius, count, total
#'   and fraction columns. The denominator includes all labelled neighbors in
#'   that shell/radius. No neighbors gives count = 0 and fraction = NA. Self is
#'   excluded by cell ID, but distinct cells at identical coordinates are kept.
#'   Unselected cells have no output rows. The source and parameters attributes
#'   describe coordinates and the calculation. There is no tissue-edge correction
#'   or dedicated plot; the table can be used with ordinary R plotting functions.
#' @export
#' @examples
#' # Constructed coordinates and labels, not a biological example.
#' counts <- matrix(1L, 2, 4, dimnames = list(c("g1", "g2"), letters[1:4]))
#' object <- SeuratObject::CreateSeuratObject(counts)
#' object$col <- c(0, 1, 3, 5)
#' object$row <- c(0, 0, 0, 0)
#' object$celltype <- c("T", "M", "M", "T")
#' profile <- SpatialNeighborhoodProfile(
#'   object, "celltype", radii = c(1, 3), cells = "a", verbose = FALSE
#' )
#' profile
SpatialNeighborhoodProfile <- function(
  object, group.by, radii, cells = NULL, sample.by = NULL, image = NULL,
  coord.cols = c("col", "row"), cumulative = FALSE, verbose = TRUE
) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object", call. = FALSE)
  validate_scalar_string(group.by, "group.by")
  if (!is.null(sample.by)) validate_scalar_string(sample.by, "sample.by")
  if (!is.numeric(radii) || !length(radii) || any(!is.finite(radii)) ||
      any(radii <= 0) || any(diff(radii) <= 0)) {
    stop("radii must be positive, finite and strictly increasing", call. = FALSE)
  }
  if (!is.logical(cumulative) || length(cumulative) != 1L || is.na(cumulative)) {
    stop("cumulative must be one nonmissing logical value", call. = FALSE)
  }
  verbose <- thisutils::get_verbose(verbose)
  meta <- object[[]]
  if (!all(c(group.by, sample.by) %in% names(meta))) {
    stop("group.by and sample.by must identify metadata columns", call. = FALSE)
  }
  resolved <- SpatialCoordinates(object, image = image, coord.cols = coord.cols, space = "raw")
  co <- resolved$data
  idx <- match(co$cell_id, rownames(meta))
  if (!nrow(co) || anyNA(idx) || anyNA(co$cell_id) || any(!nzchar(co$cell_id)) ||
      anyDuplicated(co$cell_id)) stop("invalid context cell IDs", call. = FALSE)
  labels <- as.character(meta[[group.by]][idx])
  samples <- if (is.null(sample.by)) rep("all", nrow(co)) else as.character(meta[[sample.by]][idx])
  if (anyNA(labels) || anyNA(samples) || any(!nzchar(labels)) || any(!nzchar(samples))) {
    stop("context labels and sample IDs must be nonmissing and nonempty", call. = FALSE)
  }
  if (is.null(cells)) cells <- co$cell_id
  if (!is.character(cells) || anyNA(cells) || anyDuplicated(cells) ||
      !all(cells %in% co$cell_id)) stop("cells must be unique IDs in the resolved context", call. = FALSE)
  xy <- as.matrix(co[, c("x", "y")])
  if (any(!is.finite(xy)) || min(radii) < sqrt(.Machine$double.xmin) ||
      max(radii) > sqrt(.Machine$double.xmax) / 4 ||
      max(abs(xy)) > sqrt(.Machine$double.xmax) / 4 || max(abs(xy)) / max(radii) > 1e12) {
    stop("unsupported coordinate/radius magnitude; recenter or rescale", call. = FALSE)
  }
  groups <- sort(unique(labels))
  nq <- length(cells); ng <- length(groups); nr <- length(radii)
  if (as.double(nq) * ng * nr > .Machine$integer.max) {
    stop("output too large; select fewer cells, groups or radii", call. = FALSE)
  }
  target <- match(cells, co$cell_id)
  counts <- array(0L, c(nq, ng, nr))
  for (sample in unique(samples[target])) {
    context <- which(samples == sample)
    rows <- which(samples[target] == sample)
    counts[rows, , ] <- spatial_neighborhood_profile_cpp(
      xy[context, , drop = FALSE], match(labels[context], groups),
      match(target[rows], context), radii, ng
    )
  }
  if (cumulative && nr > 1L) {
    for (b in 2:nr) counts[, , b] <- counts[, , b] + counts[, , b - 1L]
  }
  out <- expand.grid(cell_id = cells, group = groups, radius = radii,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$sample <- rep(samples[target], ng * nr)
  out$lower <- rep(if (cumulative) rep(0, nr) else c(0, head(radii, -1L)), each = nq * ng)
  out$count <- as.vector(counts)
  totals <- matrix(0, nq, nr)
  if (nq) for (b in seq_len(nr)) totals[, b] <- rowSums(matrix(counts[, , b], nq, ng))
  out$total <- as.vector(totals[rep(seq_len(nq), ng), , drop = FALSE])
  out$fraction <- out$count / out$total
  out$fraction[out$total == 0] <- NA_real_
  out <- out[, c("cell_id", "sample", "group", "lower", "radius", "count", "total", "fraction")]
  attr(out, "source") <- resolved$source
  attr(out, "parameters") <- list(group.by = group.by, sample.by = sample.by,
    radii = radii, cumulative = cumulative, context_n = nrow(co), backend = "native")
  spatial_run_receipt(
    done = "Spatial neighborhood profile completed",
    scope = sprintf("%s target cells; %s context cells; %s samples; %s distances (raw coordinate units).",
                    nq, nrow(co), length(unique(samples)), nr),
    inspect = "Returned data.frame: count, total and fraction; input object unchanged.",
    verbose = verbose
  )
  out
}
