#' Bias-corrected ATAC footprinting
#'
#' @description Compute Signac insertion bias and footprints with an
#' explicit sequence reference. No reference is downloaded. Existing cached
#' footprints are never overwritten.
#' @param object A Seurat object containing a ChromatinAssay.
#' @param genome A loaded BSgenome or indexed Rsamtools FaFile.
#' @param genome.id Verified reference label matching the assay.
#' @param bias.region One explicit chromosome interval string, for example
#' `chr1-4-1000000`. It must fit the reference. Bias estimation requires
#' adequate coverage of all hexamers and can fail on small datasets.
#' @param motifs Unique motif IDs present in the selected assay's positions.
#' @param assay Source assay; defaults to the active assay.
#' @param upstream,downstream Nonnegative flanking widths (at least 50 each).
#' @param in.peaks Restrict motif instances to peaks.
#' @param verbose Print progress.
#' @return A Seurat object with positionEnrichment and a
#' `scop_footprints` receipt in assay Misc. Input counts, cells and default assay
#' are preserved. At least two assay cells are required. A receipt is a snapshot; recalculate after changing cells,
#' peaks, fragments or reference.
#' @details Requires indexed local fragment files with complete cell/barcode
#' mappings. Bias is recalculated using the explicit region, never silently
#' reused from another reference. Estimation errors are preserved.
#' Footprints summarize Tn5 insertion patterns, not proven TF binding or
#' causality. Factors sharing a motif cannot be distinguished by this analysis.
#' @seealso [ATACFootprintPlot()]
#' @export
RunATACFootprint <- function(object, genome, genome.id, bias.region, motifs,
                             assay = NULL, upstream = 250L, downstream = 250L,
                             in.peaks = FALSE, verbose = TRUE) {
  report <- CheckATACData(object, assay = assay, genome = genome.id)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  if (!report$fragments_ready) stop("Footprinting requires accessible indexed fragments with complete cell mappings.", call. = FALSE)
  assay <- report$assay
  # Native fragment lookups cache fastmatch attributes by reference. Isolate
  # the selected assay so even those caches cannot mutate the caller.
  chromatin <- unserialize(serialize(object[[assay]], NULL))
  if (ncol(chromatin) < 2L) stop("Footprinting requires at least two assay cells.", call. = FALSE)
  atac_fp_reference(genome, genome.id, chromatin)
  for (value in list(upstream, downstream)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 50 || value %% 1 != 0) {
      stop("upstream and downstream must be integer widths >= 50 for flank normalization.", call. = FALSE)
    }
  }
  if (!is.logical(in.peaks) || length(in.peaks) != 1L || is.na(in.peaks)) stop("in.peaks must be TRUE or FALSE.", call. = FALSE)
  if (!is.character(bias.region) || length(bias.region) != 1L || is.na(bias.region) ||
    !grepl("^[^-]+-[0-9]+-[0-9]+$", bias.region)) {
    stop("bias.region must be one explicit chromosome-start-end interval string.", call. = FALSE)
  }
  interval_fields <- strsplit(bias.region, "-", fixed = TRUE)[[1L]]
  interval <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    chr = interval_fields[[1L]], start = interval_fields[[2L]], end = interval_fields[[3L]]
  ))
  chr <- as.character(GenomicRanges::seqnames(interval))
  lengths <- GenomeInfoDb::seqlengths(genome)
  if (!chr %in% names(lengths) || is.na(lengths[[chr]]) ||
    GenomicRanges::start(interval) < 4L || GenomicRanges::end(interval) > lengths[[chr]] - 3L ||
    GenomicRanges::width(interval) < 6L) {
    stop("bias.region must fit the reference with three-base edge margins.", call. = FALSE)
  }
  positions <- Signac::GetMotifData(chromatin, slot = "positions")
  if (!atac_unique_ids(motifs) || !all(motifs %in% names(positions))) stop("motifs must be unique IDs with stored positions in the selected assay.", call. = FALSE)
  existing <- SeuratObject::GetAssayData(chromatin, layer = "positionEnrichment")
  if (any(motifs %in% names(existing))) stop("Requested footprint keys already exist; use a fresh object or explicitly remove stale results first.", call. = FALSE)
  for (id in motifs) {
    regions <- positions[[id]]
    if (!length(regions) || length(unique(GenomicRanges::width(regions))) != 1L) stop("Every motif requires nonempty equal-width positions.", call. = FALSE)
    labels <- unique(as.character(GenomeInfoDb::genome(regions)))
    labels <- labels[!is.na(labels) & nzchar(labels)]
    known <- GenomeInfoDb::seqlengths(regions)
    shared <- intersect(names(known), names(lengths))
    if ((length(labels) && any(labels != genome.id)) ||
      any(known[shared] != lengths[shared], na.rm = TRUE)) {
      stop("Motif position reference labels or chromosome lengths conflict with genome.id/reference.", call. = FALSE)
    }
    if (in.peaks) regions <- IRanges::subsetByOverlaps(regions, Signac::granges(chromatin))
    if (!length(regions)) stop("No motif positions remain after in.peaks filtering.", call. = FALSE)
    rchr <- as.character(GenomicRanges::seqnames(regions))
    minus <- as.character(GenomicRanges::strand(regions)) == "-"
    left <- ifelse(minus, downstream, upstream)
    right <- ifelse(minus, upstream, downstream)
    if (any(!rchr %in% names(lengths)) || anyNA(lengths[rchr]) ||
      any(GenomicRanges::start(regions) - left < 4L) ||
      any(GenomicRanges::end(regions) + right > lengths[rchr] - 3L)) {
      stop("Motif windows must fit the declared reference, including bias margins; filter positions explicitly.", call. = FALSE)
    }
  }
  verified.fragments <- atac_indexed_fragments(chromatin, colnames(chromatin), verbose)
  # Native fragment readers can retain non-assay mappings; work on a local assay.
  fragments <- Signac::Fragments(chromatin)
  chromatin@fragments <- Filter(function(f) any(names(Signac::GetFragmentData(f, "cells")) %in% colnames(chromatin)), fragments)
  chromatin <- Signac::InsertionBias(chromatin, genome = genome, region = bias.region, verbose = verbose)
  bias <- SeuratObject::GetAssayData(chromatin, layer = "bias")
  if (!length(bias) || any(!is.finite(bias))) stop("Insertion bias is non-finite; choose a sufficiently covered reference region.", call. = FALSE)
  chromatin <- Signac::Footprint(chromatin, genome = genome,
    motif.name = motifs, key = motifs, upstream = upstream, downstream = downstream,
    compute.expected = TRUE, in.peaks = in.peaks, verbose = verbose)
  chromatin@fragments <- fragments
  receipt <- SeuratObject::Misc(chromatin, "scop_footprints")
  if (is.null(receipt)) receipt <- list()
  snapshot <- list(
    cells = as.vector(colnames(chromatin)), peaks = as.vector(rownames(chromatin)),
    genome = genome.id, bias.region = bias.region, upstream = upstream, downstream = downstream,
    compute.expected = TRUE, in.peaks = in.peaks,
    fragments = lapply(verified.fragments, function(x) x[c("path", "index", "hash", "map")]),
    backend = as.character(utils::packageVersion("Signac"))
  )
  for (id in motifs) receipt[[id]] <- snapshot
  SeuratObject::Misc(chromatin, "scop_footprints") <- receipt
  object[[assay]] <- chromatin
  object
}

atac_fp_reference <- function(genome, genome.id, chromatin) {
  if (!inherits(genome, "BSgenome") && !inherits(genome, "FaFile")) stop("genome must be a loaded BSgenome or indexed FaFile.", call. = FALSE)
  labels <- unique(as.character(GenomeInfoDb::genome(genome)))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) && any(labels != genome.id)) stop("Sequence-reference labels must match genome.id.", call. = FALSE)
  sizes <- GenomeInfoDb::seqlengths(genome)
  known <- GenomeInfoDb::seqlengths(chromatin)
  shared <- intersect(names(sizes), names(known))
  if (any(sizes[shared] != known[shared], na.rm = TRUE)) stop("Reference chromosome lengths conflict with assay.", call. = FALSE)
  ranges <- Signac::granges(chromatin)
  chr <- as.character(GenomicRanges::seqnames(ranges))
  if (any(!chr %in% names(sizes)) || anyNA(sizes[chr]) || any(GenomicRanges::end(ranges) > sizes[chr])) stop("Source peaks must fit the loaded sequence reference.", call. = FALSE)
  invisible(TRUE)
}

#' Plot bias-corrected ATAC footprints
#' @param object A Seurat object returned by RunATACFootprint.
#' @param motifs Unique footprint keys.
#' @param assay Chromatin assay.
#' @param group.by Metadata column, or NULL for identities.
#' @param normalization Signac normalization: `subtract`, `divide`, or NULL.
#' @param show.expected Show expected insertion profile.
#' @param ... Arguments passed to other functions.
#' @return A ggplot/patchwork object from Signac.
#' @export
ATACFootprintPlot <- function(object, motifs, assay = NULL, group.by = NULL,
                              normalization = "subtract", show.expected = TRUE, ...) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  assay <- report$assay
  if (!atac_unique_ids(motifs)) stop("motifs must be unique footprint keys.", call. = FALSE)
  receipt <- SeuratObject::Misc(object[[assay]], "scop_footprints")
  for (id in motifs) {
    r <- receipt[[id]]
    if (is.null(r) || !identical(r$cells, as.vector(colnames(object[[assay]]))) ||
      !identical(r$peaks, as.vector(rownames(object[[assay]]))) ||
      !identical(r$genome, report$genome)) {
      stop("Footprint receipt is missing or stale; recompute after subsetting or renaming.", call. = FALSE)
    }
  }
  if (!is.null(normalization) && (!is.character(normalization) || length(normalization) != 1L ||
    !normalization %in% c("subtract", "divide"))) {
    stop("normalization must be subtract, divide or NULL.", call. = FALSE)
  }
  cells <- colnames(object[[assay]])
  if (!is.null(group.by) && (!is.character(group.by) || length(group.by) != 1L ||
    !group.by %in% colnames(object[[]]) || anyNA(object[[]][cells, group.by]))) {
    stop("group.by must name complete metadata for assay cells.", call. = FALSE)
  }
  local <- object
  if (!identical(as.vector(colnames(object)), as.vector(cells))) {
    saved <- object[[assay]]@positionEnrichment
    local[[assay]]@positionEnrichment <- list()
    local <- subset(local, cells = cells)
    local[[assay]]@positionEnrichment <- saved
  }
  Signac::PlotFootprint(local, features = motifs, assay = assay,
    group.by = group.by, normalization = normalization, show.expected = show.expected, ...)
}
