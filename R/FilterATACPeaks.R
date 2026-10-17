#' Filter an ATAC peak universe with explicit criteria
#'
#' @md
#' @description Select peaks by raw-count support, interval width, chromosome,
#' and optional blacklist overlap, then create a clean counts-only assay.
#'
#' @param object A Seurat object with a validated `ChromatinAssay`.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param new.assay An unused assay name, default `"filtered_peaks"`.
#' @param min.cells Minimum number of cells with a nonzero raw count, inclusive.
#' @param min.counts Minimum sum of raw counts across cells, inclusive.
#' @param min.width,max.width Inclusive integer interval-width limits in bases.
#'   The default retains every positive width; `max.width` may be `Inf`.
#' @param seqlevels Optional unique chromosome names to retain. Names must
#'   occur in the source reference; chromosome styles are never guessed.
#' @param blacklist Optional `GRanges` of excluded 1-based closed intervals.
#'   Any overlap excludes a peak regardless of strand. The blacklist must
#'   declare the same genome label as the source assay. Known chromosome lengths
#'   must agree and interval ends must not exceed them. Chromosomes absent from
#'   the source assay are retained as supplied; names are never converted.
#' @param verbose Whether to print progress messages.
#'
#' @return The Seurat object with the filtered `ChromatinAssay`. All cells and
#'   retained peak order/counts are preserved. The new assay retains reference,
#'   gene annotation and fragment references, but not normalized/scaled data,
#'   variable features, motifs, links or downstream results. Native initialization
#'   puts raw counts in the data layer. Source assays, active assay and existing
#'   reductions are unchanged. `Misc(object[[new.assay]], "scop_peak_filter")`
#'   stores one row per source peak with each criterion and the final keep flag.
#'   At least two peaks and two cells must remain for native assay construction.
#' @seealso [CheckATACData()], [AnnotateATACPeaks()]
#' @export
#' @examples
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
#'   dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
#' )
#' ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
#' object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
#' filtered <- FilterATACPeaks(object, min.cells = 1, verbose = FALSE)
#' rownames(filtered[["filtered_peaks"]])
FilterATACPeaks <- function(
  object, assay = NULL, new.assay = "filtered_peaks",
  min.cells = 0L, min.counts = 0, min.width = 1L, max.width = Inf,
  seqlevels = NULL, blacklist = NULL, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message(paste(report$errors, collapse = " "), message_type = "error")
  assay <- report$assay
  atac_new_assay_name(object, new.assay)
  for (value in list(min.cells, min.counts)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value %% 1 != 0) {
      log_message("min.cells and min.counts must be non-negative integers.", message_type = "error")
    }
  }
  if (!is.numeric(min.width) || length(min.width) != 1L || !is.finite(min.width) ||
    min.width < 1 || min.width %% 1 != 0 || !is.numeric(max.width) ||
    length(max.width) != 1L || is.na(max.width) || max.width < min.width ||
    (is.finite(max.width) && max.width %% 1 != 0)) {
    log_message("Width limits must be positive ordered integers; max.width may be Inf.", message_type = "error")
  }
  chromatin <- object[[assay]]
  ranges <- Signac::granges(chromatin)
  chromosomes <- as.character(GenomicRanges::seqnames(ranges))
  if (!is.null(seqlevels) && (!atac_unique_ids(seqlevels) ||
    !all(seqlevels %in% atac_backend_fun("GenomeInfoDb", "seqlevels")(chromatin)))) {
    log_message("seqlevels must contain unique chromosome names from the source reference.", message_type = "error")
  }
  blacklisted <- rep(FALSE, length(ranges))
  if (!is.null(blacklist)) {
    atac_peak_reference(blacklist, chromatin, report$genome, "Blacklist")
    blacklisted <- !is.na(GenomicRanges::findOverlaps(ranges, blacklist,
      select = "first", ignore.strand = TRUE
    ))
  }
  counts <- SeuratObject::GetAssayData(chromatin, layer = "counts")
  detection <- Matrix::rowSums(counts > 0)
  totals <- Matrix::rowSums(counts)
  widths <- GenomicRanges::width(ranges)
  inventory <- data.frame(
    peak = rownames(counts), n_cells = detection,
    count = totals, width = widths, pass_cells = detection >= min.cells,
    pass_counts = totals >= min.counts,
    pass_width = widths >= min.width & widths <= max.width,
    pass_seqlevels = if (is.null(seqlevels)) TRUE else chromosomes %in% seqlevels,
    pass_blacklist = !blacklisted, row.names = NULL
  )
  inventory$keep <- with(inventory, pass_cells & pass_counts & pass_width & pass_seqlevels & pass_blacklist)
  keep <- inventory$keep
  if (sum(keep) < 2L || ncol(counts) < 2L) {
    log_message("Filtering must retain at least two peaks and two cells for native assay construction.", message_type = "error")
  }
  log_message("Retaining {sum(keep)} of {length(keep)} ATAC peaks.", message_type = "info", verbose = verbose)
  filtered <- Signac::CreateChromatinAssay(
    counts = counts[keep, , drop = FALSE],
    ranges = ranges[keep], genome = atac_backend_fun("GenomeInfoDb", "seqinfo")(chromatin),
    annotation = Signac::Annotation(chromatin), fragments = Signac::Fragments(chromatin),
    min.cells = 0, min.features = 0, verbose = verbose
  )
  SeuratObject::Misc(filtered, "scop_peak_filter") <- list(
    source_assay = assay,
    cells = colnames(counts), genome = report$genome, peaks = inventory,
    parameters = list(
      min.cells = min.cells, min.counts = min.counts,
      min.width = min.width, max.width = max.width, seqlevels = seqlevels
    ),
    blacklist = blacklist, backend = c(Signac = as.character(utils::packageVersion("Signac")))
  )
  object[[new.assay]] <- filtered
  object
}


atac_peak_reference <- function(regions, chromatin, genome, label) {
  if (!inherits(regions, "GRanges")) log_message(paste0(label, " must be a GRanges object."), message_type = "error")
  if (any(GenomicRanges::start(regions) < 1L |
    GenomicRanges::end(regions) < GenomicRanges::start(regions))) {
    log_message(paste0(label, " requires positive 1-based closed intervals."), message_type = "error")
  }
  get_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")
  source_lengths <- get_lengths(chromatin)
  region_lengths <- get_lengths(regions)
  shared <- intersect(names(source_lengths), names(region_lengths))
  if (any(source_lengths[shared] != region_lengths[shared], na.rm = TRUE)) {
    log_message(paste0(label, " and assay chromosome lengths conflict."), message_type = "error")
  }
  chromosomes <- as.character(GenomicRanges::seqnames(regions))
  if (any(GenomicRanges::end(regions) > source_lengths[chromosomes] |
    GenomicRanges::end(regions) > region_lengths[chromosomes], na.rm = TRUE)) {
    log_message(paste0(label, " intervals exceed known chromosome lengths."), message_type = "error")
  }
  labels <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(regions)))
  if (length(genome) != 1L || !identical(labels, genome)) {
    log_message(paste0(label, " and assay genome labels must declare the same verified reference."), message_type = "error")
  }
}
