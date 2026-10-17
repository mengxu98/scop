#' Annotate ATAC peaks with nearest genomic features
#'
#' @md
#' @description Attach verified gene annotation and record native nearest-feature
#' distances, retaining native nearest-feature hits and unmatched peaks.
#'
#' @param object A Seurat object with a validated `ChromatinAssay`.
#' @param annotation A non-empty Signac-compatible `GRanges` of 1-based closed
#'   gene/transcript intervals with gene_id, gene_name, gene_biotype, type and
#'   tx_id (or transcript_id) metadata. It must declare the source assay's
#'   genome label. Known chromosome lengths must agree and interval ends must
#'   not exceed them. Extra chromosomes are retained without name conversion.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param verbose Whether to print progress messages.
#'
#' @details Nearest distances are computed by
#'   `GenomicRanges::distanceToNearest(select = "all", ignore.strand = TRUE)`.
#'   Distance counts the bases between intervals; overlaps and adjacent intervals
#'   have distance zero. These are distances to annotation intervals, not to TSSs,
#'   promoter classifications, regulatory links, or evidence of gene regulation.
#'   Different transcripts can produce several native hits. `select = "all"`
#'   does not guarantee every mathematical tie on both sides of a peak: native
#'   GenomicRanges tie selection is preserved. No genome
#'   download or chromosome-style conversion is performed. Existing annotation
#'   must be identical after native Signac normalization (including transcript_id
#'   to tx_id); existing results are not replaced. Genome labels and known
#'   lengths are checked, not reference sequence contents. Nonmatching chromosome
#'   names remain unmatched rather than being translated.
#'
#' @return The Seurat object with native Signac gene annotation and
#'   `Misc(object[[assay]], "scop_peak_annotation")`. The receipt contains
#'   `nearest`, an `S4Vectors::DataFrame` with peak ID, annotation row index,
#'   distance and prefixed annotation metadata; `unannotated_peaks` lists peaks
#'   without a hit. Rows follow peak order and annotation index within ties.
#'   Counts, cells, other assays and reductions are unchanged.
#' @seealso [FilterATACPeaks()], [CheckATACData()]
#' @export
#' @examples
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
#'   dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
#' )
#' ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
#' object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
#' annotation <- GenomicRanges::GRanges("chr1:1-80")
#' GenomeInfoDb::genome(annotation) <- "hg38"
#' S4Vectors::mcols(annotation) <- S4Vectors::DataFrame(
#'   gene_id = "gene1", gene_name = "GENE1", gene_biotype = "protein_coding",
#'   type = "gene", tx_id = "tx1"
#' )
#' object <- AnnotateATACPeaks(object, annotation, verbose = FALSE)
#' SeuratObject::Misc(object[["peaks"]], "scop_peak_annotation")$nearest
AnnotateATACPeaks <- function(object, annotation, assay = NULL, verbose = TRUE) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message(paste(report$errors, collapse = " "), message_type = "error")
  assay <- report$assay
  chromatin <- object[[assay]]
  atac_peak_reference(annotation, chromatin, report$genome, "Gene annotation")
  if (!length(annotation) || !all(c("gene_id", "gene_name", "gene_biotype", "type") %in%
    colnames(S4Vectors::mcols(annotation))) ||
    !any(c("tx_id", "transcript_id") %in% colnames(S4Vectors::mcols(annotation)))) {
    log_message("annotation requires gene_id, gene_name, gene_biotype, type and tx_id or transcript_id metadata.", message_type = "error")
  }
  existing <- Signac::Annotation(chromatin)
  Signac::Annotation(chromatin) <- annotation
  annotation <- Signac::Annotation(chromatin)
  if (length(existing) && !identical(existing, annotation)) {
    log_message("The assay already contains different gene annotation; use a separate assay.", message_type = "error")
  }
  if (!is.null(SeuratObject::Misc(chromatin, "scop_peak_annotation"))) {
    log_message("Peak annotation results already exist on this assay; use a separate assay.", message_type = "error")
  }
  ranges <- Signac::granges(chromatin)
  hits <- GenomicRanges::distanceToNearest(ranges, annotation, select = "all", ignore.strand = TRUE)
  query <- S4Vectors::queryHits(hits)
  subject <- S4Vectors::subjectHits(hits)
  order <- order(query, subject)
  query <- query[order]
  subject <- subject[order]
  metadata <- S4Vectors::mcols(annotation)[subject, , drop = FALSE]
  colnames(metadata) <- paste0("annotation_", colnames(metadata))
  nearest <- cbind(S4Vectors::DataFrame(
    peak = rownames(chromatin)[query],
    annotation_index = subject, distance = S4Vectors::mcols(hits)$distance[order]
  ), metadata)
  SeuratObject::Misc(chromatin, "scop_peak_annotation") <- list(
    nearest = nearest,
    unannotated_peaks = rownames(chromatin)[setdiff(seq_along(ranges), query)],
    genome = report$genome, peaks = rownames(chromatin),
    parameters = list(select = "all", ignore.strand = TRUE),
    backend = c(
      GenomicRanges = as.character(utils::packageVersion("GenomicRanges")),
      Signac = as.character(utils::packageVersion("Signac"))
    )
  )
  log_message("Annotated {length(unique(query))} of {length(ranges)} peaks with nearest genomic features.",
    message_type = "info", verbose = verbose
  )
  object[[assay]] <- chromatin
  object
}
