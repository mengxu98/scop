#' Add native motif matches and sequence statistics
#'
#' @md
#' @description Match supplied position frequency/weight matrices to the exact
#' source peak universe through Signac and motifmatchr.
#' @param object A Seurat object containing a ChromatinAssay, or a `GRanges`
#'   of final peaks. The ranges-only path needs no counts or cells and can run
#'   before recount. Coordinates use 1-based closed intervals.
#' @param genome Loaded BSgenome, DNAStringSet or indexed FaFile. No download or
#'   installation is performed.
#' @param pfm A named TFBSTools PFMatrixList or PWMatrixList with unique motif IDs.
#' @param genome.id Explicit non-empty genome label matching the source assay
#'   and any known label on the loaded reference. For an unlabeled reference,
#'   this is the caller's provenance declaration, not an assembly inference.
#'   Known chromosome lengths must agree; source peaks must fit the reference.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param replace Whether to replace existing motif matches. Defaults to FALSE.
#' @param verbose Whether to display native progress messages.
#' @return The object with native Signac Motifs and RegionStats in the selected
#'   assay. Counts and cell order are unchanged. Misc(scop_atac_motifs) records
#'   ordered peaks, motif IDs, genome label and backend versions at addition time.
#'   This receipt is a snapshot, not a replacement for current native assay state.
#'   Existing sequence-statistic columns are refreshed by name; unrelated feature
#'   metadata are preserved. Explicit replacement emits the native receipt-overwrite warning.
#'   For `GRanges`, return a list with native `motifs`, aligned sequence
#'   `statistics`, the input `ranges`, and a calculation-time `receipt`.
#'   Motif row names are canonical Signac interval IDs. After recount over the
#'   same ranges, assign the native Motif with `Signac::Motifs(assay) <- motifs`
#'   and sequence statistics with `SeuratObject::AddMetaData()`; check exact
#'   peak order and reference identity before assignment. No placeholder counts
#'   are constructed and changing counts alone does not change sequence matches.
#' @seealso [RunMotifEnrichment()]
#' @examples
#' if (all(vapply(c("Biostrings", "BSgenome", "TFBSTools", "motifmatchr"),
#'   requireNamespace, logical(1),
#'   quietly = TRUE
#' ))) {
#'   reference <- Biostrings::DNAStringSet(c(chr1 = paste(rep("ACGT", 1200), collapse = "")))
#'   ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'     seqnames = "chr1", start = seq(1, by = 120, length.out = 40),
#'     end = seq(100, by = 120, length.out = 40)
#'   ))
#'   counts <- Matrix::Matrix(matrix(1, 40, 2,
#'     dimnames = list(paste0("p", seq_len(40)), c("A", "B"))
#'   ), sparse = TRUE)
#'   atac <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
#'   profile <- matrix(0L, 4, 8, dimnames = list(c("A", "C", "G", "T"), NULL))
#'   profile[cbind(rep(seq_len(4), 2), seq_len(8))] <- 100L
#'   pfm <- TFBSTools::PFMatrixList(M1 = TFBSTools::PFMatrix(
#'     ID = "M1", name = "toy_ACGT", profileMatrix = profile
#'   ))
#'   matched <- AddATACMotifs(atac, reference, pfm, genome.id = "toy", verbose = FALSE)
#'   SeuratObject::Misc(matched[["peaks"]], "scop_atac_motifs")
#' }
#' @export
AddATACMotifs <- function(
  object, genome, pfm, genome.id, assay = NULL,
  replace = FALSE, verbose = TRUE
) {
  ranges_only <- inherits(object, "GRanges")
  if (ranges_only) {
    problems <- atac_ranges_problems(object)
    if (length(problems)) log_message("{paste(problems, collapse = ' ')}", message_type = "error")
    ranges <- object
    if (!length(ranges) || !atac_unique_ids(atac_range_ids(ranges))) {
      log_message("Supply non-empty, unique peak intervals.", message_type = "error")
    }
    labels <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(ranges)))
    labels <- labels[!is.na(labels) & nzchar(labels)]
    if (length(labels) && !identical(labels, genome.id)) {
      log_message("Known peak genome labels must agree with genome.id.", message_type = "error")
    }
  } else {
    report <- CheckATACData(object, assay = assay)
    if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
    assay <- report$assay
    ranges <- Signac::granges(object[[assay]])
    if (!identical(genome.id, report$genome)) {
      log_message("Supply a genome.id matching the source assay.", message_type = "error")
    }
  }
  if (!atac_unique_ids(genome.id) || length(genome.id) != 1L ||
    !any(vapply(c("BSgenome", "DNAStringSet", "FaFile"), function(x) inherits(genome, x), logical(1)))) {
    log_message("Supply a loaded sequence reference and matching genome.id.", message_type = "error")
  }
  if (!inherits(pfm, "PFMatrixList") && !inherits(pfm, "PWMatrixList") || !atac_unique_ids(names(pfm))) {
    log_message("pfm must be a PFMatrixList or PWMatrixList with unique non-empty motif IDs.", message_type = "error")
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    log_message("replace must be TRUE or FALSE.", message_type = "error")
  }
  if (!ranges_only && !replace && !is.null(Signac::Motifs(object[[assay]]))) {
    log_message("The assay already contains motifs; set replace = TRUE explicitly.", message_type = "error")
  }
  chromosomes <- unique(as.character(GenomicRanges::seqnames(ranges)))
  reference_chromosomes <- atac_backend_fun("GenomeInfoDb", "seqlevels")(genome)
  if (!all(chromosomes %in% reference_chromosomes)) {
    log_message("The sequence reference must contain every source-peak chromosome.", message_type = "error")
  }
  get_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")
  reference_labels <- unique(atac_backend_fun("GenomeInfoDb", "genome")(genome))
  reference_labels <- reference_labels[!is.na(reference_labels) & nzchar(reference_labels)]
  if (length(reference_labels) && any(reference_labels != genome.id)) {
    log_message("Sequence-reference genome labels must agree with genome.id and the source assay.", message_type = "error")
  }
  source_lengths <- get_lengths(if (ranges_only) ranges else object[[assay]])
  shared_lengths <- get_lengths(genome)[names(source_lengths)]
  if (any(source_lengths != shared_lengths, na.rm = TRUE)) {
    log_message("Known source and sequence-reference chromosome lengths must agree.", message_type = "error")
  }
  reference_lengths <- get_lengths(genome)[as.character(GenomicRanges::seqnames(ranges))]
  if (anyNA(reference_lengths) || any(GenomicRanges::end(ranges) > reference_lengths)) {
    log_message("Source peaks must fit within the known sequence-reference lengths.", message_type = "error")
  }
  atac_backend_fun("motifmatchr", "matchMotifs")
  receipt <- list(
    genome = genome.id, sequence_class = class(genome),
    peaks = atac_range_ids(ranges), motifs = names(pfm), backend = c(
      Signac = as.character(utils::packageVersion("Signac")),
      motifmatchr = as.character(utils::packageVersion("motifmatchr"))
    )
  )
  if (ranges_only) {
    motifs <- Signac::AddMotifs(ranges, genome = genome, pfm = pfm, verbose = verbose)
    statistics <- Signac::RegionStats(ranges, genome = genome, verbose = verbose)
    rownames(statistics) <- receipt$peaks
    return(list(motifs = motifs, statistics = statistics, ranges = ranges, receipt = receipt))
  }
  result <- object
  matched <- Signac::AddMotifs(object[[assay]], genome = genome, pfm = pfm, verbose = verbose)
  statistics <- Signac::RegionStats(ranges, genome = genome, verbose = verbose)
  rownames(statistics) <- rownames(object[[assay]])
  result[[assay]] <- SeuratObject::AddMetaData(matched, metadata = statistics)
  incidence <- Signac::GetMotifData(result, assay = assay, slot = "data")
  if (!identical(as.vector(rownames(incidence)), as.vector(rownames(object[[assay]]))) ||
    !identical(as.vector(colnames(incidence)), as.vector(names(pfm))) ||
    !identical(colnames(result[[assay]]), colnames(object[[assay]]))) {
    log_message("Native motif output did not retain the source peak, motif or cell order.", message_type = "error")
  }
  SeuratObject::Misc(result[[assay]], "scop_atac_motifs") <- list(
    genome = genome.id, sequence_class = class(genome), peaks = rownames(incidence),
    motifs = colnames(incidence), backend = c(
      Signac = as.character(utils::packageVersion("Signac")),
      motifmatchr = as.character(utils::packageVersion("motifmatchr"))
    )
  )
  result
}
