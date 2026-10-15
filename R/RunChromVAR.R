#' Compute native chromVAR deviations from an explicit sequence reference
#'
#' @md
#' @description Compute native GC-corrected motif deviations and store their
#'   Z-scores in a new assay without changing source counts.
#' @param object A Seurat object with aligned peak counts and native Signac motifs.
#'   Motif IDs must survive native Seurat assay creation unchanged (for example,
#'   underscore-containing IDs are rejected rather than silently renamed).
#' @param genome A loaded BSgenome, DNAStringSet or indexed Rsamtools FaFile.
#'   Character package names are not accepted and no reference is downloaded.
#' @param genome.id Explicit non-empty reference label, identical to the source
#'   assay and any known label on the loaded reference. Unlabeled reference
#'   provenance is the caller's declaration, not an assembly inference. Known
#'   chromosome lengths must agree and source peaks must fit the reference.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param new.assay Unused name for the chromVAR Z-score assay.
#' @param background.peaks Optional integer matrix with one row per source peak,
#'   in exactly the same order and with matching row names. Values are 1-based
#'   source-peak indices. With NULL, use native GC/accessibility matching. For
#'   explicit backgrounds, the caller is responsible for appropriate matching;
#'   this function validates the current peak-index universe only.
#' @param niterations,w,bs Parameters passed unchanged to getBackgroundPeaks.
#' @param seed Seed for native background sampling; caller RNG is restored.
#'
#' @details
#' This calls chromVAR directly and retains its native deviations object and actual background
#' indices in assay Misc, while the assay data layer contains native Z-scores.
#' No peaks or cells are silently dropped, and non-finite native scores are not
#' replaced with zero. Every source peak and cell must have positive total counts.
#' At least two source-assay cells are required. Z-scores quantify motif-associated
#' accessibility deviations, not TF expression, occupancy, or causality. Explicit
#' stored sparse zeros are removed locally before native motif indexing.
#' Native computeDeviations uses the caller's BiocParallel configuration; this
#' function does not register a global parallel backend.
#'
#' @return The Seurat object with a new assay. Its data layer contains Z-scores;
#'   Misc(new.assay, "scop_chromvar") contains native deviations, GC bias,
#'   actual background peaks, ordered source identifiers and parameters. This
#'   receipt is a calculation-time snapshot: later native subset/RenameCells
#'   operations do not synchronize its nested experiments or identifiers. Source
#'   metadata and default assay are unchanged; no count QC is created from scores.
#' @examples
#' if (all(vapply(c("chromVAR", "Biostrings", "BSgenome", "BiocParallel"),
#'   requireNamespace, logical(1),
#'   quietly = TRUE
#' ))) {
#'   (function() {
#'     old <- BiocParallel::bpparam()
#'     on.exit(BiocParallel::register(old, default = TRUE))
#'     BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
#'     n <- 200L
#'     counts <- withr::with_seed(321, Matrix::Matrix(matrix(
#'       1 + stats::rpois(n * 4, 4), n, 4,
#'       dimnames = list(paste0("p", seq_len(n)), paste0("cell", 1:4))
#'     ), sparse = TRUE))
#'     ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'       seqnames = "chr1", start = (seq_len(n) - 1L) * 100L + 1L, end = seq_len(n) * 100L
#'     ))
#'     atac <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
#'     gc <- 10L + seq_len(n) %% 81L
#'     reference <- Biostrings::DNAStringSet(c(chr1 = paste0(vapply(
#'       gc,
#'       function(k) paste0(strrep("G", k), strrep("A", 100L - k)), ""
#'     ), collapse = "")))
#'     # Toy native annotations illustrate the interface, not a biological TF claim.
#'     incidence <- Matrix::Matrix(cbind(Grich = gc > 50, Arich = gc < 30), sparse = TRUE)
#'     rownames(incidence) <- rownames(atac)
#'     pwm <- stats::setNames(rep(list(matrix(1, 4, 3,
#'       dimnames = list(c("A", "C", "G", "T"), NULL)
#'     )), 2), colnames(incidence))
#'     Signac::Motifs(atac[["peaks"]]) <- Signac::CreateMotifObject(data = incidence, pwm = pwm)
#'     scored <- RunChromVAR(atac, reference, genome.id = "toy", niterations = 10, w = 0.2, bs = 10)
#'     head(SeuratObject::LayerData(scored[["chromvar"]], layer = "data"))
#'   })()
#' }
#' @export
RunChromVAR <- function(
  object, genome, genome.id, assay = NULL,
  new.assay = "chromvar", background.peaks = NULL,
  niterations = 50L, w = 0.1, bs = 50L, seed = 1L
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  atac_new_assay_name(object, new.assay)
  if (!is.character(genome.id) || length(genome.id) != 1L || is.na(genome.id) ||
    !nzchar(genome.id) || !identical(genome.id, report$genome)) {
    log_message("genome.id must declare the verified source-assay reference label.", message_type = "error")
  }
  if (!any(vapply(
    c("BSgenome", "DNAStringSet", "FaFile"),
    function(x) inherits(genome, x), logical(1)
  ))) {
    log_message("genome must be a loaded BSgenome, DNAStringSet or indexed FaFile.", message_type = "error")
  }
  add_gc <- atac_backend_fun("chromVAR", "addGCBias")
  get_background <- atac_backend_fun("chromVAR", "getBackgroundPeaks")
  compute <- atac_backend_fun("chromVAR", "computeDeviations")
  with_seed <- atac_backend_fun("withr", "with_seed")
  atac_backend_fun("withr", "local_preserve_seed")()
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
    seed < 0 || seed > .Machine$integer.max || seed %% 1 != 0) {
    log_message("seed must be one nonnegative integer within the R integer range.", message_type = "error")
  }
  counts <- SeuratObject::LayerData(object[[assay]], layer = "counts")
  ranges <- Signac::granges(object[[assay]])
  incidence <- Signac::GetMotifData(object, assay = assay, slot = "data")
  ids <- rownames(counts)
  cells <- as.vector(colnames(counts))
  # CheckATACData already verifies coordinate order. GRanges display names
  # may retain source IDs, so give the native experiment its canonical IDs.
  names(ranges) <- ids
  if (!inherits(incidence, "sparseMatrix") || !identical(as.vector(rownames(incidence)), as.vector(ids)) ||
    !atac_unique_ids(colnames(incidence)) ||
    (!inherits(incidence, "nMatrix") &&
      (!all(is.finite(incidence@x)) || !all(incidence@x %in% c(0, 1))))) {
    log_message("Motif incidence must be binary, sparse and aligned exactly with source peaks.", message_type = "error")
  }
  # chromVAR indexes sparse structure, so explicit stored zeros must be removed.
  incidence <- Matrix::drop0(incidence)
  if (ncol(counts) < 2L) {
    log_message("Native chromVAR requires at least two source-assay cells.", message_type = "error")
  }
  get_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")
  reference_labels <- unique(atac_backend_fun("GenomeInfoDb", "genome")(genome))
  reference_labels <- reference_labels[!is.na(reference_labels) & nzchar(reference_labels)]
  if (length(reference_labels) && any(reference_labels != genome.id)) {
    log_message("Sequence-reference genome labels must agree with genome.id and the source assay.", message_type = "error")
  }
  chromosomes <- as.character(GenomicRanges::seqnames(ranges))
  if (!all(chromosomes %in% atac_backend_fun("GenomeInfoDb", "seqlevels")(genome))) {
    log_message("The sequence reference must contain every source-peak chromosome.", message_type = "error")
  }
  source_lengths <- get_lengths(object[[assay]])
  reference_lengths <- get_lengths(genome)
  if (any(source_lengths != reference_lengths[names(source_lengths)], na.rm = TRUE)) {
    log_message("Known source and sequence-reference chromosome lengths must agree.", message_type = "error")
  }
  lengths <- reference_lengths[chromosomes]
  if (anyNA(lengths) || any(GenomicRanges::end(ranges) > lengths)) {
    log_message("Source peaks must fit within known sequence-reference lengths.", message_type = "error")
  }
  if (any(Matrix::rowSums(counts) <= 0) || any(Matrix::colSums(counts) <= 0)) {
    log_message("chromVAR requires positive total counts for every peak and cell; subset explicitly before analysis.", message_type = "error")
  }
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowRanges = ranges
  )
  se <- add_gc(se, genome = genome)
  bias <- SummarizedExperiment::rowRanges(se)$bias
  if (length(bias) != nrow(counts) || any(!is.finite(bias))) {
    log_message("Native GC computation produced missing/non-finite bias; verify reference sequences.", message_type = "error")
  }
  sampled <- is.null(background.peaks)
  if (sampled) {
    if (!is.numeric(niterations) || length(niterations) != 1L ||
      !is.finite(niterations) || niterations < 2 || niterations %% 1 != 0 ||
      !is.numeric(bs) || length(bs) != 1L || !is.finite(bs) || bs < 2 || bs %% 1 != 0 ||
      !is.numeric(w) || length(w) != 1L || !is.finite(w) || w <= 0) {
      log_message("niterations and bs must be integers >= 2, and w must be positive.", message_type = "error")
    }
    background.peaks <- with_seed(seed, get_background(se,
      niterations = niterations, w = w, bs = bs
    ))
    rownames(background.peaks) <- ids
  }
  if (!is.matrix(background.peaks) || !is.numeric(background.peaks) ||
    nrow(background.peaks) != length(ids) || ncol(background.peaks) < 2L ||
    !identical(as.vector(rownames(background.peaks)), as.vector(ids)) || any(!is.finite(background.peaks)) ||
    any(background.peaks %% 1 != 0) || any(background.peaks < 1) ||
    any(background.peaks > length(ids))) {
    log_message("background.peaks must contain valid 1-based indices with exact source-peak row names/order and >= 2 iterations.", message_type = "error")
  }
  deviations <- compute(se, annotations = incidence, background_peaks = background.peaks)
  scores <- chromVAR::deviationScores(deviations)
  if (!identical(as.vector(rownames(scores)), as.vector(colnames(incidence))) ||
    !identical(as.vector(colnames(scores)), cells)) {
    log_message("Native chromVAR output identifiers do not match motif and cell inputs.", message_type = "error")
  }
  activity <- SeuratObject::CreateAssayObject(data = scores)
  # Native Seurat writeback must not calculate count QC for a data-only score assay.
  activity@misc$calcN <- FALSE
  if (!identical(as.vector(rownames(activity)), as.vector(rownames(scores)))) {
    log_message("Motif IDs changed during native assay creation.", message_type = "error")
  }
  SeuratObject::Misc(activity, "scop_chromvar") <- list(
    source_assay = assay, genome = genome.id, sequence_class = class(genome),
    peaks = ids, cells = cells, motifs = colnames(incidence),
    gc_bias = stats::setNames(bias, ids), background_peaks = background.peaks,
    deviations = deviations, background_source = if (sampled) "chromVAR::getBackgroundPeaks" else "explicit",
    parameters = if (sampled) list(niterations = niterations, w = w, bs = bs, seed = seed) else NULL,
    nonfinite_scores = sum(!is.finite(scores)),
    backend = c(
      chromVAR = as.character(utils::packageVersion("chromVAR")),
      Signac = as.character(utils::packageVersion("Signac"))
    )
  )
  object[[new.assay]] <- activity
  object
}
