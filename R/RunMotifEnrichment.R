#' Test motif enrichment within an explicit peak universe
#'
#' @md
#' @description
#' Use native Signac background matching and motif enrichment with an auditable
#' foreground and population. Sequence scanning remains available directly via
#' `Signac::AddMotifs()` and does not require this function or a DAR analysis.
#'
#' @param object A Seurat object containing a validated ChromatinAssay and an
#'   aligned native Signac Motif object.
#' @param features Unique foreground peak IDs, all present in `universe`.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param universe Eligible peak IDs. With NULL, use all source-assay peaks. For
#'   DAR-derived foregrounds, explicitly supply the full tested peak universe.
#' @param background Optional explicit population peak IDs, including every
#'   foreground peak and contained in `universe`. With NULL, select reference
#'   peaks using native `Signac::MatchRegionStats()`, then include the foreground
#'   in the tested population.
#' @param n.background Maximum number of non-foreground reference peaks selected
#'   when matching. If fewer are eligible, all eligible peaks are requested.
#' @param match.features Feature-metadata columns passed to native background
#'   matching. Compute verified sequence metadata with `Signac::RegionStats()`.
#' @param seed Random seed for native matching, with the caller's RNG restored.
#' @param p.adjust.method Multiple-testing method passed to `Signac::FindMotifs()`.
#' @param verbose Whether native backends print progress.
#'
#' @details
#' The hypergeometric population includes the foreground, so the foreground is
#' a subset of the tested population. `reference` in the receipt contains only
#' the matched non-foreground controls; `background` contains the actual native
#' test population. This differs from passing a numeric background directly to
#' Signac, whose matched list excludes the query. No enrichment formula is
#' reimplemented. All motifs are tested before any downstream display filtering.
#' P-values test one-sided over-enrichment; a fold enrichment below one is not
#' an under-enrichment test. For directional DAR analyses, select up/down
#' foregrounds separately. An enrichment result is a sequence association, not
#' evidence of TF binding. Native undefined fold enrichment for a zero-hit motif
#' is retained, together with its native p-value and multiple-testing result.
#'
#' @return A list with the complete native `results` table and a `receipt` with
#'   exact foreground, eligible universe, tested population, matching parameters,
#'   motif IDs and backend versions. The input object is unchanged.
#' @seealso [RunDEtest()]
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
#'   peaks <- rownames(matched[["peaks"]])
#'   enrichment <- RunMotifEnrichment(matched, peaks[1:10],
#'     universe = peaks, background = peaks, verbose = FALSE
#'   )
#'   enrichment$results
#' }
#' @export
RunMotifEnrichment <- function(
  object, features, assay = NULL, universe = NULL,
  background = NULL, n.background = 40000L,
  match.features = c("GC.percent", "sequence.length"), seed = 1L,
  p.adjust.method = "BH", verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  peaks <- rownames(object[[assay]])
  incidence <- Signac::GetMotifData(object, assay = assay, slot = "data")
  if (!inherits(incidence, "sparseMatrix") ||
    !identical(as.vector(rownames(incidence)), as.vector(peaks)) || !atac_unique_ids(colnames(incidence)) ||
    (!inherits(incidence, "nMatrix") &&
      (!all(is.finite(incidence@x)) || !all(incidence@x %in% c(0, 1))))) {
    log_message("Motif incidence must be binary, sparse and aligned exactly with the source peaks.", message_type = "error")
  }
  if (is.null(universe)) universe <- peaks
  if (!atac_unique_ids(universe) || !length(universe) || !all(universe %in% peaks)) {
    log_message("universe must contain unique eligible source peak IDs.", message_type = "error")
  }
  if (!atac_unique_ids(features) || !length(features) || !all(features %in% universe)) {
    log_message("features must contain unique foreground IDs within universe.", message_type = "error")
  }
  if (!is.character(p.adjust.method) || length(p.adjust.method) != 1L ||
    !p.adjust.method %in% stats::p.adjust.methods) {
    log_message("p.adjust.method must name one supported stats::p.adjust method.", message_type = "error")
  }
  matched <- is.null(background)
  matching <- NULL
  if (matched) {
    if (!is.numeric(n.background) || length(n.background) != 1L ||
      !is.finite(n.background) || n.background < 1 || n.background %% 1 != 0) {
      log_message("n.background must be one positive integer.", message_type = "error")
    }
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed %% 1 != 0) {
      log_message("seed must be one nonnegative integer within the R integer range.", message_type = "error")
    }
    if (!atac_unique_ids(match.features) || !length(match.features)) {
      log_message("match.features must name unique feature-metadata columns.", message_type = "error")
    }
    reference <- setdiff(universe, features)
    if (!length(reference)) log_message("No non-foreground peaks are available for background matching.", message_type = "error")
    metadata <- object[[assay]][[]]
    if (!all(match.features %in% colnames(metadata)) ||
      !all(universe %in% rownames(metadata))) {
      log_message("Required matching metadata are missing; run native RegionStats or provide verified metadata.", message_type = "error")
    }
    metadata <- metadata[universe, match.features, drop = FALSE]
    if (!all(vapply(metadata, function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) {
      log_message("Matching metadata must be finite numeric values for every eligible peak.", message_type = "error")
    }
    n.requested <- min(n.background, length(reference))
    with_seed <- atac_backend_fun("withr", "with_seed")
    reference <- with_seed(seed, Signac::MatchRegionStats(
      meta.feature = metadata[reference, , drop = FALSE],
      query.feature = metadata[features, , drop = FALSE],
      features.match = match.features, n = n.requested, verbose = verbose
    ))
    if (!atac_unique_ids(reference) || length(reference) != n.requested ||
      !all(reference %in% setdiff(universe, features))) {
      log_message("Native background matching returned an invalid reference set.", message_type = "error")
    }
    background <- c(features, reference)
    matching <- list(
      features = match.features, seed = seed,
      maximum_reference_peaks = n.background, actual_reference_peaks = length(reference)
    )
  } else {
    if (!atac_unique_ids(background) || !length(background) ||
      !all(background %in% universe) || !all(features %in% background)) {
      log_message("Explicit background must be a unique population within universe and include every foreground peak.", message_type = "error")
    }
    reference <- setdiff(background, features)
  }
  results <- Signac::FindMotifs(
    object = object, assay = assay, features = features,
    background = background, p.adjust.method = p.adjust.method, verbose = verbose
  )
  list(results = results, receipt = list(
    assay = assay, genome = report$genome,
    foreground = features, universe = universe, reference = reference,
    background = background, motif_ids = colnames(incidence),
    matching = matching, p.adjust.method = p.adjust.method,
    population = "explicit background including the foreground",
    background_source = if (matched) "Signac::MatchRegionStats plus foreground" else "explicit",
    backend = c(Signac = as.character(utils::packageVersion("Signac")))
  ))
}
