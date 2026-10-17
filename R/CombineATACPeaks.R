#' Combine explicit peak sets into one common interval universe
#'
#' @md
#' @description
#' Apply native `GenomicRanges::reduce()` within declared groups, optionally
#' require support from several input peak sets within each group, then take
#' the union of retained group intervals. No count matrices are combined.
#'
#' @param object A uniquely named list of `GRanges`. Each name identifies one
#'   peak-calling unit, not necessarily an independent biological replicate.
#' @param groups Optional named character vector assigning every peak-set name
#'   to one support group, for example a condition. With `NULL`, all peak sets
#'   form one group. Support filtering is done within groups before their union.
#' @param min.support Minimum number of distinct input peak sets supporting a
#'   reduced interval within a group. The default, one, does not require shared
#'   peaks across samples or conditions.
#' @param min.width,max.width Inclusive width limits applied to the final union.
#' @param genome Optional verified reference label. Known labels must agree;
#'   input sets without a reference label require this explicit declaration.
#'   Absent labels are not inferred from chromosome names. Known chromosome
#'   lengths must be compatible; the reduced union is checked again after it
#'   inherits known lengths from other inputs, before width filtering.
#'
#' @details
#' Reduction uses `ignore.strand = TRUE` and `min.gapwidth = 1L`, so adjacent
#' intervals are joined. Support counts distinct input sets contributing to a
#' connected reduced interval; it does not mean every base has that support.
#' Group-specific peaks can survive without support in other conditions.
#' Final metadata records contributing sets and retained support groups, not
#' caller p-values or a newly calculated consensus significance score.
#'
#' @return A list with `ranges` in native sorted order, named by canonical peak
#'   IDs, and a `receipt` recording input sets, support groups, filtering and
#'   backend versions. Recount fragments over these ranges before comparing
#'   counts from previously different peak universes.
#' @seealso [RecountATAC()]
#' @examples
#' first <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 201), width = 100))
#' second <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(51, 401), width = 100))
#' common <- CombineATACPeaks(object = list(first = first, second = second), genome = "toy")
#' common$ranges
#' common$receipt$input_sets
#' @export
CombineATACPeaks <- function(
  object, groups = NULL, min.support = 1L,
  min.width = 1L, max.width = Inf, genome = NULL
) {
  peaks <- object
  if (!is.list(peaks) || !length(peaks) || !atac_unique_ids(names(peaks))) {
    log_message("object must be a non-empty uniquely named list of GRanges.", message_type = "error")
  }
  for (name in names(peaks)) {
    problems <- atac_ranges_problems(peaks[[name]])
    if (length(problems)) log_message("{name}: {paste(problems, collapse = ' ')}", message_type = "error")
  }
  if (!sum(lengths(peaks))) log_message("No input peaks are available to combine.", message_type = "error")
  if (!is.numeric(min.support) || length(min.support) != 1L ||
    !is.finite(min.support) || min.support < 1 || min.support %% 1 != 0) {
    log_message("min.support must be one positive integer.", message_type = "error")
  }
  if (!is.numeric(min.width) || length(min.width) != 1L ||
    !is.finite(min.width) || min.width < 1 || min.width %% 1 != 0 ||
    !is.numeric(max.width) || length(max.width) != 1L || is.na(max.width) ||
    max.width < min.width || (is.finite(max.width) && max.width %% 1 != 0)) {
    log_message("Width limits must be positive ordered integers; max.width may be Inf.", message_type = "error")
  }
  if (is.null(groups)) groups <- stats::setNames(rep("all", length(peaks)), names(peaks))
  if (!is.character(groups) || !atac_unique_ids(names(groups)) ||
    !setequal(names(groups), names(peaks)) || anyNA(groups) || any(!nzchar(groups))) {
    log_message("groups must assign every peak-set name to one non-empty support group.", message_type = "error")
  }
  groups <- groups[names(peaks)]
  get_genome <- atac_backend_fun("GenomeInfoDb", "genome")
  known <- unique(unlist(lapply(peaks, function(x) as.character(get_genome(x))), use.names = FALSE))
  known <- known[!is.na(known) & nzchar(known)]
  if (!is.null(genome) && (!is.character(genome) || length(genome) != 1L ||
    is.na(genome) || !nzchar(genome))) {
    log_message("genome must be one non-empty verified reference label.", message_type = "error")
  }
  labels <- unique(c(known, genome))
  if (length(labels) > 1L) log_message("Input peaks contain conflicting reference labels.", message_type = "error")
  unlabeled <- vapply(peaks, function(x) {
    label <- get_genome(x)[unique(as.character(GenomicRanges::seqnames(x)))]
    length(x) > 0L && any(is.na(label) | !nzchar(label))
  }, logical(1))
  if (is.null(genome) && any(unlabeled)) {
    log_message("Unlabeled input peak sets require an explicit verified genome.", message_type = "error")
  }
  # reduce is a re-exported S4 generic, not defined in the GenomicRanges
  # namespace itself. Resolve its public export rather than a private lookup.
  reduce <- GenomicRanges::reduce
  character_list <- atac_backend_fun("IRanges", "CharacterList")
  strip <- function(x) {
    S4Vectors::mcols(x) <- NULL
    names(x) <- NULL
    x
  }
  retained <- list()
  support_inventory <- list()
  for (group in unique(groups)) {
    members <- names(groups)[groups == group]
    combined <- do.call(c, unname(lapply(peaks[members], strip)))
    if (!length(combined)) next
    source <- rep(members, lengths(peaks[members]))
    reduced <- reduce(combined, ignore.strand = TRUE, min.gapwidth = 1L, with.revmap = TRUE)
    sets <- lapply(as.list(reduced$revmap), function(i) unique(source[i]))
    support <- lengths(sets)
    keep <- support >= min.support
    support_inventory[[group]] <- data.frame(
      peak = atac_range_ids(reduced),
      support = support, kept = keep, row.names = NULL
    )
    reduced <- reduced[keep]
    reduced$revmap <- NULL
    reduced$input_sets <- character_list(sets[keep])
    reduced$support_group <- rep(group, length(reduced))
    retained[[group]] <- reduced
  }
  if (!sum(lengths(retained))) log_message("No peaks remain after within-group support filtering.", message_type = "error")
  combined <- do.call(c, unname(retained))
  ranges <- reduce(combined, ignore.strand = TRUE, min.gapwidth = 1L, with.revmap = TRUE)
  problems <- atac_ranges_problems(ranges)
  if (length(problems)) log_message("{paste(problems, collapse = ' ')}", message_type = "error")
  index <- as.list(ranges$revmap)
  # Extract ordinary vectors once: repeated S4 subsetting per interval is costly
  # for full peak universes and does not change the native revmap semantics.
  combined_input_sets <- as.list(combined$input_sets)
  combined_support_group <- as.character(combined$support_group)
  ranges$input_sets <- character_list(lapply(index, function(i) {
    unique(unlist(combined_input_sets[i], use.names = FALSE))
  }))
  ranges$support_groups <- character_list(lapply(index, function(i) {
    unique(combined_support_group[i])
  }))
  ranges$revmap <- NULL
  width <- GenomicRanges::width(ranges)
  keep <- width >= min.width & width <= max.width
  discarded <- atac_range_ids(ranges)[!keep]
  ranges <- ranges[keep]
  if (!length(ranges)) log_message("No peaks remain after final-union width filtering.", message_type = "error")
  if (length(labels)) ranges <- atac_backend_fun("GenomeInfoDb", "genome<-")(ranges, value = labels)
  names(ranges) <- atac_range_ids(ranges)
  if (!atac_unique_ids(names(ranges))) log_message("Canonical peak IDs are not unique.", message_type = "error")
  list(ranges = ranges, receipt = list(
    input_sets = data.frame(
      name = names(peaks), group = unname(groups),
      n_peaks = lengths(peaks), row.names = NULL
    ),
    genome = labels, support_inventory = support_inventory,
    discarded_by_width = discarded, peak_ids = names(ranges),
    parameters = list(
      min.support = min.support, min.width = min.width,
      max.width = max.width, ignore.strand = TRUE, min.gapwidth = 1L
    ),
    support_unit = "distinct input peak sets per connected reduced interval within group",
    backend = c(GenomicRanges = as.character(utils::packageVersion("GenomicRanges")))
  ))
}
