#' Infer peak coaccessibility with Cicero
#'
#' @md
#' @description Run the official monocle3-compatible Cicero workflow on the
#' selected ATAC counts and an explicitly selected cell embedding. No Pearson
#' fallback is substituted if Cicero or its model fitting is unavailable.
#' @param object Seurat object with a ChromatinAssay.
#' @param genome.sizes Named positive integer chromosome lengths for the verified
#' reference, covering all source peaks. No reference is downloaded.
#' @param genome.id Reference label matching the selected assay.
#' @param assay ATAC assay; defaults to active assay.
#' @param reduction Existing reduction whose cells cover every assay cell.
#' @param dims Explicit dimensions used for nearest-neighbor aggregation.
#' @param k Number of cells per aggregate; at least 2 and less than cell count.
#' @param window Genomic window length, greater than 250000 because the
#' run_cicero distance constraint is fixed at 250000.
#' @param sample.num Number of windows sampled for distance estimation.
#' @param seed Local random seed. Caller RNG state is restored.
#' @param result.name Unused assay Misc key for the result.
#' @param verbose Print progress messages.
#' @details Requires the official Cicero monocle3 branch and monocle3. See the
#' [official installation instructions](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/).
#' The installed Cicero package must declare a monocle3 dependency. Binary
#' accessibility is aggregated by Cicero before size-factor
#' normalization and distance-penalized graphical-model fitting. Model fitting
#' failures (including inadequate genomic coverage) are preserved. Peak IDs and
#' cell identities are not renamed or silently dropped. Returned NA scores are
#' unresolved/discordant estimates, not zero. Associations are not
#' experimentally verified physical interactions or causal regulation.
#' @return Input Seurat object with a result list in assay Misc(result.name),
#' containing connections, aggregate membership, parameters and source
#' identifiers. Links, source counts, metadata and default assay are
#' unchanged. Receipt identifiers are calculation-time snapshots.
#' @export
RunATACCoaccessibility <- function(
  object,
  genome.sizes,
  genome.id,
  assay = NULL,
  reduction = "lsi",
  dims = 2:30,
  k = 50L,
  window = 500000L,
  sample.num = 100L,
  seed = 1L,
  result.name = "scop_coaccessibility",
  verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay, genome = genome.id)
  if (!report$valid) {
    stop(paste(report$errors, collapse = " "), call. = FALSE)
  }
  assay <- report$assay
  if (
    !atac_unique_ids(result.name) ||
      length(result.name) != 1L ||
      !is.null(SeuratObject::Misc(object[[assay]], result.name))
  ) {
    stop("result.name must be an unused assay Misc key.", call. = FALSE)
  }
  if (
    !is.numeric(genome.sizes) ||
      !atac_unique_ids(names(genome.sizes)) ||
      any(
        !is.finite(genome.sizes) |
          genome.sizes < 1 |
          genome.sizes %% 1 != 0 |
          genome.sizes > .Machine$integer.max
      )
  ) {
    stop(
      "genome.sizes must be uniquely named positive R-range integer chromosome lengths.",
      call. = FALSE
    )
  }
  ranges <- Signac::granges(object[[assay]])
  chr <- as.character(GenomicRanges::seqnames(ranges))
  known <- GenomeInfoDb::seqlengths(object[[assay]])
  shared <- intersect(names(known), names(genome.sizes))
  if (
    any(!chr %in% names(genome.sizes)) ||
      any(GenomicRanges::end(ranges) > genome.sizes[chr]) ||
      any(known[shared] != genome.sizes[shared], na.rm = TRUE)
  ) {
    stop(
      "genome.sizes must cover source peaks and agree with known assay lengths.",
      call. = FALSE
    )
  }
  # Cicero strips the literal 'chr'; disallow aliases that would collapse.
  if (anyDuplicated(gsub("chr", "", names(genome.sizes), fixed = TRUE))) {
    stop(
      "Chromosome names collide under Cicero's native chr-prefix normalization.",
      call. = FALSE
    )
  }
  cells <- as.vector(colnames(object[[assay]]))
  if (
    !is.character(reduction) ||
      length(reduction) != 1L ||
      is.na(reduction) ||
      !reduction %in% SeuratObject::Reductions(object)
  ) {
    stop("reduction must name an existing cell embedding.", call. = FALSE)
  }
  embedding <- SeuratObject::Embeddings(object[[reduction]])
  if (
    !is.numeric(dims) ||
      !length(dims) ||
      anyNA(dims) ||
      any(dims < 1 | dims %% 1 != 0 | dims > ncol(embedding)) ||
      anyDuplicated(dims)
  ) {
    stop(
      "dims must select unique existing embedding dimensions.",
      call. = FALSE
    )
  }
  if (!all(cells %in% rownames(embedding))) {
    stop("Embedding must cover every selected ATAC cell.", call. = FALSE)
  }
  coords <- embedding[cells, dims, drop = FALSE]
  if (any(!is.finite(coords))) {
    stop("Selected coordinates must be finite.", call. = FALSE)
  }
  for (value in list(k, window, sample.num, seed)) {
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value < 0 ||
        value %% 1 != 0 ||
        value > .Machine$integer.max
    ) {
      stop(
        "k, window, sample.num and seed must be nonnegative R-range integers.",
        call. = FALSE
      )
    }
  }
  if (k < 2 || k >= length(cells) || window <= 250000 || sample.num < 1) {
    stop(
      "Require 2 <= k < cell count, window > 250000, and sample.num >= 1.",
      call. = FALSE
    )
  }
  counts <- SeuratObject::GetAssayData(object[[assay]], layer = "counts")
  if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (any(Matrix::colSums(counts) == 0) || any(Matrix::rowSums(counts) == 0)) {
    stop(
      "Cicero input requires positive depth for each cell and peak; filter explicitly.",
      call. = FALSE
    )
  }
  counts <- Matrix::drop0(counts)
  counts@x[] <- 1
  check_r(c("cole-trapnell-lab/monocle3", "cicero", "withr"), verbose = FALSE)
  new_cds <- get_namespace_fun("monocle3", "new_cell_data_set")
  make_cds <- get_namespace_fun("cicero", "make_cicero_cds")
  run <- get_namespace_fun("cicero", "run_cicero")
  dependencies <- paste(
    utils::packageDescription("cicero")[c("Imports", "Depends")],
    collapse = " "
  )
  if (!grepl("monocle3", dependencies, fixed = TRUE)) {
    stop(
      "Install the official monocle3-compatible Cicero branch (cole-trapnell-lab/cicero-release@monocle3); the monocle2 implementation is not supported by this entry point.",
      call. = FALSE
    )
  }
  withr::local_preserve_seed()
  set.seed(seed)
  metadata <- data.frame(
    gene_short_name = rownames(counts),
    chr = chr,
    bp1 = GenomicRanges::start(ranges),
    bp2 = GenomicRanges::end(ranges),
    row.names = rownames(counts)
  )
  cds <- new_cds(
    counts,
    cell_metadata = data.frame(row.names = cells),
    gene_metadata = metadata
  )
  aggregated <- make_cds(
    cds,
    reduced_coordinates = coords,
    k = k,
    size_factor_normalize = TRUE,
    silent = !verbose,
    return_agg_info = TRUE
  )
  if (
    !is.list(aggregated) ||
      length(aggregated) != 2L ||
      ncol(aggregated[[1L]]) < 2L
  ) {
    stop(
      "Native Cicero returned fewer than two aggregates; use more cells or a smaller k.",
      call. = FALSE
    )
  }
  connections <- run(
    aggregated[[1L]],
    genomic_coords = data.frame(
      chromosome = names(genome.sizes),
      length = as.numeric(genome.sizes)
    ),
    window = window,
    sample_num = sample.num,
    silent = !verbose
  )
  if (
    !is.data.frame(connections) ||
      !all(c("Peak1", "Peak2", "coaccess") %in% colnames(connections))
  ) {
    stop(
      "Native Cicero returned an unexpected connection schema.",
      call. = FALSE
    )
  }
  if (
    !all(connections$Peak1 %in% rownames(counts)) ||
      !all(connections$Peak2 %in% rownames(counts))
  ) {
    stop("Native Cicero returned unknown peak identifiers.", call. = FALSE)
  }
  SeuratObject::Misc(object[[assay]], result.name) <- list(
    connections = connections,
    aggregation = aggregated[[2L]],
    cells = cells,
    peaks = as.vector(rownames(counts)),
    genome = genome.id,
    genome.sizes = genome.sizes,
    reduction = reduction,
    dims = dims,
    coordinates = coords,
    k = k,
    window = window,
    sample.num = sample.num,
    seed = seed,
    binarized = TRUE,
    backend = c(
      cicero = as.character(utils::packageVersion("cicero")),
      monocle3 = as.character(utils::packageVersion("monocle3"))
    )
  )
  object
}
