#' Link ATAC peaks to paired measured expression
#'
#' @md
#' @description Run native Signac peak-gene correlation and matched-background
#' testing on paired RNA and ATAC assays, retaining native link statistics.
#'
#' @param object A Seurat object containing paired chromatin and measured RNA
#'   assays. Identical cell names must identify the same biological cells. RNA
#'   must cover all cells of the selected ATAC assay; additional RNA-only cells
#'   are retained but excluded from this calculation. Missing pairs are rejected.
#' @param assay Chromatin assay; defaults to the active assay.
#' @param expression.assay Measured RNA assay, default `"RNA"`. Gene activity,
#'   imputed expression and chromVAR scores are not measured RNA.
#' @param peak.slot,expression.slot Single existing native layers, default
#'   `"counts"` and `"data"`. Split layers must be joined explicitly beforehand.
#' @param method Correlation method, `"pearson"` or `"spearman"`.
#' @param gene.coords Optional native gene-coordinate `GRanges` passed to
#'   `Signac::LinkPeaks()`. With `NULL`, native Signac collapses assay transcript
#'   annotation. Coordinates must declare the same genome label as the assay;
#'   known chromosome lengths must agree. Explicit coordinates require one
#'   interval per gene identity. No chromosome-style conversion is performed.
#' @param distance Maximum genomic distance to TSS in bases, default 500000.
#' @param min.distance Optional non-negative minimum distance to TSS.
#' @param min.cells Native strict detection threshold: peaks and genes must
#'   occur in more than this many cells, default 10.
#' @param genes.use Optional unique RNA feature IDs to test.
#' @param n.sample Number of matched background peaks, default 200.
#' @param pvalue.cutoff,score.cutoff Native strict p-value and absolute
#'   correlation thresholds, default 0.05. P-values are native unadjusted
#'   normal-tail values `pnorm(-abs(zscore))`, without a two-tail multiplier or
#'   multiple-testing correction; they are not adjusted significance estimates.
#' @param gene.id Whether RNA feature IDs match annotation gene_id rather than
#'   gene_name, default `FALSE`. Use `TRUE` for [ReadMultiomeData()] gene IDs.
#' @param seed Random seed used locally for native background sampling.
#' @param verbose Whether to print progress messages.
#'
#' @details
#' Run Signac RegionStats with an appropriate genome sequence before this function;
#' finite GC.percent between 0 and 100 and sequence.length equal to peak width
#' are required. These metadata checks do not verify the supplied DNA provenance.
#' Native accessibility background counts are recomputed on the paired cells in
#' a temporary assay, avoiding stale count metadata after cell subsetting.
#' Background sampling, correlation, z scores and significance filtering are
#' delegated unchanged to Signac. Its minimum support and chromosome-background
#' requirements still apply. Source counts and the active assay are preserved.
#' Existing native links or SCOP link results are never overwritten. No retained
#' links is a valid empty result, not evidence that regulation is absent.
#' Native errors for constant expression or insufficient background support are
#' preserved; no alternative correlation or significance model is substituted.
#'
#' @return The Seurat object with native `Signac::Links(object[[assay]])` and a
#'   receipt at `Misc(object[[assay]], "scop_link_peaks")`, recording source
#'   assays, layers, cells, input features, parameters, backend version and the
#'   number of retained links. Native link metadata contains gene, peak, score,
#'   zscore and pvalue for non-empty results. Links are associations, not causal
#'   or experimentally verified regulatory connections.
#' @seealso [ReadMultiomeData()], [Signac::LinkPeaks()], [Signac::RegionStats()]
#' @examples
#' if (requireNamespace("Signac", quietly = TRUE) &&
#'   requireNamespace("Biostrings", quietly = TRUE) &&
#'   requireNamespace("BSgenome", quietly = TRUE)) {
#'   local({
#'     withr::local_preserve_seed()
#'     set.seed(83)
#'     counts <- Matrix::Matrix(matrix(rpois(60 * 40, 4) + 1, 60,
#'       dimnames = list(paste0("p", 1:60), paste0("cell", 1:40))
#'     ), sparse = TRUE)
#'     ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'       seqnames = rep(c("chr1", "chr2"), each = 30),
#'       start = rep(seq(1, by = 200, length.out = 30), 2),
#'       end = rep(seq(100, by = 200, length.out = 30), 2)
#'     ))
#'     paired <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
#'     rna <- counts[c(1, 31), , drop = FALSE]
#'     rownames(rna) <- c("g1", "g2")
#'     paired[["RNA"]] <- SeuratObject::CreateAssayObject(counts = rna)
#'     genes <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'       seqnames = c("chr1", "chr2"), start = 201, end = 400, strand = "+",
#'       gene_id = c("g1", "g2"), gene_name = c("g1", "g2")
#'     ), keep.extra.columns = TRUE)
#'     GenomeInfoDb::genome(genes) <- "toy"
#'     # Small synthetic sequence solely for an executable method example.
#'     dna <- Biostrings::DNAStringSet(setNames(vapply(
#'       1:2, function(i) {
#'         paste(sample(c("A", "C", "G", "T"), 6000, replace = TRUE), collapse = "")
#'       },
#'       character(1)
#'     ), c("chr1", "chr2")))
#'     paired <- Signac::RegionStats(paired, assay = "peaks", genome = dna, verbose = FALSE)
#'     linked <- RunLinkPeaks(paired,
#'       assay = "peaks", expression.slot = "counts",
#'       gene.coords = genes, min.cells = 0, n.sample = 10, distance = 10000,
#'       pvalue.cutoff = 1, verbose = FALSE
#'     )
#'     Signac::Links(linked[["peaks"]])
#'   })
#' }
#' @export
RunLinkPeaks <- function(
  object, assay = NULL, expression.assay = "RNA",
  peak.slot = "counts", expression.slot = "data", method = "pearson",
  gene.coords = NULL, distance = 500000, min.distance = NULL, min.cells = 10L,
  genes.use = NULL, n.sample = 200L, pvalue.cutoff = 0.05, score.cutoff = 0.05,
  gene.id = FALSE, seed = 11L, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  atac_backend_fun("withr", "local_preserve_seed")()
  if (!atac_unique_ids(expression.assay) || length(expression.assay) != 1L ||
    !expression.assay %in% SeuratObject::Assays(object) || identical(expression.assay, assay) ||
    inherits(object[[expression.assay]], "ChromatinAssay")) {
    log_message("expression.assay must identify a distinct measured RNA assay.", message_type = "error")
  }
  rna <- object[[expression.assay]]
  if (!is.null(SeuratObject::Misc(rna, "scop_gene_activity")) ||
    !is.null(SeuratObject::Misc(rna, "scop_chromvar"))) {
    log_message("Gene activity or chromVAR scores cannot substitute for paired measured RNA.", message_type = "error")
  }
  chromatin <- object[[assay]]
  cells <- vapply(colnames(chromatin), identity, character(1), USE.NAMES = FALSE)
  if (!all(cells %in% colnames(rna))) {
    log_message("Measured RNA must cover every selected ATAC cell; missing pairs are not silently dropped.", message_type = "error")
  }
  for (entry in list(c(assay, peak.slot), c(expression.assay, expression.slot))) {
    if (length(entry) != 2L || anyNA(entry) ||
      !entry[[2L]] %in% SeuratObject::Layers(object[[entry[[1L]]]], search = NA)) {
      log_message("Requested layers must exist exactly; join split layers explicitly before linking peaks.", message_type = "error")
    }
    data <- SeuratObject::GetAssayData(object[[entry[[1L]]]], layer = entry[[2L]])
    if (!all(cells %in% colnames(data))) {
      log_message("Native link layers must contain every paired cell.", message_type = "error")
    }
    data <- data[, cells, drop = FALSE]
    values <- if (inherits(data, "sparseMatrix") && "x" %in% methods::slotNames(data)) data@x else as.vector(data)
    if (any(!is.finite(values))) {
      log_message("Native link layers must contain finite values and every paired cell in assay order.", message_type = "error")
    }
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
    !method %in% c("pearson", "spearman")) {
    log_message("method must be pearson or spearman.", message_type = "error")
  }
  if (!is.logical(gene.id) || length(gene.id) != 1L || is.na(gene.id)) {
    log_message("gene.id must be TRUE or FALSE.", message_type = "error")
  }
  for (value in list(distance, min.cells, n.sample, seed)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value %% 1 != 0 || value > .Machine$integer.max) {
      log_message("distance, min.cells, n.sample and seed must be non-negative R-range integers.", message_type = "error")
    }
  }
  if (length(cells) < 2L || distance < 1 || n.sample < 2 || min.cells >= length(cells) ||
    (!is.null(min.distance) && (!is.numeric(min.distance) || length(min.distance) != 1L ||
      !is.finite(min.distance) || min.distance < 0 || min.distance >= distance))) {
    log_message("Require at least two paired cells, distance > 0, n.sample >= 2, min.cells < cell count, and 0 <= min.distance < distance.", message_type = "error")
  }
  for (value in list(pvalue.cutoff, score.cutoff)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > 1) {
      log_message("pvalue.cutoff and score.cutoff must be between zero and one.", message_type = "error")
    }
  }
  if (!is.null(genes.use) && (!atac_unique_ids(genes.use) || !all(genes.use %in% rownames(rna)))) {
    log_message("genes.use must contain unique existing RNA feature IDs.", message_type = "error")
  }
  coordinates <- if (is.null(gene.coords)) Signac::Annotation(chromatin) else gene.coords
  if (!inherits(coordinates, "GRanges") || !length(coordinates)) {
    log_message("Gene coordinates must be a non-empty GRanges object.", message_type = "error")
  }
  source_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(chromatin)
  coordinate_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(coordinates)
  chromosomes <- as.character(atac_backend_fun("GenomeInfoDb", "seqnames")(coordinates))
  shared <- intersect(names(source_lengths), names(coordinate_lengths))
  if (any(source_lengths[shared] != coordinate_lengths[shared], na.rm = TRUE) ||
    any(GenomicRanges::start(coordinates) < 1L | GenomicRanges::end(coordinates) < GenomicRanges::start(coordinates)) ||
    any(GenomicRanges::end(coordinates) > source_lengths[chromosomes], na.rm = TRUE) ||
    any(GenomicRanges::end(coordinates) > coordinate_lengths[chromosomes], na.rm = TRUE)) {
    log_message("Gene coordinates conflict with known chromosome lengths or valid 1-based closed intervals.", message_type = "error")
  }
  labels <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(coordinates)))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  fields <- if (is.null(gene.coords)) c("gene_id", "gene_name", "gene_biotype") else if (gene.id) "gene_id" else "gene_name"
  if (!length(coordinates) || !all(fields %in% colnames(S4Vectors::mcols(coordinates))) ||
    length(report$genome) != 1L || !identical(labels, report$genome)) {
    log_message("Gene coordinates require native gene metadata and matching source genome labels.", message_type = "error")
  }
  gene_ids <- as.character(coordinates$gene_id)
  gene_names <- as.character(coordinates$gene_name)
  gene_key <- if (gene.id) gene_ids else gene_names
  if (anyNA(gene_key) || any(!nzchar(gene_key)) ||
    (is.null(gene.coords) && (anyNA(gene_ids) || any(!nzchar(gene_ids))))) {
    log_message("Gene coordinates need non-empty values for the selected gene identity and annotation grouping key.", message_type = "error")
  }
  if (is.null(gene.coords)) {
    strands <- as.character(GenomicRanges::strand(coordinates))
    strands[strands == "*"] <- "+"
    identities <- unique(data.frame(gene_id = gene_ids, chromosome = chromosomes, strand = strands))
    if (anyDuplicated(identities$gene_id)) {
      log_message("Each annotation gene_id must describe one contig and native strand.", message_type = "error")
    }
  } else if (anyDuplicated(gene_key[gene_key %in% (genes.use %||% rownames(rna))])) {
    log_message("Explicit gene.coords must contain one interval per selected gene identity; use gene.id = TRUE for duplicated display symbols.", message_type = "error")
  }
  feature_data <- chromatin[[]]
  if (!all(c("GC.percent", "sequence.length") %in% colnames(feature_data)) ||
    !is.numeric(feature_data$GC.percent) || !is.numeric(feature_data$sequence.length) ||
    any(!is.finite(feature_data$GC.percent)) || any(!is.finite(feature_data$sequence.length)) ||
    any(feature_data$GC.percent < 0 | feature_data$GC.percent > 100) ||
    any(feature_data$sequence.length != GenomicRanges::width(Signac::granges(chromatin)))) {
    log_message("Finite GC.percent and positive sequence.length are required; run native Signac RegionStats first.", message_type = "error")
  }
  if (!identical(
    as.vector(Signac::GRangesToString(Signac::granges(chromatin))),
    as.vector(rownames(chromatin))
  )) {
    log_message("Native LinkPeaks requires unmodified coordinate feature IDs; this assay contains sanitized contig names.", message_type = "error")
  }
  if (length(Signac::Links(chromatin)) || !is.null(SeuratObject::Misc(chromatin, "scop_link_peaks"))) {
    log_message("Native links or SCOP linking results already exist; use a separate assay.", message_type = "error")
  }
  parameters <- list(
    peak.slot = peak.slot, expression.slot = expression.slot,
    method = method, gene.coords = gene.coords, distance = distance,
    min.distance = min.distance, min.cells = min.cells, genes.use = genes.use,
    n_sample = n.sample, pvalue_cutoff = pvalue.cutoff,
    score_cutoff = score.cutoff, gene.id = gene.id
  )
  link_peaks <- atac_backend_fun("Signac", "LinkPeaks")
  with_seed <- atac_backend_fun("withr", "with_seed")
  # Restrict both native assays to the same selected cells without changing the source object.
  paired <- subset(object, cells = cells)
  SeuratObject::DefaultAssay(paired) <- expression.assay
  if (!identical(colnames(paired[[assay]]), colnames(paired[[expression.assay]]))) {
    log_message("Native paired assay selection did not produce aligned cell IDs.", message_type = "error")
  }
  # Native LinkPeaks otherwise reuses possibly stale all-cell accessibility covariates.
  paired[[assay]]@meta.features <- paired[[assay]]@meta.features[
    , colnames(paired[[assay]]@meta.features) != "count",
    drop = FALSE
  ]
  paired <- with_seed(seed, do.call(link_peaks, c(list(
    object = paired,
    peak.assay = assay, expression.assay = expression.assay, verbose = verbose
  ), parameters)))
  links <- Signac::Links(paired[[assay]])
  Signac::Links(object[[assay]]) <- links
  SeuratObject::Misc(object[[assay]], "scop_link_peaks") <- list(
    source_assay = assay,
    expression_assay = expression.assay, cells = cells, peaks = rownames(chromatin),
    genes = rownames(rna), genome = report$genome, parameters = parameters,
    seed = seed, n_links = length(links),
    background_count = "Recomputed by native Signac on the selected paired cells",
    sequence_statistics = "Caller-supplied RegionStats; sequence provenance not verified by this function",
    interpretation = "Native paired RNA-ATAC associations; not causal regulatory evidence.",
    backend = c(Signac = as.character(utils::packageVersion("Signac")))
  )
  object
}
