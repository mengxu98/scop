#' Quantify gene activity from indexed ATAC fragments
#'
#' @md
#' @description
#' Run native `Signac::GeneActivity()` with explicit annotation and counting
#' parameters, then store the sparse result in a new gene-activity assay.
#' Activity measures accessibility over gene bodies and promoters; it is not
#' measured RNA expression and does not establish regulatory links.
#'
#' @param object A Seurat object with a validated `ChromatinAssay`.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param new.assay An unused assay name. Existing assays are never overwritten.
#' @param features,extend.upstream,extend.downstream,biotypes,max.width,process_n,gene.id
#'   Parameters passed unchanged to `Signac::GeneActivity()`.
#' @param normalize Whether to populate the new assay's data layer using
#'   native Seurat LogNormalize. Raw activity counts remain in counts.
#' @param scale.factor LogNormalize scale factor, used only when normalizing.
#' @param verbose Whether native backends print progress messages.
#'
#' @details
#' Attach a verified transcript annotation with `Signac::Annotation()` before
#' calling this function. Its reference label must equal the source assay's
#' label. Required fields are gene_id, gene_name and gene_biotype; gene_id must
#' be non-missing and non-empty on every row. Repeated gene IDs may describe
#' multiple transcripts only on the same contig and native strand (unstranded
#' rows are treated as plus strand). Native Signac
#' collapses transcript coordinates to gene spans, extends promoter intervals,
#' and counts fragment-record overlaps; this wrapper does not substitute peak
#' sums. Native handling of missing chromosomes and gene filtering is retained.
#' Output genes and all source-assay cells are recorded in the assay receipt.
#' Native Signac disambiguates repeated gene symbols; use `gene.id = TRUE` for
#' stable gene identifiers. Final output IDs must be unique and survive native
#' Seurat assay creation unchanged. Native backends that drop matrix dimensions
#' for a single retained gene or cell are rejected rather than reconstructed.
#' Width filtering uses the gene span before promoter extension; native gene
#' filtering and missing-contig behavior do not establish sequencing coverage.
#'
#' Duplicate native output IDs, including distinct genes sharing an identical
#' counting window, are rejected explicitly. No private version-specific recovery,
#' synthetic IDs, row summation or silent gene dropping is performed. Select
#' unambiguous genes or revise the annotation explicitly before retrying. Gene
#' IDs alone do not resolve identical counting windows; their counts cannot
#' distinguish the associated genes.
#'
#' @return The Seurat object with a new activity assay and a source receipt in
#'   `Misc(object[[new.assay]], "scop_gene_activity")`. The default assay,
#'   source counts, original metadata columns and existing reductions are unchanged.
#'   Native count-QC columns may be added for the new assay. The receipt records
#'   active fragment resources and calculation-time IDs; later native subsetting
#'   or renaming does not synchronize this snapshot. Caller RNG state is restored.
#' @seealso [ReadATACData()], [CheckATACData()]
#' @examples
#' if (requireNamespace("Rsamtools", quietly = TRUE)) {
#'   (function() {
#'     directory <- tempfile("activity_example_")
#'     dir.create(directory)
#'     on.exit(unlink(directory, recursive = TRUE))
#'     plain <- file.path(directory, "fragments.tsv")
#'     writeLines(c(
#'       "chr1\t100\t120\tA\t5", "chr1\t300\t320\tB\t2",
#'       "chr1\t610\t630\tC\t3"
#'     ), plain)
#'     path <- paste0(plain, ".gz")
#'     Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
#'     Rsamtools::indexTabix(path, format = "bed")
#'     ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'       seqnames = "chr1", start = c(101, 301, 601), end = c(150, 350, 650)
#'     ))
#'     annotation <- ranges
#'     GenomeInfoDb::genome(annotation) <- "toy"
#'     annotation$gene_id <- paste0("GENE", 1:3)
#'     annotation$gene_name <- annotation$gene_id
#'     annotation$gene_biotype <- "protein_coding"
#'     annotation$type <- "gene"
#'     annotation$tx_id <- annotation$gene_id
#'     counts <- Matrix::sparseMatrix(
#'       i = 1:3, j = 1:3, x = 1,
#'       dimnames = list(paste0("p", 1:3), c("A", "B", "C"))
#'     )
#'     atac <- ReadATACData(counts,
#'       ranges = ranges, genome = "toy", annotation = annotation,
#'       fragments = path, fragment.cells = c(A = "A", B = "B", C = "C"), verbose = FALSE
#'     )
#'     activity <- RunGeneActivity(atac, gene.id = TRUE, extend.upstream = 0, verbose = FALSE)
#'     SeuratObject::GetAssayData(activity, assay = "ACTIVITY", layer = "counts")
#'   })()
#' }
#' @export
RunGeneActivity <- function(
  object, assay = NULL, new.assay = "ACTIVITY",
  features = NULL, extend.upstream = 2000, extend.downstream = 0,
  biotypes = "protein_coding", max.width = 500000, process_n = 2000,
  gene.id = FALSE, normalize = TRUE, scale.factor = 10000, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  atac_backend_fun("withr", "local_preserve_seed")()
  atac_new_assay_name(object, new.assay)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    log_message("normalize must be TRUE or FALSE.", message_type = "error")
  }
  if (normalize && (!is.numeric(scale.factor) || length(scale.factor) != 1L ||
    !is.finite(scale.factor) || scale.factor <= 0)) {
    log_message("scale.factor must be one finite positive number.", message_type = "error")
  }
  if (length(report$genome) != 1L) log_message("Gene activity requires one verified genome label.", message_type = "error")
  annotation <- Signac::Annotation(object[[assay]])
  if (!length(annotation) || !all(c("gene_id", "gene_name", "gene_biotype") %in% colnames(S4Vectors::mcols(annotation)))) {
    log_message("Gene activity requires transcript annotation with gene_id, gene_name and gene_biotype.", message_type = "error")
  }
  if (anyNA(annotation$gene_id) || any(!nzchar(as.character(annotation$gene_id)))) {
    log_message("Every annotation row requires a non-empty gene_id; repeated IDs may represent transcripts of the same gene.", message_type = "error")
  }
  starts <- GenomicRanges::start(annotation)
  ends <- GenomicRanges::end(annotation)
  chromosomes <- as.character(GenomicRanges::seqnames(annotation))
  strands <- as.character(GenomicRanges::strand(annotation))
  strands[strands == "*"] <- "+"
  identities <- unique(data.frame(
    gene_id = as.character(annotation$gene_id),
    chromosome = chromosomes, strand = strands, stringsAsFactors = FALSE
  ))
  if (anyDuplicated(identities$gene_id)) {
    log_message("Each gene_id must describe one contig and native strand; revise incompatible gene identity annotations explicitly.", message_type = "error")
  }
  get_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")
  annotation_lengths <- get_lengths(annotation)
  source_lengths <- get_lengths(object[[assay]])
  if (any(starts < 1L | ends < starts) || any(!nzchar(chromosomes)) ||
    any(ends > annotation_lengths[chromosomes], na.rm = TRUE) ||
    any(ends > source_lengths[chromosomes], na.rm = TRUE)) {
    log_message("Gene annotation requires positive closed intervals within known chromosome lengths.", message_type = "error")
  }
  if (any(annotation_lengths != source_lengths[names(annotation_lengths)], na.rm = TRUE)) {
    log_message("Known annotation and assay chromosome lengths must agree.", message_type = "error")
  }
  labels <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(annotation)))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!identical(labels, report$genome)) log_message("Gene annotation must declare the verified source genome.", message_type = "error")
  cells <- vapply(colnames(object[[assay]]), identity, character(1), USE.NAMES = FALSE)
  if (!report$fragments_ready) {
    log_message("Attached local indexed fragments are required: an active file or index is missing, or cell mappings are incomplete.", message_type = "error")
  }
  fragments <- Signac::Fragments(object[[assay]])[report$fragments$active]
  resources <- report$fragments[report$fragments$active, , drop = FALSE]
  resources$path <- normalizePath(resources$path, mustWork = TRUE)
  if (anyDuplicated(resources$path)) log_message("Active fragments must refer to distinct files.", message_type = "error")
  for (i in seq_along(fragments)) {
    mapping <- Signac::GetFragmentData(fragments[[i]], slot = "cells")
    mapping <- mapping[names(mapping) %in% cells]
    if (!Signac::ValidateHash(fragments[[i]], verbose = verbose)) {
      log_message("An active fragment file/index hash differs from its stored source.", message_type = "error")
    }
    if (!Signac::ValidateCells(fragments[[i]],
      cells = unname(mapping),
      tolerance = 0, max.lines = NULL, verbose = verbose
    )) {
      log_message("Selected fragment barcodes are absent from an active file.", message_type = "error")
    }
    resources$n_cells[i] <- length(mapping)
    hash <- Signac::GetFragmentData(fragments[[i]], slot = "hash")
    resources$fragment_md5[i] <- hash[[1L]]
    resources$index_md5[i] <- hash[[2L]]
  }
  counting_object <- object
  counting_object[[assay]]@fragments <- fragments
  # Native matching caches cell-name vectors in place; isolate the counting copy.
  colnames(counting_object[[assay]]) <- cells
  parameters <- list(
    features = features, extend.upstream = extend.upstream,
    extend.downstream = extend.downstream, biotypes = biotypes, max.width = max.width,
    process_n = process_n, gene.id = gene.id
  )
  counts <- do.call(Signac::GeneActivity, c(list(
    object = counting_object, assay = assay,
    verbose = verbose
  ), parameters))
  if (anyDuplicated(rownames(counts))) {
    log_message("Native GeneActivity returned duplicate gene IDs, potentially from shared counting windows. Select unambiguous genes or revise annotation explicitly; no rows were merged or discarded.", message_type = "error")
  }
  if (!inherits(counts, "sparseMatrix")) {
    log_message("Native GeneActivity did not return a sparse gene-by-cell matrix; this backend may drop dimensions for a single retained gene or cell.", message_type = "error")
  }
  if (!atac_unique_ids(rownames(counts))) {
    log_message("Native GeneActivity returned missing or non-unique gene IDs.", message_type = "error")
  }
  if (!atac_unique_ids(colnames(counts)) || !setequal(colnames(counts), cells)) {
    log_message("GeneActivity did not return every source-assay cell exactly once.", message_type = "error")
  }
  counts <- methods::as(counts[, cells, drop = FALSE], "dgCMatrix")
  problems <- atac_counts_problems(counts)
  if (length(problems)) log_message("{paste(problems, collapse = ' ')}", message_type = "error")
  activity <- SeuratObject::CreateAssayObject(counts = counts)
  if (!identical(as.vector(rownames(activity)), as.vector(rownames(counts)))) {
    log_message("Gene IDs changed during native assay creation; provide supported unique gene IDs.", message_type = "error")
  }
  if (isTRUE(normalize)) {
    activity <- Seurat::NormalizeData(activity,
      normalization.method = "LogNormalize",
      scale.factor = scale.factor, verbose = verbose
    )
  }
  SeuratObject::Misc(activity, "scop_gene_activity") <- list(
    source_assay = assay,
    genome = report$genome, annotation_rows = length(annotation), cells = cells,
    genes = rownames(counts), fragments = resources,
    parameters = parameters, normalization = if (isTRUE(normalize)) {
      list(method = "LogNormalize", scale.factor = scale.factor)
    } else {
      NULL
    },
    count_unit = "native fragment-record overlaps over gene bodies and promoters",
    backend = c(
      Signac = as.character(utils::packageVersion("Signac")),
      Seurat = as.character(utils::packageVersion("Seurat"))
    )
  )
  object[[new.assay]] <- activity
  object
}
