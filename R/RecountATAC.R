#' Recount fragments over an explicit common ATAC peak universe
#'
#' @md
#' @description
#' Use native `Signac::FeatureMatrix()` to count local indexed fragments over
#' one supplied set of genomic ranges. Rows and cells are explicitly aligned;
#' no existing peak matrices are merged, padded, normalized or binarized.
#'
#' @param object A Seurat object with a validated Signac `ChromatinAssay` and
#'   attached local fragment objects with explicit cell-to-barcode mappings.
#' @param ranges A non-empty `GRanges` in 1-based closed coordinates. Its order
#'   is the output feature order. Intervals and optional names must be unique.
#'   The source assay must declare one verified genome label, and any labels
#'   on these ranges must agree. Unlabeled ranges explicitly inherit the source
#'   label; this does not verify their external provenance. Known lengths must
#'   agree, and source lengths for matching contigs are retained in the result.
#'   Chromosome spelling does not establish assembly.
#' @param assay Source chromatin assay; defaults to the object's default assay.
#' @param cells Unique source-assay cell IDs in requested matrix-return order.
#'   The default includes every source-assay cell. An empty selection is an
#'   error. Seurat writeback uses the original object cell order; the receipt
#'   records both requested and stored orders.
#' @param new.assay Optional unused syntactically valid assay name. `NULL` returns a list; a name
#'   returns a Seurat object with a new assay, leaving the source assay and its
#'   analysis results unchanged. The native constructor requires at least two
#'   peaks and two cells for this writeback route.
#' @param process_n Number of regions per native FeatureMatrix chunk.
#' @param verbose Whether native backends print progress messages.
#'
#' @details
#' Counts are fragment-file record overlaps: a record contributes once to each
#' overlapping feature, irrespective of its fifth-column read support. A record
#' may overlap several peaks, including disjoint peaks. Counts are not split or
#' reweighted, so matrix column sums are not necessarily unique fragments in peaks.
#'
#' Requested raw barcodes and stored file/index hashes are checked with Signac's
#' native validators before counting. Different captures may reuse raw barcodes,
#' but each final cell must map to exactly one file. No barcode suffixes, genome
#' labels or chromosome styles are guessed. Only files needed by selected cells
#' are read. Features on chromosomes present in each index are counted with
#' `keep_all_features = FALSE`; ranges on absent chromosomes are explicitly
#' mapped back as structural zero rows only when declared in the source assay's
#' `Seqinfo` or observed in another selected file's index. Unknown contigs fail.
#' A Seqinfo declaration establishes a reference name, not sequencing coverage;
#' zero means no record in the supplied file, not biological inaccessibility. This
#' avoids Signac 1.17.1's multi-chunk `keep_all_features = TRUE` zeroing bug
#' without increasing `process_n`. A selected file with no matching requested
#' chromosome fails clearly. Genuine zero-overlap features and cells are kept.
#'
#' Writeback happens only after validation and counting succeed. The new assay
#' contains the recounted matrix, supplied ranges, reference labels and selected
#' fragment references. Native initialization puts raw counts in its data layer;
#' no normalized data, scaled data, variable features, gene annotation, motifs,
#' links or old peak metadata are copied from the source. Existing reductions
#' remain associated with their original assay; compute a new representation
#' explicitly. The default assay is not changed. No resources are downloaded,
#' dependencies installed or external paths repaired.
#'
#' @return With `new.assay = NULL`, a list containing sparse `counts`, aligned
#'   `ranges` named by the output peak IDs, and a `receipt` with original range
#'   names, cell/peak IDs, file and reference information,
#'   native count semantics and backend version. With `new.assay` specified,
#'   a Seurat object; its receipt is stored at
#'   `Misc(object[[new.assay]], "scop_atac_recount")`.
#' @seealso [ReadATACData()], [CheckATACData()], [RunATACQC()]
#' @examples
#' if (requireNamespace("Rsamtools", quietly = TRUE)) {
#'   (function() {
#'     directory <- tempfile("atac-recount-")
#'     dir.create(directory)
#'     on.exit(unlink(directory, recursive = TRUE), add = TRUE)
#'     path <- file.path(directory, "fragments.tsv")
#'     writeLines(c("chr1\t10\t20\tA\t5", "chr1\t110\t120\tB\t3"), path)
#'     compressed <- Rsamtools::bgzip(path, dest = paste0(path, ".gz"))
#'     Rsamtools::indexTabix(compressed, format = "bed")
#'     peaks <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
#'     counts <- Matrix::sparseMatrix(
#'       i = 1:2, j = 1:2, x = 1, dims = c(2, 2),
#'       dimnames = list(c("p1", "p2"), c("cellA", "cellB"))
#'     )
#'     atac <- ReadATACData(counts,
#'       ranges = peaks, genome = "toy",
#'       fragments = compressed, fragment.cells = c(cellA = "A", cellB = "B"),
#'       verbose = FALSE
#'     )
#'     recounted <- RecountATAC(atac, ranges = peaks, verbose = FALSE)
#'     recounted$counts
#'   })()
#' }
#' @export
RecountATAC <- function(
  object, ranges, assay = NULL, cells = NULL,
  new.assay = NULL, process_n = 2000L, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  chromatin <- object[[assay]]
  if (length(report$genome) != 1L) {
    log_message("Recounting requires one verified source-assay genome label; chromosome names do not establish a reference.", message_type = "error")
  }
  problems <- atac_ranges_problems(ranges)
  if (length(problems)) log_message("{paste(problems, collapse = ' ')}", message_type = "error")
  if (!length(ranges)) log_message("ranges must contain at least one interval.", message_type = "error")
  if (!is.null(names(ranges)) && !atac_unique_ids(names(ranges))) {
    log_message("Range names must be unique non-empty IDs.", message_type = "error")
  }
  labels <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(ranges)))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) && !identical(labels, report$genome)) {
    log_message("Recount ranges and assay genome labels must agree.", message_type = "error")
  }
  source_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(chromatin)
  range_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(ranges)
  shared_lengths <- source_lengths[names(range_lengths)]
  if (any(range_lengths != shared_lengths, na.rm = TRUE)) {
    log_message("Recount ranges and source-assay chromosome lengths conflict.", message_type = "error")
  }
  chromosomes <- as.character(GenomicRanges::seqnames(ranges))
  if (any(GenomicRanges::end(ranges) > source_lengths[chromosomes], na.rm = TRUE)) {
    log_message("Recount ranges exceed source-assay chromosome lengths.", message_type = "error")
  }
  known_lengths <- !is.na(shared_lengths)
  range_lengths[known_lengths] <- shared_lengths[known_lengths]
  ranges <- atac_backend_fun("GenomeInfoDb", "seqlengths<-")(ranges, value = range_lengths)
  if (length(report$genome) == 1L) {
    ranges <- atac_backend_fun("GenomeInfoDb", "genome<-")(ranges, value = report$genome)
  }
  input_peak_ids <- names(ranges)
  feature_ids <- atac_range_ids(ranges)
  native_ids <- Signac::GRangesToString(ranges, sep = c("-", "-"))
  names(ranges) <- feature_ids
  if (!atac_unique_ids(feature_ids)) log_message("Range IDs collide after native feature naming.", message_type = "error")
  if (is.null(cells)) cells <- colnames(chromatin)
  if (!atac_unique_ids(cells) || !all(cells %in% colnames(chromatin))) {
    log_message("cells must be unique non-empty source-assay cell IDs.", message_type = "error")
  }
  # Keep fastmatch's native lookup cache off the caller's character vector.
  cells <- vapply(cells, identity, character(1), USE.NAMES = FALSE)
  if (!is.numeric(process_n) || length(process_n) != 1L || !is.finite(process_n) ||
    process_n < 1 || process_n %% 1 != 0 || process_n > .Machine$integer.max) {
    log_message("process_n must be a positive integer.", message_type = "error")
  }
  if (!is.null(new.assay)) {
    atac_new_assay_name(object, new.assay)
    if (length(ranges) < 2L || length(cells) < 2L) {
      log_message("New-assay writeback requires at least two peaks and two cells with the native constructor.", message_type = "error")
    }
  }

  fragments <- Signac::Fragments(chromatin)
  if (!length(fragments)) log_message("Recounting requires attached local indexed fragments.", message_type = "error")
  maps <- lapply(fragments, Signac::GetFragmentData, slot = "cells")
  maps <- lapply(maps, function(x) x[names(x) %in% cells])
  active <- lengths(maps) > 0L
  fragments <- fragments[active]
  maps <- maps[active]
  if (!all(vapply(maps, function(x) atac_unique_ids(x) && atac_unique_ids(names(x)), logical(1)))) {
    log_message("Selected fragments require explicit unique cell-to-barcode mappings.", message_type = "error")
  }
  mapped <- unlist(lapply(maps, names), use.names = FALSE)
  if (anyDuplicated(mapped) || !setequal(mapped, cells)) {
    log_message("Fragment mappings must cover each requested cell exactly once.", message_type = "error")
  }
  paths <- vapply(fragments, Signac::GetFragmentData, character(1), slot = "path")
  if (!all(file.exists(paths))) log_message("Selected fragment files are missing or not local.", message_type = "error")
  paths <- normalizePath(paths, mustWork = TRUE)
  if (anyDuplicated(paths)) log_message("Selected fragment objects must refer to different files.", message_type = "error")
  indexes <- vapply(paths, function(path) {
    candidates <- paste0(path, c(".tbi", ".csi"))
    index <- candidates[file.exists(candidates)][1L]
    if (is.na(index)) log_message("Fragment index is missing for {path}.", message_type = "error")
    index
  }, character(1), USE.NAMES = FALSE)
  tabix_file <- atac_backend_fun("Rsamtools", "TabixFile")
  seqnames_tabix <- atac_backend_fun("Rsamtools", "seqnamesTabix")
  missing_chromosomes <- vector("list", length(fragments))
  available_features <- vector("list", length(fragments))
  indexed_chromosomes <- vector("list", length(fragments))
  for (i in seq_along(fragments)) {
    if (!Signac::ValidateHash(fragments[[i]], verbose = verbose)) {
      log_message("Fragment file/index hash differs from the stored source: {paths[[i]]}.", message_type = "error")
    }
    if (!Signac::ValidateCells(fragments[[i]],
      cells = unname(maps[[i]]),
      tolerance = 0, max.lines = NULL, verbose = verbose
    )) {
      log_message("Requested fragment barcodes are absent from {paths[[i]]}.", message_type = "error")
    }
    available <- seqnames_tabix(tabix_file(paths[[i]], index = indexes[[i]]))
    indexed_chromosomes[[i]] <- available
    missing_chromosomes[[i]] <- setdiff(unique(chromosomes), available)
    available_features[[i]] <- which(chromosomes %in% available)
    if (!any(chromosomes %in% available)) {
      log_message("{paste0(\"No requested chromosomes occur in selected fragment file: \", paths[[i]], \"; verify reference and chromosome names.\")}", message_type = "error")
    }
  }
  declared_contigs <- names(source_lengths)
  known_contigs <- union(declared_contigs, unlist(indexed_chromosomes, use.names = FALSE))
  unknown_contigs <- setdiff(unique(chromosomes), known_contigs)
  if (length(unknown_contigs)) {
    log_message("{paste0(\"Requested chromosomes are not declared in source Seqinfo or observed in selected fragment indexes: \", paste(unknown_contigs, collapse = \", \"), \"; verify reference contigs before recounting.\")}", message_type = "error")
  }
  parts <- lapply(seq_along(fragments), function(i) {
    features <- available_features[[i]]
    file_cells <- cells[cells %in% names(maps[[i]])]
    native <- Signac::FeatureMatrix(
      fragments = fragments[[i]], features = ranges[features],
      cells = file_cells, keep_all_features = FALSE,
      process_n = as.integer(process_n), verbose = verbose
    )
    if (!inherits(native, "sparseMatrix") || !atac_unique_ids(rownames(native)) ||
      !atac_unique_ids(colnames(native)) || !setequal(rownames(native), native_ids[features]) ||
      !setequal(colnames(native), file_cells)) {
      log_message("FeatureMatrix did not return the requested sparse feature/cell universe.", message_type = "error")
    }
    native <- methods::as(native[native_ids[features], file_cells, drop = FALSE], "dgCMatrix")
    # Only index-confirmed absent chromosomes get structural zero rows. These
    # are counts from this same fragment source, not padded pre-existing matrices.
    Matrix::sparseMatrix(
      i = features[native@i + 1L], p = native@p, x = native@x,
      dims = c(length(ranges), length(file_cells)), dimnames = list(feature_ids, file_cells)
    )
  })
  counts <- do.call(cbind, parts)[feature_ids, cells, drop = FALSE]
  counts <- methods::as(counts, "dgCMatrix")
  problems <- atac_counts_problems(counts)
  if (length(problems)) log_message("{paste(problems, collapse = ' ')}", message_type = "error")
  cells <- vapply(cells, identity, character(1), USE.NAMES = FALSE)
  genome <- unique(as.character(atac_backend_fun("GenomeInfoDb", "genome")(ranges)))
  genome <- genome[!is.na(genome) & nzchar(genome)]
  hashes <- lapply(fragments, Signac::GetFragmentData, slot = "hash")
  receipt <- list(
    source_assay = assay, new_assay = new.assay, cells = cells, requested_cells = cells,
    peak_ids = feature_ids, input_peak_ids = input_peak_ids, genome = genome,
    fragments = data.frame(
      path = paths, index = indexes, n_cells = lengths(maps),
      fragment_md5 = vapply(hashes, `[`, character(1), 1L),
      index_md5 = vapply(hashes, `[`, character(1), 2L), row.names = NULL
    ),
    missing_chromosomes = stats::setNames(missing_chromosomes, paths),
    declared_reference_contigs = declared_contigs,
    indexed_contigs = stats::setNames(indexed_chromosomes, paths),
    unmatched_cells = character(), fragment_validation = "native hash and all requested barcodes",
    count_unit = "fragment-file record overlaps per feature",
    column5_weighted = FALSE, binary = FALSE,
    multi_peak_records = "one record may contribute to more than one feature; no splitting or reweighting",
    zero_count_cells = cells[Matrix::colSums(counts) == 0],
    zero_count_peaks = feature_ids[Matrix::rowSums(counts) == 0],
    backend = c(Signac = as.character(utils::packageVersion("Signac"))),
    parameters = list(keep_all_features = FALSE, process_n = as.integer(process_n)),
    absent_chromosomes = "structural zeros require source Seqinfo declaration or a selected index observation; neither establishes biological coverage",
    cell_order = "requested order", warnings = character()
  )
  if (is.null(new.assay)) {
    return(list(counts = counts, ranges = ranges, receipt = receipt))
  }

  stored_cells <- colnames(object)[colnames(object) %in% cells]
  counts <- counts[, stored_cells, drop = FALSE]
  receipt$cells <- stored_cells
  receipt$zero_count_cells <- stored_cells[Matrix::colSums(counts) == 0]
  receipt$cell_order <- "original Seurat object order on writeback"
  new_chromatin <- Signac::CreateChromatinAssay(
    counts = counts, ranges = ranges,
    fragments = fragments, genome = atac_backend_fun("GenomeInfoDb", "seqinfo")(ranges),
    min.cells = 0, min.features = 0, verbose = verbose
  )
  if (!identical(as.character(rownames(new_chromatin)), as.character(feature_ids)) ||
    !identical(as.character(colnames(new_chromatin)), as.character(stored_cells))) {
    log_message("Native assay construction changed the requested feature/cell order.", message_type = "error")
  }
  SeuratObject::Misc(new_chromatin, slot = "scop_atac_recount") <- receipt
  object[[new.assay]] <- new_chromatin
  if (!identical(as.character(rownames(object[[new.assay]])), as.character(feature_ids)) ||
    !identical(as.character(colnames(object[[new.assay]])), as.character(stored_cells)) ||
    !identical(unname(SeuratObject::GetAssayData(object, assay = new.assay, layer = "counts")), unname(counts))) {
    log_message("Native writeback changed the recounted matrix; no result returned.", message_type = "error")
  }
  object
}
