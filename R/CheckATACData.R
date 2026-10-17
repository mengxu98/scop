#' Check an ATAC object without changing it
#'
#' @md
#' @description
#' Report raw-count, peak-coordinate, cell/metadata and reference alignment for
#' an existing Seurat `ChromatinAssay`. Fragment path/index availability is
#' reported separately, so an otherwise valid matrix can be analysed without
#' fragments. This check does not install packages, retrieve resources, modify
#' the object, or inspect dependencies of downstream optional analyses.
#'
#' @param object A Seurat object with a Signac `ChromatinAssay`.
#' @param assay Assay to check; defaults to the object's default assay.
#' @param genome Optional expected genome label. Unknown references are
#'   reported, never inferred from chromosome spelling.
#' @param sample.by Optional metadata column containing sample IDs. This checks
#'   complete labels, not whether samples are independent biological replicates.
#'
#' @return A list with `valid`, `errors`, `warnings`, `assay`, `n_cells`,
#'   `n_peaks`, `genome`, `samples`, a `fragments` resource table, and
#'   `fragments_ready`. `valid` describes the matrix contract. Fragment readiness
#'   checks local file/index existence and cell mapping for resources with cells
#'   in the selected assay. The resource table retains all attached files; its
#'   `active` column identifies those required for this assay. Empty mappings
#'   left by native subsetting do not block readiness. It does not hash files,
#'   scan their contents, or guarantee fragment-dependent computation will pass.
#'   Use Signac's `ValidateFragments()` for content validation when required.
#' @seealso [ReadATACData()], [RunStandardWorkflow()]
#' @export
#' @examples
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 2), j = c(1, 2), x = c(2, 3), dims = c(2, 2),
#'   dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
#' )
#' ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
#' object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
#' report <- CheckATACData(object)
#' report$valid
#' report$fragments_ready
CheckATACData <- function(object, assay = NULL, genome = NULL, sample.by = NULL) {
  report <- list(
    valid = FALSE, errors = character(), warnings = character(), assay = assay,
    n_cells = NA_integer_, n_peaks = NA_integer_, genome = character(),
    samples = NULL,
    fragments = data.frame(
      path = character(), index = character(),
      file_exists = logical(), index_exists = logical(), n_cells = integer(),
      active = logical()
    ),
    fragments_ready = FALSE
  )
  if (!inherits(object, "Seurat")) {
    report$errors <- "object must be a Seurat object."
    return(report)
  }
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  report$assay <- assay
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) ||
    !assay %in% SeuratObject::Assays(object) || !inherits(object[[assay]], "ChromatinAssay")) {
    report$errors <- "assay must identify a Signac ChromatinAssay."
    return(report)
  }
  check_r("Signac", verbose = FALSE)

  chromatin <- object[[assay]]
  counts <- SeuratObject::GetAssayData(chromatin, layer = "counts")
  ranges <- getExportedValue("Signac", "granges")(chromatin)
  report$n_cells <- ncol(counts)
  report$n_peaks <- nrow(counts)
  report$errors <- c(atac_counts_problems(counts), atac_ranges_problems(ranges))
  if (length(ranges) != nrow(counts) ||
    !identical(as.vector(atac_range_ids(ranges)), as.vector(rownames(counts)))) {
    report$errors <- c(report$errors, "Matrix peak IDs and genomic intervals are not aligned.")
  }
  cells <- colnames(counts)
  metadata <- object[[]]
  if (!atac_unique_ids(rownames(metadata)) || !all(cells %in% rownames(metadata))) {
    report$errors <- c(report$errors, "Object metadata does not cover the assay cell IDs uniquely.")
  } else if (!is.null(sample.by)) {
    if (!is.character(sample.by) || length(sample.by) != 1L || is.na(sample.by) ||
      !sample.by %in% colnames(metadata)) {
      report$errors <- c(report$errors, "sample.by must name an existing metadata column.")
    } else {
      samples <- as.character(metadata[cells, sample.by])
      if (anyNA(samples) || any(!nzchar(samples))) {
        report$errors <- c(report$errors, "Sample labels must be present for every assay cell.")
      } else {
        report$samples <- table(samples)
      }
    }
  }
  get_genome <- atac_backend_fun("GenomeInfoDb", "genome")
  labels <- unique(c(as.character(get_genome(ranges)), as.character(get_genome(chromatin))))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  report$genome <- labels
  annotation <- get_namespace_fun("Signac", "Annotation")(chromatin)
  annotation_labels <- if (length(annotation)) unique(as.character(get_genome(annotation))) else character()
  annotation_labels <- annotation_labels[!is.na(annotation_labels) & nzchar(annotation_labels)]
  if (length(unique(c(labels, annotation_labels))) > 1L) {
    report$errors <- c(report$errors, "Assay seqinfo, ranges and annotation contain conflicting genome labels.")
  }
  assay_lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(chromatin)
  range_chromosomes <- as.character(GenomicRanges::seqnames(ranges))
  if (any(GenomicRanges::end(ranges) > assay_lengths[range_chromosomes], na.rm = TRUE)) {
    report$errors <- c(report$errors, "Ranges exceed chromosome lengths in the assay seqinfo.")
  }
  if (!is.null(genome)) {
    if (!is.character(genome) || length(genome) != 1L || is.na(genome) || !nzchar(genome)) {
      report$errors <- c(report$errors, "genome must be one non-empty expected reference label.")
    } else if (length(labels) != 1L || !identical(labels, genome)) {
      report$errors <- c(report$errors, "The expected genome is missing or differs from the assay reference.")
    }
  } else if (!length(labels)) {
    report$warnings <- c(report$warnings, "The assay genome is unknown; supply verified reference information before sequence analysis.")
  }
  fragment_objects <- get_namespace_fun("Signac", "Fragments")(chromatin)
  if (length(fragment_objects)) {
    report$fragments <- do.call(rbind, lapply(fragment_objects, function(fragment) {
      path <- get_namespace_fun("Signac", "GetFragmentData")(fragment, slot = "path")
      candidates <- paste0(path, c(".tbi", ".csi"))
      index <- candidates[file.exists(candidates)][1L]
      if (is.na(index)) index <- candidates[1L]
      data.frame(
        path = path, index = index,
        file_exists = file.exists(path), index_exists = file.exists(index),
        n_cells = length(get_namespace_fun("Signac", "GetFragmentData")(fragment, slot = "cells")),
        active = any(names(get_namespace_fun("Signac", "GetFragmentData")(fragment, slot = "cells")) %in% cells)
      )
    }))
    mappings <- lapply(fragment_objects, get_namespace_fun("Signac", "GetFragmentData"), slot = "cells")
    mappings <- lapply(mappings, function(x) x[names(x) %in% cells])
    mappings <- mappings[report$fragments$active]
    mapped_cells <- unlist(lapply(mappings, names), use.names = FALSE)
    mapping_ok <- all(vapply(mappings, function(x) {
      atac_unique_ids(x) && atac_unique_ids(names(x))
    }, logical(1))) && !anyDuplicated(mapped_cells) && all(cells %in% mapped_cells)
    resources_ready <- all((report$fragments$file_exists & report$fragments$index_exists)[report$fragments$active])
    report$fragments_ready <- resources_ready && mapping_ok
    if (!mapping_ok) {
      report$warnings <- c(
        report$warnings,
        "Fragment mappings require unique non-empty barcodes per file and unique coverage of assay cells."
      )
    }
    if (!resources_ready) {
      report$warnings <- c(report$warnings, "External fragment files or indexes are missing; matrix-only analysis remains available.")
    }
  } else {
    report$warnings <- c(report$warnings, "No fragments are attached; only matrix-supported analyses are available.")
  }
  report$valid <- !length(report$errors)
  report$fragments_ready <- report$valid && report$fragments_ready
  report
}
