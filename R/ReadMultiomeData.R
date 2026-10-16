#' Read paired 10x RNA and ATAC counts
#'
#' @md
#' @description Import both measured modalities from one local 10x ARC count
#' matrix without dropping cells, normalizing counts, or inferring RNA from ATAC.
#'
#' Existing compatible Seurat objects can enter the analysis workflow directly;
#' this constructor is only needed when importing counts or files.
#'
#' @param object Path to one 10x ARC H5 file or one matrix directory containing
#'   `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`. The input must
#'   contain exactly the `Gene Expression` and `Peaks` modalities. H5 input
#'   requires at least two features per modality and two cells because the
#'   native Seurat H5 reader otherwise drops matrix dimensions.
#' @param metadata,cell.names,sample,genome,annotation,fragments,fragment.cells,project,verbose
#'   Passed to [ReadATACData()]. Metadata and sample vectors use original
#'   barcodes; fragment mappings use final cell names.
#' @param assay Name of the chromatin assay, default `"peaks"`.
#' @param rna.assay Name of the measured RNA assay, default `"RNA"`.
#'
#' @details
#' Gene IDs, rather than potentially duplicated gene symbols, identify RNA
#' features. Both modalities must contain the same unique barcodes. Their counts
#' are explicitly aligned to the ATAC barcode order. Peak coordinates in 10x
#' files are 0-based half-open and are converted to 1-based closed coordinates.
#' No barcode suffixes are removed and no absent modality is zero-filled.
#' H5 input requires the optional hdf5r package. Fragment files remain external
#' references, as described in [ReadATACData()].
#'
#' @return A Seurat object with a `ChromatinAssay` and a measured RNA `Assay`,
#'   both containing raw counts for the same ordered cells. The active assay is
#'   `assay`. Provenance is stored in `Misc(object, "scop_multiome_input")` and
#'   in the native ATAC import receipt. Native assay initialization also copies
#'   raw counts to the data layer; this is not normalized expression.
#' @seealso [ReadATACData()], [CheckATACData()]
#' @export
#' @examples
#' \dontrun{
#' multiome <- ReadMultiomeData("filtered_feature_bc_matrix.h5", genome = "hg38")
#' }
ReadMultiomeData <- function(
  object, metadata = NULL, cell.names = NULL,
  sample = NULL, genome = NULL, annotation = NULL, fragments = NULL,
  fragment.cells = NULL, assay = "peaks", rna.assay = "RNA",
  project = "Multiome", verbose = TRUE
) {
  if (!is.character(object) || length(object) != 1L || is.na(object) || !file.exists(object)) {
    log_message("object must identify one existing local 10x ARC file or directory.", message_type = "error")
  }
  if (!atac_unique_ids(rna.assay) || length(rna.assay) != 1L ||
    !identical(make.names(rna.assay), rna.assay) || identical(rna.assay, assay)) {
    log_message("rna.assay must be one syntactically valid name distinct from assay.", message_type = "error")
  }
  path <- normalizePath(object, mustWork = TRUE)
  log_message("Reading paired 10x RNA and ATAC counts...", message_type = "running", verbose = verbose)
  if (dir.exists(path)) {
    files <- file.path(path, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
    if (!all(file.exists(files))) {
      log_message("A 10x ARC directory requires matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz.", message_type = "error")
    }
    features <- utils::read.delim(files[[2L]], header = FALSE, colClasses = "character")
    if (ncol(features) < 3L || !atac_unique_ids(features[[1L]]) ||
      !setequal(unique(features[[3L]]), c("Gene Expression", "Peaks"))) {
      log_message("features.tsv.gz must contain unique feature IDs and exactly Gene Expression and Peaks modalities.", message_type = "error")
    }
    matrices <- Seurat::Read10X(
      data.dir = unname(path), gene.column = 1L,
      unique.features = TRUE, strip.suffix = FALSE
    )
    source <- "10x_arc_matrix_market"
  } else {
    atac_backend_fun("hdf5r", "h5version")
    matrices <- Seurat::Read10X_h5(filename = path, use.names = FALSE, unique.features = FALSE)
    files <- path
    source <- "10x_arc_h5"
  }
  if (!is.list(matrices) || !setequal(names(matrices), c("Gene Expression", "Peaks")) ||
    length(matrices) != 2L) {
    log_message("The source must contain exactly one Gene Expression and one Peaks matrix.", message_type = "error")
  }
  peaks <- matrices[["Peaks"]]
  rna <- matrices[["Gene Expression"]]
  if (source == "10x_arc_h5" && (is.null(dim(peaks)) || is.null(dim(rna)))) {
    log_message("The native Seurat H5 reader dropped matrix dimensions; H5 input requires at least two features per modality and two cells.", message_type = "error")
  }
  for (counts in list(peaks, rna)) {
    problems <- atac_counts_problems(counts)
    if (length(problems)) log_message(paste(problems, collapse = " "), message_type = "error")
  }
  if (!setequal(colnames(peaks), colnames(rna))) {
    log_message("RNA and ATAC barcodes must match exactly; cells cannot be silently intersected or filled.", message_type = "error")
  }
  rna <- rna[, colnames(peaks), drop = FALSE]
  ids <- rownames(peaks)
  fields <- regmatches(ids, regexec("^(.+):([0-9]+)-([0-9]+)$", ids))
  if (any(lengths(fields) != 4L)) {
    log_message("10x peak IDs must use chromosome:start-end BED coordinates.", message_type = "error")
  }
  ranges <- atac_bed_ranges(as.data.frame(do.call(rbind, lapply(fields, `[`, 2:4)),
    stringsAsFactors = FALSE
  ))
  names(ranges) <- ids
  result <- ReadATACData(peaks,
    ranges = ranges, metadata = metadata,
    cell.names = cell.names, sample = sample, genome = genome, annotation = annotation,
    fragments = fragments, fragment.cells = fragment.cells, assay = assay,
    project = project, verbose = verbose
  )
  atac_new_assay_name(result, rna.assay)
  colnames(rna) <- colnames(result[[assay]])
  rna_object <- SeuratObject::CreateAssayObject(counts = rna)
  if (!identical(as.vector(rownames(rna_object)), as.vector(rownames(rna)))) {
    log_message("Native RNA assay creation changed gene IDs; provide IDs supported by Seurat.", message_type = "error")
  }
  result[[rna.assay]] <- rna_object
  receipt <- SeuratObject::Misc(result[[assay]], "scop_atac_input")
  receipt$source <- source
  receipt$files <- files
  receipt$input_coordinate_system <- "0-based half-open"
  result[[assay]]@misc$scop_atac_input <- receipt
  SeuratObject::Misc(result, "scop_multiome_input") <- list(
    source = source,
    files = files, peak_assay = assay, expression_assay = rna.assay,
    cells = colnames(result[[assay]]), input_barcodes = colnames(peaks),
    genes = rownames(rna), peaks = rownames(result[[assay]]),
    expression_source = "measured Gene Expression counts from the same 10x ARC matrix",
    backend = c(
      Seurat = as.character(utils::packageVersion("Seurat")),
      Signac = as.character(utils::packageVersion("Signac"))
    )
  )
  result
}
