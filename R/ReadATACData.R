#' Read raw ATAC counts into a ChromatinAssay
#'
#' @md
#' @description
#' Import one explicitly shared peak universe from a sparse matrix and genomic
#' ranges, a 10x ATAC/ARC H5 file, or Matrix Market/peaks BED/barcodes files.
#' SCOP checks and aligns identifiers; Signac creates the native assay. No cells
#' or peaks are filtered, and no resources are downloaded or installed.
#'
#' @param object A numeric sparse peak-by-cell matrix, a local 10x ATAC/ARC
#'   H5 path, or a local Matrix Market path. Separate matrices with different peak
#'   universes must be recounted before combining; this function never fills
#'   unmeasured peaks with zeros.
#'   The supported Signac constructor requires at least two cells and two peaks.
#' @param ranges A `GRanges` in 1-based closed coordinates, required with
#'   a sparse matrix. Named ranges are aligned to the matrix row names; unnamed ranges
#'   explicitly use matrix row order. Duplicate genomic intervals are rejected.
#' @param metadata Optional data frame with unique row names matching the input
#'   matrix cell IDs exactly. Rows are reordered by ID, never by position.
#'   Supplied `nCount_<assay>`/`nFeature_<assay>` columns must agree with the
#'   raw matrix; conflicting QC metadata is rejected.
#' @param peaks.file,barcodes.file Three-column peaks BED and one-column
#'   barcode paths, both required when `object` is a Matrix Market path.
#'   Gzip compression is supported. Without companion files a path is read as
#'   H5, importing only the explicitly declared `Peaks` modality from its
#'   standard feature table; H5 requires optional `hdf5r`.
#'   BED and 10x H5 peak IDs use 0-based half-open coordinates, converted to
#'   1-based closed `GRanges` (start + 1, end unchanged).
#' @param cell.names Optional unique final cell IDs, in input column order.
#'   For example, explicitly prefix barcodes from separate captures. Suffixes
#'   are never removed or guessed.
#' @param sample Optional single sample ID, or a named character vector keyed
#'   by input cell IDs, stored as `sample` metadata. A sample label does not
#'   establish biological replication; retain donor/condition/capture metadata.
#' @param genome Optional genome label such as `"hg38"`, checked against known
#'   labels on `ranges` and `annotation`. This sets a label, not a reference
#'   download. Supply seqlengths in `ranges` for chromosome-bound checks.
#' @param annotation Optional Signac-compatible gene annotation `GRanges`.
#' @param fragments Optional local fragment path(s). Files must be bgzip
#'   compressed and indexed. Signac validates the paths and barcode membership.
#' @param fragment.cells Explicit fragment barcode mapping: a named character
#'   vector for one file, or a list of vectors in `fragments` order. Vector
#'   names are final object cell IDs; values are exact barcodes in that file.
#'   Multiple files require disjoint mappings covering every imported cell.
#'   With one file and no mapping, input matrix barcodes are used unchanged.
#' @param assay Name of the new chromatin assay.
#' @param project Project name passed to `CreateSeuratObject()`.
#' @param verbose Whether Signac prints progress messages.
#'
#' @details
#' Signac uses genomic coordinates as feature names. Original matrix feature
#' IDs are retained in assay feature metadata `input_peak_id`; input barcodes
#' are retained as cell metadata `input_barcode`. Import provenance is stored
#' next to the data in `Misc(object[[assay]], "scop_atac_input")`.
#'
#' Fragment files remain external references. `saveRDS()` does not copy them.
#' After `readRDS()`, use [CheckATACData()] to inspect their availability and
#' export its `fragments` table as a resource manifest. Missing external files
#' do not invalidate otherwise valid matrix-only analysis. Update moved paths
#' explicitly using Signac's `UpdatePath()`; paths are never repaired silently.
#'
#' @return A Seurat object containing a Signac `ChromatinAssay` with raw counts.
#' @seealso [CheckATACData()], [RunStandardWorkflow()]
#' @export
#' @examples
#' peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'   seqnames = c("chr1", "chr1", "chr2"),
#'   start = c(1L, 101L, 201L), end = c(50L, 150L, 250L)
#' ))
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 3, 2, 3), j = c(1, 1, 2, 2), x = c(2, 1, 3, 1),
#'   dims = c(3, 2), dimnames = list(c("p1", "p2", "p3"), c("AA-1", "BB-1"))
#' )
#' atac <- ReadATACData(counts, ranges = peaks, genome = "hg38", sample = "S1")
#' CheckATACData(atac)$valid
ReadATACData <- function(
  object,
  ranges = NULL,
  metadata = NULL,
  peaks.file = NULL,
  barcodes.file = NULL,
  cell.names = NULL,
  sample = NULL,
  genome = NULL,
  annotation = NULL,
  fragments = NULL,
  fragment.cells = NULL,
  assay = "peaks",
  project = "ATAC",
  verbose = TRUE
) {
  check_r("Signac", verbose = FALSE)

  counts <- object
  h5 <- matrix.file <- NULL
  if (is.character(object)) {
    if (!is.null(peaks.file) || !is.null(barcodes.file)) {
      matrix.file <- object
    } else {
      h5 <- object
    }
  }
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
    log_message("assay must be one non-empty name.", message_type = "error")
  }
  source <- "sparse_matrix"
  input_files <- character()
  if (!is.null(h5)) {
    if (!is.null(ranges) || !is.null(peaks.file) || !is.null(barcodes.file)) {
      log_message("h5 already defines peaks and barcodes; do not supply ranges or companion files.", message_type = "error")
    }
    h5 <- atac_input_path(h5, "h5")
    atac_backend_fun("hdf5r", "h5version")
    h5_file <- hdf5r::H5File$new(h5, mode = "r")
    feature_types <- tryCatch(
      {
        if (!h5_file$exists("matrix/features/feature_type")) {
          log_message("H5 requires a standard matrix/features/feature_type table declaring Peaks.", message_type = "error")
        }
        h5_file[["matrix/features/feature_type"]][]
      },
      finally = h5_file$close_all()
    )
    if (!"Peaks" %in% feature_types) {
      log_message("H5 must explicitly declare a Peaks modality; coordinate-like IDs alone are insufficient.", message_type = "error")
    }
    counts <- Seurat::Read10X_h5(h5, use.names = FALSE, unique.features = FALSE)
    if (is.list(counts)) {
      if (!"Peaks" %in% names(counts)) {
        log_message("H5 must contain one Peaks matrix; multi-genome lists are not supported.", message_type = "error")
      }
      counts <- counts[["Peaks"]]
    }
    ids <- rownames(counts)
    matched <- regmatches(ids, regexec("^(.+):([0-9]+)-([0-9]+)$", ids))
    if (length(matched) == 0L || any(lengths(matched) != 4L)) {
      log_message("H5 peak IDs must use 10x BED coordinates: chromosome:start-end.", message_type = "error")
    }
    bed <- as.data.frame(do.call(rbind, lapply(matched, `[`, 2:4)), stringsAsFactors = FALSE)
    ranges <- atac_bed_ranges(bed)
    names(ranges) <- ids
    source <- "10x_h5"
    input_files <- c(h5 = h5)
  } else if (!is.null(matrix.file)) {
    if (!is.null(ranges) || is.null(peaks.file) || is.null(barcodes.file)) {
      log_message("A Matrix Market path requires peaks.file and barcodes.file, without ranges.", message_type = "error")
    }
    input_files <- c(
      matrix = atac_input_path(matrix.file, "matrix.file"),
      peaks = atac_input_path(peaks.file, "peaks.file"),
      barcodes = atac_input_path(barcodes.file, "barcodes.file")
    )
    con <- gzfile(input_files[["matrix"]], open = "rt")
    counts <- tryCatch(Matrix::readMM(con), finally = close(con))
    bed <- utils::read.delim(input_files[["peaks"]], header = FALSE, colClasses = "character")
    barcodes <- utils::read.delim(input_files[["barcodes"]], header = FALSE, colClasses = "character")
    if (ncol(bed) != 3L || ncol(barcodes) != 1L ||
      nrow(bed) != nrow(counts) || nrow(barcodes) != ncol(counts)) {
      log_message("Matrix dimensions must match three-column peaks BED and one-column barcodes.", message_type = "error")
    }
    ranges <- atac_bed_ranges(bed)
    rownames(counts) <- paste0(bed[[1]], ":", bed[[2]], "-", bed[[3]])
    colnames(counts) <- barcodes[[1]]
    names(ranges) <- rownames(counts)
    source <- "matrix_market_bed"
  } else if (!is.null(peaks.file) || !is.null(barcodes.file)) {
    log_message("peaks.file and barcodes.file are only used with a Matrix Market path.", message_type = "error")
  }

  problems <- atac_counts_problems(counts)
  if (length(problems)) log_message(paste(problems, collapse = " "), message_type = "error")
  if (any(dim(counts) < 2L)) {
    log_message("ReadATACData requires at least two cells and two peaks with the supported Signac constructor.", message_type = "error")
  }
  problems <- atac_ranges_problems(ranges)
  if (length(problems)) log_message(paste(problems, collapse = " "), message_type = "error")
  if (length(ranges) != nrow(counts)) {
    log_message("ranges must have one interval per matrix row.", message_type = "error")
  }
  input_peak_ids <- rownames(counts)
  if (!is.null(names(ranges))) {
    if (!atac_unique_ids(names(ranges)) || !setequal(names(ranges), input_peak_ids)) {
      log_message("Named ranges must match matrix peak IDs exactly.", message_type = "error")
    }
    ranges <- ranges[match(input_peak_ids, names(ranges))]
  }
  ranges <- atac_set_genome(ranges, genome, annotation)
  peak_ids <- atac_range_ids(ranges)
  if (anyDuplicated(peak_ids)) {
    log_message("Genomic intervals collide after native Seurat feature-name conversion.", message_type = "error")
  }
  input_cells <- as.vector(colnames(counts))
  if (is.null(cell.names)) cell.names <- input_cells
  if (!atac_unique_ids(cell.names) || length(cell.names) != length(input_cells)) {
    log_message("cell.names must contain one unique non-empty ID per matrix column.", message_type = "error")
  }
  if (is.null(metadata)) metadata <- data.frame(row.names = input_cells)
  if (!is.data.frame(metadata) || !atac_unique_ids(rownames(metadata)) ||
    !setequal(rownames(metadata), input_cells) ||
    (ncol(metadata) > 0L && !atac_unique_ids(colnames(metadata)))) {
    log_message("metadata must be a data frame with unique rows matching input cell IDs exactly.", message_type = "error")
  }
  metadata <- metadata[input_cells, , drop = FALSE]
  qc_columns <- paste0(c("nCount_", "nFeature_"), assay)
  for (column in intersect(qc_columns, colnames(metadata))) {
    expected <- if (identical(column, qc_columns[[1L]])) {
      Matrix::colSums(counts)
    } else {
      Matrix::colSums(counts > 0)
    }
    if (!is.numeric(metadata[[column]]) ||
      !identical(as.numeric(metadata[[column]]), as.numeric(expected))) {
      log_message("Metadata {column} conflicts with the imported raw count matrix.", message_type = "error")
    }
  }
  if ("input_barcode" %in% names(metadata) &&
    !identical(as.character(metadata$input_barcode), input_cells)) {
    log_message("metadata input_barcode conflicts with matrix cell IDs.", message_type = "error")
  }
  metadata$input_barcode <- input_cells
  if (!is.null(sample)) {
    if (!is.character(sample) || anyNA(sample) || any(!nzchar(sample))) {
      log_message("sample must contain non-empty character labels.", message_type = "error")
    }
    if (length(sample) == 1L && is.null(names(sample))) {
      sample <- rep(sample, length(input_cells))
    } else {
      if (!atac_unique_ids(names(sample)) || !setequal(names(sample), input_cells)) {
        log_message("A per-cell sample vector must be named by all input cell IDs.", message_type = "error")
      }
      sample <- unname(sample[input_cells])
    }
    if ("sample" %in% names(metadata) && !identical(as.character(metadata$sample), sample)) {
      log_message("sample conflicts with metadata$sample.", message_type = "error")
    }
    metadata$sample <- sample
  }
  rownames(metadata) <- cell.names
  colnames(counts) <- cell.names
  fragment_objects <- atac_import_fragments(fragments, fragment.cells, input_cells, cell.names, verbose)
  chromatin <- get_namespace_fun("Signac", "CreateChromatinAssay")(
    counts = counts, ranges = ranges,
    genome = atac_backend_fun("GenomeInfoDb", "seqinfo")(ranges),
    annotation = annotation, fragments = fragment_objects,
    min.cells = 0, min.features = 0, verbose = verbose
  )
  # Native fragment lookup can attach separate fastmatch cache pointers to
  # otherwise identical ID vectors. Compare the ordered identifiers, not caches.
  if (!identical(as.vector(rownames(chromatin)), as.vector(peak_ids)) ||
    !identical(as.vector(colnames(chromatin)), as.vector(cell.names))) {
    log_message("Signac changed the expected feature or cell alignment during import.", message_type = "error")
  }
  chromatin <- SeuratObject::AddMetaData(chromatin, data.frame(
    input_peak_id = input_peak_ids, row.names = rownames(chromatin)
  ))
  SeuratObject::Misc(chromatin, "scop_atac_input") <- list(
    source = source, files = input_files, assay = assay,
    coordinate_system = "1-based closed",
    input_coordinate_system = if (source == "sparse_matrix") "1-based closed" else "0-based half-open",
    backend_versions = vapply(
      c("Signac", "Seurat", "SeuratObject"),
      function(x) as.character(utils::packageVersion(x)), character(1)
    )
  )
  result <- SeuratObject::CreateSeuratObject(
    counts = chromatin, assay = assay,
    meta.data = metadata, project = project
  )
  # Fragment barcode lookup may attach a fastmatch cache to the shared ID vector.
  rownames(metadata) <- as.vector(rownames(metadata))
  metadata$input_barcode <- as.vector(metadata$input_barcode)
  if (!all(colnames(metadata) %in% colnames(result[[]]))) {
    log_message("Native object creation omitted supplied metadata.", message_type = "error")
  }
  observed_metadata <- result[[]][, colnames(metadata), drop = FALSE]
  rownames(observed_metadata) <- as.vector(rownames(observed_metadata))
  observed_metadata$input_barcode <- as.vector(observed_metadata$input_barcode)
  if (!identical(observed_metadata, metadata)) {
    log_message("Native object creation changed or omitted supplied metadata.", message_type = "error")
  }
  result
}

atac_backend_fun <- function(package, fun) {
  available <- check_r(package, install = FALSE, verbose = FALSE)
  if (!isTRUE(all(unlist(available, use.names = FALSE)))) {
    log_message(sprintf("Package '%s' is required; install it separately before ATAC analysis.", package), message_type = "error")
  }
  resolved <- getExportedValue(package, fun)
  if (!is.function(resolved)) {
    log_message(sprintf("Function '%s' is not available in package '%s'.", fun, package), message_type = "error")
  }
  resolved
}

atac_unique_ids <- function(x) {
  is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
}

atac_input_path <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
    !file.exists(path) || dir.exists(path)) {
    log_message(sprintf("%s must be an existing local file.", label), message_type = "error")
  }
  normalizePath(path, mustWork = TRUE)
}

atac_counts_problems <- function(counts) {
  if (!inherits(counts, "sparseMatrix") || !"x" %in% methods::slotNames(counts) ||
    !is.numeric(counts@x)) {
    return("counts must be a numeric sparse matrix of raw counts.")
  }
  errors <- character()
  if (any(dim(counts) == 0L)) errors <- c(errors, "counts must contain cells and peaks.")
  if (!atac_unique_ids(rownames(counts))) errors <- c(errors, "Peak IDs must be unique and non-empty.")
  if (!atac_unique_ids(colnames(counts))) errors <- c(errors, "Cell IDs must be unique and non-empty.")
  if (any(!is.finite(counts@x) | counts@x < 0 | counts@x %% 1 != 0)) {
    errors <- c(errors, "counts must contain finite, non-negative integers; normalized values are not raw counts.")
  }
  errors
}

atac_ranges_problems <- function(ranges) {
  if (!inherits(ranges, "GRanges")) {
    return("ranges must be a GRanges object.")
  }
  start <- GenomicRanges::start(ranges)
  end <- GenomicRanges::end(ranges)
  seqnames <- as.character(atac_backend_fun("GenomeInfoDb", "seqnames")(ranges))
  lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(ranges)[seqnames]
  errors <- character()
  if (any(start < 1L | end < start) || any(!nzchar(seqnames))) {
    errors <- c(errors, "Ranges require positive 1-based starts and non-empty closed intervals.")
  }
  if (any(end > lengths, na.rm = TRUE)) errors <- c(errors, "Ranges exceed known chromosome lengths.")
  if (anyDuplicated(paste(seqnames, start, end, sep = ":"))) {
    errors <- c(errors, "Duplicate genomic intervals are not allowed.")
  }
  errors
}

atac_range_ids <- function(ranges) {
  gsub("[_|]", "-", get_namespace_fun("Signac", "GRangesToString")(ranges, sep = c("-", "-")))
}

atac_bed_ranges <- function(bed) {
  start <- suppressWarnings(as.numeric(bed[[2]]))
  end <- suppressWarnings(as.numeric(bed[[3]]))
  if (anyNA(bed[[1]]) || any(!nzchar(bed[[1]])) ||
    any(!is.finite(start) | !is.finite(end) | start < 0 | end <= start |
      start %% 1 != 0 | end %% 1 != 0 | end > .Machine$integer.max)) {
    log_message("BED requires chromosome, non-negative integer start, and end greater than start.", message_type = "error")
  }
  atac_backend_fun("GenomicRanges", "makeGRangesFromDataFrame")(
    data.frame(seqnames = bed[[1]], start = start, end = end), starts.in.df.are.0based = TRUE
  )
}

atac_set_genome <- function(ranges, genome, annotation) {
  get_genome <- atac_backend_fun("GenomeInfoDb", "genome")
  known <- unique(as.character(get_genome(ranges)))
  if (!is.null(annotation)) {
    if (!inherits(annotation, "GRanges")) log_message("annotation must be a GRanges object.", message_type = "error")
    known <- unique(c(known, as.character(get_genome(annotation))))
  }
  known <- known[!is.na(known) & nzchar(known)]
  if (length(known) > 1L || (!is.null(genome) &&
    (!is.character(genome) || length(genome) != 1L || is.na(genome) || !nzchar(genome) ||
      (length(known) && !identical(known, genome))))) {
    log_message("Genome labels on ranges, annotation, and genome must agree.", message_type = "error")
  }
  if (!is.null(genome)) {
    ranges <- atac_backend_fun("GenomeInfoDb", "genome<-")(ranges, value = genome)
  }
  ranges
}

atac_import_fragments <- function(paths, mappings, input_cells, cells, verbose) {
  if (is.null(paths)) {
    if (!is.null(mappings)) log_message("fragment.cells requires fragments.", message_type = "error")
    return(NULL)
  }
  if (!is.character(paths) || !length(paths) || anyNA(paths) || anyDuplicated(paths)) {
    log_message("fragments must be unique local fragment file paths.", message_type = "error")
  }
  paths <- vapply(paths, atac_input_path, character(1), label = "fragments")
  if (anyDuplicated(paths)) {
    log_message("fragments must identify different files after resolving local paths.", message_type = "error")
  }
  if (is.null(mappings) && length(paths) == 1L) mappings <- stats::setNames(input_cells, cells)
  if (is.character(mappings) && length(paths) == 1L) mappings <- list(mappings)
  if (!is.list(mappings) || length(mappings) != length(paths) ||
    !all(vapply(mappings, function(x) atac_unique_ids(x) && atac_unique_ids(names(x)), logical(1)))) {
    log_message("Provide one explicit named barcode mapping per fragment file.", message_type = "error")
  }
  mapped_cells <- unlist(lapply(mappings, names), use.names = FALSE)
  if (anyDuplicated(mapped_cells) || !setequal(mapped_cells, cells)) {
    log_message("Fragment mappings must cover every object cell exactly once across files.", message_type = "error")
  }
  lapply(seq_along(paths), function(i) {
    get_namespace_fun("Signac", "CreateFragmentObject")(
      path = paths[[i]], cells = mappings[[i]], validate.fragments = TRUE,
      tolerance = 0, max.lines = NULL, verbose = verbose
    )
  })
}

atac_new_assay_name <- function(object, name) {
  if (!atac_unique_ids(name) || length(name) != 1L ||
    !identical(make.names(name), name) || name %in% c(names(object), colnames(object[[]])) ||
    any(paste0(c("nCount_", "nFeature_"), name) %in% colnames(object[[]]))) {
    log_message("new.assay must be a syntactically valid unused name without assay, object or QC metadata collisions.", message_type = "error")
  }
  invisible(name)
}


atac_qc_ranges <- function(regions, chromatin, label) {
  if (!inherits(regions, "GRanges")) log_message(paste0(label, " must be a GRanges object."), message_type = "error")
  if (any(GenomicRanges::start(regions) < 1L |
    GenomicRanges::end(regions) < GenomicRanges::start(regions))) {
    log_message(paste0(label, " requires positive 1-based closed intervals."), message_type = "error")
  }
  get_genome <- atac_backend_fun("GenomeInfoDb", "genome")
  labels <- unique(c(
    as.character(get_genome(chromatin)),
    as.character(get_genome(Signac::granges(chromatin))), as.character(get_genome(regions))
  ))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) > 1L) log_message(paste0(label, " and assay genome labels must agree."), message_type = "error")
}
