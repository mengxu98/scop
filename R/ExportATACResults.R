#' Export an ATAC object and explicit result/resource tables
#'
#' @md
#' @description Save the complete object, supplied result tables and external
#' resource inventory in a new directory with file checksums.
#'
#' @param object A validated Seurat ChromatinAssay object.
#' @param outdir A new output directory. Existing paths are never overwritten.
#' @param assay Source chromatin assay; defaults to the active assay.
#' @param tables Uniquely named list of data frames, such as complete DAR tests,
#'   design/contrast tables, motif enrichment or peak annotation. Names must be
#'   safe filenames composed of letters, numbers, dots, underscores or hyphens,
#'   unique even when compared case-insensitively.
#' @param counts Whether to additionally write the sparse peak counts as Matrix
#'   Market. The full object RDS always retains its assays and analysis results.
#' @param compress Compression passed to saveRDS for the complete object.
#'
#' @details
#' The output contains object.rds, cells.csv, peaks.bed (0-based half-open),
#' fragments.csv, fragment_cells.csv, check.rds, sessionInfo.txt, supplied result
#' tables and a file manifest with MD5 checksums. Fragments and reference genomes
#' remain external; no automatic copying, downloading or path repair occurs.
#' The BED has three columns and, with counts.mtx and barcodes.tsv when
#' `counts = TRUE`, can be read by [ReadATACData()] when the selected assay
#' satisfies its native minimum of two peaks and two cells. Smaller subsets can
#' still be restored from object.rds. cells.csv and barcodes.tsv
#' contain the selected assay cells in matrix column order. object.rds preserves
#' all assays and cells. CSV files preserve values, not R column classes.
#' The `fragment` column joins fragment_cells.csv to fragments.csv. The resource
#' table retains inactive resources, with `active` marking current assay use;
#' `n_cells` is the attached mapping size, not the exported cell count.
#' Missing fragment files do not prevent a valid matrix-only export. Stored
#' Signac resource hashes are labeled as such, and resource availability is only
#' a snapshot. After moving resources, explicitly update paths with native Signac
#' UpdatePath and rerun CheckATACData before fragment-dependent computations.
#' Files are staged in the destination parent and moved to the requested new
#' directory only after all writes succeed. This does not modify the object.
#'
#' @return Invisibly, a list with the normalized directory, export file manifest
#'   and validation report. Tables are exported as supplied, without filtering.
#' @export
#' @examples
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
#'   dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
#' )
#' ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
#' object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
#' destination <- tempfile("atac-export-")
#' result <- ExportATACResults(object, destination, counts = TRUE)
#' result$manifest$file
#' restored <- readRDS(file.path(destination, "object.rds"))
#' CheckATACData(restored)$valid
#' matrix_only <- ReadATACData(file.path(destination, "counts.mtx"),
#'   peaks.file = file.path(destination, "peaks.bed"),
#'   barcodes.file = file.path(destination, "barcodes.tsv"), genome = "hg38"
#' )
#' CheckATACData(matrix_only)$valid
#' unlink(destination, recursive = TRUE)
ExportATACResults <- function(
  object, outdir, assay = NULL, tables = list(),
  counts = FALSE, compress = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message(paste(report$errors, collapse = " "), message_type = "error")
  assay <- report$assay
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) {
    log_message("outdir must name one new output directory.", message_type = "error")
  }
  if (atac_export_path_exists(outdir)) log_message("outdir already exists; select a new directory.", message_type = "error")
  parent <- normalizePath(dirname(outdir), mustWork = TRUE)
  outdir <- file.path(parent, basename(outdir))
  if (!is.logical(counts) || length(counts) != 1L || is.na(counts)) {
    log_message("counts must be TRUE or FALSE.", message_type = "error")
  }
  if (!is.list(tables) || (length(tables) && (!atac_unique_ids(names(tables)) || anyDuplicated(tolower(names(tables))) ||
    any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", names(tables))) ||
    !all(vapply(tables, is.data.frame, logical(1)))))) {
    log_message("tables must be a uniquely and safely named list of data frames.", message_type = "error")
  }
  cells <- as.vector(colnames(object[[assay]]))
  ranges <- Signac::granges(object[[assay]])
  fragments <- Signac::Fragments(object[[assay]])
  resources <- report$fragments
  resources$fragment <- seq_len(nrow(resources))
  mapping <- data.frame(fragment = integer(), cell = character(), barcode = character())
  if (length(fragments)) {
    resources$stored_signac_fragment_hash <- vapply(
      fragments,
      function(x) Signac::GetFragmentData(x, "hash")[[1L]], character(1)
    )
    resources$stored_signac_index_hash <- vapply(
      fragments,
      function(x) Signac::GetFragmentData(x, "hash")[[2L]], character(1)
    )
    mapping <- do.call(rbind, lapply(seq_along(fragments), function(i) {
      map <- Signac::GetFragmentData(fragments[[i]], "cells")
      map <- map[names(map) %in% cells]
      data.frame(
        fragment = rep(i, length(map)), cell = as.character(names(map)),
        barcode = unname(map), row.names = NULL
      )
    }))
  }
  staging <- tempfile(pattern = paste0(".", basename(outdir), "-"), tmpdir = parent)
  if (!dir.create(staging)) log_message("Cannot create export staging directory.", message_type = "error")
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  saveRDS(object, file.path(staging, "object.rds"), compress = compress)
  saveRDS(report, file.path(staging, "check.rds"))
  utils::write.csv(object[[]][cells, , drop = FALSE], file.path(staging, "cells.csv"))
  bed <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(ranges)),
    start = GenomicRanges::start(ranges) - 1L, end = GenomicRanges::end(ranges)
  )
  utils::write.table(bed, file.path(staging, "peaks.bed"),
    sep = "\t",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  utils::write.csv(resources, file.path(staging, "fragments.csv"), row.names = FALSE)
  utils::write.csv(mapping, file.path(staging, "fragment_cells.csv"), row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(staging, "sessionInfo.txt"))
  writeLines(
    c(
      "External fragment and reference files are not embedded by this export.",
      "peaks.bed uses 0-based half-open coordinates; object.rds ranges are 1-based closed.",
      "Stored Signac hashes are source fingerprints; existence fields are an export-time snapshot.",
      "After readRDS, use CheckATACData. Explicitly use Signac::UpdatePath for moved fragments.",
      "Result tables retain the complete rows supplied by the caller; no significance filtering is applied."
    ),
    file.path(staging, "README.txt")
  )
  if (counts) {
    Matrix::writeMM(SeuratObject::LayerData(object[[assay]], layer = "counts"), file.path(staging, "counts.mtx"))
    writeLines(cells, file.path(staging, "barcodes.tsv"))
  }
  if (length(tables)) {
    dir.create(file.path(staging, "tables"))
    for (name in names(tables)) {
      utils::write.csv(tables[[name]],
        file.path(staging, "tables", paste0(name, ".csv")),
        row.names = TRUE
      )
    }
  }
  files <- list.files(staging, recursive = TRUE, full.names = TRUE)
  manifest <- data.frame(
    file = substring(files, nchar(staging) + 2L),
    bytes = file.info(files)$size, md5 = unname(tools::md5sum(files))
  )
  utils::write.csv(manifest, file.path(staging, "manifest.csv"), row.names = FALSE)
  if (atac_export_path_exists(outdir) || !file.rename(staging, outdir)) {
    log_message("Could not finalize export into the new destination directory.", message_type = "error")
  }
  invisible(list(directory = outdir, manifest = manifest, check = report))
}

atac_export_path_exists <- function(path) {
  link <- Sys.readlink(path)
  file.exists(path) || dir.exists(path) || (!is.na(link) && nzchar(link))
}
