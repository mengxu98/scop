# Exact record-overlap fraction, using native indexed IO and interval operations.
atac_fragment_frip <- function(object, assay, fragment_counts, chunk.size) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = \" \")}", message_type = "error")
  if (is.null(fragment_counts)) log_message("Requested FRiP requires available fragments.", message_type = "error")
  if (!is.numeric(chunk.size) || length(chunk.size) != 1L || !is.finite(chunk.size) ||
    chunk.size < 1 || chunk.size > .Machine$integer.max || chunk.size %% 1 != 0) {
    log_message("frip.chunk.size must be one positive integer within the R integer range.", message_type = "error")
  }
  cells <- as.vector(colnames(object[[assay]]))
  validated <- atac_validate_fragment_inputs(object, assay, cells)
  peaks <- Signac::granges(object[[assay]])
  regions <- GenomicRanges::reduce(peaks, ignore.strand = TRUE)
  numerator <- denominator <- stats::setNames(numeric(length(cells)), cells)
  for (fragment in Signac::Fragments(object[[assay]])) {
    map <- Signac::GetFragmentData(fragment, "cells")
    map <- map[names(map) %in% cells]
    if (!length(map)) next
    file <- Signac::GetFragmentData(fragment, "path")
    resource <- match(normalizePath(file), validated$resources$path)
    handle <- Rsamtools::TabixFile(file,
      index = validated$resources$index[resource],
      yieldSize = as.integer(chunk.size)
    )
    open(handle)
    tryCatch(
      {
        repeat {
          lines <- Rsamtools::scanTabix(handle)[[1L]]
          if (!length(lines)) break
          fields <- data.table::fread(
            text = paste(lines, collapse = "\n"), header = FALSE,
            sep = "\t", colClasses = "character", showProgress = FALSE, data.table = FALSE
          )
          if (!ncol(fields) %in% c(5L, 6L)) log_message("Fragment records require five or six tab-separated fields.", message_type = "error")
          index <- match(fields[[4L]], unname(map))
          keep <- !is.na(index)
          if (!any(keep)) next
          fields <- fields[keep, , drop = FALSE]
          cell_index <- match(names(map)[index[keep]], cells)
          start <- suppressWarnings(as.numeric(fields[[2L]]))
          end <- suppressWarnings(as.numeric(fields[[3L]]))
          if (any(!is.finite(start)) || any(!is.finite(end)) || any(start < 0) ||
            any(end <= start) || any(start %% 1 != 0) || any(end %% 1 != 0)) {
            log_message("Selected fragment records contain invalid BED coordinates.", message_type = "error")
          }
          in_peak_contigs <- fields[[1L]] %in% GenomeInfoDb::seqlevels(regions)
          hit <- rep(FALSE, nrow(fields))
          # Each source line is represented once, preserving duplicate-line units.
          if (any(in_peak_contigs)) {
            records <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
              seqnames = fields[[1L]][in_peak_contigs],
              start = start[in_peak_contigs] + 1, end = end[in_peak_contigs]
            ))
            hit[in_peak_contigs] <- !is.na(GenomicRanges::findOverlaps(
              records, regions,
              ignore.strand = TRUE, select = "first"
            ))
          }
          denominator <- denominator + tabulate(cell_index, nbins = length(cells))
          numerator <- numerator + tabulate(cell_index[hit], nbins = length(cells))
        }
      },
      finally = close(handle)
    )
  }
  native_total <- fragment_counts[cells, "total.fragments"]
  if (!isTRUE(all.equal(unname(denominator), unname(native_total), check.attributes = FALSE))) {
    log_message("Sequential record totals disagree with native Signac CountFragments; FRiP was not stored.", message_type = "error")
  }
  fraction <- numerator / denominator
  metadata <- data.frame(frip = fraction, pct_reads_in_peaks = fraction, row.names = cells)
  list(metadata = metadata, receipt = list(
    assay = assay, cells = cells,
    peak_ids = atac_range_ids(peaks), union_ranges = regions,
    numerator = numerator, denominator = denominator,
    count_unit = "fragment-file records, including repeated identical lines, without read-support weighting",
    numerator_rule = "one count if record overlaps any interval in current peak union",
    denominator_rule = "all records of the same selected cells and files, across all indexed contigs",
    fragments = validated$resources, chunk_size = chunk.size,
    backend = c(
      Signac = as.character(utils::packageVersion("Signac")),
      Rsamtools = as.character(utils::packageVersion("Rsamtools")),
      GenomicRanges = as.character(utils::packageVersion("GenomicRanges"))
    )
  ))
}

# Native validators are deliberately separate: ValidateFragments forwards its
# extra arguments to ValidateHash, which does not accept barcode-scan options.
atac_validate_fragment_inputs <- function(object, assay, cells, verbose = FALSE) {
  fragments <- Signac::Fragments(object[[assay]])
  if (!length(fragments)) log_message("Attached local indexed fragments are required.", message_type = "error")
  maps <- lapply(fragments, Signac::GetFragmentData, slot = "cells")
  maps <- lapply(maps, function(x) x[names(x) %in% cells])
  active <- lengths(maps) > 0L
  fragments <- fragments[active]
  maps <- maps[active]
  if (!all(vapply(maps, function(x) atac_unique_ids(x) && atac_unique_ids(names(x)), logical(1)))) {
    log_message("Fragments require explicit unique cell-to-barcode mappings.", message_type = "error")
  }
  mapped <- unlist(lapply(maps, names), use.names = FALSE)
  if (anyDuplicated(mapped) || !setequal(mapped, cells)) {
    log_message("Fragment mappings must cover every selected cell exactly once.", message_type = "error")
  }
  paths <- vapply(fragments, Signac::GetFragmentData, character(1), slot = "path")
  if (!all(file.exists(paths))) log_message("Selected fragment files are missing or not local.", message_type = "error")
  paths <- normalizePath(paths, mustWork = TRUE)
  if (anyDuplicated(paths)) log_message("Fragment objects must refer to distinct files.", message_type = "error")
  indexes <- vapply(paths, function(path) {
    present <- paste0(path, c(".tbi", ".csi"))
    present <- present[file.exists(present)]
    if (!length(present)) log_message("{paste0(\"Fragment index is missing: \", path)}", message_type = "error")
    present[[1L]]
  }, character(1), USE.NAMES = FALSE)
  for (i in seq_along(fragments)) {
    if (!Signac::ValidateHash(fragments[[i]], verbose = verbose)) {
      log_message("{paste0(\"Fragment file/index hash differs from the stored source: \", paths[[i]])}", message_type = "error")
    }
    if (!Signac::ValidateCells(fragments[[i]],
      cells = unname(maps[[i]]),
      tolerance = 0, max.lines = NULL, verbose = verbose
    )) {
      log_message("{paste0(\"Selected fragment barcodes are absent: \", paths[[i]])}", message_type = "error")
    }
  }
  hashes <- lapply(fragments, Signac::GetFragmentData, slot = "hash")
  list(
    fragments = fragments, maps = maps,
    resources = data.frame(
      path = paths, index = indexes, n_cells = lengths(maps),
      fragment_md5 = vapply(hashes, `[`, character(1), 1L),
      index_md5 = vapply(hashes, `[`, character(1), 2L), row.names = NULL
    )
  )
}
