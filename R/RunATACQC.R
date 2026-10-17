#' @title Run scATAC quality control metrics
#'
#' @description
#' Calculate common scATAC QC metrics and optionally filter cells by thresholds.
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @param tss.positions TSS positions passed to `Signac::TSSEnrichment`. An
#'   explicit value requests recalculation, even when an imported TSS score exists.
#' @param blacklist A `GRanges` blacklist used to compute `blacklist_ratio`.
#' @param fast Whether to use the fast mode in `Signac::TSSEnrichment`.
#' @param min_pct_reads_in_peaks,min_TSS_enrichment,max_nucleosome_signal,max_blacklist_ratio
#' Optional thresholds used for filtering cells. Fraction thresholds use 0--1,
#' not percentages. A requested metric must be present and finite for every cell;
#' missing or failed calculations never silently disable a requested filter.
#' @param frip Explicitly recalculate a fragment-record fraction overlapping the
#'   current peak union. Each record counts once, even across several peaks.
#'   Duplicate input rows remain separate records; column-five support is not
#'   used as a weight. Stores `frip` and the legacy 0--1 alias.
#' @param frip.chunk.size Maximum fragment records read per sequential block.
#' @param frip.columns Optional named character vector with `in_peaks` and
#'   `total` metadata column names. Supplying these explicitly asserts that both
#'   count the same read/fragment unit over the same files, cells and filtering
#'   rules, and that the numerator refers to the current peak universe. Legacy
#'   column names are not guessed. Without either this argument or `frip = TRUE`,
#'   an existing `pct_reads_in_peaks` column is preserved, but a new FRiP is not
#'   inferred from matrix column sums.
#'
#' @details
#' `total.fragments` counts fragment-file records (`frequency_count`), not
#' column-five read support (`reads_count`). Each file's explicit Signac barcode
#' mapping is applied separately, so captures may reuse the same raw barcode.
#' Resources without selected-assay cells are ignored during calculation but
#' remain attached to the returned assay.
#' `blacklist_ratio` is the fraction of peak-matrix counts assigned to peaks
#' overlapping any blacklist interval. Each peak contributes once; it is not a
#' fragment-level blacklist fraction. Existing QC columns are imported values;
#' callers must verify their source, units and relevance after changing peaks.
#' Explicit FRiP source columns must not use names overwritten by fragment
#' counting (`frequency_count`, `reads_count`, `mononucleosomal`,
#' `nucleosome_free`, `total.fragments`); retain upstream values under distinct
#' names when attaching fragments. TSS enrichment uses Signac's R implementation
#' and does not select an external backend merely because a newer Signac release
#' exports it.
#' Missing optional inputs leave metrics unavailable. This function does not
#' perform doublet detection or establish that all fragment QC is complete.
#'
#' @return A `Seurat` object, retaining input cell order or the ordered subset
#'   passing all requested thresholds. Available metadata include fragment
#'   counts, `nucleosome_signal`, `nucleosome_percentile`, `TSS.enrichment`,
#'   `pct_reads_in_peaks` and `blacklist_ratio`. Direct fragment FRiP also writes
#'   `frip` and a `scop_frip` receipt in the assay's `misc` slot, recording the
#'   peak union, mapped cells, source hashes and per-cell numerator/denominator.
#'   Explicit source-column recalculation also synchronizes `frip` and clears any
#'   previous direct-fragment receipt. Metrics are calculated only for cells in
#'   the selected assay; other cells and their existing metadata are preserved.
#'   Filtering removes only failed cells from the selected assay; at least one
#'   cell in that assay must pass the requested thresholds.
#'   Zero denominators produce `NA`; requested filters reject unavailable or
#'   non-finite metrics. Fraction metrics use 0--1, not percentages.
#' @seealso [ReadATACData()], [CheckATACData()]
#' @export
#' @examples
#' counts <- Matrix::sparseMatrix(
#'   i = c(1, 2, 1, 2), j = c(1, 1, 2, 2), x = c(2, 1, 1, 3),
#'   dimnames = list(c("chr1-1-10", "chr1-21-30"), c("cell1", "cell2"))
#' )
#' ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
#'   seqnames = "chr1", start = c(1L, 21L), end = c(10L, 30L)
#' ))
#' object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
#' object$source_in_peaks <- c(3, 4)
#' object$source_total <- c(5, 8)
#' object <- RunATACQC(object,
#'   frip.columns = c(in_peaks = "source_in_peaks", total = "source_total"),
#'   min_pct_reads_in_peaks = 0.5, verbose = FALSE
#' )
#' object$pct_reads_in_peaks
RunATACQC <- function(
  object,
  assay = NULL,
  tss.positions = NULL,
  blacklist = NULL,
  fast = TRUE,
  min_pct_reads_in_peaks = NULL,
  min_TSS_enrichment = NULL,
  max_nucleosome_signal = NULL,
  max_blacklist_ratio = NULL,
  verbose = TRUE,
  srt = NULL,
  frip.columns = NULL,
  frip = FALSE,
  frip.chunk.size = 100000L
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    log_message(
      "{.arg assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }

  get_fragments <- atac_backend_fun("Signac", "Fragments")
  original <- NULL
  assay_cells <- colnames(srt[[assay]])
  if (!setequal(assay_cells, colnames(srt))) {
    original <- srt
    srt <- subset(srt, cells = assay_cells)
  }

  if (!is.logical(frip) || length(frip) != 1L || is.na(frip)) {
    log_message("frip must be TRUE or FALSE.", message_type = "error")
  }
  if (frip && !is.null(frip.columns)) {
    log_message("Choose either direct fragment FRiP or explicit source columns.", message_type = "error")
  }

  if (!is.null(frip.columns) && length(get_fragments(srt[[assay]])) &&
    any(frip.columns %in% c(
      "frequency_count", "reads_count", "mononucleosomal",
      "nucleosome_free", "total.fragments"
    ))) {
    log_message("frip.columns cannot use fields overwritten by fragment counting; copy verified source counts to distinct metadata columns.", message_type = "error")
  }

  log_message("Calculating ATAC QC metrics...", verbose = verbose)

  fragment_counts <- atac_frag_counts(
    srt = srt,
    assay = assay,
    verbose = verbose
  )
  if (!is.null(fragment_counts)) {
    srt <- SeuratObject::AddMetaData(srt, metadata = fragment_counts)
  }

  if (!"nucleosome_signal" %in% colnames(srt@meta.data)) {
    nucleosome_md <- atac_nuc_meta(fragment_counts)
    if (!is.null(nucleosome_md)) {
      srt <- SeuratObject::AddMetaData(srt, metadata = nucleosome_md)
    } else {
      log_message("Nucleosome signal unavailable: no fragment counts available",
        message_type = "warning", verbose = verbose
      )
    }
  }

  if (!is.null(tss.positions)) atac_qc_ranges(tss.positions, srt[[assay]], "TSS positions")
  if (!is.null(tss.positions) && is.null(fragment_counts)) {
    log_message("Explicit tss.positions requires available fragment counts.", message_type = "error")
  }
  if ((!"TSS.enrichment" %in% colnames(srt@meta.data) || !is.null(tss.positions)) &&
    !is.null(fragment_counts)) {
    srt <- tryCatch(
      {
        fragments <- get_namespace_fun("Signac", "Fragments")(srt[[assay]])
        tss_object <- srt
        tss_object[[assay]]@fragments <- Filter(function(fragment) {
          any(names(get_namespace_fun("Signac", "GetFragmentData")(fragment, "cells")) %in% colnames(srt[[assay]]))
        }, fragments)
        tss_object <- suppressWarnings(
          get_namespace_fun("Signac", "TSSEnrichment")(
            object = tss_object,
            assay = assay,
            tss.positions = tss.positions,
            fast = fast,
            verbose = verbose
          )
        )
        tss_object[[assay]]@fragments <- fragments
        tss_object
      },
      error = function(error) {
        if (!is.null(tss.positions)) stop(error)
        log_message(
          "Skip TSS enrichment: {.val {conditionMessage(error)}}",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )
  }

  if (frip) {
    direct_frip <- atac_fragment_frip(srt, assay, fragment_counts, frip.chunk.size)
    srt <- SeuratObject::AddMetaData(srt, metadata = direct_frip$metadata)
    SeuratObject::Misc(srt[[assay]], "scop_frip") <- direct_frip$receipt
  } else if (!is.null(frip.columns)) {
    frip_md <- atac_frip_meta(srt, frip.columns = frip.columns)
    srt <- SeuratObject::AddMetaData(srt, metadata = frip_md)
    srt[[assay]]@misc$scop_frip <- NULL
  } else if (!"pct_reads_in_peaks" %in% colnames(srt@meta.data)) {
    log_message(
      "FRiP unavailable: supply explicit frip.columns with matching count units and the current peak universe",
      message_type = "warning", verbose = verbose
    )
  }

  if ("pct_reads_in_peaks" %in% colnames(srt@meta.data)) {
    frip <- srt$pct_reads_in_peaks
    if (!is.numeric(frip) || any(!is.finite(frip) & !is.na(frip)) ||
      any(frip < 0 | frip > 1, na.rm = TRUE)) {
      log_message("pct_reads_in_peaks must contain fractions in 0--1; percentages are not converted automatically.", message_type = "error")
    }
  }

  if (!is.null(blacklist)) {
    problems <- CheckATACData(srt, assay = assay)$errors
    if (length(problems)) log_message("{paste(problems, collapse = \" \")}", message_type = "error")
    atac_qc_ranges(blacklist, srt[[assay]], "Blacklist")
    counts <- SeuratObject::GetAssayData(srt, assay = assay, layer = "counts")
    overlaps <- !is.na(GenomicRanges::findOverlaps(atac_backend_fun("Signac", "granges")(srt[[assay]]),
      blacklist,
      ignore.strand = TRUE, select = "first"
    ))
    total <- Matrix::colSums(counts)
    fraction <- Matrix::colSums(counts[overlaps, , drop = FALSE]) / total
    fraction[total == 0] <- NA_real_
    srt <- SeuratObject::AddMetaData(srt, metadata = fraction, col.name = "blacklist_ratio")
  }

  thresholds <- list(
    pct_reads_in_peaks = min_pct_reads_in_peaks,
    TSS.enrichment = min_TSS_enrichment, nucleosome_signal = max_nucleosome_signal,
    blacklist_ratio = max_blacklist_ratio
  )
  keep <- atac_qc_keep(srt, assay, thresholds)
  if (!any(keep)) log_message("No cells satisfy the requested ATAC QC thresholds.", message_type = "error")
  if (!is.null(original)) {
    metadata <- original[[]]
    qc_metadata <- srt[[]]
    for (column in setdiff(colnames(qc_metadata), colnames(metadata))) {
      metadata[[column]] <- NA_real_
    }
    metadata[rownames(qc_metadata), colnames(qc_metadata)] <- qc_metadata
    original <- SeuratObject::AddMetaData(original, metadata = metadata)
    original[[assay]] <- srt[[assay]]
    full_keep <- stats::setNames(rep(TRUE, ncol(original)), colnames(original))
    full_keep[names(keep)] <- keep
    keep <- full_keep
    srt <- original
  }

  if (any(!keep)) {
    log_message(
      "Filtering {.val {sum(!keep)}} cells using ATAC QC thresholds",
      verbose = verbose
    )
    srt <- subset(srt, cells = colnames(srt)[keep])
  }

  log_message(
    "Available ATAC QC metrics calculated",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )
  srt
}

atac_frag_counts <- function(srt, assay, verbose = TRUE) {
  fragments <- get_namespace_fun("Signac", "Fragments")(srt[[assay]])
  if (!length(fragments)) {
    return(NULL)
  }
  cells <- colnames(srt[[assay]])
  maps <- lapply(fragments, get_namespace_fun("Signac", "GetFragmentData"), slot = "cells")
  maps <- lapply(maps, function(x) x[names(x) %in% cells])
  active <- lengths(maps) > 0L
  fragments <- fragments[active]
  maps <- maps[active]
  if (!all(vapply(maps, function(x) atac_unique_ids(x) && atac_unique_ids(names(x)), logical(1)))) {
    log_message("Fragment counts require explicit unique cell-to-barcode mappings.", message_type = "error")
  }
  paths <- vapply(fragments, get_namespace_fun("Signac", "GetFragmentData"), character(1), slot = "path")
  if (!all(file.exists(paths))) {
    log_message("Fragment counts unavailable: local fragment files are missing",
      message_type = "warning", verbose = verbose
    )
    return(NULL)
  }
  paths <- normalizePath(paths, mustWork = TRUE)
  if (anyDuplicated(paths)) {
    log_message("Fragment objects must refer to different files, not aliases of the same path.", message_type = "error")
  }
  mapped <- unlist(lapply(maps, names), use.names = FALSE)
  if (anyDuplicated(mapped) || !setequal(mapped, cells)) {
    log_message("Fragment mappings must cover each assay cell exactly once.", message_type = "error")
  }
  counts <- lapply(seq_along(fragments), function(i) {
    map <- maps[[i]]
    if (!length(map)) {
      return(NULL)
    }
    result <- get_namespace_fun("Signac", "CountFragments")(
      fragments = paths[[i]], cells = unname(map),
      max_lines = NULL, verbose = verbose
    )
    if (anyDuplicated(result$CB) || !all(unname(map) %in% result$CB)) {
      log_message("{paste0(\"Fragment counts are missing mapped barcodes in \", paths[[i]], \".\")}", message_type = "error")
    }
    result <- result[match(unname(map), result$CB), , drop = FALSE]
    rownames(result) <- names(map)
    result$CB <- NULL
    result
  })
  counts <- do.call(rbind, counts)[cells, , drop = FALSE]
  counts$total.fragments <- counts$frequency_count
  counts
}

atac_nuc_meta <- function(fragment_counts) {
  if (is.null(fragment_counts)) {
    return(NULL)
  }
  required_cols <- c("mononucleosomal", "nucleosome_free")
  if (!all(required_cols %in% colnames(fragment_counts))) {
    return(NULL)
  }

  nucleosome_signal <- ifelse(
    fragment_counts$nucleosome_free > 0,
    fragment_counts$mononucleosomal / fragment_counts$nucleosome_free,
    NA_real_
  )
  md <- data.frame(
    nucleosome_signal = nucleosome_signal,
    row.names = rownames(fragment_counts)
  )
  finite_signal <- md$nucleosome_signal[is.finite(md$nucleosome_signal)]
  md$nucleosome_percentile <- NA_real_
  if (length(finite_signal) > 0) {
    percentile <- stats::ecdf(finite_signal)
    md$nucleosome_percentile[is.finite(md$nucleosome_signal)] <- round(
      percentile(md$nucleosome_signal[is.finite(md$nucleosome_signal)]),
      2
    )
  }
  md
}

atac_frip_meta <- function(srt, frip.columns) {
  if (!is.character(frip.columns) || length(frip.columns) != 2L ||
    !setequal(names(frip.columns), c("in_peaks", "total")) || anyNA(frip.columns) ||
    !all(frip.columns %in% colnames(srt@meta.data))) {
    log_message("frip.columns must name existing in_peaks and total metadata columns.", message_type = "error")
  }
  numerator <- srt[[frip.columns[["in_peaks"]], drop = TRUE]]
  denominator <- srt[[frip.columns[["total"]], drop = TRUE]]
  if (!is.numeric(numerator) || !is.numeric(denominator) ||
    any(!is.finite(numerator)) || any(!is.finite(denominator)) ||
    any(numerator < 0) || any(denominator < 0) || any(numerator > denominator)) {
    log_message("FRiP requires finite non-negative matched counts with in_peaks <= total.", message_type = "error")
  }
  frip <- ifelse(denominator > 0, numerator / denominator, NA_real_)
  data.frame(frip = frip, pct_reads_in_peaks = frip, row.names = colnames(srt))
}


atac_qc_keep <- function(object, assay, thresholds) {
  cells <- colnames(object[[assay]])
  keep <- stats::setNames(rep(TRUE, ncol(object)), colnames(object))
  for (metric in names(thresholds)) {
    threshold <- thresholds[[metric]]
    if (is.null(threshold)) next
    if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 ||
      (metric %in% c("pct_reads_in_peaks", "blacklist_ratio") && threshold > 1)) {
      log_message("{paste0(\"Invalid threshold for \", metric, \"; fraction thresholds use 0--1.\")}", message_type = "error")
    }
    if (!metric %in% colnames(object@meta.data)) {
      log_message("{paste0(\"Cannot filter on unavailable metric: \", metric, \".\")}", message_type = "error")
    }
    values <- object[[metric, drop = TRUE]][cells]
    if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0) ||
      (metric %in% c("pct_reads_in_peaks", "blacklist_ratio") && any(values > 1))) {
      log_message("{paste0(\"Cannot filter on missing or invalid values in \", metric, \"; fraction values must use 0--1.\")}", message_type = "error")
    }
    keep[cells] <- keep[cells] & if (metric %in% c("pct_reads_in_peaks", "TSS.enrichment")) {
      values >= threshold
    } else {
      values <= threshold
    }
  }

  keep
}
