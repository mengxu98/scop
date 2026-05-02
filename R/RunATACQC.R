#' @title Run scATAC quality control metrics
#'
#' @description
#' Calculate common scATAC QC metrics and optionally filter cells by thresholds.
#'
#' @md
#' @inheritParams standard_scop
#' @param tss.positions TSS positions passed to `Signac::TSSEnrichment`.
#' @param blacklist A `GRanges` blacklist used to compute `blacklist_ratio`.
#' @param fast Whether to use the fast mode in `Signac::TSSEnrichment`.
#' @param min_pct_reads_in_peaks,min_TSS_enrichment,max_nucleosome_signal,max_blacklist_ratio
#' Optional thresholds used for filtering cells.
#'
#' @return A `Seurat` object with QC metadata added.
#' @export
#' @examples
#' \donttest{
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub <- RunATACQC(
#'   pbmcmultiome_sub,
#'   assay = "peaks",
#'   fast = TRUE
#' )
#' }
RunATACQC <- function(
  srt,
  assay = NULL,
  tss.positions = NULL,
  blacklist = NULL,
  fast = TRUE,
  min_pct_reads_in_peaks = NULL,
  min_TSS_enrichment = NULL,
  max_nucleosome_signal = NULL,
  max_blacklist_ratio = NULL,
  verbose = TRUE
) {
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
      srt <- tryCatch(
        {
          suppressWarnings(Signac::NucleosomeSignal(object = srt, assay = assay))
        },
        error = function(error) {
          log_message(
            "Skip nucleosome signal: {.val {conditionMessage(error)}}",
            message_type = "warning",
            verbose = verbose
          )
          srt
        }
      )
    }
  }

  used_modern_qc <- FALSE
  if (
    "ATACqc" %in% getNamespaceExports("Signac") &&
      is.null(tss.positions) &&
      atac_has_annotation(srt[[assay]])
  ) {
    atacqc_result <- tryCatch(
      {
        list(
          object = Signac::ATACqc(
            object = srt,
            assay = assay,
            verbose = verbose
          ),
          success = TRUE
        )
      },
      error = function(error) {
        log_message(
          "Skip Signac::ATACqc(): {.val {conditionMessage(error)}}",
          message_type = "warning",
          verbose = verbose
        )
        list(object = srt, success = FALSE)
      }
    )
    srt <- atacqc_result$object
    used_modern_qc <- atacqc_result$success
  }

  if (!"TSS.enrichment" %in% colnames(srt@meta.data) && !used_modern_qc) {
    srt <- tryCatch(
      {
        suppressWarnings(
          Signac::TSSEnrichment(
            object = srt,
            assay = assay,
            tss.positions = tss.positions,
            fast = fast
          )
        )
      },
      error = function(error) {
        log_message(
          "Skip TSS enrichment: {.val {conditionMessage(error)}}",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )
  }

  total_fragments_col <- atac_total_fragments_col(srt, assay = assay)
  frip_md <- atac_frip_meta(srt, assay = assay, total_fragments_col = total_fragments_col)
  if (!is.null(frip_md)) {
    srt <- SeuratObject::AddMetaData(srt, metadata = frip_md)
  } else if (!is.null(total_fragments_col)) {
    srt <- tryCatch(
      {
        Signac::FRiP(
          object = srt,
          assay = assay,
          total.fragments = total_fragments_col,
          col.name = "pct_reads_in_peaks"
        )
      },
      error = function(error) {
        log_message(
          "Skip FRiP calculation: {.val {conditionMessage(error)}}",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )
  } else {
    log_message(
      "Skip FRiP calculation: no total fragment count column or local fragments available",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (!is.null(blacklist)) {
    srt <- tryCatch(
      {
        Signac::FractionCountsInRegion(
          object = srt,
          assay = assay,
          regions = blacklist,
          col.name = "blacklist_ratio"
        )
      },
      error = function(error) {
        log_message(
          "Skip blacklist ratio: {.val {conditionMessage(error)}}",
          message_type = "warning",
          verbose = verbose
        )
        srt
      }
    )
  } else if (!"blacklist_ratio" %in% colnames(srt@meta.data)) {
    srt$blacklist_ratio <- NA_real_
  }

  keep <- rep(TRUE, ncol(srt))
  if (!is.null(min_pct_reads_in_peaks) && "pct_reads_in_peaks" %in% colnames(srt@meta.data)) {
    keep <- keep & srt$pct_reads_in_peaks >= min_pct_reads_in_peaks
  }
  if (!is.null(min_TSS_enrichment) && "TSS.enrichment" %in% colnames(srt@meta.data)) {
    keep <- keep & srt$TSS.enrichment >= min_TSS_enrichment
  }
  if (!is.null(max_nucleosome_signal) && "nucleosome_signal" %in% colnames(srt@meta.data)) {
    keep <- keep & srt$nucleosome_signal <= max_nucleosome_signal
  }
  if (!is.null(max_blacklist_ratio) && "blacklist_ratio" %in% colnames(srt@meta.data)) {
    keep <- keep & (is.na(srt$blacklist_ratio) | srt$blacklist_ratio <= max_blacklist_ratio)
  }

  if (any(!keep)) {
    log_message(
      "Filtering {.val {sum(!keep)}} cells using ATAC QC thresholds",
      verbose = verbose
    )
    srt <- subset(srt, cells = colnames(srt)[keep])
  }

  log_message(
    "ATAC QC completed",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )
  srt
}

atac_frag_counts <- function(srt, assay, verbose = TRUE) {
  fragment_paths <- atac_frag_paths(srt = srt, assay = assay)
  if (length(fragment_paths) == 0) {
    return(NULL)
  }
  if (any(grepl("^(https?|ftp)://", fragment_paths, ignore.case = TRUE))) {
    return(NULL)
  }
  if (!all(file.exists(fragment_paths))) {
    return(NULL)
  }

  counts <- tryCatch(
    {
      Signac::CountFragments(
        fragments = as.list(fragment_paths),
        cells = colnames(srt),
        verbose = verbose
      )
    },
    error = function(...) NULL
  )
  if (is.null(counts) || nrow(counts) == 0) {
    return(NULL)
  }

  rownames(counts) <- counts$CB
  counts$CB <- NULL
  counts <- counts[colnames(srt), , drop = FALSE]
  if ("reads_count" %in% colnames(counts) && !"total.fragments" %in% colnames(counts)) {
    counts$total.fragments <- counts$reads_count
  }
  counts
}

atac_frag_paths <- function(srt, assay) {
  frags <- tryCatch(
    Signac::Fragments(srt[[assay]]),
    error = function(...) list()
  )
  if (length(frags) == 0) {
    return(character(0))
  }
  vapply(
    frags,
    FUN = function(fragment) Signac::GetFragmentData(fragment, slot = "path"),
    FUN.VALUE = character(1)
  )
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

atac_has_annotation <- function(chrom_assay) {
  annotation <- tryCatch(
    Signac::Annotation(chrom_assay),
    error = function(...) NULL
  )
  !is.null(annotation) && length(annotation) > 0
}

atac_total_fragments_col <- function(srt, assay) {
  candidates <- c(
    "passed_filters",
    paste0(assay, "_passed_filters"),
    "total.fragments",
    "total_fragments",
    paste0(assay, "_total.fragments"),
    paste0(assay, "_total_fragments")
  )
  candidates <- candidates[candidates %in% colnames(srt@meta.data)]
  if (length(candidates) > 0) {
    return(candidates[1])
  }
  NULL
}

atac_frip_meta <- function(srt, assay, total_fragments_col = NULL) {
  peak_fragment_candidates <- c(
    "peak_region_fragments",
    paste0(assay, "_peak_region_fragments")
  )
  peak_fragment_candidates <- peak_fragment_candidates[
    peak_fragment_candidates %in% colnames(srt@meta.data)
  ]
  if (length(peak_fragment_candidates) == 0 || is.null(total_fragments_col)) {
    return(NULL)
  }

  total_fragments <- srt[[total_fragments_col, drop = TRUE]]
  peak_fragments <- srt[[peak_fragment_candidates[1], drop = TRUE]]
  frip <- ifelse(total_fragments > 0, peak_fragments / total_fragments, NA_real_)
  names(frip) <- colnames(srt)
  data.frame(
    pct_reads_in_peaks = frip,
    row.names = colnames(srt)
  )
}
