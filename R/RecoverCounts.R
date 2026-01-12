#' @title Attempt to recover raw counts from the normalized matrix
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @param trans The transformation function to applied when data is presumed to be log-normalized.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#' However, due to decimal point preservation during normalization,
#' the calculated nCount is usually a floating point number close to the integer.
#' The tolerance is its difference from the integer.
#' Default is `0.1`
#' @param sf Set the scaling factor manually.
#'
#' @seealso
#' [CheckDataType]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' raw_counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#'
#' # Normalized the data
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#'
#' # Now replace counts with the log-normalized data matrix
#' data <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "data"
#' )
#' new_pancreas_sub <- SeuratObject::SetAssayData(
#'   object = pancreas_sub,
#'   layer = "counts",
#'   new.data = data,
#'   assay = "RNA"
#' )
#' # Recover the counts and compare with the raw counts matrix
#' pancreas_sub <- RecoverCounts(new_pancreas_sub)
#' new_counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#' identical(raw_counts, new_counts)
RecoverCounts <- function(
    srt,
    assay = NULL,
    trans = c("expm1", "exp", "none"),
    min_count = c(1, 2, 3),
    tolerance = 0.1,
    sf = NULL,
    verbose = TRUE) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  counts <- GetAssayData5(
    srt,
    layer = "counts",
    assay = assay
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- SeuratObject::as.sparse(
      counts[seq_len(nrow(counts)), , drop = FALSE]
    )
  }

  status <- CheckDataType(counts)
  if (status == "raw_counts") {
    log_message(
      "The data is raw counts",
      verbose = verbose
    )
    return(srt)
  }

  if (status == "log_normalized_counts") {
    log_message(
      "The data is presumed to be log-normalized",
      verbose = verbose
    )
    trans <- match.arg(trans)
    if (trans %in% c("expm1", "exp")) {
      log_message(
        "Perform {.val {trans}} on the raw data",
        verbose = verbose
      )
      counts <- do.call(trans, list(counts))
    }
  }
  if (status == "raw_normalized_counts") {
    log_message(
      "The data is presumed to be normalized without log transformation",
      verbose = verbose
    )
  }
  if (is.null(sf)) {
    sf <- unique(round(Matrix::colSums(counts)))
  }

  if (length(sf) == 1) {
    counts <- counts / sf
    elements <- split(
      counts@x,
      rep(seq_len(ncol(counts)), diff(counts@p))
    )
    min_norm <- sapply(elements, min)
    nCount <- NULL
    for (m in min_count) {
      if (is.null(nCount)) {
        presumed_nCount <- m / min_norm
        diff_value <- abs(presumed_nCount - round(presumed_nCount))
        if (max(diff_value, na.rm = TRUE) < tolerance) {
          nCount <- round(presumed_nCount)
        }
      }
    }

    if (is.null(nCount)) {
      top10_cells <- utils::head(colnames(counts)[diff_value < tolerance], 10)
      log_message(
        "The presumed nCount of some cells is not valid: {.val {top10_cells}}",
        message_type = "warning",
        verbose = verbose
      )
      return(srt)
    }
    counts@x <- round(counts@x * rep(nCount, diff(counts@p)))
    srt <- SeuratObject::SetAssayData(
      srt,
      new.data = counts,
      assay = assay,
      layer = "counts"
    )

    nCount <- stats::setNames(
      rep(nCount, length.out = ncol(srt)),
      colnames(srt)
    )
    srt[[paste0("nCount_", assay)]] <- nCount
  } else {
    log_message(
      "Scale factor is not unique. No changes to be made",
      message_type = "warning",
      verbose = verbose
    )
  }

  return(srt)
}
