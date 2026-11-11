#' Plot motif sequence logos from a selected ATAC assay
#' @param object A Seurat object with motifs.
#' @param motifs Unique motif IDs (not ambiguous common names).
#' @param assay Chromatin assay.
#' @param use.names Use motif names in plot labels.
#' @param ... Arguments passed to other functions.
#' @return A motif-logo plot. Motif similarity is not TF binding evidence.
#' @details Requires the optional ggseqlogo package.
#' @export
ATACMotifPlot <- function(object, motifs, assay = NULL, use.names = TRUE, ...) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  assay <- report$assay
  pwm <- Signac::GetMotifData(object, assay = assay, slot = "pwm")
  if (!atac_unique_ids(motifs) || !all(motifs %in% names(pwm))) stop("motifs must be unique motif IDs in the selected assay PWM list.", call. = FALSE)
  check_r("ggseqlogo", verbose = FALSE)
  Signac::MotifPlot(object, motifs = motifs, assay = assay, use.names = use.names, ...)
}
