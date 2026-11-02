#' Plot motif sequence logos from a selected ATAC assay
#' @param object A Seurat object with native motifs.
#' @param motifs Unique motif IDs (not ambiguous common names).
#' @param assay Chromatin assay.
#' @param use.names Use native motif names in plot labels.
#' @param ... Arguments passed to Signac::MotifPlot and ggseqlogo.
#' @return Native motif-logo plot. Motif similarity is not TF binding evidence.
#' @export
ATACMotifPlot <- function(object, motifs, assay = NULL, use.names = TRUE, ...) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  assay <- report$assay
  pwm <- Signac::GetMotifData(object, assay = assay, slot = "pwm")
  if (!atac_unique_ids(motifs) || !all(motifs %in% names(pwm))) stop("motifs must be unique motif IDs in the selected assay PWM list.", call. = FALSE)
  atac_backend_fun("ggseqlogo", "ggseqlogo")
  atac_backend_fun("Signac", "MotifPlot")(object, motifs = motifs, assay = assay, use.names = use.names, ...)
}
