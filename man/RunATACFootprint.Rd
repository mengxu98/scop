% Generated documentation for ATAC extensions
\name{RunATACFootprint}
\alias{RunATACFootprint}
\title{Bias-corrected ATAC footprinting}
\usage{
RunATACFootprint(object, genome, genome.id, bias.region, motifs,
  assay = NULL, upstream = 250L, downstream = 250L,
  in.peaks = FALSE, verbose = TRUE)
}
\arguments{
\item{object}{Seurat object with a ChromatinAssay.}
\item{genome}{Loaded BSgenome or indexed FaFile; no reference is downloaded.}
\item{genome.id}{Verified reference label matching the assay.}
\item{bias.region}{Explicit chromosome-start-end interval with three-base boundary margins.}
\item{motifs}{Unique motif IDs with stored positions.}
\item{assay}{Chromatin assay.}
\item{upstream, downstream}{Integer flanking widths, at least 50 bases.}
\item{in.peaks}{Restrict motif positions to peaks.}
\item{verbose}{Print progress.}
}
\value{Seurat object with native bias and footprints, plus assay Misc scop_footprints receipts.}
\description{Bias-corrected ATAC footprinting.}
\details{Native Signac bias is recalculated on an explicit region. Existing footprint keys are not overwritten. Adequate hexamer coverage is required; native failures are preserved. Footprints are insertion patterns, not proof of TF binding. Recompute after changing reference, cells, peaks or fragments.}
