% Generated documentation for ATAC extensions
\name{ATACMotifPlot}
\alias{ATACMotifPlot}
\title{Plot motif sequence logos}
\usage{
ATACMotifPlot(object, motifs, assay = NULL, use.names = TRUE, ...)
}
\arguments{
\item{object}{Seurat object with motifs.}
\item{motifs}{Unique motif IDs in the selected assay PWM list.}
\item{assay}{Chromatin assay.}
\item{use.names}{Use native motif display names.}
\item{...}{Arguments passed to Signac::MotifPlot and ggseqlogo.}
}
\value{Native motif sequence logo plot.}
\description{Plot motif sequence logos.}
\details{Uses exact assay-specific IDs rather than ambiguous common names. Requires ggseqlogo.}
