% Generated documentation for ATAC extensions
\name{ATACFootprintPlot}
\alias{ATACFootprintPlot}
\title{Plot bias-corrected footprints}
\usage{
ATACFootprintPlot(object, motifs, assay = NULL, group.by = NULL,
  normalization = "subtract", show.expected = TRUE, ...)
}
\arguments{
\item{object}{Seurat object returned by RunATACFootprint.}
\item{motifs}{Unique stored footprint IDs.}
\item{assay}{Chromatin assay.}
\item{group.by}{Complete cell metadata column, or identities when NULL.}
\item{normalization}{Native subtract, divide, or NULL normalization.}
\item{show.expected}{Show expected insertion profile.}
\item{...}{Additional Signac::PlotFootprint arguments.}
}
\value{Native ggplot or patchwork object.}
\description{Plot bias-corrected footprints.}
\details{Missing or stale cell/peak receipts are rejected. This plot does not establish TF occupancy.}
