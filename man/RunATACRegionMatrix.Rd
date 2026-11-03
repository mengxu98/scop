\name{RunATACRegionMatrix}
\alias{RunATACRegionMatrix}
\title{Count strand-oriented insertion endpoints around ATAC regions}
\usage{
RunATACRegionMatrix(object, regions, assay = NULL, group.by = NULL,
  cells = NULL, upstream = 3000L, downstream = 3000L, bin.size = 50L,
  max.entries = 5e7, verbose = TRUE)
}
\arguments{
\item{object}{Seurat object with a validated ChromatinAssay, local indexed fragments and explicit cell-to-barcode mappings.}
\item{regions}{Non-empty GRanges in 1-based closed coordinates. The center is start + floor((width - 1)/2). Minus-strand regions are reversed. Genomic intervals and names, if supplied, must be unique.}
\item{assay}{Chromatin assay; NULL selects the active assay.}
\item{group.by}{Metadata grouping column; NULL uses active identities.}
\item{cells}{Unique source-assay cells; NULL includes all.}
\item{upstream,downstream}{Non-negative integer flank widths. Windows must fit inside known chromosome boundaries.}
\item{bin.size}{Positive integer positional bin width. A shorter last bin is retained.}
\item{max.entries}{Dense-output entry guard across all groups. Errors rather than sampling when exceeded; transient single-query memory is not capped.}
\item{verbose}{Whether to print progress messages.}
}
\description{Build group-specific region-by-position matrices from indexed fragments without subsampling. This explicit SCOP counter is not Signac RegionMatrix. Repeated endpoints are counted with multiplicity.}
\details{A BED fragment [start,end) contributes endpoints start + 1 and end. Only endpoints inside each window count. Column 5 is not a weight; no extra Tn5 shift, bias correction or depth normalization is applied. Overlapping windows can count the same endpoint. The source object is unchanged. Group signal is descriptive, not a replicate-level statistical test or evidence of binding.}
\value{A scop_atac_regions list with raw matrices, regions, bin positions and widths, selected group cell counts and a receipt.}
\seealso{\code{\link{ATACRegionHeatmap}}, \code{\link{ATACRegionPlot}}}
