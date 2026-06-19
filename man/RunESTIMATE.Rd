% Generated manually for ESTIMATE integration.
\name{RunESTIMATE}
\alias{RunESTIMATE}
\title{Run ESTIMATE tumor microenvironment scoring}
\usage{
RunESTIMATE(
  object = NULL,
  count_matrix = NULL,
  assay = NULL,
  layer = "data",
  bulk_assay = "counts",
  sample.by = NULL,
  group.by = NULL,
  aggregate_fun = c("mean", "sum"),
  platform = c("rnaseq", "affymetrix", "agilent", "illumina"),
  filter_common_genes = TRUE,
  min_sig_genes = 10,
  purity = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{object}{Optional expression matrix, \code{SummarizedExperiment}, or \code{Seurat} object.}
\item{count_matrix}{Optional genes-by-samples expression matrix.}
\item{assay}{Assay used for \code{Seurat} input.}
\item{layer}{Assay layer used for \code{Seurat} input.}
\item{bulk_assay}{Assay name used for \code{SummarizedExperiment} input.}
\item{sample.by}{Metadata column used to aggregate \code{Seurat} cells into pseudo-bulk samples.}
\item{group.by}{Optional grouping metadata column. If \code{sample.by} is absent, this is also used for pseudo-bulk aggregation.}
\item{aggregate_fun}{Pseudo-bulk aggregation function for \code{Seurat} input.}
\item{platform}{Platform label stored in result metadata.}
\item{filter_common_genes}{Whether to restrict expression to the ESTIMATE common-gene universe.}
\item{min_sig_genes}{Minimum observed stromal and immune signature genes required after filtering.}
\item{purity}{Whether to calculate ESTIMATE tumor purity.}
\item{verbose}{Whether to print progress messages.}
}
\value{
A result bundle for matrix input, or the modified \code{SummarizedExperiment} or \code{Seurat} object.
}
\description{
Compute ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity scores from bulk or pseudo-bulk expression data without requiring the external \code{estimate} package.
}
\details{
\code{ESTIMATEScore} is \code{StromalScore + ImmuneScore}. \code{TumorPurity} is calculated as \code{cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)}. The purity formula was calibrated in the original ESTIMATE work primarily for Affymetrix data.
}
\references{
Yoshihara et al. (2013) \doi{10.1038/ncomms3612}.
Barbie et al. (2009) \doi{10.1038/nature08460}.
}
