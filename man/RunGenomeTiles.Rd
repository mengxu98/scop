\name{RunGenomeTiles}
\alias{RunGenomeTiles}
\title{Count fragments in fixed-width genome tiles}
\usage{
RunGenomeTiles(object, assay = NULL, seqlengths = NULL, chromosomes = NULL,
  bin.size = 5000L, cells = NULL, new.assay = NULL, process_n = 2000L,
  max.tiles = 2000000L, verbose = TRUE)
}
\arguments{
\item{object}{Seurat object with validated ChromatinAssay and fragments.}
\item{assay}{Source assay; NULL selects active assay.}
\item{seqlengths}{Named positive integer reference lengths, or NULL to use known assay lengths. Supplied lengths must agree with known source lengths.}
\item{chromosomes}{Optional unique contigs in output order.}
\item{bin.size}{Positive integer tile width; last tile per chromosome is truncated.}
\item{cells}{Source-assay cells; see RecountATAC.}
\item{new.assay}{Unused assay name for writeback; NULL returns a list.}
\item{process_n}{Native FeatureMatrix chunk size.}
\item{max.tiles}{Positive upper bound on tile count before allocation.}
\item{verbose}{Whether to print progress messages.}
}
\description{Generate non-overlapping genome tiles and delegate sparse fragment counting to RecountATAC. No peak universe is required to define tiles; the source assay provides fragment mappings and reference provenance.}
\value{RecountATAC's return value with tile parameters in its receipt. Counts are overlapping fragment records, not insertion endpoints or unique fragments across bins. A fragment crossing a boundary can count in both bins.}
\seealso{\code{\link{RecountATAC}}}
