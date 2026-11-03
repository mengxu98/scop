% Generated documentation for ATAC extensions
\name{RunATACCoaccessibility}
\alias{RunATACCoaccessibility}
\title{Infer coaccessibility using native Cicero}
\usage{
RunATACCoaccessibility(object, genome.sizes, genome.id, assay = NULL,
  reduction = "lsi", dims = 2:30, k = 50L, window = 500000L,
  sample.num = 100L, seed = 1L,
  result.name = "scop_coaccessibility", verbose = TRUE)
}
\arguments{
\item{object}{Seurat object with ATAC counts.}
\item{genome.sizes}{Named positive integer chromosome lengths for the explicit reference.}
\item{genome.id}{Verified reference label matching assay.}
\item{assay}{Chromatin assay.}
\item{reduction}{Existing embedding covering all assay cells.}
\item{dims}{Unique embedding dimensions.}
\item{k}{Cells per aggregate, at least 2 and less than cell count.}
\item{window}{Native window size greater than 250000.}
\item{sample.num}{Number of windows for distance estimation.}
\item{seed}{Local seed; caller RNG restored.}
\item{result.name}{Unused assay Misc key.}
\item{verbose}{Print native progress.}
}
\value{Seurat object with native connection table, aggregate membership and provenance in assay Misc. Native Links and source assays are unchanged.}
\description{Infer coaccessibility using native Cicero.}
\details{Requires official monocle3-compatible Cicero and monocle3. Native binary-accessibility aggregation, size-factor normalization and distance-penalized graphical modelling are used, without a correlation fallback. Native NA scores and model errors are preserved. Small contigs or inadequate peaks may not support estimation. These associations do not establish physical interactions or causal regulation.}
