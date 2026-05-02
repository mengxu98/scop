\name{RunBulk}
\alias{RunBulk}
\title{Unified bulk analysis entrypoint}
\usage{
RunBulk(
  srt = NULL,
  bulk_se = NULL,
  method = "de_edgeR_qlf",
  sample.by = NULL,
  condition.by = NULL,
  group.by = NULL,
  ref_srt = NULL,
  celltype.by = NULL,
  assay = NULL,
  layer = "counts",
  bulk_assay = "counts",
  ref_assay = NULL,
  ref_layer = "counts",
  condition1 = NULL,
  condition2 = NULL,
  de_markers_type = "single",
  markers_type = NULL,
  run_enrichment = FALSE,
  run_gsea = FALSE,
  enrichment_args = list(),
  gsea_args = list(),
  de_args = list(),
  deconv_args = list(),
  csde_args = list(),
  method_args = list(),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{srt}{A Seurat object used for Seurat-centered bulk strategy workflows.}
\item{bulk_se}{A SummarizedExperiment object used for true bulk input.}
\item{method}{Character vector of methods to run. Supported values are \code{de_limma_voom}, \code{de_edgeR_qlf}, \code{de_DESeq2}, \code{de_dream}, \code{deconv_MuSiC}, \code{deconv_BisqueRNA}, \code{deconv_BayesPrism}, and \code{csde_TOAST}.}
\item{sample.by}{Metadata column containing sample IDs. Required for Seurat input.}
\item{condition.by}{Metadata column containing condition labels. Required for DE and CSDE methods.}
\item{group.by}{Optional metadata column used to stratify DE by subgroup.}
\item{ref_srt}{Optional Seurat reference for deconvolution and CSDE.}
\item{celltype.by}{Metadata column defining reference cell types.}
\item{assay, layer, bulk_assay, ref_assay, ref_layer}{Assay and layer controls for expression extraction.}
\item{condition1, condition2}{Condition labels for pairwise comparisons.}
\item{de_markers_type, markers_type}{DE comparison mode: \code{single}, \code{all}, or \code{paired}. \code{markers_type} is an alias.}
\item{run_enrichment, run_gsea}{Whether to run \code{\link{RunEnrichment}} or \code{\link{RunGSEA}} from bulk DE results.}
\item{enrichment_args, gsea_args}{Named lists of additional parameters forwarded to \code{\link{RunEnrichment}} and \code{\link{RunGSEA}}.}
\item{de_args, deconv_args, csde_args}{Module-level method parameters. Current deconvolution wrappers use \code{backend = "internal"} by default, and current \code{csde_TOAST} uses \code{backend = "limma_interaction"} by default. The computational engine is recorded in each method result bundle.}
\item{method_args}{Named list of method-specific parameters. Names should match RunBulk method names and override broader argument lists.}
\item{verbose}{Whether to print progress messages.}
\item{...}{Additional parameters forwarded to method functions.}
}
\value{
Returns the input object with results stored in \code{@tools[["Bulk"]]} for Seurat input or \code{metadata(bulk_se)[["Bulk"]]} for SummarizedExperiment input. Standard tables are available at \code{Bulk$results$de}, \code{Bulk$results$deconv}, and \code{Bulk$results$csde}. When \code{run_enrichment} or \code{run_gsea} succeeds in Seurat mode, results are also stored in \code{@tools[["Enrichment_Bulk_wilcox"]]} and \code{@tools[["GSEA_Bulk_wilcox"]]} for \code{\link{EnrichmentPlot}} and \code{\link{GSEAPlot}}.
}
\description{
\code{RunBulk()} provides one concise entrypoint for basic bulk DE, deconvolution, and CSDE workflows while keeping results in standard tables that can be passed to existing scop visualization and enrichment functions.
}
\examples{
\dontrun{
srt <- RunBulk(
  srt = srt,
  method = "de_edgeR_qlf",
  sample.by = "sample_id",
  condition.by = "condition",
  method_args = list(
    de_edgeR_qlf = list(logfc.threshold = 0.25, p.adjust.method = "BH")
  )
)
DEtestPlot(srt, group.by = "Bulk")
srt <- RunEnrichment(srt, group.by = "Bulk")
srt <- RunGSEA(srt, group.by = "Bulk")
GroupHeatmap(srt, features = srt@tools$Bulk$results$de$gene)
}
}
