\name{scFEAVolcanoPlot}
\alias{scFEAVolcanoPlot}
\title{Plot scFEA flux Cohen's d volcano plots}
\usage{
scFEAVolcanoPlot(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAflux",
  layer = "data",
  label_by = c("reaction", "module", "module_reaction"),
  pathways = NULL,
  p_adj_cutoff = 0.001,
  cohen_cutoff = 0.2,
  combine = TRUE,
  width = 12,
  height = 10.4,
  verbose = TRUE
)
}
\arguments{
\item{srt}{A Seurat object returned by \code{\link[=RunscFEA]{RunscFEA()}}.}
\item{group.by}{Metadata column defining groups. If \code{ident.1} and \code{ident.2} are both \code{NULL}, each group is compared against all remaining cells.}
\item{ident.1, ident.2}{Group names to compare. Cohen's d is \code{mean(ident.1) - mean(ident.2)} divided by pooled standard deviation. If both are \code{NULL}, one-vs-rest contrasts are generated automatically.}
\item{assay}{Flux assay name.}
\item{layer}{Flux assay layer.}
\item{label_by}{Label mode for significant pathway modules.}
\item{pathways}{Optional \code{SM_anno} classes to plot.}
\item{p_adj_cutoff}{Adjusted p-value cutoff for highlighting, labels, and the horizontal threshold line.}
\item{cohen_cutoff}{Absolute Cohen's d cutoff for highlighting, labels, and the vertical threshold lines.}
\item{combine}{Whether to combine pathway plots into 2 x 2 paged panels. If \code{FALSE}, returns one plot per pathway.}
\item{width, height}{Suggested export size. Stored as attributes on the returned plot object.}
\item{verbose}{Whether to print messages.}
}
\value{
For a single contrast, a paged list of \pkg{ggplot2} plots when \code{combine = TRUE}, otherwise a named list of pathway \pkg{ggplot2} plots. Statistics are stored in the \code{"data"} attribute. For automatic one-vs-rest mode, a named list of single-contrast results is returned and combined statistics are stored in the outer \code{"data"} attribute.
}
\description{
Compares scFEA module flux between two Seurat metadata groups, computes Wilcoxon adjusted p-values and Cohen's d, and draws one volcano plot per M168 pathway class. All modules are shown as grey background points; modules in the current pathway are colored red when Cohen's d is positive and blue when negative; threshold-passing modules in the current pathway are overlaid as larger green labeled points.
}
