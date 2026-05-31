\name{ScFEABalanceBarPlot}
\alias{ScFEABalanceBarPlot}
\title{Plot scFEA metabolite balance changes}
\usage{
ScFEABalanceBarPlot(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAbalance",
  layer = "data",
  top_n = NULL,
  p_adj_cutoff = 0.01,
  title = NULL
)
}
\arguments{
\item{srt}{A Seurat object returned by \code{\link[=RunScFEA]{RunScFEA()}}.}
\item{group.by}{Metadata column defining groups. If \code{ident.1} and \code{ident.2} are both \code{NULL}, each group is compared against all remaining cells.}
\item{ident.1, ident.2}{Group names to compare. Difference is \code{mean(ident.1) - mean(ident.2)}. If both are \code{NULL}, one-vs-rest contrasts are generated automatically.}
\item{assay}{Balance assay name.}
\item{layer}{Balance assay layer.}
\item{top_n}{If \code{NULL}, plot all metabolites. If positive, plot the top \code{top_n} increased and top \code{top_n} decreased metabolites after the adjusted-p-value filter.}
\item{p_adj_cutoff}{Adjusted p-value cutoff used when \code{top_n} is not \code{NULL}.}
\item{title}{Optional plot title.}
}
\value{
For a single contrast, a list with \code{plot} and plotted \code{data}. For automatic one-vs-rest mode, a named list of single-contrast results is returned and combined statistics are stored in the outer \code{"data"} attribute.
}
\description{
Computes scFEA balance mean differences between two Seurat metadata groups and draws a metabolite bar plot.
}
