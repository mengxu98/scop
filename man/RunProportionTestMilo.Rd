% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/RunProportionTest.R
\name{RunProportionTestMilo}
\alias{RunProportionTestMilo}
\title{Milo differential abundance wrapper}
\usage{
RunProportionTestMilo(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  milo_k = 20L,
  milo_d = 30L,
  n_bootstrap = 500,
  seed = 11,
  verbose = TRUE
)
}
\arguments{
\item{srt}{A Seurat object.}
\item{group.by}{Cell-type metadata column.}
\item{split.by}{Condition metadata column.}
\item{sample.by}{Sample metadata column.}
\item{comparison}{Optional comparisons to perform.}
\item{milo_k}{k used for Milo graph building.}
\item{milo_d}{Number of dimensions used by Milo.}
\item{n_bootstrap}{Bootstrap iterations for summary/fallback.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print messages.}
}
\description{
Method-specific Milo wrapper used by \link{RunProportionTest}.
Stores both group-level and neighborhood-level result lists. When \pkg{miloR}
is available, graph-ready neighborhood node/edge data are also stored under
\code{details$milo_graph_data} for visualization.
}
\seealso{\link{RunProportionTest}}

