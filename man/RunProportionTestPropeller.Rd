% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/RunProportionTest.R
\name{RunProportionTestPropeller}
\alias{RunProportionTestPropeller}
\title{Propeller differential abundance wrapper}
\usage{
RunProportionTestPropeller(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  n_bootstrap = 1000,
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
\item{n_bootstrap}{Bootstrap iterations for confidence intervals.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print messages.}
}
\description{
Method-specific propeller-style wrapper used by \link{RunProportionTest}.
}
\seealso{\link{RunProportionTest}}
