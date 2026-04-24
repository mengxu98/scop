% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/RunProportionTest.R
\name{RunProportionTestPermutation}
\alias{RunProportionTestPermutation}
\title{Permutation-based proportion test}
\usage{
RunProportionTestPermutation(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  n_permutations = 1000,
  include_all_cells = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{srt}{A Seurat object.}
\item{group.by}{Cell-type metadata column.}
\item{split.by}{Condition metadata column.}
\item{comparison}{Optional comparisons to perform.}
\item{n_permutations}{Number of permutations.}
\item{include_all_cells}{Whether to include all cell types in grid completion.}
\item{verbose}{Whether to print messages.}
}
\description{
Method-specific permutation implementation used by \link{RunProportionTest}.
This is a permutation-based statistical test for compositional shift
(not a separate biological model). Bootstrap confidence intervals are included
as uncertainty estimates for observed log2 fold-differences.
}
\seealso{\link{RunProportionTest}}
