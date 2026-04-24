% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/RunProportionTest.R
\name{RunProportionTestScCODA}
\alias{RunProportionTestScCODA}
\title{scCODA differential abundance wrapper}
\usage{
RunProportionTestScCODA(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  reference_cell_type = NULL,
  credible_effect_threshold = 0.95,
  n_mcmc_samples = 20000L,
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
\item{reference_cell_type}{Optional scCODA reference cell type.}
\item{credible_effect_threshold}{Inclusion probability threshold for credible effects.}
\item{n_mcmc_samples}{Requested MCMC sample count.}
\item{n_bootstrap}{Fallback bootstrap iterations.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print messages.}
}
\description{
Method-specific scCODA wrapper used by \link{RunProportionTest}.
Calls Python helper via reticulate and falls back to sample-level summaries when unavailable.
}
\seealso{\link{RunProportionTest}}
