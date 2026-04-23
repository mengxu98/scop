% Generated manually because roxygen2 is unavailable in this environment
% Please edit documentation in R/RunProportionTest.R
\name{RunProportionTest}
\alias{RunProportionTest}
\title{Proportion Test}
\usage{
RunProportionTest(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  proportion_method,
  sample.by = NULL,
  n_permutations = 1000,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  include_all_cells = FALSE,
  seed = 11,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{srt}{A Seurat object.}

\item{group.by}{Name of a metadata column for cell type/group labels.}

\item{split.by}{Name of a metadata column defining conditions to compare.}

\item{comparison}{Optional comparisons to perform. Supports
\code{list(c("A", "B"))} or character values like \code{"A_vs_B"}.}

\item{proportion_method}{Differential abundance method. Canonical values are
\code{"permutation"}, \code{"milo"}, \code{"sccoda"}, and
\code{"propeller"}. Alias values (for example \code{"permutation_test"} or
\code{"perm"}) are accepted and normalized to \code{"permutation"}.}

\item{sample.by}{Metadata column for biological sample IDs. Required for
\code{"milo"}, \code{"sccoda"}, and \code{"propeller"}.}

\item{n_permutations}{Number of permutations for permutation mode.}

\item{FDR_threshold}{FDR value cutoff for significance.}

\item{log2FD_threshold}{Absolute value of log2FD cutoff for significance.}

\item{include_all_cells}{Whether to include all cell types in permutation complete grid.}

\item{seed}{Random seed.}

\item{verbose}{Whether to print messages.}

\item{...}{Additional arguments passed to the selected method function.}
}
\description{
\link{RunProportionTest} performs differential abundance testing for cell proportions
through a unified dispatcher.
}
\examples{
data(pancreas_sub)
pancreas_sub <- RunProportionTest(
  pancreas_sub,
  group.by = "CellType",
  split.by = "Phase",
  proportion_method = "permutation",
  comparison = list(c("G2M", "G1"))
)

ProportionTestPlot(pancreas_sub)
}
\seealso{
\link{RunProportionTestPermutation},
\link{RunProportionTestMilo},
\link{RunProportionTestScCODA},
\link{RunProportionTestPropeller},
\link{ProportionTestPlot}
}

