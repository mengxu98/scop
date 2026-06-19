% Generated manually for ESTIMATE integration.
\name{EstimateScorePlot}
\alias{EstimateScorePlot}
\title{ESTIMATE score plots}
\usage{
EstimateScorePlot(
  object = NULL,
  score.data = NULL,
  group.by = NULL,
  group.data = NULL,
  plot_type = c("violin", "box", "heatmap", "cor"),
  scores = c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
  add_stat = TRUE,
  ...
)
}
\arguments{
\item{object}{Optional \code{RunESTIMATE()} bundle, \code{SummarizedExperiment}, or \code{Seurat} object containing ESTIMATE results.}
\item{score.data}{Optional score matrix or data frame with samples in rows.}
\item{group.by}{Optional grouping column for grouped violin and box plots.}
\item{group.data}{Optional named vector or data frame containing sample groups.}
\item{plot_type}{Plot type.}
\item{scores}{ESTIMATE score columns to plot.}
\item{add_stat}{Whether to add group comparison labels to violin or box plots.}
\item{...}{Additional plotting arguments.}
}
\value{
A \code{ggplot} object.
}
\description{
Visualize ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity scores from a \code{RunESTIMATE()} result.
}
