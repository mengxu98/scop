\name{ATACRegionHeatmap}
\alias{ATACRegionHeatmap}
\alias{ATACRegionPlot}
\title{Plot positional ATAC region signal}
\usage{
ATACRegionHeatmap(x, normalize = TRUE, per.base = TRUE, groups = NULL,
  log1p = FALSE, order.rows = FALSE, low = "white", high = "#2166AC")
ATACRegionPlot(x, normalize = TRUE, per.base = TRUE, groups = NULL)
}
\arguments{
\item{x}{Output of RunATACRegionMatrix.}
\item{normalize}{Divide by selected cells per group.}
\item{per.base}{Divide by actual bin widths, including the shorter last bin.}
\item{groups}{Unique group names in display order; NULL includes all.}
\item{log1p}{Apply log1p to displayed heatmap values.}
\item{order.rows}{Order all panels by total signal across selected groups.}
\item{low,high}{Heatmap color endpoints.}
}
\description{Display raw or explicitly normalized endpoint signal from ATAC regions.}
\value{A ggplot with summarized plotting data. The line plot averages across regions, not across independent biological replicates; it is not a confidence interval or enrichment test.}
\seealso{\code{\link{RunATACRegionMatrix}}}
