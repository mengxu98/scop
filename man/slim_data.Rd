% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/themes.R
\name{slim_data}
\alias{slim_data}
\alias{slim_data.ggplot}
\alias{slim_data.patchwork}
\title{Drop unused data from the plot to reduce the object size}
\usage{
slim_data(p)

\method{slim_data}{ggplot}(p)

\method{slim_data}{patchwork}(p)
}
\arguments{
\item{p}{A \code{patchwork} object.}
}
\description{
Drop unused data from the plot to reduce the object size
}
\examples{
library(ggplot2)
p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
  geom_point()
object.size(p)
colnames(p$data)

p_slim <- slim_data(p)
object.size(p_slim)
colnames(p_slim$data)

}
