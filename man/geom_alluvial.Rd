% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggsankey.R
\name{geom_alluvial}
\alias{geom_alluvial}
\title{geom_alluvial}
\usage{
geom_alluvial(
  mapping = NULL,
  data = NULL,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  space = 0,
  width = 0.1,
  smooth = 8,
  inherit.aes = TRUE,
  ...
)
}
\arguments{
\item{mapping}{provide you own mapping. both x and y need to be numeric.}

\item{data}{provide you own data}

\item{position}{change position}

\item{na.rm}{remove missing values}

\item{show.legend}{show legend in plot}

\item{space}{space between nodes in the y-direction}

\item{width}{width of nodes}

\item{smooth}{how much smooth should the curve have? More means steeper curve.}

\item{inherit.aes}{should the geom inherits aestethics}

\item{...}{other arguments to be passed to the geo}
}
\value{
ggplot layer
}
\description{
Creates an alluvial plot which visualize flows between nodes.
Each observation needs to have a `x` aesthetic as well as a `next_x` column which declares where that observation should flow.
Also each observation should have a `node` and a `next_node` aesthetic which provide information about which group in the y-direction.
}
