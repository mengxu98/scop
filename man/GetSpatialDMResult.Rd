\name{GetSpatialDMResult}
\alias{GetSpatialDMResult}
\title{Access stored SpatialDM results}
\usage{GetSpatialDMResult(object, result.name = NULL,
  type = c("global", "local", "weights"), pair = NULL)}
\arguments{
  \item{object}{A Seurat object with SpatialDM results.}
  \item{result.name}{Stored result name.}
  \item{type}{Result payload to return.}
  \item{pair}{Optional LR interaction name for local results.}
}
\value{The requested global table, local matrices, or spatial weight matrices.}
\seealso{\link{RunSpatialDM}, \link{SpatialDMPlot}}
