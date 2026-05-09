% Please edit documentation in R/DataConvert.R
\name{loom_to_adata}
\alias{loom_to_adata}
\title{Read a \code{.loom} file as an AnnData object}
\usage{
loom_to_adata(path, verbose = TRUE, ...)
}
\arguments{
\item{path}{Path to a \code{.loom} file (passed to \code{scanpy.read_loom()}).}

\item{verbose}{Whether to print messages.}

\item{...}{Additional arguments passed to \code{scanpy.read_loom()}.}
}
\value{
A Python \code{anndata.AnnData} object.
}
\description{
Read a \code{.loom} file as an AnnData object
}
\details{
This is a Python-backed wrapper and requires \code{reticulate} plus a Python
environment with \code{scanpy} and \code{loompy} available. It is independent from
\code{\link[=loom_to_srt]{loom_to_srt()}}, which reads loom files directly in R without initializing
Python.
}
\examples{
\dontrun{
adata <- loom_to_adata("path/to/data.loom")
adata
}
}
\seealso{
\code{\link[=loom_to_srt]{loom_to_srt()}}, \code{\link[=adata_to_srt]{adata_to_srt()}}, \code{\link[=srt_to_adata]{srt_to_adata()}}
}
