% Please edit documentation in R/DataConvert.R
\name{loom_to_srt}
\alias{loom_to_srt}
\title{Read a \code{.loom} file and convert to a \code{Seurat}}
\usage{
loom_to_srt(
  path,
  layers = c("spliced", "unspliced"),
  verbose = TRUE,
  chunk_rows = 1000
)
}
\arguments{
\item{path}{Path to a \code{.loom} file.}

\item{layers}{Character vector of loom layers to import as additional Seurat
assays. Missing layers are skipped with a warning. Default is
\code{c("spliced", "unspliced")}.}

\item{verbose}{Whether to print messages.}

\item{chunk_rows}{Number of feature rows to read from each matrix dataset per
chunk. Larger values can be faster but use more memory.}
}
\value{
A \code{Seurat} object.
}
\description{
Read a \code{.loom} file and convert to a \code{Seurat}
}
\details{
This function reads loom/HDF5 files directly in R using \code{rhdf5}; it does not
call \code{reticulate}, \code{scanpy}, or \code{\link[=loom_to_adata]{loom_to_adata()}}. The loom \code{/matrix} dataset
is imported as the \code{RNA} assay. Requested \code{/layers/*} datasets are imported
as additional assays, which is useful for velocity-style loom files with
\code{spliced} and \code{unspliced} layers.
}
\examples{
\dontrun{
srt <- loom_to_srt("path/to/data.loom")
srt
}
}
\seealso{
\code{\link[=loom_to_adata]{loom_to_adata()}}, \code{\link[=adata_to_srt]{adata_to_srt()}}, \code{\link[=srt_to_adata]{srt_to_adata()}}
}
