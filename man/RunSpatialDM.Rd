\name{RunSpatialDM}
\alias{RunSpatialDM}
\title{Run official SpatialDM spatial ligand-receptor association analysis}
\usage{
RunSpatialDM(srt, species = c("human", "mouse"), lr.database = NULL,
  assay = NULL, layer = "data", counts.layer = "counts", image = NULL,
  coord.cols = c("col", "row"), coordinate.unit = NULL,
  method = c("z-score", "permutation"), run.local = TRUE, l = NULL,
  eff_dist = NULL, cutoff = 0.1, n_neighbors = NULL,
  n_neighbor_layers = 6, single_cell = FALSE,
  complex.mean = c("algebra", "geometric"), min_cell = 3, n_perm = 1000,
  nproc = 1, seed = 1, global.fdr = TRUE, global.threshold = 0.1,
  local.fdr = FALSE, local.threshold = 0.1, result.name = "default",
  envname = "spatialdm_env", overwrite = FALSE, verbose = TRUE)
}
\arguments{
  \item{srt}{A Seurat spatial object.}
  \item{species}{Built-in SpatialDM species.}
  \item{lr.database}{Optional custom LR data frame with \code{ligand}, \code{receptor}, and \code{annotation}.}
  \item{assay,layer,counts.layer}{Assay and normalized/count expression layers.}
  \item{image,coord.cols}{Spatial image or metadata coordinate columns.}
  \item{coordinate.unit}{Documented unit for raw coordinates.}
  \item{method}{Global/local test mode.}
  \item{run.local}{Whether to run local spot statistics.}
  \item{l,eff_dist}{RBF scale parameters; one must be supplied.}
  \item{cutoff,n_neighbors,n_neighbor_layers,single_cell}{SpatialDM weight parameters.}
  \item{complex.mean,min_cell,n_perm,nproc,seed}{Backend and reproducibility parameters.}
  \item{global.fdr,global.threshold,local.fdr,local.threshold}{Selection rules.}
  \item{result.name,envname,overwrite,verbose}{Result, environment, overwrite, and logging controls.}
}
\value{The input Seurat object with a schema-v1 SpatialDM result bundle in \code{srt@tools$SpatialDM}.}
\description{Runs the pinned optional Python SpatialDM backend using raw analysis coordinates and normalized/log-transformed expression.}
\references{Li et al. (2023), Nature Communications, 14, 3995.}
\seealso{\link{SpatialDMPlot}, \link{GetSpatialDMResult}}
