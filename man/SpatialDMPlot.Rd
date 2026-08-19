\name{SpatialDMPlot}
\alias{SpatialDMPlot}
\title{Plot stored SpatialDM results}
\usage{SpatialDMPlot(object, result.name = NULL,
  plot_type = c("weights", "global", "local"), pair = NULL,
  spot = NULL, signaling = c("secreted", "contact"), highlight = NULL,
  palette = "RdBu", palcolor = NULL, theme_use = "theme_scop",
  theme_args = list(), ...)}
\arguments{
  \item{object}{A Seurat object with SpatialDM results.}
  \item{result.name}{Stored result name.}
  \item{plot_type}{Plot type.}
  \item{pair}{LR interaction name for local plotting.}
  \item{spot}{Spot identifier for weight plotting.}
  \item{signaling}{Stored weight matrix to use.}
  \item{highlight}{LR interactions to label in the global plot.}
  \item{palette,palcolor,theme_use,theme_args}{SCOP plotting controls.}
  \item{...}{Additional arguments passed to SCOP spatial plotting functions.}
}
\value{A ggplot or patchwork object.}
\seealso{\link{RunSpatialDM}, \link{GetSpatialDMResult}}
