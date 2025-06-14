% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunHarmony2.R
\name{RunHarmony2}
\alias{RunHarmony2}
\alias{RunHarmony2.Seurat}
\title{Run Harmony algorithm}
\usage{
RunHarmony2(object, ...)

\method{RunHarmony2}{Seurat}(
  object,
  group.by.vars,
  reduction = "pca",
  dims.use = 1:30,
  project.dim = TRUE,
  reduction.name = "Harmony",
  reduction.key = "Harmony_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
}
\arguments{
\item{object}{A Seurat object.}

\item{...}{Additional arguments to be passed to the \link[harmony:RunHarmony]{harmony::RunHarmony} function.}

\item{group.by.vars}{A character vector specifying the batch variable name.}

\item{reduction}{A character string specifying the reduction to be used. Default is "pca".}

\item{dims.use}{An integer vector specifying the dimensions to be used. Default is 1:30.}

\item{project.dim}{A logical value indicating whether to project dimension reduction loadings. Default is TRUE.}

\item{reduction.name}{A character string specifying the name of the reduction to be stored in the Seurat object. Default is "Harmony".}

\item{reduction.key}{A character string specifying the prefix for the column names of the Harmony embeddings. Default is "Harmony_".}

\item{verbose}{A logical value indicating whether to print verbose output. Default is TRUE.}

\item{seed.use}{An integer specifying the random seed to be used. Default is 11.}
}
\description{
This is a modified version of harmony::RunHarmony specifically designed for compatibility with RunSymphonyMap.
}
\examples{
panc8_sub <- RunHarmony2(
  panc8_sub,
  group.by.vars = "tech",
  reduction = "pca"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "pca"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "Harmony"
)
}
