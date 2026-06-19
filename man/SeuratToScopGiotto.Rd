\name{SeuratToScopGiotto}
\alias{SeuratToScopGiotto}
\alias{CreateScopGiotto}
\alias{scop_giotto}
\alias{print.scop_giotto}
\alias{plot.scop_giotto}
\alias{GiottoPlot.scop_giotto}
\alias{RunGiottoWorkflow}
\alias{GiottoPreprocess}
\alias{GiottoReduce}
\alias{GiottoCluster}
\alias{GiottoSpatialNetwork}
\alias{GiottoSpatialGenes}
\alias{GiottoSpatialModules}
\alias{GiottoCellProximity}
\alias{GiottoHMRF}
\alias{AddGiottoToSeurat}
\title{scop-managed Giotto object workflow}
\usage{
SeuratToScopGiotto(
  srt,
  assay = NULL,
  layer = "counts",
  sct.assay = "SCT",
  use_sct = c("auto", "none", "normalized"),
  image = NULL,
  coord.cols = c("x", "y"),
  features = NULL,
  conversion_params = list(),
  use_official = TRUE,
  verbose = TRUE,
  seed = 11
)

CreateScopGiotto(...)

scop_giotto(
  giotto,
  source = list(),
  results = list(),
  active = NULL,
  history = list(),
  parameters = list()
)

RunGiottoWorkflow(
  x,
  steps = c("basic", "full"),
  group.by = NULL,
  verbose = TRUE,
  seed = 11,
  ...
)

GiottoPreprocess(
  x,
  filter_params = list(),
  norm_params = list(),
  stat_params = list(),
  hvf_params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoReduce(
  x,
  reduction = c("pca", "umap"),
  dims = 1:20,
  name = NULL,
  features = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoCluster(
  x,
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  network_name = "scop_NN",
  cluster_name = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoSpatialNetwork(
  x,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  params = list(),
  verbose = TRUE
)

GiottoSpatialGenes(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  top_n = 100,
  params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoSpatialModules(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  detect_params = list(),
  cluster_params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoCellProximity(
  x,
  group.by,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  number_of_simulations = 1000,
  adjust_method = "fdr",
  params = list(),
  verbose = TRUE,
  seed = 11
)

GiottoHMRF(
  x,
  spatial_genes = NULL,
  network_name = "Delaunay_full",
  k = 20,
  betas = c(0, 10, 20),
  hmrf_name = "scop_HMRF",
  params = list(),
  verbose = TRUE,
  seed = 11
)

AddGiottoToSeurat(
  srt,
  x,
  result = c("cluster", "hmrf"),
  name = NULL,
  tool_name = "Giotto",
  store_result = TRUE
)
}
\arguments{
\item{srt}{A Seurat object.}
\item{assay, layer}{Assay and layer used for Giotto input. Raw counts are recommended.}
\item{sct.assay}{Name of the SCT assay.}
\item{use_sct}{How SCT data are handled. The default keeps counts as the main Giotto input.}
\item{image, coord.cols}{Spatial image name or metadata coordinate columns.}
\item{features}{Features used for conversion or downstream Giotto methods.}
\item{conversion_params, filter_params, norm_params, stat_params, hvf_params, params, detect_params, cluster_params}{Named lists passed to Giotto functions.}
\item{use_official}{Whether to try \code{Giotto::seuratToGiottoV5()} before using the scop fallback converter.}
\item{verbose}{Whether to print progress messages.}
\item{seed}{Random seed passed to Giotto methods.}
\item{giotto}{A Giotto object.}
\item{source, results, active, history, parameters}{Components of a \code{scop_giotto} object.}
\item{x}{A \code{scop_giotto} object, or a Seurat object for \code{RunGiottoWorkflow()}.}
\item{steps}{Workflow preset: \code{"basic"} or \code{"full"}.}
\item{group.by}{Metadata group column for cell proximity enrichment.}
\item{reduction, dims, name, method, k, resolution, network_name, cluster_name, network_method, bin_method, top_n, cor_method, number_of_simulations, adjust_method, spatial_genes, betas, hmrf_name}{Giotto analysis parameters.}
\item{result}{Result to bridge or plot.}
\item{tool_name}{Name used when explicitly storing the Giotto object in Seurat tools.}
\item{store_result}{Whether \code{AddGiottoToSeurat()} stores the full \code{scop_giotto} object in \code{srt@tools}.}
\item{...}{Additional arguments passed through to the underlying converter or method.}
}
\description{
These functions implement a standalone Giotto object flow inside scop. Seurat and SCT data are used as inputs, but Giotto results are kept in a \code{scop_giotto} object and are not written back to Seurat unless \code{AddGiottoToSeurat()} is called explicitly.
}
\value{
A \code{scop_giotto} object for conversion and Giotto analysis functions; a Seurat object for \code{AddGiottoToSeurat()}.
}
\examples{
\dontrun{
# Convert Seurat/SCT data into a standalone Giotto workflow object.
g <- SeuratToScopGiotto(
  spatial_seurat,
  assay = "RNA",
  layer = "counts",
  use_sct = "normalized"
)

# Run the basic Giotto pipeline and plot directly from the Giotto object.
g <- RunGiottoWorkflow(g, steps = "basic")
GiottoPlot(g, plot_type = "cluster")
GiottoPlot(g, plot_type = "network")
GiottoPlot(g, plot_type = "dim")

# Add selected Giotto analyses without writing to the Seurat object.
g <- GiottoSpatialGenes(g, top_n = 50)
g <- GiottoCellProximity(g, group.by = "celltype", number_of_simulations = 100)
GiottoPlot(g, plot_type = "spatial_genes", top_n = 20)
GiottoPlot(g, plot_type = "cell_proximity")

# Export a result back to Seurat only when explicitly requested.
spatial_seurat <- AddGiottoToSeurat(
  spatial_seurat,
  g,
  result = "cluster",
  name = "Giotto_cluster"
)
}
}
