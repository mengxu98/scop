% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunSCVELO.R
\name{RunSCVELO}
\alias{RunSCVELO}
\title{Run scVelo workflow}
\usage{
RunSCVELO(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group_by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  mode = "stochastic",
  fitting_by = "stochastic",
  magic_impute = FALSE,
  knn = 5,
  t = 2,
  min_shared_counts = 30,
  n_pcs = 30,
  n_neighbors = 30,
  filter_genes = TRUE,
  min_counts = 3,
  min_counts_u = 3,
  normalize_per_cell = TRUE,
  log_transform = TRUE,
  use_raw = FALSE,
  diff_kinetics = FALSE,
  stream_smooth = NULL,
  stream_density = 2,
  arrow_length = 5,
  arrow_size = 5,
  arrow_density = 0.5,
  denoise = FALSE,
  denoise_topn = 3,
  kinetics = FALSE,
  kinetics_topn = 100,
  calculate_velocity_genes = FALSE,
  compute_velocity_confidence = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_paga = TRUE,
  top_n = 6,
  n_jobs = 1,
  palette = "Paired",
  palcolor = NULL,
  show_plot = TRUE,
  save = FALSE,
  dpi = 300,
  dirpath = "./",
  fileprefix = "",
  return_seurat = !is.null(srt)
)
}
\arguments{
\item{srt}{A Seurat object.}

\item{assay_x}{Assay to convert as the main data matrix (X) in the anndata object.}

\item{layer_x}{Layer name for assay_x in the Seurat object.}

\item{assay_y}{Assays to convert as layers in the anndata object.}

\item{layer_y}{Layer names for the assay_y in the Seurat object.}

\item{adata}{An anndata object.}

\item{group_by}{Variable to use for grouping cells in the Seurat object.}

\item{linear_reduction}{Linear reduction method to use, e.g., "PCA".}

\item{nonlinear_reduction}{Non-linear reduction method to use, e.g., "UMAP".}

\item{basis}{The basis to use for reduction, e.g., "UMAP".}

\item{mode}{Velocity estimation models to use. Can be a vector containing "deterministic", "stochastic", and/or "dynamical".}

\item{fitting_by}{Method used to fit gene velocities for dynamical modeling, e.g., "stochastic".}

\item{magic_impute}{Flag indicating whether to perform magic imputation.}

\item{knn}{The number of nearest neighbors for magic.MAGIC.}

\item{t}{power to which the diffusion operator is powered for magic.MAGIC.}

\item{min_shared_counts}{Minimum number of counts (both unspliced and spliced) required for a gene.}

\item{n_pcs}{Number of principal components (PCs) used for velocity estimation.}

\item{n_neighbors}{Number of nearest neighbors used for velocity estimation.}

\item{filter_genes}{Whether to filter genes based on minimum counts.}

\item{min_counts}{Minimum counts for gene filtering.}

\item{min_counts_u}{Minimum unspliced counts for gene filtering.}

\item{normalize_per_cell}{Whether to normalize counts per cell.}

\item{log_transform}{Whether to apply log transformation.}

\item{use_raw}{Whether to use raw data for dynamical modeling.}

\item{diff_kinetics}{Whether to use differential kinetics.}

\item{stream_smooth}{Multiplication factor for scale in Gaussian kernel around grid point.}

\item{stream_density}{Controls the closeness of streamlines.
When density = 2 (default), the domain is divided into a 60x60 grid,
whereas density linearly scales this grid.
Each cell in the grid can have, at most, one traversing streamline.}

\item{arrow_length}{Length of arrows.}

\item{arrow_size}{Size of arrows.}

\item{arrow_density}{Amount of velocities to show.}

\item{denoise}{Boolean flag indicating whether to denoise.}

\item{denoise_topn}{Number of genes with highest likelihood selected to infer velocity directions.}

\item{kinetics}{Boolean flag indicating whether to estimate RNA kinetics.}

\item{kinetics_topn}{Number of genes with highest likelihood selected to infer velocity directions.}

\item{calculate_velocity_genes}{Boolean flag indicating whether to calculate velocity genes.}

\item{compute_velocity_confidence}{Whether to compute velocity confidence metrics.}

\item{compute_terminal_states}{Whether to compute terminal states (root and end points).}

\item{compute_pseudotime}{Whether to compute velocity pseudotime.}

\item{compute_paga}{Whether to compute PAGA (Partition-based graph abstraction).}

\item{top_n}{The number of top features to plot.}

\item{n_jobs}{The number of parallel jobs to run.}

\item{palette}{The palette to use for coloring cells.}

\item{palcolor}{A vector of colors to use as the palette.}

\item{show_plot}{Whether to show the plot.}

\item{save}{Whether to save the plots.}

\item{dpi}{The DPI (dots per inch) for saving the plot.}

\item{dirpath}{The directory to save the plots.}

\item{fileprefix}{The file prefix to use for the plots.}

\item{return_seurat}{Whether to return a Seurat object instead of an anndata object. Default is TRUE.}
}
\description{
scVelo is a scalable toolkit for RNA velocity analysis in single cells.
This function runs an enhanced scVelo workflow on a Seurat object with improved
error handling, version compatibility, and modular design.
}
\examples{
\dontrun{
data(pancreas_sub)
pancreas_sub <- RunSCVELO(
  srt = pancreas_sub,
  assay_x = "RNA",
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)
head(pancreas_sub[[]])
names(pancreas_sub@assays)

FeatureDimPlot(
  pancreas_sub,
  c("stochastic_length", "stochastic_confidence")
)
FeatureDimPlot(
  pancreas_sub,
  "stochastic_pseudotime"
)
VelocityPlot(
  pancreas_sub,
  reduction = "UMAP",
  plot_type = "stream"
)
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = NA,
  velocity = "stochastic"
)

pancreas_sub <- RunSCVELO(
  pancreas_sub,
  assay_x = "RNA",
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  mode = c("deterministic", "stochastic"),
  filter_genes = TRUE,
  min_counts = 5,
  compute_velocity_confidence = TRUE,
  compute_terminal_states = TRUE,
  compute_pseudotime = TRUE,
  compute_paga = TRUE
)
}
}
\seealso{
\link{srt_to_adata} \link{VelocityPlot} \link{CellDimPlot} \link{RunPAGA}
}
