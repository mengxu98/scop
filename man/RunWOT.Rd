% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunWOT.R
\name{RunWOT}
\alias{RunWOT}
\title{Run WOT analysis}
\usage{
RunWOT(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group_by = NULL,
  time_field = "Time",
  growth_iters = 3L,
  tmap_out = "tmaps/tmap_out",
  time_from = NULL,
  time_to = NULL,
  get_coupling = FALSE,
  recalculate = FALSE,
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

\item{time_field}{A character string specifying the column name in `adata.obs` or `srt@meta.data` that contains the time information.}

\item{growth_iters}{An integer specifying the number of growth iterations to perform during the OT Model computation. Default is 3.}

\item{tmap_out}{A character string specifying the path to store the computed transport maps.}

\item{time_from}{A numeric value specifying the starting time point for trajectory and fate analysis.}

\item{time_to}{A numeric value specifying the ending time point for trajectory and fate analysis. If not provided, only trajectory and fate analysis for the specified `time_from` will be performed.}

\item{get_coupling}{A logical value indicating whether to compute and store the coupling matrix between the specified `time_from` and `time_to`. Default is FALSE.}

\item{recalculate}{A logical value indicating whether to recalculate the transport maps even if they already exist at the specified `tmap_out` location. Default is FALSE.}

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
Run WOT analysis
}
\examples{
\dontrun{
data(pancreas_sub)
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

print(range(pancreas_sub$Lineage1, na.rm = TRUE))

pancreas_sub <- RunWOT(
  pancreas_sub,
  group_by = "SubCellType",
  time_field = "Lineage1",
  time_from = min(pancreas_sub$Lineage1, na.rm = TRUE),
  time_to = max(pancreas_sub$Lineage1, na.rm = TRUE),
  get_coupling = TRUE,
  tmap_out = "tmaps/lineage_tmap"
)

pancreas_sub$Custom_Time <- sample(
  1:10,
  ncol(pancreas_sub),
  replace = TRUE
)
pancreas_sub <- RunWOT(
  pancreas_sub,
  group_by = "CellType",
  time_field = "Custom_Time",
  time_from = 1,
  time_to = 10,
  tmap_out = "tmaps/custom_tmap"
)
}
}
\seealso{
\code{\link{srt_to_adata}}
}
