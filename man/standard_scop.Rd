% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_scop.R
\name{standard_scop}
\alias{standard_scop}
\title{Standard scop}
\usage{
standard_scop(
  srt,
  prefix = "Standard",
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_method = "vst",
  nHVF = 2000,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  seed = 11
)
}
\arguments{
\item{srt}{A Seurat object.}

\item{prefix}{A prefix to add to the names of intermediate objects created by the function (default is "Standard").}

\item{assay}{The name of the assay to use for the analysis. If NULL, the default assay of the Seurat object will be used.}

\item{do_normalization}{A logical value indicating whether to perform normalization. If NULL, normalization will be performed if the specified assay does not have scaled data.}

\item{normalization_method}{The method to use for normalization. Options are "LogNormalize", "SCT", or "TFIDF" (default is "LogNormalize").}

\item{do_HVF_finding}{A logical value indicating whether to perform high variable feature finding. If TRUE, the function will force to find the highly variable features (HVF) using the specified HVF method.}

\item{HVF_method}{The method to use for finding highly variable features. Options are "vst", "mvp" or "disp" (default is "vst").}

\item{nHVF}{The number of highly variable features to select. If NULL, all highly variable features will be used.}

\item{HVF}{A vector of feature names to use as highly variable features. If NULL, the function will use the highly variable features identified by the HVF method.}

\item{do_scaling}{A logical value indicating whether to perform scaling. If TRUE, the function will force to scale the data using the ScaleData function.}

\item{vars_to_regress}{A vector of feature names to use as regressors in the scaling step. If NULL, no regressors will be used.}

\item{regression_model}{The regression model to use for scaling. Options are "linear", "poisson", or "negativebinomial" (default is "linear").}

\item{linear_reduction}{The linear dimensionality reduction method to use. Options are "pca", "svd", "ica", "nmf", "mds", or "glmpca" (default is "pca").}

\item{linear_reduction_dims}{The number of dimensions to keep after linear dimensionality reduction (default is 50).}

\item{linear_reduction_dims_use}{The dimensions to use for downstream analysis. If NULL, all dimensions will be used.}

\item{linear_reduction_params}{A list of parameters to pass to the linear dimensionality reduction method.}

\item{force_linear_reduction}{A logical value indicating whether to force linear dimensionality reduction even if the specified reduction is already present in the Seurat object.}

\item{nonlinear_reduction}{The nonlinear dimensionality reduction method to use. Options are "umap","umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", or "fr" (default is "umap").}

\item{nonlinear_reduction_dims}{The number of dimensions to keep after nonlinear dimensionality reduction. If a vector is provided, different numbers of dimensions can be specified for each method (default is c(2, 3)).}

\item{nonlinear_reduction_params}{A list of parameters to pass to the nonlinear dimensionality reduction method.}

\item{force_nonlinear_reduction}{A logical value indicating whether to force nonlinear dimensionality reduction even if the specified reduction is already present in the Seurat object.}

\item{neighbor_metric}{The distance metric to use for finding neighbors. Options are "euclidean", "cosine", "manhattan", or "hamming" (default is "euclidean").}

\item{neighbor_k}{The number of nearest neighbors to use for finding neighbors (default is 20).}

\item{cluster_algorithm}{The clustering algorithm to use. Options are "louvain", "slm", or "leiden" (default is "louvain").}

\item{cluster_resolution}{The resolution parameter to use for clustering. Larger values result in fewer clusters (default is 0.6).}

\item{seed}{The random seed to use for reproducibility (default is 11).}
}
\value{
A \code{Seurat} object.
}
\description{
This function performs a standard single-cell analysis workflow.
}
\examples{
library(Matrix)
PrepareEnv()
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType"
)

# Use a combination of different linear
# or non-linear dimension reduction methods
linear_reductions <- c(
  "pca", "nmf", "mds", "glmpca"
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  linear_reduction = linear_reductions,
  nonlinear_reduction = "umap"
)
plist1 <- lapply(
  linear_reductions, function(lr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standard", lr, "UMAP2D"
      ),
      xlab = "", ylab = "", title = lr,
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plotlist = plist1)

nonlinear_reductions <- c(
  "umap", "tsne", "dm", "phate",
  "pacmap", "trimap", "largevis", "fr"
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  linear_reduction = "pca",
  nonlinear_reduction = nonlinear_reductions
)
plist2 <- lapply(
  nonlinear_reductions, function(nr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standardpca", toupper(nr), "2D"
      ),
      xlab = "", ylab = "", title = nr,
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plotlist = plist2)
}
\seealso{
\code{\link{integration_scop}}
}
