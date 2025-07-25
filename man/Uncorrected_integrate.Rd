% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{Uncorrected_integrate}
\alias{Uncorrected_integrate}
\title{Uncorrected_integrate}
\usage{
Uncorrected_integrate(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
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
\item{srt_merge}{A merged Seurat object that includes the batch information.}

\item{batch}{A character string specifying the batch variable name.}

\item{append}{Logical, if TRUE, the integrated data will be appended to the original Seurat object (srt_merge).}

\item{srt_list}{A list of Seurat objects to be checked and preprocessed.}

\item{assay}{The name of the assay to be used for downstream analysis.}

\item{do_normalization}{A logical value indicating whether data normalization should be performed.}

\item{normalization_method}{The normalization method to be used. Possible values are "LogNormalize", "SCT", and "TFIDF". Default is "LogNormalize".}

\item{do_HVF_finding}{A logical value indicating whether highly variable feature (HVF) finding should be performed. Default is TRUE.}

\item{HVF_source}{The source of highly variable features. Possible values are "global" and "separate". Default is "separate".}

\item{HVF_method}{The method for selecting highly variable features. Default is "vst".}

\item{nHVF}{The number of highly variable features to select. Default is 2000.}

\item{HVF_min_intersection}{The feature needs to be present in batches for a minimum number of times in order to be considered as highly variable. The default value is 1.}

\item{HVF}{A vector of highly variable features. Default is NULL.}

\item{do_scaling}{A logical value indicating whether to perform scaling. If TRUE, the function will force to scale the data using the ScaleData function.}

\item{vars_to_regress}{A vector of variable names to include as additional regression variables. Default is NULL.}

\item{regression_model}{The regression model to use for scaling. Options are "linear", "poisson", or "negativebinomial" (default is "linear").}

\item{scale_within_batch}{Whether to scale data within each batch.
Only valid when the \code{integration_method} is one of \code{"Uncorrected"}, \code{"Seurat"}, \code{"MNN"}, \code{"Harmony"}, \code{"BBKNN"}, \code{"CSS"}, \code{"ComBat"}.}

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

\item{seed}{An integer specifying the random seed for reproducibility. Default is 11.}
}
\description{
Uncorrected_integrate
}
