# Standard workflow for scop

This function performs a standard single-cell or spot-level spatial
analysis workflow.

## Usage

``` r
standard_scop(
  srt,
  prefix = "Standard",
  workflow = c("single_cell", "spatial"),
  assay = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  do_spot_qc = TRUE,
  spot_qc_params = list(),
  do_spatial_variable_features = TRUE,
  spatial_variable_features_params = list(),
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = NULL,
  bayesspace_params = list(),
  reference = NULL,
  reference_label = NULL,
  reference_assay = NULL,
  do_deconvolution = !is.null(reference),
  deconvolution_method = "RCTD",
  deconvolution_params = list(),
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
  verbose = TRUE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- prefix:

  A prefix to add to the names of intermediate objects created by the
  function. Default is `"Standard"`.

- workflow:

  Workflow to run. `"single_cell"` keeps the original standard workflow.
  `"spatial"` runs a spot-level spatial workflow that wraps the original
  workflow with spot QC, spatial variable features, optional spatial
  clustering, and optional RCTD deconvolution.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- do_spot_qc:

  Whether to run
  [`RunSpotQC()`](https://mengxu98.github.io/scop/reference/RunSpotQC.md)
  in the spatial workflow.

- spot_qc_params:

  Named list of additional arguments passed to
  [`RunSpotQC()`](https://mengxu98.github.io/scop/reference/RunSpotQC.md).

- do_spatial_variable_features:

  Whether to run
  [`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md)
  in the spatial workflow.

- spatial_variable_features_params:

  Named list of additional arguments passed to
  [`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md).

- do_spatial_cluster:

  Whether to run spatial-aware clustering in the spatial workflow.

- spatial_cluster_method:

  Spatial clustering method. Only `"BayesSpace"` is supported in this
  workflow.

- spatial_q:

  Number of spatial clusters for
  [`RunBayesSpace()`](https://mengxu98.github.io/scop/reference/RunBayesSpace.md).
  If `NULL`, the number of ordinary spot clusters is used.

- bayesspace_params:

  Named list of additional arguments passed to
  [`RunBayesSpace()`](https://mengxu98.github.io/scop/reference/RunBayesSpace.md).

- reference:

  Optional single-cell reference used for spatial deconvolution.

- reference_label:

  Metadata column in `reference` containing cell type labels.

- reference_assay:

  Assay used in `reference` for deconvolution.

- do_deconvolution:

  Whether to run deconvolution in the spatial workflow. If `NULL`,
  deconvolution is run only when `reference` is provided.

- deconvolution_method:

  Deconvolution method. Only `"RCTD"` is supported in this workflow.

- deconvolution_params:

  Named list of additional arguments passed to
  [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md).

- do_normalization:

  Whether to perform normalization. If `NULL`, normalization will be
  performed if the specified assay does not have scaled data.

- normalization_method:

  The method to use for normalization. Options are `"LogNormalize"`,
  `"SCT"`, `"TFIDF"`, or `"scran"`. When `"SCT"` is used on an RNA
  assay, downstream reductions and clustering are run on the generated
  `"SCT"` assay. Default is `"LogNormalize"`.

- do_HVF_finding:

  Whether to perform high variable feature finding. If `TRUE`, the
  function will force to find the highly variable features (HVF) using
  the specified HVF method.

- HVF_method:

  The method to use for finding highly variable features. Options are
  `"vst"`, `"mvp"`, `"disp"`, or `"scran"`. Default is `"vst"`.

- nHVF:

  The number of highly variable features to select. If NULL, all highly
  variable features will be used. Default is `2000`.

- HVF:

  A vector of feature names to use as highly variable features. If NULL,
  the function will use the highly variable features identified by the
  HVF method.

- do_scaling:

  Whether to perform scaling. If `TRUE`, the function will force to
  scale the data using the
  [Seurat::ScaleData](https://satijalab.org/seurat/reference/ScaleData.html)
  function.

- vars_to_regress:

  A vector of feature names to use as regressors in the scaling step. If
  NULL, no regressors will be used.

- regression_model:

  The regression model to use for scaling. Options are `"linear"`,
  `"poisson"`, or `"negativebinomial"`. Default is `"linear"`.

- linear_reduction:

  The linear dimensionality reduction method to use. Options are
  `"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, or `"glmpca"`. Default is
  `"pca"`.

- linear_reduction_dims:

  The number of dimensions to keep after linear dimensionality
  reduction. Default is `50`.

- linear_reduction_dims_use:

  The dimensions to use for downstream analysis. If `NULL`, estimated
  dimensions stored in the linear reduction will be used when available;
  otherwise, the first up to `50` dimensions will be used as a fallback.

- linear_reduction_params:

  A list of parameters to pass to the linear dimensionality reduction
  method.

- force_linear_reduction:

  Whether to force linear dimensionality reduction even if the specified
  reduction is already present in the Seurat object.

- nonlinear_reduction:

  The nonlinear dimensionality reduction method to use. Options are
  `"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`,
  `"trimap"`, `"largevis"`, or `"fr"`. Default is `"umap"`.

- nonlinear_reduction_dims:

  The number of dimensions to keep after nonlinear dimensionality
  reduction. If a vector is provided, different numbers of dimensions
  can be specified for each method. Default is `c(2, 3)`.

- nonlinear_reduction_params:

  A list of parameters to pass to the nonlinear dimensionality reduction
  method.

- force_nonlinear_reduction:

  Whether to force nonlinear dimensionality reduction even if the
  specified reduction is already present in the Seurat object. Default
  is `TRUE`.

- neighbor_metric:

  The distance metric to use for finding neighbors. Options are
  `"euclidean"`, `"cosine"`, `"manhattan"`, or `"hamming"`. Default is
  `"euclidean"`.

- neighbor_k:

  The number of nearest neighbors to use for finding neighbors. Default
  is `20`.

- cluster_algorithm:

  The clustering algorithm to use. Options are `"louvain"`, `"slm"`, or
  `"leiden"`. Default is `"louvain"`.

- cluster_resolution:

  The resolution parameter to use for clustering. Larger values result
  in fewer clusters. Default is `0.6`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- ...:

  Additional parameters to pass to the dimensionality reduction methods.

## Value

A `Seurat` object.

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType"
)

# Use a combination of different linear
# or nonlinear dimension reduction methods
linear_reductions <- c(
  "pca", "nmf", "mds"
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
      xlab = "", ylab = "",
      title = paste0(lr, "_umap"),
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist1)

nonlinear_reductions <- c(
  "umap", "tsne", "fr"
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
        "Standardpca", nr, "2D"
      ),
      xlab = "", ylab = "",
      title = paste0("pca_", nr),
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist2)

if (requireNamespace("scran", quietly = TRUE)) {
  data(pancreas_sub)
  pancreas_scran <- pancreas_sub[, 1:80]
  pancreas_scran <- standard_scop(
    pancreas_scran,
    assay = "RNA",
    do_normalization = TRUE,
    normalization_method = "scran",
    do_HVF_finding = TRUE,
    HVF_method = "scran",
    nHVF = 100,
    linear_reduction_dims = 10,
    linear_reduction_dims_use = 1:5,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = 2
  )
  CellDimPlot(
    pancreas_scran,
    reduction = "StandardUMAP2D",
    group.by = "Standardclusters"
  )
}

if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- standard_scop(
  visium_human_pancreas_sub,
  workflow = "spatial",
  assay = "Spatial",
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  do_deconvolution = FALSE,
  deconvolution_method = "RCTD",
  linear_reduction_dims = 10,
  linear_reduction_dims_use = 1:5,
  nonlinear_reduction_dims = 2,
  spatial_variable_features_params = list(nfeatures = 50)
)
SpatialSpotPlot(spatial, group.by = "SpotQC")
SpatialSpotPlot(spatial, group.by = "Standardclusters")
SpatialSpotPlot(
  spatial,
  features = spatial@misc[["SpatialVariableFeatures"]][1:2]
)

spatial_bayes <- standard_scop(
  visium_human_pancreas_sub,
  workflow = "spatial",
  assay = "Spatial",
  do_spatial_cluster = TRUE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = 3,
  do_deconvolution = FALSE,
  deconvolution_method = "RCTD",
  bayesspace_params = list(
    n.PCs = 5,
    n.HVGs = 200,
    store_sce = FALSE,
    spatial_cluster_params = list(
      nrep = 200,
      burn.in = 50,
      thin = 10,
      save.chain = FALSE
    )
  )
)
SpatialSpotPlot(spatial_bayes, group.by = "BayesSpace_cluster")

data(panc8_sub)
spatial_rctd <- standard_scop(
  visium_human_pancreas_sub,
  workflow = "spatial",
  assay = "Spatial",
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  do_deconvolution = TRUE,
  deconvolution_method = "RCTD",
  reference = panc8_sub,
  reference_assay = "RNA",
  reference_label = "celltype",
  deconvolution_params = list(
    max_cores = 1,
    min_cells = 25
  )
)
SpatialSpotPlot(spatial_rctd, group.by = "RCTD_dominant_type")

rctd_cols <- grep(
  "^RCTD_prop_",
  colnames(spatial_rctd@meta.data),
  value = TRUE
)
SpatialSpotPlot(
  spatial_rctd,
  group.by = rctd_cols[1:min(4, length(rctd_cols))]
)
SpatialSpotPlot(
  spatial_rctd,
  group.by = "RCTD_dominant_type",
  plot_type = "pie"
)
} # }
```
