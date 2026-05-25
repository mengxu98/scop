# The WNN integration function

The WNN integration function

## Usage

``` r
WNN_integrate(
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
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt_merge:

  A merged \`Seurat\` object that includes the batch information.

- batch:

  A character string specifying the batch variable name.

- append:

  Whether the integrated data will be appended to the original Seurat
  object (`srt_merge`). Default is `TRUE`.

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- do_normalization:

  Whether data normalization should be performed. Default is `TRUE`.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, and `"scran"`. Default is
  `"LogNormalize"`.

- do_HVF_finding:

  Whether to perform high variable feature finding. If `TRUE`, the
  function will force to find the highly variable features (HVF) using
  the specified HVF method.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`. Default is `"separate"`.

- HVF_method:

  The method to use for finding highly variable features. Options are
  `"vst"`, `"mvp"`, `"disp"`, or `"scran"`. Default is `"vst"`.

- nHVF:

  The number of highly variable features to select. If NULL, all highly
  variable features will be used. Default is `2000`.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable. Default is `1`.

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

  A vector of variable names to include as additional regression
  variables. Default is `NULL`.

- regression_model:

  The regression model to use for scaling. Options are `"linear"`,
  `"poisson"`, or `"negativebinomial"`. Default is `"linear"`.

- scale_within_batch:

  Whether to scale data within each batch. Only valid when the
  `integration_method` is one of `"Uncorrected"`, `"Seurat"`, `"MNN"`,
  `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.

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

## Examples

``` r
data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
pbmcmultiome_sub <- WNN_integrate(
  srt_merge = pbmcmultiome_sub,
  batch = "batch",
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10
)
#> ℹ [2026-05-25 08:46:49] Start standard processing workflow...
#> ℹ [2026-05-25 08:46:49] Auto preprocess assays: "RNA" and "peaks"
#> ℹ [2026-05-25 08:46:49] Start standard processing workflow...
#> ℹ [2026-05-25 08:46:49] Checking a list of <Seurat>...
#> ! [2026-05-25 08:46:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 08:46:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 08:46:51] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> Warning: pseudoinverse used at -2.3979
#> Warning: neighborhood radius 0.30103
#> Warning: reciprocal condition number  9.9917e-16
#> ℹ [2026-05-25 08:46:52] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 08:46:52] Number of available HVF: 2000
#> ℹ [2026-05-25 08:46:52] Finished check
#> ℹ [2026-05-25 08:46:52] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 08:46:52] Perform pca linear dimension reduction
#> ℹ [2026-05-25 08:46:53] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 08:46:53] Reorder clusters...
#> ℹ [2026-05-25 08:46:53] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 08:46:53] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 08:46:53] Perform umap nonlinear dimension reduction using RNApca (1:10)
#> ℹ [2026-05-25 08:46:57] Perform umap nonlinear dimension reduction using RNApca (1:10)
#> ✔ [2026-05-25 08:47:02] Standard processing workflow completed
#> ℹ [2026-05-25 08:47:02] Start standard processing workflow...
#> ℹ [2026-05-25 08:47:02] Checking a list of <Seurat>...
#> ! [2026-05-25 08:47:02] Data 1/1 of the `srt_list` is "raw_counts"
#> ℹ [2026-05-25 08:47:02] Perform `RunTFIDF()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 08:47:02] Perform `FindTopFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 08:47:02] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 08:47:02] Number of available HVF: 11413
#> ℹ [2026-05-25 08:47:03] Finished check
#> ℹ [2026-05-25 08:47:03] `normalization_method` is TFIDF. Use lsi workflow
#> ℹ [2026-05-25 08:47:03] Perform svd linear dimension reduction
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-05-25 08:47:04] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 08:47:04] Reorder clusters...
#> ℹ [2026-05-25 08:47:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 08:47:04] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 08:47:04] Perform umap nonlinear dimension reduction using ATACsvd (1:10)
#> ℹ [2026-05-25 08:47:09] Perform umap nonlinear dimension reduction using ATACsvd (1:10)
#> ✔ [2026-05-25 08:47:13] Standard processing workflow completed
#> ℹ [2026-05-25 08:47:13] Adjust neighbor k from 20 to 20 for small-sample WNN graph construction
#> ℹ [2026-05-25 08:47:13] Adjust WNN knn.range to 80 for small-sample graph construction
#> ℹ [2026-05-25 08:47:13] Perform WNN integration using RNApca and ATAClsi
#> Calculating cell-specific modality weights
#> Finding 20 nearest neighbors for each modality.
#> Calculating kernel bandwidths
#> Finding multimodal neighbors
#> Constructing multimodal KNN graph
#> Constructing multimodal SNN graph
#> ℹ [2026-05-25 08:47:15] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-25 08:47:15] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-25 08:47:15] Reorder clusters...
#> ℹ [2026-05-25 08:47:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 08:47:15] Perform umap nonlinear dimension reduction using WNN
#> ℹ [2026-05-25 08:47:20] Perform umap nonlinear dimension reduction using WNN
#> Warning: Key ‘RNApcaUMAP2D_’ taken, using ‘rnaumap2d_’ instead
#> Warning: Key ‘RNApcaUMAP3D_’ taken, using ‘rnaumap3d_’ instead
```
