# The integration workflow

Integrate single-cell RNA-seq data using various integration methods.

## Usage

``` r
integration_scop(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_method = c("Uncorrected", "Seurat", "CCA", "RPCA", "scVI", "scVI5", "MNN",
    "fastMNN", "fastMNN5", "Harmony", "Harmony5", "Scanorama", "BBKNN", "CSS",
    "Coralysis", "LIGER", "Conos", "ComBat"),
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
  seed = 11,
  verbose = TRUE,
  ...
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
  will be used.

- integration_method:

  A character vector specifying the integration method to use. Supported
  methods are: `"Uncorrected"`, `"Seurat"`, `"CCA"`, `"RPCA"`, `"scVI"`,
  `"scVI5"`, `"MNN"`, `"fastMNN"`, `"fastMNN5"`, `"Harmony"`,
  `"Harmony5"`, `"Scanorama"`, `"BBKNN"`, `"CSS"`, `"Coralysis"`,
  `"LIGER"`, `"Conos"`, `"ComBat"`. Default is `"Uncorrected"`.

- do_normalization:

  Whether data normalization should be performed. Default is `TRUE`.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, and `"TFIDF"`. Default is `"LogNormalize"`.

- do_HVF_finding:

  Whether to perform high variable feature finding. If `TRUE`, the
  function will force to find the highly variable features (HVF) using
  the specified HVF method.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`. Default is `"separate"`.

- HVF_method:

  The method to use for finding highly variable features. Options are
  `"vst"`, `"mvp"`, or `"disp"`. Default is `"vst"`.

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

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to the integration method functions.

## Value

A `Seurat` object.

## See also

[Seurat_integrate](https://mengxu98.github.io/scop/reference/Seurat_integrate.md),
[scVI_integrate](https://mengxu98.github.io/scop/reference/scVI_integrate.md),
[MNN_integrate](https://mengxu98.github.io/scop/reference/MNN_integrate.md),
[fastMNN_integrate](https://mengxu98.github.io/scop/reference/fastMNN_integrate.md),
[Harmony_integrate](https://mengxu98.github.io/scop/reference/Harmony_integrate.md),
[Scanorama_integrate](https://mengxu98.github.io/scop/reference/Scanorama_integrate.md),
[BBKNN_integrate](https://mengxu98.github.io/scop/reference/BBKNN_integrate.md),
[CSS_integrate](https://mengxu98.github.io/scop/reference/CSS_integrate.md),
[Coralysis_integrate](https://mengxu98.github.io/scop/reference/Coralysis_integrate.md),
[LIGER_integrate](https://mengxu98.github.io/scop/reference/LIGER_integrate.md),
[Conos_integrate](https://mengxu98.github.io/scop/reference/Conos_integrate.md),
[ComBat_integrate](https://mengxu98.github.io/scop/reference/ComBat_integrate.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "LIGER"
)
#> ◌ [2026-04-06 05:20:39] Run integration workflow...
#> ℹ [2026-04-06 05:20:39] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-06 05:20:40] Checking a list of <Seurat>...
#> ! [2026-04-06 05:20:40] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:42] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-04-06 05:20:42] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:42] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:44] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-04-06 05:20:44] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:44] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:46] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-04-06 05:20:46] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:46] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:48] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-04-06 05:20:48] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:50] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-06 05:20:51] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:20:51] Number of available HVF: 2000
#> ℹ [2026-04-06 05:20:51] Finished check
#> Warning: Layer ‘ligerScaleData’ is empty
#> ℹ [2026-04-06 05:20:56] Prepare rliger layer "ligerScaleData" ...
#> ℹ [2026-04-06 05:20:57] Perform LIGER integration
#> ℹ [2026-04-06 05:21:09] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-06 05:21:09] Reorder clusters...
#> ℹ [2026-04-06 05:21:11] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:21:12] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ℹ [2026-04-06 05:21:16] Perform umap nonlinear dimension reduction using LIGER (1:20)
#> ✔ [2026-04-06 05:21:25] LIGER integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)


integration_methods <- c(
  "Uncorrected", "Seurat", "CCA", "RPCA", "scVI", "scVI5",
  "MNN", "fastMNN", "fastMNN5", "Harmony", "Harmony5",
  "Scanorama", "BBKNN", "CSS", "Coralysis", "LIGER", "Conos", "ComBat"
)
p_list <- list()
for (method in integration_methods) {
  panc8_sub <- integration_scop(
    panc8_sub,
    batch = "tech",
    integration_method = method,
    linear_reduction_dims_use = 1:50,
    nonlinear_reduction = "umap"
  )
  p_list[[method]] <- CellDimPlot(
    panc8_sub,
    group.by = c("tech", "celltype"),
    reduction = paste0(method, "UMAP2D"),
    xlab = "", ylab = "",
    title = method,
    legend.position = "none",
    theme_use = "theme_blank"
  )
}
#> ◌ [2026-04-06 05:21:25] Run integration workflow...
#> ℹ [2026-04-06 05:21:25] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-06 05:21:26] Checking a list of <Seurat>...
#> ℹ [2026-04-06 05:21:26] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:26] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:27] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:27] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:27] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:27] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:28] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:28] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:28] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:28] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:29] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:21:29] Number of available HVF: 2000
#> ℹ [2026-04-06 05:21:29] Finished check
#> ℹ [2026-04-06 05:21:35] Perform Uncorrected integration
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-04-06 05:21:35] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 05:21:36] Perform "pca" linear dimension reduction
#> ℹ [2026-04-06 05:21:38] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-06 05:21:38] Reorder clusters...
#> ℹ [2026-04-06 05:21:39] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:21:39] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-04-06 05:21:44] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-04-06 05:21:51] Uncorrected integration completed
#> ◌ [2026-04-06 05:21:51] Run integration workflow...
#> ℹ [2026-04-06 05:21:51] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-06 05:21:52] Checking a list of <Seurat>...
#> ℹ [2026-04-06 05:21:52] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:52] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:53] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:53] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:53] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:53] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:54] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:54] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:55] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:21:55] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-06 05:21:55] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:21:55] Number of available HVF: 2000
#> ℹ [2026-04-06 05:21:56] Finished check
#> ℹ [2026-04-06 05:22:15] Perform FindIntegrationAnchors
#> ℹ [2026-04-06 05:23:01] Perform Seurat integration
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Warning: Different cells in new layer data than already exists for scale.data
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Warning: Different cells in new layer data than already exists for scale.data
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Warning: Different cells in new layer data than already exists for scale.data
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> ℹ [2026-04-06 05:23:43] Perform ScaleData on `srt_integrated`
#> ℹ [2026-04-06 05:23:44] Perform "pca" linear dimension reduction
#> ℹ [2026-04-06 05:23:45] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-06 05:23:45] Reorder clusters...
#> ℹ [2026-04-06 05:23:45] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:23:45] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ℹ [2026-04-06 05:23:50] Perform umap nonlinear dimension reduction using Seuratpca (1:50)
#> ✔ [2026-04-06 05:24:02] Seurat integration completed
#> ◌ [2026-04-06 05:24:02] Run integration workflow...
#> ! [2026-04-06 05:24:03] Data is "unknown". Will perform `Seurat::NormalizeData()`
#> ℹ [2026-04-06 05:24:03] Perform `Seurat::NormalizeData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-06 05:24:05] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-04-06 05:24:06] Number of available HVF: 2000
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-04-06 05:24:07] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-06 05:24:08] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-04-06 05:24:09] Perform Seurat v5 integration with `CCAIntegration()`
#> ℹ [2026-04-06 05:24:30] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-06 05:24:30] Reorder clusters...
#> ℹ [2026-04-06 05:24:30] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:24:30] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ℹ [2026-04-06 05:24:36] Perform umap nonlinear dimension reduction using CCA (1:50)
#> ✔ [2026-04-06 05:24:46] CCA integration completed
#> ◌ [2026-04-06 05:24:46] Run integration workflow...
#> ! [2026-04-06 05:24:48] Data is "unknown". Will perform `Seurat::NormalizeData()`
#> ℹ [2026-04-06 05:24:48] Perform `Seurat::NormalizeData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-06 05:24:49] Perform `Seurat::FindVariableFeatures()` per batch (`HVF_source = 'separate'`)
#> ℹ [2026-04-06 05:24:51] Number of available HVF: 2000
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-04-06 05:24:52] Perform `Seurat::ScaleData()` on split layers for Seurat v5 integration
#> ℹ [2026-04-06 05:24:52] Perform PCA on split layers before `Seurat::IntegrateLayers()`
#> ℹ [2026-04-06 05:24:54] Perform Seurat v5 integration with `RPCAIntegration()`
#> Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 10 globals exported for future expression (‘FUN()’) is 4.41 GiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option The three largest globals are ‘FUN’ (4.37 GiB of class ‘function’), ‘object.list’ (35.49 MiB of class ‘list’) and ‘NNHelper’ (9.67 KiB of class ‘function’)

nonlinear_reductions <- c(
  "umap", "tsne", "dm", "phate",
  "pacmap", "trimap", "largevis", "fr"
)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat",
  linear_reduction_dims_use = 1:50,
  nonlinear_reduction = nonlinear_reductions
)
#> ◌ [2026-04-06 05:24:56] Run integration workflow...
#> ℹ [2026-04-06 05:24:56] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-06 05:24:58] Checking a list of <Seurat>...
#> ℹ [2026-04-06 05:24:59] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:24:59] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-04-06 05:24:59] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:24:59] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-04-06 05:25:00] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:25:00] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-04-06 05:25:00] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:25:00] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-04-06 05:25:01] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 05:25:01] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-06 05:25:01] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:25:02] Number of available HVF: 2000
#> ℹ [2026-04-06 05:25:02] Finished check
#> ℹ [2026-04-06 05:26:50] Perform FindIntegrationAnchors
#> Error in getGlobalsAndPackages(expr, envir = envir, globals = globals): The total size of the 3 globals exported for future expression (‘FUN()’) is 605.91 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (605.88 MiB of class ‘function’), ‘anchor.features’ (25.71 KiB of class ‘character’) and ‘ScaleData’ (406 bytes of class ‘function’)
for (nr in nonlinear_reductions) {
  print(
    CellDimPlot(
      panc8_sub,
      group.by = c("tech", "celltype"),
      reduction = paste0("Seurat", nr, "2D"),
      xlab = "", ylab = "", title = nr,
      legend.position = "none", theme_use = "theme_blank"
    )
  )
}







```
