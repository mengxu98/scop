# The integration_scop function

Integrate single-cell RNA-seq data using various integration methods.

## Usage

``` r
integration_scop(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_method = c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony",
    "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat"),
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
  methods are: `"Uncorrected"`, `"Seurat"`, `"scVI"`, `"MNN"`,
  `"fastMNN"`, `"Harmony"`, `"Scanorama"`, `"BBKNN"`, `"CSS"`,
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
  [Seurat::ScaleData](https://rdrr.io/pkg/Seurat/man/ScaleData.html)
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

  The dimensions to use for downstream analysis. If `NULL`, all
  dimensions will be used.

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
[LIGER_integrate](https://mengxu98.github.io/scop/reference/LIGER_integrate.md),
[Conos_integrate](https://mengxu98.github.io/scop/reference/Conos_integrate.md),
[ComBat_integrate](https://mengxu98.github.io/scop/reference/ComBat_integrate.md),
[standard_scop](https://mengxu98.github.io/scop/reference/standard_scop.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-01-29 13:45:23] Run Uncorrected integration...
#> ℹ [2026-01-29 13:45:23] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-29 13:45:23] Checking a list of <Seurat>...
#> ! [2026-01-29 13:45:23] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:45:24] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:25] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ! [2026-01-29 13:45:26] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:45:26] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 2/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:27] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ! [2026-01-29 13:45:28] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:45:28] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 3/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:29] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ! [2026-01-29 13:45:30] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:45:30] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 4/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:31] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ! [2026-01-29 13:45:32] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:45:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:33] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:34] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:45:34] Number of available HVF: 2000
#> ℹ [2026-01-29 13:45:35] Finished check
#> ℹ [2026-01-29 13:45:37] Perform Uncorrected integration
#> ℹ [2026-01-29 13:45:41] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:45:41] Perform linear dimension reduction("pca")
#> ℹ [2026-01-29 13:45:43] Perform Seurat::FindClusters ("louvain")
#> ℹ [2026-01-29 13:45:43] Reorder clusters...
#> ℹ [2026-01-29 13:45:43] Perform nonlinear dimension reduction ("umap")
#> ℹ [2026-01-29 13:45:43] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ℹ [2026-01-29 13:45:48] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ✔ [2026-01-29 13:45:54] Run Uncorrected integration done
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)


panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected",
  HVF_min_intersection = 5
)
#> ◌ [2026-01-29 13:45:55] Run Uncorrected integration...
#> ℹ [2026-01-29 13:45:55] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-29 13:45:56] Checking a list of <Seurat>...
#> ℹ [2026-01-29 13:45:56] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:45:56] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:56] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:45:57] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:57] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:45:57] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:58] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:45:58] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:58] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:45:58] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-29 13:45:58] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:45:59] Number of available HVF: 270
#> ℹ [2026-01-29 13:45:59] Finished check
#> ℹ [2026-01-29 13:46:04] Perform Uncorrected integration
#> ℹ [2026-01-29 13:46:04] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:46:04] Perform linear dimension reduction("pca")
#> ℹ [2026-01-29 13:46:06] Perform Seurat::FindClusters ("louvain")
#> ℹ [2026-01-29 13:46:06] Reorder clusters...
#> ℹ [2026-01-29 13:46:06] Perform nonlinear dimension reduction ("umap")
#> ℹ [2026-01-29 13:46:06] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-12) as input
#> ℹ [2026-01-29 13:46:11] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-12) as input
#> ✔ [2026-01-29 13:46:18] Run Uncorrected integration done
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)


panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected",
  HVF_min_intersection = 5,
  scale_within_batch = TRUE
)
#> ◌ [2026-01-29 13:46:19] Run Uncorrected integration...
#> ℹ [2026-01-29 13:46:19] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-29 13:46:20] Checking a list of <Seurat>...
#> ℹ [2026-01-29 13:46:20] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:46:20] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ℹ [2026-01-29 13:46:20] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:46:20] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ℹ [2026-01-29 13:46:21] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:46:21] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ℹ [2026-01-29 13:46:21] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:46:21] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ℹ [2026-01-29 13:46:22] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-29 13:46:22] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-29 13:46:22] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:46:23] Number of available HVF: 270
#> ℹ [2026-01-29 13:46:23] Finished check
#> ℹ [2026-01-29 13:46:42] Perform Uncorrected integration
#> ℹ [2026-01-29 13:46:42] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:46:43] Perform linear dimension reduction("pca")
#> ℹ [2026-01-29 13:46:44] Perform Seurat::FindClusters ("louvain")
#> ℹ [2026-01-29 13:46:44] Reorder clusters...
#> ℹ [2026-01-29 13:46:44] Perform nonlinear dimension reduction ("umap")
#> ℹ [2026-01-29 13:46:44] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-13) as input
#> ℹ [2026-01-29 13:46:49] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-13) as input
#> ✔ [2026-01-29 13:47:00] Run Uncorrected integration done
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)


if (FALSE) { # \dontrun{
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat"
)
CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))

panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat",
  FindIntegrationAnchors_params = list(reduction = "rpca")
)
CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))

integration_methods <- c(
  "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony",
  "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat"
)
for (method in integration_methods) {
  panc8_sub <- integration_scop(
    panc8_sub,
    batch = "tech",
    integration_method = method,
    linear_reduction_dims_use = 1:50,
    nonlinear_reduction = "umap"
  )
  print(
    CellDimPlot(panc8_sub,
      group.by = c("tech", "celltype"),
      reduction = paste0(method, "UMAP2D"),
      xlab = "", ylab = "", title = method,
      legend.position = "none", theme_use = "theme_blank"
    )
  )
}

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
} # }
```
