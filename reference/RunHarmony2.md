# Run Harmony algorithm

This is a modified version of harmony::RunHarmony specifically designed
for compatibility with
[RunSymphonyMap](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md).

## Usage

``` r
RunHarmony2(object, ...)

# S3 method for class 'Seurat'
RunHarmony2(
  object,
  group.by.vars,
  assay = NULL,
  reduction = "pca",
  dims.use = 1:30,
  project.dim = TRUE,
  reduction.name = "Harmony",
  reduction.key = "Harmony_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- ...:

  Additional arguments to be passed to harmony::RunHarmony.

- group.by.vars:

  The batch variable name.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- reduction:

  Which dimensionality reduction to use. Default is `"pca"`.

- dims.use:

  The dimensions to be used. Default is `1:30`.

- project.dim:

  Whether to project dimension reduction loadings. Default is `TRUE`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"Harmony"`.

- reduction.key:

  The prefix for the column names of the Harmony embeddings. Default is
  `"Harmony_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-04-02 16:42:08] Start standard processing workflow...
#> ℹ [2026-04-02 16:42:08] Checking a list of <Seurat>...
#> ! [2026-04-02 16:42:08] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:42:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:42:10] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:42:11] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:42:11] Number of available HVF: 2000
#> ℹ [2026-04-02 16:42:11] Finished check
#> ℹ [2026-04-02 16:42:12] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:42:12] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:42:16] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:42:16] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:42:16] Reorder clusters...
#> ℹ [2026-04-02 16:42:16] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:42:16] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:42:16] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:42:16] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:42:16] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:42:20] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:42:23] Standard processing workflow completed
panc8_sub <- RunHarmony2(
  panc8_sub,
  group.by.vars = "tech",
  reduction = "pca"
)
#> Error in loadNamespace(x): there is no package called ‘harmony’

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


panc8_sub <- standard_scop(
  panc8_sub,
  prefix = "Harmony",
  linear_reduction = "Harmony"
)
#> ℹ [2026-04-02 16:42:34] Start standard processing workflow...
#> Error in standard_scop(panc8_sub, prefix = "Harmony", linear_reduction = "Harmony"): `linear_reduction` must be one of: "pca", "svd", "ica", "nmf", "mds",
#> and "glmpca"

CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "StandardpcaUMAP2D"
)


CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "HarmonyUMAP2D"
)
```
