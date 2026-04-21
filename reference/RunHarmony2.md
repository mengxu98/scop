# Run Harmony algorithm

This is a modified version of
[harmony::RunHarmony](https://rdrr.io/pkg/harmony/man/RunHarmony.html)
specifically designed for compatibility with
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

  Additional arguments to be passed to
  [harmony::RunHarmony](https://rdrr.io/pkg/harmony/man/RunHarmony.html).

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
#> ℹ [2026-04-21 07:40:16] Start standard processing workflow...
#> ℹ [2026-04-21 07:40:17] Checking a list of <Seurat>...
#> ! [2026-04-21 07:40:17] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-21 07:40:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:40:19] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:40:20] Use the separate HVF from `srt_list`
#> ℹ [2026-04-21 07:40:20] Number of available HVF: 2000
#> ℹ [2026-04-21 07:40:20] Finished check
#> ℹ [2026-04-21 07:40:20] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-21 07:40:21] Perform pca linear dimension reduction
#> ℹ [2026-04-21 07:40:22] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-21 07:40:22] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-21 07:40:22] Reorder clusters...
#> ℹ [2026-04-21 07:40:23] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-21 07:40:23] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-21 07:40:23] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-21 07:40:28] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-21 07:40:32] Standard processing workflow completed
panc8_sub <- RunHarmony2(
  panc8_sub,
  group.by.vars = "tech",
  reduction = "pca"
)
#> Transposing data matrix
#> Initializing state using k-means centroids initialization
#> Harmony 1/10
#> Harmony 2/10
#> Harmony 3/10
#> Harmony 4/10
#> Harmony 5/10
#> Harmony 6/10
#> Harmony 7/10
#> Harmony converged after 7 iterations

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
#> ℹ [2026-04-21 07:40:43] Start standard processing workflow...
#> ℹ [2026-04-21 07:40:43] Checking a list of <Seurat>...
#> ℹ [2026-04-21 07:40:44] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-04-21 07:40:44] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:40:44] Use the separate HVF from `srt_list`
#> ℹ [2026-04-21 07:40:44] Number of available HVF: 2000
#> ℹ [2026-04-21 07:40:45] Finished check
#> ℹ [2026-04-21 07:40:45] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-21 07:40:45] Perform Harmony linear dimension reduction
#> ℹ [2026-04-21 07:40:45] `linear_reduction` Harmony is already existed. Skip calculation
#> ℹ [2026-04-21 07:40:45] Use stored estimated dimensions 1:20 for HarmonyHarmony
#> ℹ [2026-04-21 07:40:46] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-21 07:40:46] Reorder clusters...
#> ℹ [2026-04-21 07:40:46] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-21 07:40:46] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-21 07:40:46] Perform umap nonlinear dimension reduction using HarmonyHarmony (1:20)
#> ℹ [2026-04-21 07:40:51] Perform umap nonlinear dimension reduction using HarmonyHarmony (1:20)
#> ✔ [2026-04-21 07:40:56] Standard processing workflow completed

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
