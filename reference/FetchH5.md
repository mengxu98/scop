# Fetch data from the hdf5 file and returns a Seurat object

Fetch data from the hdf5 file and returns a Seurat object

## Usage

``` r
FetchH5(
  data_file,
  meta_file,
  name = NULL,
  features = NULL,
  layer = NULL,
  assay = NULL,
  metanames = NULL,
  reduction = NULL
)
```

## Arguments

- data_file:

  The path to the hdf5 file containing the data.

- meta_file:

  The path to the hdf5 file containing the metadata.

- name:

  The name of the dataset in the hdf5 file. If not specified, the
  function will attempt to find the shared group name in both files.

- features:

  The names of the genes or features to fetch. If specified, only these
  features will be fetched.

- layer:

  The layer for the counts in the hdf5 file. If not specified, the first
  layer will be used.

- assay:

  The name of the assay to use. If not specified, the default assay in
  the hdf5 file will be used.

- metanames:

  The names of the metadata columns to fetch.

- reduction:

  The name of the reduction to fetch.

## Value

A Seurat object with the fetched data.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[RunSCExplorer](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-12 03:56:44] Start standard processing workflow...
#> ℹ [2026-05-12 03:56:45] Checking a list of <Seurat>...
#> ! [2026-05-12 03:56:45] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 03:56:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 03:56:46] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 03:56:47] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 03:56:47] Number of available HVF: 2000
#> ℹ [2026-05-12 03:56:47] Finished check
#> ℹ [2026-05-12 03:56:47] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 03:56:48] Perform pca linear dimension reduction
#> ℹ [2026-05-12 03:56:48] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-12 03:56:48] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-12 03:56:49] Reorder clusters...
#> ℹ [2026-05-12 03:56:49] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-12 03:56:49] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-12 03:56:49] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-12 03:56:53] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-12 03:56:56] Standard processing workflow completed
PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
#> ℹ [2026-05-12 03:56:56] Create SCExplorer base directory: ./SCExplorer
#> ℹ [2026-05-12 03:56:56] Set the project name of each <Seurat> to their dataset name
#> ℹ [2026-05-12 03:56:56] Prepare data for object: "SeuratProject"
#> ℹ [2026-05-12 03:56:56] Write the expression matrix to: ./SCExplorer/data.hdf5
#> ℹ [2026-05-12 03:56:59] Write the meta information to: ./SCExplorer/meta.hdf5
srt <- FetchH5(
  data_file = "./SCExplorer/data.hdf5",
  meta_file = "./SCExplorer/meta.hdf5",
  features = c("Ins1", "Ghrl"),
  metanames = c("SubCellType", "Phase"),
  reduction = "UMAP"
)
CellDimPlot(
  srt,
  group.by = c("SubCellType", "Phase"),
  reduction = "UMAP"
)

FeatureDimPlot(
  srt,
  features = c("Ins1", "Ghrl"),
  reduction = "UMAP"
)
```
