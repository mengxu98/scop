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
  reduction = NULL,
  verbose = TRUE
)
```

## Arguments

- data_file:

  The path to the hdf5 file containing the data.

- meta_file:

  The path to the hdf5 file containing the metadata.

- name:

  Dataset in the hdf5 file. If not specified, the function will attempt
  to find the shared group name in both files.

- features:

  The names of the genes or features to fetch. If specified, only these
  features will be fetched.

- layer:

  The layer for the counts in the hdf5 file. If not specified, the first
  layer will be used.

- assay:

  Assay to use. If not specified, the default assay in the hdf5 file
  will be used.

- metanames:

  The names of the metadata columns to fetch.

- reduction:

  Reduction to fetch.

- verbose:

  Whether to print the message. Default is `TRUE`.

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
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 21:43:23] Start standard processing workflow...
#> ℹ [2026-09-20 21:43:23] Checking a list of <Seurat>...
#> ! [2026-09-20 21:43:23] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 21:43:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:43:24] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:43:24] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 21:43:24] Number of available HVF: 2000
#> ℹ [2026-09-20 21:43:24] Finished check
#> ℹ [2026-09-20 21:43:24] Perform `ScaleData()`
#> ℹ [2026-09-20 21:43:24] Perform pca linear dimension reduction
#> ℹ [2026-09-20 21:43:24] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 21:43:24] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 21:43:24] Reorder clusters...
#> ℹ [2026-09-20 21:43:25] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 21:43:25] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 21:43:30] Standard processing workflow completed
PrepareSCExplorer(pancreas_sub, base_dir = tempdir())
#> ℹ [2026-09-20 21:43:30] Set the project name of each <Seurat> to their dataset name
#> ℹ [2026-09-20 21:43:30] Prepare data for object: "SeuratProject"
#> ℹ [2026-09-20 21:43:30] Write the expression matrix to: /tmp/Rtmpi4FlDH/data.hdf5
#> ℹ [2026-09-20 21:43:33] Write the meta information to: /tmp/Rtmpi4FlDH/meta.hdf5
srt <- FetchH5(
  data_file = paste0(tempdir(), "/data.hdf5"),
  meta_file = paste0(tempdir(), "/meta.hdf5"),
  features = c("Ins1", "Ghrl"),
  metanames = c("SubCellType", "Phase"),
  reduction = "UMAP"
)
```
