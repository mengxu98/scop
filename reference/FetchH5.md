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
#> ℹ [2026-09-27 21:49:08] Start standard processing workflow...
#> ℹ [2026-09-27 21:49:08] Checking a list of <Seurat>...
#> ! [2026-09-27 21:49:08] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 21:49:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 21:49:08] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 21:49:08] Use the separate HVF from `srt_list`
#> ℹ [2026-09-27 21:49:08] Number of available HVF: 2000
#> ℹ [2026-09-27 21:49:08] Finished check
#> ℹ [2026-09-27 21:49:08] Perform `ScaleData()`
#> ℹ [2026-09-27 21:49:08] Perform pca linear dimension reduction
#> ℹ [2026-09-27 21:49:08] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-27 21:49:09] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-27 21:49:09] Reorder clusters...
#> ℹ [2026-09-27 21:49:09] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-27 21:49:09] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-27 21:49:15] Standard processing workflow completed
PrepareSCExplorer(pancreas_sub, base_dir = tempdir())
#> ℹ [2026-09-27 21:49:15] Set the project name of each <Seurat> to their dataset name
#> ℹ [2026-09-27 21:49:15] Prepare data for object: "SeuratProject"
#> ℹ [2026-09-27 21:49:15] Write the expression matrix to: /tmp/RtmpNIRZ7c/data.hdf5
#> ℹ [2026-09-27 21:49:16] Write the meta information to: /tmp/RtmpNIRZ7c/meta.hdf5
srt <- FetchH5(
  data_file = paste0(tempdir(), "/data.hdf5"),
  meta_file = paste0(tempdir(), "/meta.hdf5"),
  features = c("Ins1", "Ghrl"),
  metanames = c("SubCellType", "Phase"),
  reduction = "UMAP"
)
```
