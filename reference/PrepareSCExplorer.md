# Prepare Seurat objects for the SCExplorer

Prepares one or multiple `Seurat` objects for the SCExplorer app. It
takes a `Seurat` object or a list of `Seurat` objects as input and
outputs two hdf5 files: one for the data and one for the metadata.

## Usage

``` r
PrepareSCExplorer(
  object,
  base_dir = "SCExplorer",
  data_file = "data.hdf5",
  meta_file = "meta.hdf5",
  assays = "RNA",
  layers = c("counts", "data"),
  ignore_nlevel = 100,
  write_tools = FALSE,
  write_misc = FALSE,
  compression_level = 6,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object or a list of `Seurat` objects.

- base_dir:

  The base directory where the SCExplorer hdf5 files will be written.

- data_file:

  Path to the output data file. If not provided, the file will be named
  `"data.hdf5"` in the current directory.

- meta_file:

  Path to the output meta file. If not provided, the file will be named
  `"meta.hdf5"` in the current directory.

- assays:

  The assays to include in the data file.

- layers:

  The layers to include in the data file.

- ignore_nlevel:

  The number of levels above which a metadata field will be ignored.

- write_tools:

  Whether to write the tools information to the meta file.

- write_misc:

  Whether to write the miscellaneous information to the meta file.

- compression_level:

  Compression level for the HDF5 dataset.

- overwrite:

  Whether to overwrite existing data in the data file.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md),
[RunSCExplorer](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 21:59:06] Start standard processing workflow...
#> ℹ [2026-09-20 21:59:06] Checking a list of <Seurat>...
#> ! [2026-09-20 21:59:06] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 21:59:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:59:06] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:59:06] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 21:59:06] Number of available HVF: 2000
#> ℹ [2026-09-20 21:59:06] Finished check
#> ℹ [2026-09-20 21:59:06] Perform `ScaleData()`
#> ℹ [2026-09-20 21:59:06] Perform pca linear dimension reduction
#> ℹ [2026-09-20 21:59:07] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 21:59:07] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 21:59:07] Reorder clusters...
#> ℹ [2026-09-20 21:59:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 21:59:07] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 21:59:14] Standard processing workflow completed
PrepareSCExplorer(pancreas_sub, base_dir = tempdir())
#> ℹ [2026-09-20 21:59:14] Set the project name of each <Seurat> to their dataset name
#> ℹ [2026-09-20 21:59:14] Prepare data for object: "SeuratProject"
#> ℹ [2026-09-20 21:59:14] Write the expression matrix to: /tmp/Rtmpi4FlDH/data.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/RNA/counts" already exists in /tmp/Rtmpi4FlDH/data.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/RNA/data" already exists in /tmp/Rtmpi4FlDH/data.hdf5
#> ℹ [2026-09-20 21:59:14] Write the meta information to: /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/orig.ident" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nCount_RNA" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nFeature_RNA" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/S_score" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/G2M_score" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nCount_spliced" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nFeature_spliced" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nCount_unspliced" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/nFeature_unspliced" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/CellType" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/SubCellType" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/Phase" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/Standardpca_SNN_res.0.6" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/ident" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/Standardpcaclusters" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata/Standardclusters" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/metadata.stat" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/reductions/Standardpca" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/reductions/StandardpcaUMAP2D" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group "/SeuratProject/reductions/StandardUMAP2D" already exists in /tmp/Rtmpi4FlDH/meta.hdf5
#> ℹ [2026-09-20 21:59:14] Group /SeuratProject/reductions.stat already exists in the /tmp/Rtmpi4FlDH/meta.hdf5
```
