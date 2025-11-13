# Fetch data from the hdf5 file

This function fetches data from an hdf5 file. It can fetch gene
expression data, metadata, and reduction data from the specified file
and returns a Seurat object.

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
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
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
} # }
```
