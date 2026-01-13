# Create HDF5 data file from Seurat object

Create HDF5 data file from Seurat object

## Usage

``` r
CreateDataFile(
  srt,
  data_file,
  name = NULL,
  assays = "RNA",
  layers = "data",
  compression_level = 6,
  overwrite = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- data_file:

  Path to the output data file. If not provided, the file will be named
  `"data.hdf5"` in the current directory.

- name:

  Name of the dataset. If not provided, the name will default to the
  Seurat object's project name.

- assays:

  The assays to include in the data file. Default is `"RNA"`.

- layers:

  The layers to include in the data file. Default is `"data"`.

- compression_level:

  Compression level for the HDF5 dataset. Default is `6`.

- overwrite:

  Whether to overwrite existing data in the data file. Default is
  `TRUE`.

## See also

[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md),
[RunSCExplorer](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
