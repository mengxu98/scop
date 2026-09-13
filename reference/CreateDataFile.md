# Create HDF5 data file from Seurat object

Create HDF5 data file from Seurat object

## Usage

``` r
CreateDataFile(
  object,
  data_file,
  name = NULL,
  assays = "RNA",
  layers = "data",
  compression_level = 6,
  overwrite = TRUE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- data_file:

  Path to the output data file. If not provided, the file will be named
  `"data.hdf5"` in the current directory.

- name:

  Name of the dataset. If not provided, the name will default to the
  Seurat object's project name.

- assays:

  The assays to include in the data file.

- layers:

  The layers to include in the data file.

- compression_level:

  Compression level for the HDF5 dataset.

- overwrite:

  Whether to overwrite existing data in the data file.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## See also

[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md),
[RunSCExplorer](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
