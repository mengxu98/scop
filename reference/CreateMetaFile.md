# Create Meta File

Creates a meta file in HDF5 format from a Seurat object.

## Usage

``` r
CreateMetaFile(
  srt,
  meta_file,
  name = NULL,
  write_tools = FALSE,
  write_misc = FALSE,
  ignore_nlevel = 100,
  compression_level = 6,
  overwrite = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- meta_file:

  Path to the output meta file. If not provided, the file will be named
  `"meta.hdf5"` in the current directory.

- name:

  Name of the dataset. If not provided, the name will default to the
  Seurat object's project name.

- write_tools:

  Whether to write the tools information to the meta file. Default is
  `FALSE`.

- write_misc:

  Whether to write the miscellaneous information to the meta file.
  Default is `FALSE`.

- ignore_nlevel:

  The number of levels above which a metadata field will be ignored.
  Default is `100`.

- compression_level:

  The level of compression for the meta file. Default is `6`.

- overwrite:

  Whether to overwrite existing metadata and reductions in the meta
  file. Default is `TRUE`.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md),
[RunSCExplorer](https://mengxu98.github.io/scop/reference/RunSCExplorer.md)
