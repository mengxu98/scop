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
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
} # }
```
