# Check and preprocess a list of `Seurat` objects

Checks and preprocesses a list of `Seurat` objects. It performs various
checks on the input, including verification of input types, assay type
consistency, feature name consistency, and batch column consistency. It
also performs data normalization and variable feature finding based on
the specified parameters. Finally, it prepares the data for integration
analysis based on the highly variable features.

## Usage

``` r
CheckDataList(
  srt_list,
  batch,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  vars_to_regress = NULL,
  cores = 1L,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- batch:

  Batch variable name.

- assay:

  Assay to use. `NULL` uses the default assay.

- do_normalization:

  Whether data normalization should be performed.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, and `"scran"`.

- do_HVF_finding, HVF_method, nHVF, HVF:

  Highly variable features. `HVF_method` is `"vst"`, `"mvp"`, `"disp"`,
  or `"scran"`.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables.

- cores:

  Number of CPU cores.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

## Value

A list containing the preprocessed `Seurat` objects, the highly variable
features, the assay name, and the type of assay.
