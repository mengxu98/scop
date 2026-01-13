# Check and preprocess a list of `Seurat` objects

This function checks and preprocesses a list of `Seurat` objects. It
performs various checks on the input, including verification of input
types, assay type consistency, feature name consistency, and batch
column consistency. It also performs data normalization and variable
feature finding based on the specified parameters. Finally, it prepares
the data for integration analysis based on the highly variable features.

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
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt_list:

  A list of `Seurat` objects to be checked and preprocessed.

- batch:

  A character string specifying the batch variable name.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- do_normalization:

  Whether data normalization should be performed. Default is `TRUE`.

- normalization_method:

  The normalization method to be used. Possible values are
  `"LogNormalize"`, `"SCT"`, and `"TFIDF"`. Default is `"LogNormalize"`.

- do_HVF_finding:

  Whether highly variable feature (HVF) finding should be performed.
  Default is `TRUE`.

- HVF_source:

  The source of highly variable features. Possible values are `"global"`
  and `"separate"`. Default is `"separate"`.

- HVF_method:

  The method for selecting highly variable features. Default is `"vst"`.

- nHVF:

  The number of highly variable features to select. Default is `2000`.

- HVF_min_intersection:

  The feature needs to be present in batches for a minimum number of
  times in order to be considered as highly variable. Default is `1`.

- HVF:

  A vector of highly variable features. Default is `NULL`.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables. Default is `NULL`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A list containing the preprocessed `Seurat` objects, the highly variable
features, the assay name, and the type of assay.
