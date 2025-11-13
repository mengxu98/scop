# Check and preprocess a merged seurat object

This function checks and preprocesses a merged seurat object.

## Usage

``` r
CheckDataMerge(
  srt_merge,
  batch = NULL,
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

- srt_merge:

  A merged \`Seurat\` object that includes the batch information.

- batch:

  A character string specifying the batch variable name.

- assay:

  The name of the assay to be used for downstream analysis.

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
  times in order to be considered as highly variable. The default value
  is `1`.

- HVF:

  A vector of highly variable features. Default is `NULL`.

- vars_to_regress:

  A vector of variable names to include as additional regression
  variables. Default is `NULL`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  An integer specifying the random seed for reproducibility. Default is
  `11`.

## See also

\[CheckDataList\]
