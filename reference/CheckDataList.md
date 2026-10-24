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

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

## Value

A list containing the preprocessed `Seurat` objects, the highly variable
features, the assay name, and the type of assay.

## Details

For chromatin assays with `HVF_source = "separate"`, native
[`Seurat::SelectIntegrationFeatures()`](https://satijalab.org/seurat/reference/SelectIntegrationFeatures.html)
ranks peaks by the number of batches selecting them, then by their
median variable-feature rank. The `nHVF` cap is also limited to the
smallest per-batch variable-feature count; `HVF_min_intersection` is
applied to the ranked selection. This preserves native within-batch peak
frequency ranking instead of sorting tied peak IDs.

## Examples

``` r
counts <- matrix(rep(6:1, 4),
  nrow = 6,
  dimnames = list(
    paste0("chr1-", seq(1, 501, 100), "-", seq(50, 550, 100)),
    paste0("cell", 1:4)
  )
)
peaks <- Signac::CreateChromatinAssay(counts,
  min.cells = 0,
  min.features = 0, verbose = FALSE
)
object <- SeuratObject::CreateSeuratObject(peaks, assay = "peaks")
checked <- CheckDataList(list(object),
  batch = "sample", assay = "peaks",
  do_normalization = FALSE, normalization_method = "TFIDF", nHVF = 3,
  verbose = FALSE
)
checked$HVF
#> [1] "chr1-1-50"    "chr1-101-150" "chr1-201-250"
```
