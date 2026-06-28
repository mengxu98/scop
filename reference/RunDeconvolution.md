# Run bulk or pseudobulk deconvolution

Estimate cell-type proportions from a bulk-like expression matrix stored
in a `SummarizedExperiment` object, using a `Seurat` reference.

## Usage

``` r
RunDeconvolution(object, ...)

# S3 method for class 'SummarizedExperiment'
RunDeconvolution(
  object,
  reference = NULL,
  method = c("MuSiC", "BisqueRNA", "BayesPrism", "CIBERSORT"),
  group.by = NULL,
  sample.by = NULL,
  cellstate.by = NULL,
  bulk_assay = "counts",
  ref_assay = NULL,
  ref_layer = "counts",
  backend = c("cpp", "r"),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing bulk-like counts.

- ...:

  Additional parameters forwarded to the internal deconvolution backend.

- reference:

  A `Seurat` reference object used to build cell-type profiles. Not
  required for `"CIBERSORT"`.

- method:

  Deconvolution method. One of `"MuSiC"`, `"BisqueRNA"`, `"BayesPrism"`,
  or `"CIBERSORT"`.

- group.by:

  Metadata column in `reference` defining reference cell types.

- sample.by:

  Metadata column in `reference` defining biological sample / donor IDs.
  Used by the `r` backends of `MuSiC` and `BisqueRNA`. If `NULL`, SCOP
  will try to infer a suitable column automatically.

- cellstate.by:

  Metadata column in `reference` defining cell states for the `r`
  backend of `BayesPrism`. If `NULL`, `group.by` is reused.

- bulk_assay:

  Assay name in `object` used as the bulk counts matrix.

- ref_assay:

  Assay name in `reference` used for the reference profiles.

- ref_layer:

  Layer name in `reference` used for reference counts.

- backend:

  Deconvolution engine backend. `"r"` uses the original method package
  implementation. `"cpp"` is reserved for native SCOP implementations
  when available.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `SummarizedExperiment` object with results stored in
`S4Vectors::metadata(object)[["Deconvolution"]]`.

## See also

[DeconvolutionPlot](https://mengxu98.github.io/scop/reference/DeconvolutionPlot.md)

## Examples

``` r
data(islet_bulk)
islet_bulk <- RunDeconvolution(
  islet_bulk,
  method = "CIBERSORT",
  backend = "cpp",
  perm = 0
)
#> ℹ [2026-06-28 21:03:33] Use 400 shared genes for CIBERSORT
DeconvolutionPlot(islet_bulk, plot_type = "bar")


DeconvolutionPlot(
  islet_bulk,
  plot_type = "heatmap",
  sample_annotation = "condition",
  sample_split = "condition"
)


DeconvolutionPlot(islet_bulk, plot_type = "box")
```
