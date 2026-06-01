# Run bulk or pseudobulk deconvolution

Estimate cell-type proportions from a bulk-like expression matrix stored
in a `SummarizedExperiment` object, using a `Seurat` reference.

## Usage

``` r
RunDeconvolution(object, ...)

# S3 method for class 'SummarizedExperiment'
RunDeconvolution(
  object,
  reference,
  method = c("MuSiC", "BisqueRNA", "BayesPrism"),
  group.by,
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

  A `Seurat` reference object used to build cell-type profiles.

- method:

  Deconvolution method. One of `"MuSiC"`, `"BisqueRNA"`, or
  `"BayesPrism"`.

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
data(panc8_sub)
islet_bulk <- RunDeconvolution(
  islet_bulk,
  reference = panc8_sub,
  method = "MuSiC",
  group.by = "celltype"
)
#> Error in build_context(mode = "pure_bulk", bulk_se = object, bulk_assay = bulk_assay): could not find function "build_context"
DeconvolutionPlot(islet_bulk, plot_type = "bar")
#> Error in resolve_deconvolution_result(object = object, res = res): Cannot find deconvolution results in
#> `S4Vectors::metadata(object)[['Deconvolution']]()`
ht <- DeconvolutionPlot(
  islet_bulk,
  plot_type = "heatmap",
  sample_annotation = "condition",
  sample_split = "condition"
)
#> Error in resolve_deconvolution_result(object = object, res = res): Cannot find deconvolution results in
#> `S4Vectors::metadata(object)[['Deconvolution']]()`
ComplexHeatmap::draw(ht)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'draw': object 'ht' not found
DeconvolutionPlot(islet_bulk, plot_type = "box")
#> Error in resolve_deconvolution_result(object = object, res = res): Cannot find deconvolution results in
#> `S4Vectors::metadata(object)[['Deconvolution']]()`
```
