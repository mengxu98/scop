# Run CIBERSORT deconvolution

Estimate immune cell proportions from a bulk expression matrix using the
external `CIBERSORT` package. `sig_matrix = "LM22"` downloads the LM22
signature matrix from `mengxu98/datasets` and caches it locally.

## Usage

``` r
RunCIBERSORT(
  object = NULL,
  count_matrix = NULL,
  sig_matrix = "LM22",
  bulk_assay = "counts",
  perm = 100,
  QN = TRUE,
  absolute = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Optional `SummarizedExperiment` object or expression matrix. When a
  `SummarizedExperiment` is provided, results are stored in
  `metadata(object)[["Deconvolution"]]`.

- count_matrix:

  Optional expression matrix with genes in rows and samples in columns.
  Used when `object` is not provided as a matrix.

- sig_matrix:

  Signature matrix, local file path, or `"LM22"`.

- bulk_assay:

  Assay name in `object` used as the bulk counts matrix.

- perm:

  Number of CIBERSORT permutations.

- QN:

  Whether CIBERSORT should use quantile normalization.

- absolute:

  Passed to CIBERSORT when supported by the installed package.

- verbose:

  Whether to print messages.

- ...:

  Additional parameters forwarded to
  [`CIBERSORT::cibersort()`](https://rdrr.io/pkg/CIBERSORT/man/CIBERSORT.html).

## Value

A deconvolution result bundle for matrix input, or the modified
`SummarizedExperiment` object for `SummarizedExperiment` input.

## Examples

``` r
data(islet_bulk)

if (FALSE) {
# Run CIBERSORT
islet_bulk <- RunCIBERSORT(
  object = islet_bulk,
  sig_matrix = "LM22",
  bulk_assay = "counts",
  perm = 100,
  QN = TRUE
)

# Immune abundance stacked bar plot
p1 <- ImmuneAbundancePlot(
  object = islet_bulk,
  plot_type = "bar",
  group.by = "condition"
)
p1

# Immune cell correlation heatmap
p2 <- ImmuneAbundancePlot(
  object = islet_bulk,
  plot_type = "cor"
)
p2

# Gene-immune correlation butterfly plot
p3 <- GeneImmuneCorPlot(
  object = islet_bulk,
  features = rownames(SummarizedExperiment::assay(islet_bulk, "counts"))[1:3]
)
p3
}
```
