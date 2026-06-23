# Run ESTIMATE tumor microenvironment scoring

Compute ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity
scores from bulk or pseudo-bulk expression data. The implementation uses
the original stromal and immune ESTIMATE signatures with an in-package
scoring routine, so it does not require the external `estimate` package.

## Usage

``` r
RunESTIMATE(
  object = NULL,
  count_matrix = NULL,
  assay = NULL,
  layer = "data",
  bulk_assay = "counts",
  sample.by = NULL,
  group.by = NULL,
  aggregate_fun = c("mean", "sum"),
  platform = c("rnaseq", "affymetrix", "agilent", "illumina"),
  filter_common_genes = TRUE,
  min_sig_genes = 10,
  purity = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  Optional expression matrix, `SummarizedExperiment`, or `Seurat`
  object.

- count_matrix:

  Optional expression matrix with genes in rows and samples in columns.
  Used when `object` is not provided as a matrix.

- assay:

  Assay used for `Seurat` input.

- layer:

  Assay layer used for `Seurat` input.

- bulk_assay:

  Assay name used for `SummarizedExperiment` input.

- sample.by:

  Metadata column used to aggregate `Seurat` cells into pseudo-bulk
  samples.

- group.by:

  Optional metadata column used for grouping pseudo-bulk samples. If
  `sample.by` is not supplied, `group.by` is used as the aggregation
  variable.

- aggregate_fun:

  Pseudo-bulk aggregation function for `Seurat` input.

- platform:

  Platform label stored in the result metadata.

- filter_common_genes:

  Whether to restrict the expression matrix to the ESTIMATE common-gene
  universe before scoring.

- min_sig_genes:

  Minimum number of observed stromal and immune signature genes required
  after filtering.

- purity:

  Whether to calculate ESTIMATE tumor purity. The purity formula was
  calibrated in the original ESTIMATE work primarily for Affymetrix
  data; for RNA-seq it is best interpreted cautiously.

- verbose:

  Whether to print progress messages.

## Value

A result bundle for matrix input, the modified `SummarizedExperiment`
for `SummarizedExperiment` input, or the modified `Seurat` object for
`Seurat` input.

## References

Yoshihara et al. (2013) <doi:10.1038/ncomms3612>. Barbie et al. (2009)
<doi:10.1038/nature08460>.
