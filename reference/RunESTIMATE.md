# Run ESTIMATE tumor microenvironment scoring

Compute ESTIMATE stromal, immune, combined ESTIMATE, and tumor-purity
scores from bulk or pseudo-bulk expression data without requiring the
external `estimate` package.

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

  Optional genes-by-samples expression matrix.

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

  Optional grouping metadata column. If `sample.by` is absent, this is
  also used for pseudo-bulk aggregation.

- aggregate_fun:

  Pseudo-bulk aggregation function for `Seurat` input.

- platform:

  Platform label stored in result metadata.

- filter_common_genes:

  Whether to restrict expression to the ESTIMATE common-gene universe.

- min_sig_genes:

  Minimum observed stromal and immune signature genes required after
  filtering.

- purity:

  Whether to calculate ESTIMATE tumor purity.

- verbose:

  Whether to print progress messages.

## Value

A result bundle for matrix input, or the modified `SummarizedExperiment`
or `Seurat` object.

## Details

`ESTIMATEScore` is `StromalScore + ImmuneScore`. `TumorPurity` is
calculated as `cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)`. The
purity formula was calibrated in the original ESTIMATE work primarily
for Affymetrix data.

## References

Yoshihara et al. (2013)
[doi:10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612) . Barbie et
al. (2009)
[doi:10.1038/nature08460](https://doi.org/10.1038/nature08460) .
