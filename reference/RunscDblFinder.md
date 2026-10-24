# Run doublet-calling with scDblFinder

Run doublet-calling with scDblFinder

## Usage

``` r
RunscDblFinder(
  object,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE,
  srt = NULL,
  sample.by = NULL,
  seed = 1L
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay used for doublet calling.

- db_rate:

  Expected doublet rate.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Additional arguments to be passed to
  [`scDblFinder::scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.html).

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- sample.by:

  For ATAC, required metadata column identifying captures or libraries
  in which doublets can form. Biological conditions are not capture IDs.
  Native scDblFinder runs independently per capture
  (`multiSampleMode = "split"`). RNA behavior is unchanged.

- seed:

  Random seed for the ATAC path, restoring the caller's RNG.

## Value

A Seurat object with `db.scDblFinder_score` and `db.scDblFinder_class`
updated only for selected-assay cells. For ATAC,
`Misc(assay, "scop_atac_doublets")` records input cells, capture
assignments, parameters, native output and backend version; RNA-only
cells are preserved.

## Details

For ChromatinAssay input, use all raw peaks with native feature
aggregation (`aggregateFeatures = TRUE`, `nfeatures = 25`,
`processing = "normFeatures"`). `nfeatures` and other native model
parameters can be supplied in `...`. With omitted `db_rate`, the ATAC
path uses the backend's per-capture expected rate; explicitly supply the
appropriate rate for the capture chemistry. Scores/classes are
predictions, not ground-truth doublet labels. No cells are removed.
Existing upstream doublet labels can instead be used directly to select
cells; running a second detector is not required.
