# Plot ATAC fragment lengths

Plot fragment lengths in an explicitly selected genomic region using
[`Signac::FragmentHistogram()`](https://stuartlab.org/signac/reference/FragmentHistogram.html).
This is a regional diagnostic, not a whole-library fragment-length
estimate.

## Usage

``` r
ATACFragmentPlot(
  object,
  region,
  assay = NULL,
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object containing a ChromatinAssay.

- region:

  Genomic region(s), as a GRanges or strings accepted by Signac.
  Required explicitly to avoid silently assuming chromosome naming or
  genome.

- assay:

  Name of the ChromatinAssay; defaults to the active assay.

- group.by:

  One metadata column used for grouping, or NULL for identities.

- cells:

  Cell names in the selected assay; NULL uses all assay cells.

- log.scale:

  Whether to use a logarithmic count axis.

- ...:

  Arguments passed to other functions.

## Value

A ggplot object. The input object is not modified.

## Details

Requires local bgzip fragment files with a `.tbi` or `.csi` index and
explicit cell/barcode mappings. The histogram counts fragment records,
not column-five read multiplicities, and displays lengths between 0 and
800 bp by default. As in Signac, a fragment returned for several query
intervals can contribute more than once; use a single region for an
unambiguous regional histogram.

## See also

[`RunATACQC()`](https://mengxu98.github.io/scop/reference/RunATACQC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ATACFragmentPlot(object, region = "chr1-1-2000000", assay = "peaks")
} # }
```
