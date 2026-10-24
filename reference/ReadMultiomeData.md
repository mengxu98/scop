# Read paired 10x RNA and ATAC counts

Import both measured modalities from one local 10x ARC count matrix
without dropping cells, normalizing counts, or inferring RNA from ATAC.

## Usage

``` r
ReadMultiomeData(
  object,
  metadata = NULL,
  cell.names = NULL,
  sample = NULL,
  genome = NULL,
  annotation = NULL,
  fragments = NULL,
  fragment.cells = NULL,
  assay = "peaks",
  rna.assay = "RNA",
  project = "Multiome",
  verbose = TRUE
)
```

## Arguments

- object:

  Path to one 10x ARC H5 file or one matrix directory containing
  `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`. The input
  must contain exactly the `Gene Expression` and `Peaks` modalities. H5
  input requires at least two features per modality and two cells
  because the native Seurat H5 reader otherwise drops matrix dimensions.

- metadata, cell.names, sample, genome, annotation, fragments,
  fragment.cells, project, verbose:

  Passed to
  [`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md).
  Metadata and sample vectors use original barcodes; fragment mappings
  use final cell names.

- assay:

  Name of the chromatin assay, default `"peaks"`.

- rna.assay:

  Name of the measured RNA assay, default `"RNA"`.

## Value

A Seurat object with a `ChromatinAssay` and a measured RNA `Assay`, both
containing raw counts for the same ordered cells. The active assay is
`assay`. Provenance is stored in `Misc(object, "scop_multiome_input")`
and in the native ATAC import receipt. Native assay initialization also
copies raw counts to the data layer; this is not normalized expression.

## Details

Gene IDs, rather than potentially duplicated gene symbols, identify RNA
features. Both modalities must contain the same unique barcodes. Their
counts are explicitly aligned to the ATAC barcode order. Peak
coordinates in 10x files are 0-based half-open and are converted to
1-based closed coordinates. No barcode suffixes are removed and no
absent modality is zero-filled. H5 input requires the optional hdf5r
package. Fragment files remain external references, as described in
[`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md).

## See also

[`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md),
[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
multiome <- ReadMultiomeData("filtered_feature_bc_matrix.h5", genome = "hg38")
} # }
```
