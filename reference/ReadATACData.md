# Read raw ATAC counts into a ChromatinAssay

Import one explicitly shared peak universe from a sparse matrix and
genomic ranges, a 10x ATAC/ARC H5 file, or Matrix Market/peaks
BED/barcodes files. SCOP checks and aligns identifiers; Signac creates
the native assay. No cells or peaks are filtered, and no resources are
downloaded or installed.

## Usage

``` r
ReadATACData(
  object,
  ranges = NULL,
  metadata = NULL,
  peaks.file = NULL,
  barcodes.file = NULL,
  cell.names = NULL,
  sample = NULL,
  genome = NULL,
  annotation = NULL,
  fragments = NULL,
  fragment.cells = NULL,
  assay = "peaks",
  project = "ATAC",
  verbose = TRUE
)
```

## Arguments

- object:

  A numeric sparse peak-by-cell matrix, a local 10x ATAC/ARC H5 path, or
  a local Matrix Market path. Separate matrices with different peak
  universes must be recounted before combining; this function never
  fills unmeasured peaks with zeros. The supported Signac constructor
  requires at least two cells and two peaks.

- ranges:

  A `GRanges` in 1-based closed coordinates, required with a sparse
  matrix. Named ranges are aligned to the matrix row names; unnamed
  ranges explicitly use matrix row order. Duplicate genomic intervals
  are rejected.

- metadata:

  Optional data frame with unique row names matching the input matrix
  cell IDs exactly. Rows are reordered by ID, never by position.
  Supplied `nCount_<assay>`/`nFeature_<assay>` columns must agree with
  the raw matrix; conflicting QC metadata is rejected.

- peaks.file, barcodes.file:

  Three-column peaks BED and one-column barcode paths, both required
  when `object` is a Matrix Market path. Gzip compression is supported.
  Without companion files a path is read as H5, importing only the
  explicitly declared `Peaks` modality from its standard feature table;
  H5 requires optional `hdf5r`. BED and 10x H5 peak IDs use 0-based
  half-open coordinates, converted to 1-based closed `GRanges` (start +
  1, end unchanged).

- cell.names:

  Optional unique final cell IDs, in input column order. For example,
  explicitly prefix barcodes from separate captures. Suffixes are never
  removed or guessed.

- sample:

  Optional single sample ID, or a named character vector keyed by input
  cell IDs, stored as `sample` metadata. A sample label does not
  establish biological replication; retain donor/condition/capture
  metadata.

- genome:

  Optional genome label such as `"hg38"`, checked against known labels
  on `ranges` and `annotation`. This sets a label, not a reference
  download. Supply seqlengths in `ranges` for chromosome-bound checks.

- annotation:

  Optional Signac-compatible gene annotation `GRanges`.

- fragments:

  Optional local fragment path(s). Files must be bgzip compressed and
  indexed. Signac validates the paths and barcode membership.

- fragment.cells:

  Explicit fragment barcode mapping: a named character vector for one
  file, or a list of vectors in `fragments` order. Vector names are
  final object cell IDs; values are exact barcodes in that file.
  Multiple files require disjoint mappings covering every imported cell.
  With one file and no mapping, input matrix barcodes are used
  unchanged.

- assay:

  Name of the new chromatin assay.

- project:

  Project name passed to `CreateSeuratObject()`.

- verbose:

  Whether Signac prints progress messages.

## Value

A Seurat object containing a Signac `ChromatinAssay` with raw counts.

## Details

Signac uses genomic coordinates as feature names. Original matrix
feature IDs are retained in assay feature metadata `input_peak_id`;
input barcodes are retained as cell metadata `input_barcode`. Import
provenance is stored next to the data in
`Misc(object[[assay]], "scop_atac_input")`.

Fragment files remain external references.
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) does not copy them.
After [`readRDS()`](https://rdrr.io/r/base/readRDS.html), use
[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md)
to inspect their availability and export its `fragments` table as a
resource manifest. Missing external files do not invalidate otherwise
valid matrix-only analysis. Update moved paths explicitly using Signac's
`UpdatePath()`; paths are never repaired silently.

## See also

[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md),
[`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)

## Examples

``` r
peaks <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
  seqnames = c("chr1", "chr1", "chr2"),
  start = c(1L, 101L, 201L), end = c(50L, 150L, 250L)
))
counts <- Matrix::sparseMatrix(
  i = c(1, 3, 2, 3), j = c(1, 1, 2, 2), x = c(2, 1, 3, 1),
  dims = c(3, 2), dimnames = list(c("p1", "p2", "p3"), c("AA-1", "BB-1"))
)
atac <- ReadATACData(counts, ranges = peaks, genome = "hg38", sample = "S1")
CheckATACData(atac)$valid
#> [1] TRUE
```
