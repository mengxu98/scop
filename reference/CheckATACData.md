# Check an ATAC object without changing it

Report raw-count, peak-coordinate, cell/metadata and reference alignment
for an existing Seurat `ChromatinAssay`. Fragment path/index
availability is reported separately, so an otherwise valid matrix can be
analysed without fragments. This check does not install packages,
retrieve resources, modify the object, or inspect dependencies of
downstream optional analyses.

## Usage

``` r
CheckATACData(object, assay = NULL, genome = NULL, sample.by = NULL)
```

## Arguments

- object:

  A Seurat object with a Signac `ChromatinAssay`.

- assay:

  Assay to check; defaults to the object's default assay.

- genome:

  Optional expected genome label. Unknown references are reported, never
  inferred from chromosome spelling.

- sample.by:

  Optional metadata column containing sample IDs. This checks complete
  labels, not whether samples are independent biological replicates.

## Value

A list with `valid`, `errors`, `warnings`, `assay`, `n_cells`,
`n_peaks`, `genome`, `samples`, a `fragments` resource table, and
`fragments_ready`. `valid` describes the matrix contract. Fragment
readiness checks local file/index existence and cell mapping for
resources with cells in the selected assay. The resource table retains
all attached files; its `active` column identifies those required for
this assay. Empty mappings left by native subsetting do not block
readiness. It does not hash files, scan their contents, or guarantee
fragment-dependent computation will pass. Use Signac's
`ValidateFragments()` for content validation when required.

## See also

[`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md),
[`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md)

## Examples

``` r
counts <- Matrix::sparseMatrix(
  i = c(1, 2), j = c(1, 2), x = c(2, 3), dims = c(2, 2),
  dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
)
ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
report <- CheckATACData(object)
report$valid
#> [1] TRUE
report$fragments_ready
#> [1] FALSE
```
