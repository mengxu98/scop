# Annotate ATAC peaks with nearest genomic features

Attach verified gene annotation and record native nearest-feature
distances, retaining native nearest-feature hits and unmatched peaks.

## Usage

``` r
AnnotateATACPeaks(object, annotation, assay = NULL, verbose = TRUE)
```

## Arguments

- object:

  A Seurat object with a validated `ChromatinAssay`.

- annotation:

  A non-empty Signac-compatible `GRanges` of 1-based closed
  gene/transcript intervals with gene_id, gene_name, gene_biotype, type
  and tx_id (or transcript_id) metadata. It must declare the source
  assay's genome label. Known chromosome lengths must agree and interval
  ends must not exceed them. Extra chromosomes are retained without name
  conversion.

- assay:

  Source chromatin assay; defaults to the active assay.

- verbose:

  Whether to print progress messages.

## Value

The Seurat object with native Signac gene annotation and
`Misc(object[[assay]], "scop_peak_annotation")`. The receipt contains
`nearest`, an
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with peak ID, annotation row index, distance and prefixed annotation
metadata; `unannotated_peaks` lists peaks without a hit. Rows follow
peak order and annotation index within ties. Counts, cells, other assays
and reductions are unchanged.

## Details

Nearest distances are computed by
`GenomicRanges::distanceToNearest(select = "all", ignore.strand = TRUE)`.
Distance counts the bases between intervals; overlaps and adjacent
intervals have distance zero. These are distances to annotation
intervals, not to TSSs, promoter classifications, regulatory links, or
evidence of gene regulation. Different transcripts can produce several
native hits. `select = "all"` does not guarantee every mathematical tie
on both sides of a peak: native GenomicRanges tie selection is
preserved. No genome download or chromosome-style conversion is
performed. Existing annotation must be identical after native Signac
normalization (including transcript_id to tx_id); existing results are
not replaced. Genome labels and known lengths are checked, not reference
sequence contents. Nonmatching chromosome names remain unmatched rather
than being translated.

## See also

[`FilterATACPeaks()`](https://mengxu98.github.io/scop/reference/FilterATACPeaks.md),
[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md)

## Examples

``` r
counts <- Matrix::sparseMatrix(
  i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
  dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
)
ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
annotation <- GenomicRanges::GRanges("chr1:1-80")
GenomeInfoDb::genome(annotation) <- "hg38"
S4Vectors::mcols(annotation) <- S4Vectors::DataFrame(
  gene_id = "gene1", gene_name = "GENE1", gene_biotype = "protein_coding",
  type = "gene", tx_id = "tx1"
)
object <- AnnotateATACPeaks(object, annotation, verbose = FALSE)
SeuratObject::Misc(object[["peaks"]], "scop_peak_annotation")$nearest
#> DataFrame with 2 rows and 8 columns
#>           peak annotation_index  distance annotation_gene_id
#>    <character>        <integer> <integer>        <character>
#> 1    chr1-1-50                1         0              gene1
#> 2 chr1-101-150                1        20              gene1
#>   annotation_gene_name annotation_gene_biotype annotation_type annotation_tx_id
#>            <character>             <character>     <character>      <character>
#> 1                GENE1          protein_coding            gene              tx1
#> 2                GENE1          protein_coding            gene              tx1
```
