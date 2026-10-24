# Link ATAC peaks to paired measured expression

Run native Signac peak-gene correlation and matched-background testing
on paired RNA and ATAC assays, retaining native link statistics.

## Usage

``` r
RunLinkPeaks(
  object,
  assay = NULL,
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10L,
  genes.use = NULL,
  n.sample = 200L,
  pvalue.cutoff = 0.05,
  score.cutoff = 0.05,
  gene.id = FALSE,
  seed = 11L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing paired chromatin and measured RNA assays.
  Identical cell names must identify the same biological cells. RNA must
  cover all cells of the selected ATAC assay; additional RNA-only cells
  are retained but excluded from this calculation. Missing pairs are
  rejected.

- assay:

  Chromatin assay; defaults to the active assay.

- expression.assay:

  Measured RNA assay, default `"RNA"`. Gene activity, imputed expression
  and chromVAR scores are not measured RNA.

- peak.slot, expression.slot:

  Single existing native layers, default `"counts"` and `"data"`. Split
  layers must be joined explicitly beforehand.

- method:

  Correlation method, `"pearson"` or `"spearman"`.

- gene.coords:

  Optional native gene-coordinate `GRanges` passed to
  [`Signac::LinkPeaks()`](https://stuartlab.org/signac/reference/LinkPeaks.html).
  With `NULL`, native Signac collapses assay transcript annotation.
  Coordinates must declare the same genome label as the assay; known
  chromosome lengths must agree. Explicit coordinates require one
  interval per gene identity. No chromosome-style conversion is
  performed.

- distance:

  Maximum genomic distance to TSS in bases, default 500000.

- min.distance:

  Optional non-negative minimum distance to TSS.

- min.cells:

  Native strict detection threshold: peaks and genes must occur in more
  than this many cells, default 10.

- genes.use:

  Optional unique RNA feature IDs to test.

- n.sample:

  Number of matched background peaks, default 200.

- pvalue.cutoff, score.cutoff:

  Native strict p-value and absolute correlation thresholds, default
  0.05. P-values are native unadjusted normal-tail values
  `pnorm(-abs(zscore))`, without a two-tail multiplier or
  multiple-testing correction; they are not adjusted significance
  estimates.

- gene.id:

  Whether RNA feature IDs match annotation gene_id rather than
  gene_name, default `FALSE`. Use `TRUE` for
  [`ReadMultiomeData()`](https://mengxu98.github.io/scop/reference/ReadMultiomeData.md)
  gene IDs.

- seed:

  Random seed used locally for native background sampling.

- verbose:

  Whether to print progress messages.

## Value

The Seurat object with native `Signac::Links(object[[assay]])` and a
receipt at `Misc(object[[assay]], "scop_link_peaks")`, recording source
assays, layers, cells, input features, parameters, backend version and
the number of retained links. Native link metadata contains gene, peak,
score, zscore and pvalue for non-empty results. Links are associations,
not causal or experimentally verified regulatory connections.

## Details

Run Signac RegionStats with an appropriate genome sequence before this
function; finite GC.percent between 0 and 100 and sequence.length equal
to peak width are required. These metadata checks do not verify the
supplied DNA provenance. Native accessibility background counts are
recomputed on the paired cells in a temporary assay, avoiding stale
count metadata after cell subsetting. Background sampling, correlation,
z scores and significance filtering are delegated unchanged to Signac.
Its minimum support and chromosome-background requirements still apply.
Source counts and the active assay are preserved. Existing native links
or SCOP link results are never overwritten. No retained links is a valid
empty result, not evidence that regulation is absent. Native errors for
constant expression or insufficient background support are preserved; no
alternative correlation or significance model is substituted.

## See also

[`ReadMultiomeData()`](https://mengxu98.github.io/scop/reference/ReadMultiomeData.md),
[`Signac::LinkPeaks()`](https://stuartlab.org/signac/reference/LinkPeaks.html),
[`Signac::RegionStats()`](https://stuartlab.org/signac/reference/RegionStats.html)

## Examples

``` r
if (requireNamespace("Signac", quietly = TRUE) &&
  requireNamespace("Biostrings", quietly = TRUE) &&
  requireNamespace("BSgenome", quietly = TRUE)) {
  local({
    withr::local_preserve_seed()
    set.seed(83)
    counts <- Matrix::Matrix(matrix(rpois(60 * 40, 4) + 1, 60,
      dimnames = list(paste0("p", 1:60), paste0("cell", 1:40))
    ), sparse = TRUE)
    ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
      seqnames = rep(c("chr1", "chr2"), each = 30),
      start = rep(seq(1, by = 200, length.out = 30), 2),
      end = rep(seq(100, by = 200, length.out = 30), 2)
    ))
    paired <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
    rna <- counts[c(1, 31), , drop = FALSE]
    rownames(rna) <- c("g1", "g2")
    paired[["RNA"]] <- SeuratObject::CreateAssayObject(counts = rna)
    genes <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
      seqnames = c("chr1", "chr2"), start = 201, end = 400, strand = "+",
      gene_id = c("g1", "g2"), gene_name = c("g1", "g2")
    ), keep.extra.columns = TRUE)
    GenomeInfoDb::genome(genes) <- "toy"
    # Small synthetic sequence solely for an executable method example.
    dna <- Biostrings::DNAStringSet(setNames(vapply(
      1:2, function(i) {
        paste(sample(c("A", "C", "G", "T"), 6000, replace = TRUE), collapse = "")
      },
      character(1)
    ), c("chr1", "chr2")))
    paired <- Signac::RegionStats(paired, assay = "peaks", genome = dna, verbose = FALSE)
    linked <- RunLinkPeaks(paired,
      assay = "peaks", expression.slot = "counts",
      gene.coords = genes, min.cells = 0, n.sample = 10, distance = 10000,
      pvalue.cutoff = 1, verbose = FALSE
    )
    Signac::Links(linked[["peaks"]])
  })
}
#> GRanges object with 49 ranges and 5 metadata columns:
#>        seqnames    ranges strand |      score        gene           peak
#>           <Rle> <IRanges>  <Rle> |  <numeric> <character>    <character>
#>    [1]     chr1    51-201      * |   1.000000          g1     chr1-1-100
#>    [2]     chr1   201-251      * |  -0.117405          g1   chr1-201-300
#>    [3]     chr1   201-451      * |   0.353564          g1   chr1-401-500
#>    [4]     chr1   201-651      * |  -0.133679          g1   chr1-601-700
#>    [5]     chr1   201-851      * |   0.133121          g1   chr1-801-900
#>    ...      ...       ...    ... .        ...         ...            ...
#>   [45]     chr2  201-4651      * |  0.1088618          g2 chr2-4601-4700
#>   [46]     chr2  201-4851      * | -0.0558111          g2 chr2-4801-4900
#>   [47]     chr2  201-5251      * | -0.0578105          g2 chr2-5201-5300
#>   [48]     chr2  201-5451      * | -0.1845793          g2 chr2-5401-5500
#>   [49]     chr2  201-5651      * |  0.0928417          g2 chr2-5601-5700
#>           zscore      pvalue
#>        <numeric>   <numeric>
#>    [1]  9.673625 1.95151e-22
#>    [2] -0.707541 2.39615e-01
#>    [3]  2.929233 1.69900e-03
#>    [4] -1.822162 3.42152e-02
#>    [5]  0.655145 2.56187e-01
#>    ...       ...         ...
#>   [45]  0.489887    0.312107
#>   [46] -0.511184    0.304611
#>   [47] -0.474493    0.317574
#>   [48] -0.869799    0.192205
#>   [49]  0.170479    0.432317
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
