# Test motif enrichment within an explicit peak universe

Use native Signac background matching and motif enrichment with an
auditable foreground and population. Sequence scanning remains available
directly via
[`Signac::AddMotifs()`](https://stuartlab.org/signac/reference/AddMotifs.html)
and does not require this function or a DAR analysis.

## Usage

``` r
RunMotifEnrichment(
  object,
  features,
  assay = NULL,
  universe = NULL,
  background = NULL,
  n.background = 40000L,
  match.features = c("GC.percent", "sequence.length"),
  seed = 1L,
  p.adjust.method = "BH",
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing a validated ChromatinAssay and an aligned
  native Signac Motif object.

- features:

  Unique foreground peak IDs, all present in `universe`.

- assay:

  Source chromatin assay; defaults to the active assay.

- universe:

  Eligible peak IDs. With NULL, use all source-assay peaks. For
  DAR-derived foregrounds, explicitly supply the full tested peak
  universe.

- background:

  Optional explicit population peak IDs, including every foreground peak
  and contained in `universe`. With NULL, select reference peaks using
  native
  [`Signac::MatchRegionStats()`](https://stuartlab.org/signac/reference/MatchRegionStats.html),
  then include the foreground in the tested population.

- n.background:

  Maximum number of non-foreground reference peaks selected when
  matching. If fewer are eligible, all eligible peaks are requested.

- match.features:

  Feature-metadata columns passed to native background matching. Compute
  verified sequence metadata with
  [`Signac::RegionStats()`](https://stuartlab.org/signac/reference/RegionStats.html).

- seed:

  Random seed for native matching, with the caller's RNG restored.

- p.adjust.method:

  Multiple-testing method passed to
  [`Signac::FindMotifs()`](https://stuartlab.org/signac/reference/FindMotifs.html).

- verbose:

  Whether native backends print progress.

## Value

A list with the complete native `results` table and a `receipt` with
exact foreground, eligible universe, tested population, matching
parameters, motif IDs and backend versions. The input object is
unchanged.

## Details

The hypergeometric population includes the foreground, so the foreground
is a subset of the tested population. `reference` in the receipt
contains only the matched non-foreground controls; `background` contains
the actual native test population. This differs from passing a numeric
background directly to Signac, whose matched list excludes the query. No
enrichment formula is reimplemented. All motifs are tested before any
downstream display filtering. P-values test one-sided over-enrichment; a
fold enrichment below one is not an under-enrichment test. For
directional DAR analyses, select up/down foregrounds separately. An
enrichment result is a sequence association, not evidence of TF binding.
Native undefined fold enrichment for a zero-hit motif is retained,
together with its native p-value and multiple-testing result.

## See also

[`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)

## Examples

``` r
if (all(vapply(c("Biostrings", "BSgenome", "TFBSTools", "motifmatchr"),
  requireNamespace, logical(1),
  quietly = TRUE
))) {
  reference <- Biostrings::DNAStringSet(c(chr1 = paste(rep("ACGT", 1200), collapse = "")))
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1", start = seq(1, by = 120, length.out = 40),
    end = seq(100, by = 120, length.out = 40)
  ))
  counts <- Matrix::Matrix(matrix(1, 40, 2,
    dimnames = list(paste0("p", seq_len(40)), c("A", "B"))
  ), sparse = TRUE)
  atac <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  profile <- matrix(0L, 4, 8, dimnames = list(c("A", "C", "G", "T"), NULL))
  profile[cbind(rep(seq_len(4), 2), seq_len(8))] <- 100L
  pfm <- TFBSTools::PFMatrixList(M1 = TFBSTools::PFMatrix(
    ID = "M1", name = "toy_ACGT", profileMatrix = profile
  ))
  matched <- AddATACMotifs(atac, reference, pfm, genome.id = "toy", verbose = FALSE)
  peaks <- rownames(matched[["peaks"]])
  enrichment <- RunMotifEnrichment(matched, peaks[1:10],
    universe = peaks, background = peaks, verbose = FALSE
  )
  enrichment$results
}
```
