# Add native motif matches and sequence statistics

Match supplied position frequency/weight matrices to the exact source
peak universe through Signac and motifmatchr.

## Usage

``` r
AddATACMotifs(
  object,
  genome,
  pfm,
  genome.id,
  assay = NULL,
  replace = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing a ChromatinAssay, or a `GRanges` of final
  peaks. The ranges-only path needs no counts or cells and can run
  before recount. Coordinates use 1-based closed intervals.

- genome:

  Loaded BSgenome, DNAStringSet or indexed FaFile. No download or
  installation is performed.

- pfm:

  A named TFBSTools PFMatrixList or PWMatrixList with unique motif IDs.

- genome.id:

  Explicit non-empty genome label matching the source assay and any
  known label on the loaded reference. For an unlabeled reference, this
  is the caller's provenance declaration, not an assembly inference.
  Known chromosome lengths must agree; source peaks must fit the
  reference.

- assay:

  Source chromatin assay; defaults to the active assay.

- replace:

  Whether to replace existing motif matches. Defaults to FALSE.

- verbose:

  Whether to display native progress messages.

## Value

The object with native Signac Motifs and RegionStats in the selected
assay. Counts and cell order are unchanged. Misc(scop_atac_motifs)
records ordered peaks, motif IDs, genome label and backend versions at
addition time. This receipt is a snapshot, not a replacement for current
native assay state. Existing sequence-statistic columns are refreshed by
name; unrelated feature metadata are preserved. Explicit replacement
emits the native receipt-overwrite warning. For `GRanges`, return a list
with native `motifs`, aligned sequence `statistics`, the input `ranges`,
and a calculation-time `receipt`. Motif row names are canonical Signac
interval IDs. After recount over the same ranges, assign the native
Motif with `Signac::Motifs(assay) <- motifs` and sequence statistics
with
[`SeuratObject::AddMetaData()`](https://satijalab.github.io/seurat-object/reference/AddMetaData.html);
check exact peak order and reference identity before assignment. No
placeholder counts are constructed and changing counts alone does not
change sequence matches.

## See also

[`RunMotifEnrichment()`](https://mengxu98.github.io/scop/reference/RunMotifEnrichment.md)

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
  SeuratObject::Misc(matched[["peaks"]], "scop_atac_motifs")
}
```
