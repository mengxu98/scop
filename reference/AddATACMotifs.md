# Add motif matches and sequence statistics

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

  Whether to display progress messages.

## Value

The object with Signac Motifs and RegionStats in the selected assay.
Counts and cell order are unchanged. Misc(scop_atac_motifs) records
ordered peaks, motif IDs, genome label and dependency versions at
addition time. This receipt is a snapshot, not a replacement for current
assay state. Existing sequence-statistic columns are refreshed by name;
unrelated feature metadata are preserved. Explicit replacement emits the
receipt-overwrite warning. For `GRanges`, return a list with `motifs`,
aligned sequence `statistics`, the input `ranges`, and a
calculation-time `receipt`. Motif row names are canonical Signac
interval IDs. After recount over the same ranges, assign the Motif with
`Signac::Motifs(assay) <- motifs` and sequence statistics with
[`SeuratObject::AddMetaData()`](https://satijalab.github.io/seurat-object/reference/AddMetaData.html);
check exact peak order and reference identity before assignment. No
placeholder counts are constructed and changing counts alone does not
change sequence matches.

## See also

[`RunMotifEnrichment()`](https://mengxu98.github.io/scop/reference/RunMotifEnrichment.md)

## Examples

``` r
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
#> $genome
#> [1] "toy"
#> 
#> $sequence_class
#> [1] "DNAStringSet"
#> attr(,"package")
#> [1] "Biostrings"
#> 
#> $peaks
#>  [1] "chr1-1-100"     "chr1-121-220"   "chr1-241-340"   "chr1-361-460"  
#>  [5] "chr1-481-580"   "chr1-601-700"   "chr1-721-820"   "chr1-841-940"  
#>  [9] "chr1-961-1060"  "chr1-1081-1180" "chr1-1201-1300" "chr1-1321-1420"
#> [13] "chr1-1441-1540" "chr1-1561-1660" "chr1-1681-1780" "chr1-1801-1900"
#> [17] "chr1-1921-2020" "chr1-2041-2140" "chr1-2161-2260" "chr1-2281-2380"
#> [21] "chr1-2401-2500" "chr1-2521-2620" "chr1-2641-2740" "chr1-2761-2860"
#> [25] "chr1-2881-2980" "chr1-3001-3100" "chr1-3121-3220" "chr1-3241-3340"
#> [29] "chr1-3361-3460" "chr1-3481-3580" "chr1-3601-3700" "chr1-3721-3820"
#> [33] "chr1-3841-3940" "chr1-3961-4060" "chr1-4081-4180" "chr1-4201-4300"
#> [37] "chr1-4321-4420" "chr1-4441-4540" "chr1-4561-4660" "chr1-4681-4780"
#> 
#> $motifs
#> [1] "M1"
#> 
#> $backend
#>      Signac motifmatchr 
#>    "1.17.1"    "1.34.0" 
#> 
```
