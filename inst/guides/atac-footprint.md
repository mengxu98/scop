# ATAC footprint

Requires indexed mapped fragments, an explicit loaded reference and bias interval. Computes native corrected footprints; native errors are preserved. Existing caches are not overwritten. Native synthetic parity is implementation acceptance, not biological evidence.

Interfaces: `RunATACFootprint`, `ATACFootprintPlot`. See each help page for arguments.

```r
# reference is a loaded BSgenome or indexed FaFile matching the assay.
footprinted <- RunATACFootprint(atac, assay = "peaks",
  genome = reference, genome.id = "hg38",
  bias.region = "chr1-4-20000000", motifs = "MA0139.1")
ATACFootprintPlot(footprinted, assay = "peaks",
  motifs = "MA0139.1", group.by = "celltype")
```

Choose a bias interval with enough fragment coverage for native hexamer
estimation. A small interval or sparse sample can fail; no uncorrected
fallback is substituted. Recompute after changing fragments, cells, peaks
or reference. A footprint does not establish physical TF binding.
