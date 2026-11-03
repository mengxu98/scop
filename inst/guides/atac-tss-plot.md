# TSS enrichment curves

A per-cell `TSS.enrichment` number is not enough to reconstruct a curve. First
compute the positional matrix, using TSS coordinates from the correct reference.

```r
library(scop)
# `atac` contains indexed fragments; `tss` is the matched-genome TSS GRanges.
# Explicit tss.positions forces recalculation, even when scalar scores exist.
atac <- RunATACQC(atac, assay = "peaks",
  tss.positions = tss, fast = FALSE)
p <- ATACTSSPlot(atac, assay = "peaks", group.by = "celltype")
print(p)
ggplot2::ggsave("tss-enrichment.pdf", p, width = 8, height = 5)
# Existing group labels can restrict the display:
p_group <- ATACTSSPlot(atac, assay = "peaks", group.by = "celltype",
  idents = "B cell")
```

`Signac::TSSEnrichment(atac, assay = "peaks", tss.positions = tss,
fast = FALSE)` is an alternative preprocessing entry point. The plot never
silently calculates enrichment or turns scalar scores into artificial profiles.
It aligns the positional matrix to selected assay cells, supports a one-cell
selection, and excludes RNA-only cells from partial-assay multiome summaries.
Changing reference coordinates or preprocessing requires recomputation upstream.

This wrapper retains the upstream lifecycle warning: Signac 1.17 deprecates
`TSSPlot()`. It does not invent a replacement numerical estimator. The focused
tests compare positional aggregation with Signac and test missing/invalid data;
these software checks do not validate a biological reference or sample.
