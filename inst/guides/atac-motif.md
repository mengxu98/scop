# ATAC motif

Uses explicit motif IDs in the selected assay. ggseqlogo is optional. Motif similarity is not proof of TF binding.

Interfaces: `ATACMotifPlot`. See each help page for arguments.

```r
# Use an existing motif ID from the selected assay PWM list.
ATACMotifPlot(atac, assay = "peaks", motifs = "MA0139.1")
```
