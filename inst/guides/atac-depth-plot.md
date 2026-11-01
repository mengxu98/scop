# ATAC reduction-depth correlations

Inspect LSI dimensions against actual count totals without trusting potentially
stale `nCount_peaks` metadata.

```r
library(scop)
# `atac` contains a peaks ChromatinAssay and an lsi reduction.
p <- ATACDepthCorPlot(atac, assay = "peaks", reduction = "lsi", n = 20)
print(p)
ggplot2::ggsave("lsi-depth.pdf", p, width = 7, height = 4)
correlations <- ATACDepthCorPlot(atac, assay = "peaks", reduction = "lsi",
  n = NULL, return.data = TRUE)
write.csv(correlations, "lsi-depth-correlations.csv", row.names = FALSE)
# Use another depth definition only when its provenance and units are verified:
fragment_cor <- ATACDepthCorPlot(atac, assay = "peaks", reduction = "lsi",
  depth.column = "total.fragments", method = "spearman", return.data = TRUE)
```

The default recomputes column sums of the selected assay's counts layer. These
are peak-matrix totals, not necessarily total fragments or sequencing reads.
Explicit metadata columns are user-selected alternatives; they are checked for
finite, non-negative values, not independently authenticated provenance.

Only selected-assay cells are considered and every selected cell must occur in
the reduction. The function does not silently drop missing cells. At least three
cells and variable depth are required. Constant components produce NA;
`n` is a maximum and NULL selects all components. SVD component signs are
arbitrary, so assess the magnitude as well as the sign of correlations. This
plot supports a decision about dimensions; it does not automatically exclude
LSI1 or establish a biological/technical explanation.

Focused tests cover native-style correlations, explicit cell ordering,
RNA-only cells, stale metadata, a single component, and invalid inputs.
