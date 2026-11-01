# Regional ATAC fragment-length plots

This feature is independent of TSS/depth plotting. It uses Signac's native
fragment extraction, including each capture's object-cell to raw-barcode map.
It requires real local bgzip fragment files plus a `.tbi` or `.csi` index.

```r
library(scop)
# `atac` is an existing Seurat object with an indexed fragment-backed peaks assay.
p <- ATACFragmentPlot(atac, assay = "peaks",
  region = "chr1-1-2000000", group.by = "sample")
print(p)
ggplot2::ggsave("fragment-lengths.pdf", p, width = 8, height = 5)
# Optional cell restriction, using object cell IDs, not raw barcodes:
p_subset <- ATACFragmentPlot(atac, assay = "peaks",
  region = "chr1-1-2000000", cells = colnames(atac[["peaks"]])[1:20])
```

Choose a region matching the actual reference genome and chromosome labels.
The example region is not an automatic genome choice. The histogram is a
regional fragment-record distribution, not a full-library estimate or a
read-multiplicity-weighted count. Native display limits are 0–800 bp. A fragment
returned by multiple query intervals can be counted more than once; use a single
region for an unambiguous diagnostic. Missing mappings/files and empty results
produce actionable errors. Plotting does not modify the input object, including
barcode hash-cache attributes.

The focused test file constructs indexed synthetic fragments, compares native
Signac results, and exercises reused raw barcodes in different captures. These
are software tests, not claims about biological sample quality.
