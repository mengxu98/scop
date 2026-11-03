# ATAC coaccessibility

Requires official monocle3-compatible Cicero. No Pearson fallback. Results and native NA scores remain separate from peak-gene Links. Window must exceed 250000. Opt-in native test: SCOP_RUN_NATIVE_CICERO=true.

Interfaces: `RunATACCoaccessibility`. See each help page for arguments.

```r
# chromosome_lengths is an explicit named vector for the verified assembly.
coaccessible <- RunATACCoaccessibility(atac, assay = "peaks",
  genome.sizes = chromosome_lengths, genome.id = "hg38",
  reduction = "lsi", dims = 2:20, k = 50,
  window = 500000, sample.num = 100, seed = 1)
result <- SeuratObject::Misc(coaccessible[["peaks"]], "scop_coaccessibility")
head(result$connections)
```

The function never installs packages or downloads references. Inadequate
genomic coverage can prevent native distance estimation. Associations
are not experimentally validated chromatin contacts or causal regulation.
