# Infer coaccessibility using native Cicero

Infer coaccessibility using native Cicero.

## Usage

``` r
RunATACCoaccessibility(object, genome.sizes, genome.id, assay = NULL,
  reduction = "lsi", dims = 2:30, k = 50L, window = 500000L,
  sample.num = 100L, seed = 1L,
  result.name = "scop_coaccessibility", verbose = TRUE)
```

## Arguments

- object:

  Seurat object with ATAC counts.

- genome.sizes:

  Named positive integer chromosome lengths for the explicit reference.

- genome.id:

  Verified reference label matching assay.

- assay:

  Chromatin assay.

- reduction:

  Existing embedding covering all assay cells.

- dims:

  Unique embedding dimensions.

- k:

  Cells per aggregate, at least 2 and less than cell count.

- window:

  Native window size greater than 250000.

- sample.num:

  Number of windows for distance estimation.

- seed:

  Local seed; caller RNG restored.

- result.name:

  Unused assay Misc key.

- verbose:

  Print native progress.

## Value

Seurat object with native connection table, aggregate membership and
provenance in assay Misc. Native Links and source assays are unchanged.

## Details

Requires official monocle3-compatible Cicero and monocle3. Native
binary-accessibility aggregation, size-factor normalization and
distance-penalized graphical modelling are used, without a correlation
fallback. Native NA scores and model errors are preserved. Small contigs
or inadequate peaks may not support estimation. These associations do
not establish physical interactions or causal regulation.
