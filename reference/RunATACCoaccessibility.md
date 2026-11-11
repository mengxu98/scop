# Infer peak coaccessibility with Cicero

Run the official monocle3-compatible Cicero workflow on the selected
ATAC counts and an explicitly selected cell embedding. No Pearson
fallback is substituted if Cicero or its model fitting is unavailable.

## Usage

``` r
RunATACCoaccessibility(
  object,
  genome.sizes,
  genome.id,
  assay = NULL,
  reduction = "lsi",
  dims = 2:30,
  k = 50L,
  window = 500000L,
  sample.num = 100L,
  seed = 1L,
  result.name = "scop_coaccessibility",
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat object with a ChromatinAssay.

- genome.sizes:

  Named positive integer chromosome lengths for the verified reference,
  covering all source peaks. No reference is downloaded.

- genome.id:

  Reference label matching the selected assay.

- assay:

  ATAC assay; defaults to active assay.

- reduction:

  Existing reduction whose cells cover every assay cell.

- dims:

  Explicit dimensions used for nearest-neighbor aggregation.

- k:

  Number of cells per aggregate; at least 2 and less than cell count.

- window:

  Genomic window length, greater than 250000 because the run_cicero
  distance constraint is fixed at 250000.

- sample.num:

  Number of windows sampled for distance estimation.

- seed:

  Local random seed. Caller RNG state is restored.

- result.name:

  Unused assay Misc key for the result.

- verbose:

  Print progress messages.

## Value

Input Seurat object with a result list in assay Misc(result.name),
containing connections, aggregate membership, parameters and source
identifiers. Links, source counts, metadata and default assay are
unchanged. Receipt identifiers are calculation-time snapshots.

## Details

Requires the official Cicero monocle3 branch and monocle3. See the
[official installation
instructions](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/).
The installed Cicero package must declare a monocle3 dependency. Binary
accessibility is aggregated by Cicero before size-factor normalization
and distance-penalized graphical-model fitting. Model fitting failures
(including inadequate genomic coverage) are preserved. Peak IDs and cell
identities are not renamed or silently dropped. Returned NA scores are
unresolved/discordant estimates, not zero. Associations are not
experimentally verified physical interactions or causal regulation.
