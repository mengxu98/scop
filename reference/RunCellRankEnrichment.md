# Enrich CellRank trend modules

Run over-representation analysis for each stored CellRank trend module
and keep complete tables and an execution manifest in the Seurat object.

## Usage

``` r
RunCellRankEnrichment(
  srt,
  lineage,
  db = c("MSigDB_C2", "GO_BP", "GO_CC", "GO_MF", "MSigDB_H"),
  species = "Mus_musculus",
  universe = NULL,
  minGSSize = 10L,
  maxGSSize = 500L,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  p_adjust_method = "BH",
  show_category = 8L,
  output_dir = NULL,
  continue_on_error = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by \[RunCellRankTrends\].

- lineage:

  CellRank lineage whose trend modules should be enriched.

- db:

  Annotation databases accepted by \[PrepareDB\].

- species:

  Species passed to \[PrepareDB\].

- universe:

  Background genes. \`NULL\` uses genes tested for CellRank lineage
  drivers.

- minGSSize:

  Minimum gene-set size.

- maxGSSize:

  Maximum gene-set size.

- pvalue_cutoff:

  Nominal enrichment cutoff.

- qvalue_cutoff:

  Adjusted enrichment cutoff.

- p_adjust_method:

  Multiple-testing method.

- show_category:

  Number of terms shown in each dot plot.

- output_dir:

  Optional directory for CSV/PDF exports.

- continue_on_error:

  Whether to keep other modules when one enrichment call fails. Failures
  are recorded in the manifest.

- verbose:

  Whether to print progress messages.

## Value

The Seurat object with results under
\`srt@tools\$CellRank\$enrichment\[\[lineage\]\]\`.
