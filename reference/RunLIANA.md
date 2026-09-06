# Run LIANA cell-cell communication analysis

Run LIANA cell-cell communication analysis

## Usage

``` r
RunLIANA(
  srt,
  group.by,
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = NULL,
  assay = NULL,
  min_cells = 5,
  return_all = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE,
  species = c("human", "mouse"),
  consensus = c("auto", "rank", "aggregate", "none"),
  consensus_args = list(),
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- group.by:

  Metadata column defining cell groups. Passed to `liana::liana_wrap()`
  as `idents_col`.

- method:

  LIANA methods to run. Defaults to LIANA's internal methods.

- resource:

  LIANA ligand-receptor resource(s). If `NULL`, `"Consensus"` is used
  for human data and `"MouseConsensus"` for mouse data. This is a
  ligand-receptor resource choice, not a multi-method result consensus.

- assay:

  Assay used by LIANA. If `NULL`, LIANA uses the default assay.

- min_cells:

  Minimum cells per identity retained by LIANA.

- return_all:

  Whether LIANA should return all possible interactions.

- backend:

  Backend used only for result post-processing and unified CCC table
  aggregation. Upstream LIANA inference is unchanged.

- verbose:

  Whether to print the message. Default is `TRUE`.

- species:

  Species used to select the default ligand-receptor resource.

- consensus:

  Multi-method result aggregation. `"auto"` uses `rank_aggregate()` when
  at least two methods are requested and does not invent a consensus for
  a single method. `"rank"` returns magnitude and specificity consensus
  ranks, `"aggregate"` uses `liana_aggregate()`, and `"none"` keeps only
  method-specific results.

- consensus_args:

  Named arguments passed to the selected LIANA aggregation function.

- ...:

  Additional arguments passed to `liana::liana_wrap()`.

## Value

A `Seurat` object with LIANA results stored in `srt@tools[["LIANA"]]`.
