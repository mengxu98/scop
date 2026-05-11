# Convert CCC results to OmicVerse communication AnnData

Convert CCC results to OmicVerse communication AnnData

## Usage

``` r
ccc_to_adata(
  srt = NULL,
  method = NULL,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  thresh = 0.05,
  liana_res = NULL,
  score_key = "score",
  pvalue_key = "pvalue",
  inverse_score = FALSE,
  inverse_pvalue = FALSE,
  sample_col = NULL,
  h5ad_path = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- method:

  Communication result type to use.

- condition:

  Result name or comparison name.

- dataset:

  Dataset index or name.

- slot.name:

  CellChat slot name.

- signaling:

  Signaling pathway to focus on.

- pairLR.use:

  Specific ligand-receptor pair(s) to keep.

- sender.use:

  Sender cell types to keep.

- receiver.use:

  Receiver cell types to keep.

- ligand.use:

  Ligands to keep.

- receptor.use:

  Receptors to keep.

- interaction.use:

  Interaction names to keep.

- thresh:

  Significance threshold used when extracting communication results.

- liana_res:

  Optional precomputed LIANA-like data frame. If supplied, `srt` is not
  required.

- score_key:

  Column in `liana_res` used for `layers[["means"]]`.

- pvalue_key:

  Column in `liana_res` used for `layers[["pvalues"]]`.

- inverse_score:

  Whether smaller `score_key` values should be converted to larger
  communication strengths. Useful for rank-like metrics.

- inverse_pvalue:

  Whether smaller `pvalue_key` values should be inverted before writing
  to `layers[["pvalues"]]`.

- sample_col:

  Optional column used as sample/context/dataset key. If `NULL`, the
  first available column among `"sample"`, `"context"`, `"condition"`,
  and `"dataset"` is used.

- h5ad_path:

  Optional output path. If provided, the AnnData object is written to
  this file before being returned.

- verbose:

  Whether to print progress messages.

## Value

A Python `anndata.AnnData` communication object compatible with
`ov.pl.ccc_*` plotting functions.
