# Run LIANA cell-cell communication analysis

Run LIANA cell-cell communication analysis

## Usage

``` r
RunLIANA(
  srt,
  group.by,
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus",
  assay = NULL,
  min_cells = 5,
  return_all = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Metadata column defining cell groups. Passed to `liana::liana_wrap()`
  as `idents_col`.

- method:

  LIANA methods to run. Defaults to LIANA's internal methods.

- resource:

  LIANA ligand-receptor resource(s). Default is `"Consensus"`.

- assay:

  Assay used by LIANA. If `NULL`, LIANA uses the default assay.

- min_cells:

  Minimum cells per identity retained by LIANA.

- return_all:

  Whether LIANA should return all possible interactions.

- backend:

  Backend used for scop post-processing and unified CCC table
  aggregation. Upstream LIANA inference is unchanged.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to `liana::liana_wrap()`.

## Value

A `Seurat` object with LIANA results stored in `srt@tools[["LIANA"]]`.
