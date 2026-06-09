# Run common cell-cell communication analyses

Run common cell-cell communication analyses

## Usage

``` r
RunCCC(
  srt,
  group.by,
  methods = c("CellChat", "CellphoneDB", "LIANA"),
  method_params = list(),
  backend = c("cpp", "r"),
  skip_failed = FALSE,
  rebuild_unified = TRUE,
  thresh = 0.05,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- methods:

  Cell-cell communication methods to run. Currently supports
  `"CellChat"`, `"CellphoneDB"`, and `"LIANA"` in the unified scheduler.
  NicheNet and MultiNicheNet require explicit receiver/sender/contrast
  design arguments and should be called through
  [`RunNichenetr()`](https://mengxu98.github.io/scop/reference/RunNichenetr.md)
  or
  [`RunMultiNichenetr()`](https://mengxu98.github.io/scop/reference/RunMultiNichenetr.md).

- method_params:

  Named list of method-specific arguments passed to the corresponding
  wrapper. For example, use `method_params$CellphoneDB$pvalue` for
  CellphoneDB-specific parameters.

- backend:

  Backend used for scop post-processing and unified CCC table
  aggregation. The upstream CellChat, CellphoneDB, and LIANA inference
  logic is unchanged.

- skip_failed:

  Whether to keep running remaining methods if one method fails.

- rebuild_unified:

  Whether to rebuild `srt@tools[["CCC"]]` from the completed methods
  after all requested methods finish.

- thresh:

  Significance threshold used when rebuilding unified CCC tables and
  passed to
  [`RunCellChat()`](https://mengxu98.github.io/scop/reference/RunCellChat.md)
  unless overridden in `method_params$CellChat`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with method-specific results and a unified
`srt@tools[["CCC"]]` bundle.
