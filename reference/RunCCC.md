# Run common cell-cell communication analyses

Run common cell-cell communication analyses

## Usage

``` r
RunCCC(
  srt,
  group.by = NULL,
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

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- methods:

  Supported cell-cell communication methods to run. The default core
  methods are `"CellChat"`, `"CellphoneDB"`, and `"LIANA"`. NicheNet,
  MultiNicheNet, SpatialCellChat, SpatialDM, and MDIC3 can be selected
  when their design-specific arguments are supplied through
  `method_params`. LIANA's default internal method set includes its
  CellPhoneDB scorer, so standalone CellphoneDB and LIANA consensus
  results are not statistically independent evidence.

- method_params:

  Named list of method-specific arguments passed to the corresponding
  wrapper. For example, use `method_params$CellphoneDB$pvalue` for
  CellphoneDB-specific parameters.

- backend:

  Backend used only for result post-processing and unified CCC table
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
