# Run a Giotto workflow

Run basic or full Giotto analysis on a \`giotto2\` object. Seurat input
is converted first with \[SeuratToScopGiotto()\].

## Usage

``` r
RunGiottoWorkflow(
  x,
  steps = c("basic", "full"),
  group.by = NULL,
  return_seurat = inherits(x, "Seurat"),
  store_results = TRUE,
  tool_name = "Giotto",
  verbose = TRUE,
  seed = 11,
  ...
)
```

## Arguments

- x:

  A \`giotto2\` or Seurat object.

- steps:

  \`"basic"\` runs preprocessing, PCA/UMAP, nearest-network clustering,
  and spatial network construction. \`"full"\` additionally runs spatial
  genes, spatial modules, optional cell proximity, and HMRF.

- group.by:

  Metadata column used for cell proximity enrichment.

- return_seurat:

  Whether to return a Seurat object when \`x\` is Seurat. If \`FALSE\`,
  returns the internal \`giotto2\` workflow object.

- store_results:

  Whether to store the internal Giotto workflow object in
  \`srt@tools\[\[tool_name\]\]\` when returning Seurat.

- tool_name:

  Name used to store the Giotto workflow object in \`srt@tools\`.

- ...:

  Passed to \[SeuratToScopGiotto()\] when \`x\` is Seurat.

## Value

A Seurat object by default for Seurat input, otherwise a \`giotto2\`
workflow object.
