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

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

- ...:

  Passed to \[SeuratToScopGiotto()\] when \`x\` is Seurat.

## Value

A Seurat object by default for Seurat input, otherwise a \`giotto2\`
workflow object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
g <- structure(
  list(
    source = list(
      cells = colnames(spatial),
      features = rownames(spatial),
      coordinates = data.frame(
        cell_ID = colnames(spatial),
        sdimx = spatial$x,
        sdimy = spatial$y
      )
    ),
    results = list(
      cluster = list(
        table = data.frame(
          cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
          row.names = colnames(spatial)
        )
      )
    ),
    active = "cluster"
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "cluster")


if (
  isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE))
) {
g <- RunGiottoWorkflow(
  spatial,
  steps = "basic",
  assay = "Spatial",
  layer = "counts",
  coord.cols = c("x", "y"),
  return_seurat = FALSE,
  verbose = FALSE
)
}
#> Error in check_r("giotto-suite/Giotto", verbose = FALSE): could not find function "check_r"
```
