# Run CellChat analysis

Run CellChat analysis

## Usage

``` r
RunCellChat(
  srt,
  group.by,
  species = c("Homo_sapiens", "Mus_musculus", "zebrafish"),
  split.by = NULL,
  annotation_selected = NULL,
  group_column = NULL,
  group_cmp = NULL,
  thresh = 0.05,
  min.cells = 10,
  assay = NULL,
  layer = "data",
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- species:

  The species of the data, either `"Homo_sapiens"`, `"Mus_musculus"`, or
  `"zebrafish"`.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- annotation_selected:

  A vector of cell annotations of interest for running the `CellChat`
  analysis. If not provided, all cell types will be considered.

- group_column:

  Name of the metadata column in the `Seurat` object that defines
  conditions or groups.

- group_cmp:

  A list of pairwise condition comparisons for differential `CellChat`
  analysis.

- thresh:

  The threshold for computing centrality scores. Default is `0.05`.

- min.cells:

  the minmum number of expressed cells required for the genes that are
  considered for cell-cell communication analysis. Default is `10`.

- assay:

  Which assay to use. If `NULL`, the default assay of the `Seurat`
  object will be used.

- layer:

  The layer to use for the expression data. Default is `"data"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with `CellChat` results stored in
`srt@tools[["CellChat"]]`.

## References

[CellChat](https://github.com/jinworks/CellChat)

## See also

[CCCHeatmap](https://mengxu98.github.io/scop/reference/CCCHeatmap.md),
[CCCStatPlot](https://mengxu98.github.io/scop/reference/CCCStatPlot.md),
[CCCNetworkPlot](https://mengxu98.github.io/scop/reference/CCCNetworkPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)

CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  plot_type = "bipartite"
)

CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  plot_type = "heatmap"
)

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  plot_type = "violin",
  top_n = 50
)
} # }
```
