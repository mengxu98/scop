# Projection Plot

This function generates a projection plot, which can be used to compare
two groups of cells in a dimensionality reduction space.

## Usage

``` r
ProjectionPlot(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_reduction = "ref.embeddings",
  ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||% NULL,
  query_param = list(palette = "Set1", cells.highlight = TRUE),
  ref_param = list(palette = "Chinese"),
  xlim = NULL,
  ylim = NULL,
  pt.size = 0.8,
  stroke.highlight = 0.5,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  A Seurat object or count matrix representing the reference object. If
  provided, the similarities will be calculated between cells from the
  query and reference objects. If not provided, the similarities will be
  calculated within the query object.

- query_group:

  The grouping variable for the query group cells.

- ref_group:

  The grouping variable for the reference group cells.

- query_reduction:

  The name of the reduction in the query group cells.

- ref_reduction:

  The name of the reduction in the reference group cells.

- query_param:

  A list of parameters for customizing the query group plot. Available
  parameters: palette (color palette for groups) and cells.highlight
  (whether to highlight cells).

- ref_param:

  A list of parameters for customizing the reference group plot.
  Available parameters: palette (color palette for groups) and
  cells.highlight (whether to highlight cells).

- xlim:

  The x-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- ylim:

  The y-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- pt.size:

  The size of the points in the plot.

- stroke.highlight:

  The size of the stroke highlight for cells.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
CellDimPlot(
  srt_ref,
  group.by = c("celltype", "tech")
)

# Projection
srt_query <- RunKNNMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_umap = "UncorrectedUMAP2D"
)
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
```
