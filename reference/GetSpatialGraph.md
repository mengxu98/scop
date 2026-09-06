# Read or convert a stored spatial graph

Read a graph from a \`SpatialNetwork\` result without modifying the
Seurat object. Graphs can be returned as their complete list
representation, a sparse matrix, or a Seurat \`Graph\` object.

## Usage

``` r
GetSpatialGraph(
  object = NULL,
  res = NULL,
  graph.name = NULL,
  format = c("list", "sparse", "seurat"),
  value = c("weight", "distance")
)
```

## Arguments

- object:

  Optional \`Seurat\` object containing \`SpatialNetwork\` results.

- res:

  Optional \`SpatialNetwork\` result list.

- graph.name:

  Stored graph name. The active graph is used when \`NULL\`.

- format:

  Output representation.

- value:

  Edge value used for matrix conversions.

## Value

A graph list, \`dgCMatrix\`, or Seurat \`Graph\` object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- RunSpatialNetwork(
  visium_human_pancreas_sub,
  k = 6,
  coord.cols = c("x", "y"),
  verbose = FALSE
)
graph <- GetSpatialGraph(spatial, format = "list")
graph_sparse <- GetSpatialGraph(spatial, format = "sparse")
c(nodes = nrow(graph$nodes), edges = nrow(graph$edges))
#> nodes edges 
#>  1986  6358 
graph_sparse[
  seq_len(min(3, nrow(graph_sparse))),
  seq_len(min(3, ncol(graph_sparse)))
]
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>                    TGGTATCGGTCTGTAT-1 ATTATCTCGACAGATC-1 TGAGATCAAATACTCA-1
#> TGGTATCGGTCTGTAT-1                  .                  1                  .
#> ATTATCTCGACAGATC-1                  1                  .                  1
#> TGAGATCAAATACTCA-1                  .                  1                  .
```
