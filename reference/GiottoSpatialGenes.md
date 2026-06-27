# Run Giotto spatial gene detection

Run Giotto spatial gene detection

## Usage

``` r
GiottoSpatialGenes(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  top_n = 100,
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- features:

  Features to test.

- network_method:

  Spatial network method.

- network_name:

  Name for the Giotto spatial network.

- bin_method:

  Binarization method passed to \`Giotto::binSpect()\`.

- top_n:

  Number of top spatial genes to store.

- params:

  Additional parameters passed to \`Giotto::binSpect()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:80],
  features = rownames(visium_human_pancreas_sub)[1:200]
)
#> Warning: Not validating Centroids objects
#> Warning: Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating Seurat objects
coords <- data.frame(
  cell_ID = colnames(spatial),
  sdimx = spatial$x,
  sdimy = spatial$y,
  row.names = colnames(spatial)
)
spatial_gene_table <- data.frame(
  feat_ID = rownames(spatial)[1:8],
  spatGeneRank = seq_len(8),
  adj.p.value = seq(0.001, 0.04, length.out = 8)
)
g <- structure(
  list(
    source = list(cells = colnames(spatial), coordinates = coords),
    results = list(
      spatial_genes = list(
        table = spatial_gene_table,
        top_features = spatial_gene_table$feat_ID
      )
    )
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "spatial_genes", top_n = 6)


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoPreprocess(g)
  g <- GiottoSpatialNetwork(g)
  g <- GiottoSpatialGenes(g, features = rownames(spatial)[1:50], top_n = 10)
}
```
