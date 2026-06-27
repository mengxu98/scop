# Run Giotto dimensional reduction

Run Giotto dimensional reduction

## Usage

``` r
GiottoReduce(
  x,
  reduction = c("pca", "umap"),
  dims = 1:20,
  name = NULL,
  features = NULL,
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- reduction:

  Dimensional reduction to run.

- dims:

  Dimensions to use.

- name:

  Name for the Giotto reduction.

- features:

  Features used for the reduction.

- params:

  Additional parameters passed to the Giotto reduction function.

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
embedding <- cbind(
  UMAP_1 = as.numeric(scale(spatial$x)),
  UMAP_2 = as.numeric(scale(spatial$y))
)
rownames(embedding) <- colnames(spatial)
g <- structure(
  list(
    giotto = list(umap = embedding),
    source = list(cells = colnames(spatial), features = rownames(spatial)),
    results = list(
      cluster = list(
        table = data.frame(
          cell = colnames(spatial),
          cluster = spatial$coda_label,
          row.names = colnames(spatial)
        )
      )
    ),
    parameters = list(umap_name = "umap")
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "dim")


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoPreprocess(g)
  g <- GiottoReduce(g, reduction = "pca", dims = 1:10)
  g <- GiottoReduce(g, reduction = "umap", dims = 1:10)
}
```
