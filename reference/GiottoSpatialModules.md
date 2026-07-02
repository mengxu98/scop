# Run Giotto spatial co-expression modules

Run Giotto spatial co-expression modules

## Usage

``` r
GiottoSpatialModules(
  x,
  features = NULL,
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  detect_params = list(),
  cluster_params = list(),
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

- cor_method:

  Correlation method used by Giotto.

- k:

  Number of spatial co-expression modules.

- detect_params:

  Additional parameters passed to \`Giotto::detectSpatialCorFeats()\`.

- cluster_params:

  Additional parameters passed to \`Giotto::clusterSpatialCorFeats()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
features <- rownames(spatial)[1:6]
module_table <- expand.grid(
  feat_ID = features,
  variable = paste0("module_", 1:3),
  stringsAsFactors = FALSE
)
module_table$spat_cor <- seq(-0.7, 0.8, length.out = nrow(module_table))
g <- structure(
  list(
    source = list(cells = colnames(spatial), features = features),
    results = list(
      spatial_modules = list(
        module_tables = list(result.cor_DT = module_table),
        features = features
      )
    )
  ),
  class = c("giotto2", "list")
)
GiottoPlot(g, plot_type = "spatial_modules", top_n = 6)


if (
  isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE))
) {
  g <- SeuratToScopGiotto(spatial, coord.cols = c("x", "y"))
  g <- GiottoPreprocess(g)
  g <- GiottoSpatialNetwork(g)
  g <- GiottoSpatialModules(g, features = rownames(spatial)[1:50], k = 3)
}
#> Error in check_r("giotto-suite/Giotto", verbose = FALSE): could not find function "check_r"
```
