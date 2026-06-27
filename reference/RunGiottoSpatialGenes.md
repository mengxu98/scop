# Run Giotto spatial gene detection

Use Giotto `binSpect()` as a temporary backend for spatially variable
gene detection. The complete Giotto object and result tables are
returned as a standalone result; the input `Seurat` object is not
modified.

## Usage

``` r
RunGiottoSpatialGenes(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  set_variable_features = FALSE,
  top_n = 100,
  tool_name = "GiottoSpatialGenes",
  store_giotto = TRUE,
  conversion_params = list(),
  network_params = list(),
  binSpect_params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the expression matrix.

- features:

  Features to test with `Giotto::binSpect()`. If `NULL`, current
  variable features are used, falling back to all assay features.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- network_method:

  Spatial network method passed to `Giotto::createSpatialNetwork()`.

- network_name:

  Name for the Giotto spatial network.

- bin_method:

  Binarization method passed to `Giotto::binSpect()`.

- set_variable_features:

  Deprecated compatibility argument. Seurat variable features are never
  modified by this function.

- top_n:

  Number of top genes to store.

- tool_name:

  Result name recorded in returned parameters. This function does not
  write to `srt@tools`.

- store_giotto:

  Deprecated compatibility argument. The complete Giotto object is
  always returned in the `giotto` element.

- conversion_params:

  Additional parameters passed to `Giotto::createGiottoObject()`.

- network_params:

  Additional parameters passed to `Giotto::createSpatialNetwork()`.

- binSpect_params:

  Additional parameters passed to `Giotto::binSpect()`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A `giotto2_result` list containing the full Giotto object, spatial gene
table, top features, raw Giotto result, parameters, features, and cells.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:120],
  features = rownames(visium_human_pancreas_sub)[1:400]
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
spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
giotto_genes <- list(
  results = data.frame(
    feat_ID = rownames(spatial)[1:6],
    spatGeneRank = c(40, 35, 28, 20, 16, 10)
  ),
  top_features = rownames(spatial)[1:4],
  parameters = list(assay = "Spatial", layer = "data", coord.cols = c("x", "y"))
)
class(giotto_genes) <- c("giotto2_spatial_genes", "giotto2_result", "list")

head(giotto_genes$results)
#>   feat_ID spatGeneRank
#> 1  TMSB4X           40
#> 2     UBC           35
#> 3     GCG           28
#> 4    ACTB           20
#> 5  COL3A1           16
#> 6  COL1A1           10
GiottoPlot(giotto_genes, plot_type = "ranking", top_n = 6)

GiottoPlot(
  giotto_genes,
  srt = spatial,
  plot_type = "feature",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
spatial <- Seurat::FindVariableFeatures(
  spatial,
  assay = "Spatial",
  nfeatures = 300,
  verbose = FALSE
)
giotto_genes <- RunGiottoSpatialGenes(
  spatial,
  assay = "Spatial",
  layer = "data",
  features = Seurat::VariableFeatures(spatial, assay = "Spatial"),
  coord.cols = c("x", "y"),
  top_n = 50
)
}
```
