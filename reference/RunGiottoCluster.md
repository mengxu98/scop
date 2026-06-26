# Run Giotto nearest-network clustering

Run Giotto as a temporary backend for nearest-network clustering and
return the complete Giotto object together with extracted cluster
results. The input `Seurat` object is not modified.

## Usage

``` r
RunGiottoCluster(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  method = c("leiden", "louvain"),
  dims = 1:20,
  k = 20,
  resolution = 1,
  cluster_colname = "Giotto_cluster",
  tool_name = "GiottoCluster",
  store_giotto = TRUE,
  conversion_params = list(),
  preprocess_params = list(),
  network_params = list(),
  cluster_params = list(),
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

  Features used for PCA and clustering. If `NULL`, current variable
  features are used, falling back to all assay features.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- method:

  Giotto clustering method.

- dims:

  Dimensions used to build the Giotto nearest-neighbor network.

- k:

  Number of nearest neighbors used by Giotto.

- resolution:

  Resolution passed to Giotto clustering.

- cluster_colname:

  Result column name recorded in returned parameters. This function does
  not write to `srt@meta.data`.

- tool_name:

  Result name recorded in returned parameters. This function does not
  write to `srt@tools`.

- store_giotto:

  Deprecated compatibility argument. The complete Giotto object is
  always returned in the `giotto` element.

- conversion_params:

  Additional parameters passed to `Giotto::createGiottoObject()`.

- preprocess_params:

  Additional parameters passed to `Giotto::runPCA()`.

- network_params:

  Additional parameters passed to `Giotto::createNearestNetwork()`.

- cluster_params:

  Additional parameters passed to Giotto clustering.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A `giotto2_result` list containing the full Giotto object, cluster
assignments, Giotto metadata, parameters, features, and cells.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
spatial <- Seurat::FindVariableFeatures(
  spatial,
  assay = "Spatial",
  nfeatures = 500,
  verbose = FALSE
)

giotto_clusters <- RunGiottoCluster(
  spatial,
  assay = "Spatial",
  layer = "data",
  dims = 1:10,
  k = 8,
  resolution = 0.4,
  coord.cols = c("col", "row")
)

head(giotto_clusters$clusters)
GiottoPlot(
  giotto_clusters,
  srt = spatial,
  overlay_image = FALSE,
  coord.cols = c("col", "row")
)
} # }
```
