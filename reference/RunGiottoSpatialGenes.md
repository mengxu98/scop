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
