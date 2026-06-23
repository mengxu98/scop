# Run Giotto spatial co-expression modules

Use Giotto `detectSpatialCorFeats()` and `clusterSpatialCorFeats()` as a
temporary backend for feature-level spatial co-expression modules. The
complete Giotto object and module results are returned as a standalone
result; the input `Seurat` object is not modified.

## Usage

``` r
RunGiottoSpatialModules(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  tool_name = "GiottoSpatialModules",
  store_giotto = TRUE,
  conversion_params = list(),
  network_params = list(),
  detect_params = list(),
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

  Features to test for spatial co-expression modules. If `NULL`, current
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

- cor_method:

  Correlation method used by Giotto.

- k:

  Number of feature modules passed to Giotto clustering.

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

- detect_params:

  Additional parameters passed to `Giotto::detectSpatialCorFeats()`.

- cluster_params:

  Additional parameters passed to `Giotto::clusterSpatialCorFeats()`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A `giotto2_result` list containing the full Giotto object, spatial
correlation object, module object, extracted module tables, parameters,
features, and cells.
