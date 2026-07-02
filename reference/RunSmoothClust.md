# Run smoothclust spatial domain clustering

Smooth expression across spatial neighborhoods with the optional
`smoothclust` package, then cluster the smoothed profiles into spatial
domains with PCA and k-means.

## Usage

``` r
RunSmoothClust(
  srt,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("col", "row"),
  features = NULL,
  nfeatures = 2000,
  min_spots = 5,
  smooth_method = c("uniform", "kernel", "knn"),
  bandwidth = 0.05,
  k = 18,
  truncate = 0.05,
  n_threads = 1,
  n_clusters,
  n_pcs = 15,
  center = TRUE,
  scale = TRUE,
  nstart = 10,
  iter.max = 100,
  algorithm = "Hartigan-Wong",
  cluster_colname = "SmoothClust_cluster",
  tool_name = "SmoothClust",
  store_results = TRUE,
  store_smoothed = FALSE,
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay used for expression. If `NULL`, the default assay is used.

- layer:

  Assay layer used for expression values.

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- coord.cols:

  Metadata coordinate columns used when no Seurat image is available.

- features:

  Features to use. If `NULL`, current variable features are used; if no
  variable features are present, the top `nfeatures` by variance are
  used.

- nfeatures:

  Number of variance-ranked features to use when `features = NULL` and
  no variable features are present.

- min_spots:

  Minimum number of spots with non-zero expression required for a
  feature to be used.

- smooth_method:

  Smoothing method passed to
  [`smoothclust::smoothclust()`](https://rdrr.io/pkg/smoothclust/man/smoothclust.html).

- bandwidth, k, truncate, n_threads:

  Smoothing parameters passed to
  [`smoothclust::smoothclust()`](https://rdrr.io/pkg/smoothclust/man/smoothclust.html).

- n_clusters:

  Number of spatial domains for k-means clustering. This must be
  supplied explicitly.

- n_pcs:

  Number of principal components used for k-means.

- center, scale:

  Whether to center and scale features before PCA.

- nstart, iter.max, algorithm:

  Parameters passed to
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

- cluster_colname:

  Metadata column used for smoothclust clusters.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed results in `srt@tools`.

- store_smoothed:

  Whether to store the smoothed expression matrix in
  `srt@tools[[tool_name]]`. This can be large.

- seed:

  Random seed used for k-means.

- verbose:

  Whether to print progress messages.

- ...:

  Additional arguments passed to
  [`smoothclust::smoothclust()`](https://rdrr.io/pkg/smoothclust/man/smoothclust.html).

## Value

A `Seurat` object with smoothclust clusters in metadata. When
`store_results = TRUE`, detailed outputs are stored in
`srt@tools[[tool_name]]`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial$SmoothClust_cluster <- factor(
  paste0("SmoothClust", (seq_len(ncol(spatial)) - 1) %% 3 + 1)
)

SpatialSpotPlot(
  spatial,
  group.by = "SmoothClust_cluster",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  isTRUE(check_r("lmweber/smoothclust", verbose = FALSE))
) {
spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
spatial <- Seurat::FindVariableFeatures(
  spatial,
  assay = "Spatial",
  nfeatures = 200,
  verbose = FALSE
)

spatial <- RunSmoothClust(
  spatial,
  assay = "Spatial",
  n_clusters = 3,
  smooth_method = "knn",
  coord.cols = c("x", "y"),
  k = 6,
  verbose = FALSE
)

table(spatial$SmoothClust_cluster)
}
#> Error in check_r("lmweber/smoothclust", verbose = FALSE): could not find function "check_r"
```
