# Run BayesSpace spatial clustering

Run BayesSpace spatial clustering

## Usage

``` r
RunBayesSpace(
  srt,
  q,
  assay = NULL,
  platform = c("Visium", "VisiumHD", "ST"),
  image = NULL,
  use_reduction = NULL,
  dims = 1:15,
  preprocess = TRUE,
  n.PCs = 15,
  n.HVGs = 2000,
  skip.PCA = !is.null(use_reduction),
  spatial_preprocess_params = list(),
  spatial_cluster_params = list(),
  cluster_colname = "BayesSpace_cluster",
  init_colname = "BayesSpace_init",
  store_sce = TRUE,
  verbose = TRUE,
  coord.cols = c("col", "row")
)
```

## Arguments

- srt:

  A `Seurat` object.

- q:

  Number of BayesSpace clusters.

- assay:

  Assay to use. `NULL` uses the default assay.

- platform:

  Spatial sequencing platform.

- image:

  Name of the Seurat spatial image used to recover spot coordinates when
  they are not already present in metadata. For regular Visium data with
  only pixel `x`/`y` coordinates, BayesSpace array coordinates are
  inferred from the spatial grid.

- use_reduction:

  Optional Seurat reduction to pass to BayesSpace as PCA.

- dims:

  Dimensions from `use_reduction` to use.

- preprocess:

  Whether to run
  [`BayesSpace::spatialPreprocess()`](https://rdrr.io/pkg/BayesSpace/man/spatialPreprocess.html).

- n.PCs, n.HVGs:

  Parameters passed to `spatialPreprocess()`.

- skip.PCA:

  Whether to skip PCA inside `spatialPreprocess()`.

- spatial_preprocess_params:

  Additional parameters passed to
  [`BayesSpace::spatialPreprocess()`](https://rdrr.io/pkg/BayesSpace/man/spatialPreprocess.html).

- spatial_cluster_params:

  Additional parameters passed to
  [`BayesSpace::spatialCluster()`](https://rdrr.io/pkg/BayesSpace/man/spatialCluster.html).

- cluster_colname:

  Metadata column used for BayesSpace clusters.

- init_colname:

  Metadata column used for BayesSpace initial clusters.

- store_sce:

  Whether to store the BayesSpace `SingleCellExperiment` in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- coord.cols:

  Two metadata columns containing raw x/y coordinates when no spatial
  image is available.

## Value

A `Seurat` object with BayesSpace clusters in metadata and raw results
in `srt@tools[["BayesSpace"]]`.

## Examples

``` r
data(visium_human_pancreas_results_sub)
spatial <- RunBayesSpace(
  visium_human_pancreas_results_sub,
  q = 3,
  n.PCs = 5,
  n.HVGs = 200,
  store_sce = FALSE,
  spatial_cluster_params = list(nrep = 100, burn.in = 10, thin = 5),
  coord.cols = c("x", "y"),
  verbose = FALSE
)
#> Neighbors were identified for 132 out of 400 spots.
#> Fitting model...
#> Calculating labels using iterations 11 through 100.
SpatialSpotPlot(
  spatial,
  group.by = "BayesSpace_cluster",
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  pt.size = 1.5
)
```
