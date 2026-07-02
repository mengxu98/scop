# Run mistyR multiview spatial modeling

Build a small `mistyR` view composition from a spatial `Seurat` object,
train MISTy models, collect results, and store a SCOP-style result
bundle in `srt@tools`. The intraview is always created from the selected
assay layer; optional juxtaview and paraview components describe local
and broader spatial context. `mistyR` is an optional Bioconductor
dependency installable with `BiocManager::install("mistyR")`.

## Usage

``` r
RunMistyR(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  views = "para",
  para_l = 10,
  para_zoi = 0,
  para_family = c("gaussian", "exponential", "linear", "constant"),
  para_approx = 1,
  para_nn = NULL,
  juxta_neighbor_thr = 15,
  view_cached = FALSE,
  results_folder = NULL,
  seed = 42,
  target_subset = NULL,
  bypass_intra = FALSE,
  cv_folds = 10,
  model_cached = FALSE,
  append = FALSE,
  tool_name = "MistyR",
  store_results = TRUE,
  store_views = FALSE,
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

- features:

  Features used by MISTy. If `NULL`, variable features are used when
  available; otherwise all assay features are used.

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- coord.cols:

  Metadata coordinate columns used when no Seurat image coordinates are
  available.

- views:

  Spatial views to add besides the required intraview. One or both of
  `"para"` and `"juxta"`.

- para_l, para_zoi, para_family, para_approx, para_nn:

  Parameters passed to `mistyR::add_paraview()`.

- juxta_neighbor_thr:

  Neighbor threshold passed to `mistyR::add_juxtaview()`.

- view_cached:

  Whether generated mistyR views should use cache.

- results_folder:

  Folder passed to `mistyR::run_misty()`. If `NULL`, a temporary folder
  is used.

- seed, target_subset, bypass_intra, cv_folds, model_cached, append:

  Parameters passed to `mistyR::run_misty()`.

- tool_name:

  Name used to store results in `srt@tools`.

- store_results:

  Whether to store results in `srt@tools`.

- store_views:

  Whether to store the mistyR view composition in `srt@tools`. This can
  be large.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional named arguments passed to `mistyR::run_misty()`.

## Value

A `Seurat` object with results stored in `srt@tools[[tool_name]]` when
`store_results = TRUE`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)

if (
  isTRUE(check_r("mistyR", verbose = FALSE))
) {
  spatial <- RunMistyR(
    spatial,
    assay = "Spatial",
    layer = "data",
    features = rownames(spatial)[1:10],
    coord.cols = c("x", "y"),
    views = "para",
    para_l = 5,
    cv_folds = 3,
    verbose = FALSE
  )
  spatial@tools$MistyR$summary
}
#> Error in check_r("mistyR", verbose = FALSE): could not find function "check_r"
```
