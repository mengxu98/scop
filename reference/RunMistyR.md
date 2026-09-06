# Run mistyR multiview spatial modeling

Build a small `mistyR` view composition from a spatial `Seurat` object,
train MISTy models, collect results, and store a standardized result
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
  coordinate_space = c("raw", "legacy_display"),
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

  Name of the Seurat spatial image. Required when multiple images are
  present; a single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used when no Seurat image coordinates are
  available.

- views:

  Spatial views to add besides the required intraview. One or both of
  `"para"` and `"juxta"`.

- para_l, para_zoi, para_family, para_approx, para_nn:

  Parameters passed to `mistyR::add_paraview()`. `para_l` and `para_zoi`
  use the selected coordinate units; `para_nn` is a unitless neighbor
  count.

- juxta_neighbor_thr:

  Neighbor threshold passed to `mistyR::add_juxtaview()`, expressed in
  the selected coordinate units.

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

- coordinate_space:

  Coordinate system used to build MISTy views. The default is raw
  acquisition coordinates; `"legacy_display"` remains an explicit
  compatibility option.

- ...:

  Additional named arguments passed to `mistyR::run_misty()`.

## Value

A `Seurat` object with results stored in `srt@tools[[tool_name]]` when
`store_results = TRUE`.

## Examples

``` r
data(visium_human_pancreas_results_sub)
misty <- visium_human_pancreas_results_sub@tools$MistyR
misty$summary
#> $views
#> [1] "intraview"      "misty.uniqueid" "paraview.5"    
#> 
#> $n_targets
#> [1] 5
#> 
#> $n_improvement_records
#> [1] 40
#> 
#> $n_contribution_records
#> [1] 30
#> 
#> $importance_views
#> [1] "view"       "Predictor"  "Target"     "Importance" "nsamples"  
#> 
MistyRPlot(res = misty, type = "improvements", top_n = 10)


if (FALSE) { # \dontrun{
check_r("saezlab/mistyR", verbose = FALSE)
spatial <- RunMistyR(
  visium_human_pancreas_results_sub,
  assay = "Spatial", features = rownames(visium_human_pancreas_results_sub)[1:5],
  coord.cols = c("x", "y"), views = "para", para_l = 5,
  cv_folds = 3, verbose = FALSE
)
} # }
```
