# Run RCTD spatial deconvolution

Estimate spot-level cell-type proportions using RCTD and an annotated
single-cell reference.

## Usage

``` r
RunRCTD(
  object,
  reference,
  reference_label = "celltype",
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  rctd_mode = c("full", "multi", "doublet"),
  max_cores = 1,
  min_cells = 25,
  prefix = "RCTD",
  store_results = TRUE,
  round_counts = TRUE,
  create_rctd_params = list(),
  run_rctd_params = list(),
  verbose = TRUE,
  ...,
  coordinate_space = c("raw", "legacy_display"),
  tool_name = "RCTD",
  srt = NULL
)
```

## Arguments

- object:

  Spatial `Seurat` object used as the RCTD query.

- reference:

  Reference `Seurat` object containing annotated single cells.

- reference_label:

  Metadata column in `reference` with cell type labels.

- assay:

  Assay used in `srt`. If `NULL`, the default assay is used.

- reference_assay:

  Assay used in `reference`.

- layer, reference_layer:

  Assay layers used for spatial and reference raw counts.

- features:

  Features used for RCTD. If `NULL`, shared features are used.

- image:

  Name of the Seurat spatial image used to recover coordinates when
  `coord.cols` are not available.

- coord.cols:

  Metadata coordinate columns used when no image coordinate source is
  requested or available.

- rctd_mode:

  RCTD mode passed to `spacexr`. `"full"` is the default for Visium spot
  deconvolution.

- max_cores:

  Number of cores passed to `spacexr`.

- min_cells:

  Minimum number of reference cells required for each cell type. Old
  `spacexr` RCTD requires at least 25 cells per type.

- prefix:

  Prefix for metadata columns.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- round_counts:

  Whether to round non-integer values before passing data to `spacexr`.
  RCTD requires integer count matrices; rounding is a compatibility
  option and does not recover original counts.

- create_rctd_params:

  Additional parameters passed to `spacexr::createRctd()` or
  [`spacexr::create.RCTD()`](https://rdrr.io/pkg/spacexr/man/create.RCTD.html).

- run_rctd_params:

  Additional parameters passed to `spacexr::runRctd()` or
  [`spacexr::run.RCTD()`](https://rdrr.io/pkg/spacexr/man/run.RCTD.html).

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to the RCTD run step.

- coordinate_space:

  Coordinate space used for distance-sensitive input. The default is raw
  acquisition coordinates, so backend distances use raw coordinate
  units. Use `"legacy_display"` explicitly to reproduce the
  display-scaled coordinates used before scop 0.9.0.

- tool_name:

  Name used to store the plain result bundle in `srt@tools`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with RCTD proportion columns in metadata and dominant
cell type summaries. When `store_results = TRUE`, detailed results are
also stored in `srt@tools[[tool_name]]`.

## Examples

``` r
data(visium_human_pancreas_sub)
data(panc8_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
# Results are conditional on the three reference cell types used here.
reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
  c("ductal", "alpha", "beta")]
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
features_use <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(spatial)
), 300)
spatial <- RunRCTD(
  object = spatial,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = features_use,
  coord.cols = c("x", "y"),
  rctd_mode = "full",
  max_cores = 1,
  min_cells = 25,
  prefix = "RCTD",
  verbose = FALSE
)
#> Begin: process_cell_type_info
#> process_cell_type_info: number of cells in reference: 694
#> process_cell_type_info: number of genes in reference: 97
#> 
#> ductal   beta  alpha 
#>    211    304    179 
#> End: process_cell_type_info
#> create.RCTD: getting regression differentially expressed genes: 
#> get_de_genes: ductal found DE genes: 24
#> get_de_genes: beta found DE genes: 7
#> get_de_genes: alpha found DE genes: 51
#> get_de_genes: total DE genes: 82
#> create.RCTD: getting platform effect normalization differentially expressed genes: 
#> get_de_genes: ductal found DE genes: 27
#> get_de_genes: beta found DE genes: 8
#> get_de_genes: alpha found DE genes: 55
#> get_de_genes: total DE genes: 90
#> fitBulk: decomposing bulk
#> chooseSigma: using initial Q_mat with sigma =  1
#> Likelihood value: 1065.34481069009
#> Sigma value:  0.84
#> Likelihood value: 1043.70504794368
#> Sigma value:  0.69
#> Likelihood value: 1027.37241216243
#> Sigma value:  0.61
#> Likelihood value: 1021.04732656451
#> Sigma value:  0.53
#> Likelihood value: 1017.22091138995
#> Sigma value:  0.49
#> Likelihood value: 1016.52797026522
#> Sigma value:  0.48
#> Likelihood value: 1016.50517429608
#> Sigma value:  0.48
SpatialDeconvolutionPlot(
  spatial,
  tool_name = "RCTD",
  plot_type = "dominant",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
```
