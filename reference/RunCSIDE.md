# Run C-SIDE spatial differential expression

Run `spacexr` C-SIDE after RCTD to test cell type-specific spatial or
condition-aware differential expression. C-SIDE stores effect and
significance results, not spot-level cell-type proportions, so it has no
recommended proportion plot. Inspect its stored tables in
`srt@tools[[tool_name]]`. The metadata columns `<prefix>_n_sig` and
`<prefix>_mode` repeat run-level summaries across spots. The former
counts significant result records, not unique genes or spot-level
effects. The example is a non-executing template because C-SIDE requires
a completed optional RCTD result.

## Usage

``` r
RunCSIDE(
  srt,
  rctd_result = NULL,
  explanatory.variable = NULL,
  group.by = NULL,
  condition.by = NULL,
  design = NULL,
  region_list = NULL,
  barcodes = NULL,
  mode = c("auto", "single", "regions", "general", "intercept"),
  assay = NULL,
  layer = "counts",
  celltypes = NULL,
  features = NULL,
  prefix = "CSIDE",
  tool_name = "CSIDE",
  store_results = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  Spatial `Seurat` object used as the RCTD query.

- rctd_result:

  Optional old-api `spacexr` RCTD object. If `NULL`,
  `srt@tools[["RCTD"]]$object` from
  [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md) is
  used.

- explanatory.variable:

  Named numeric vector used by `spacexr::run.CSIDE.single()`. Names must
  match spatial spot names.

- group.by:

  Metadata column used to build C-SIDE regions when `region_list` is not
  supplied.

- condition.by:

  Binary metadata column converted to a 0/1 explanatory variable for
  `spacexr::run.CSIDE.single()`.

- design:

  Numeric design matrix used by `spacexr::run.CSIDE()`. Row names must
  match `barcodes` or spatial spot names.

- region_list:

  Named list of barcode vectors used by `spacexr::run.CSIDE.regions()`.

- barcodes:

  Barcodes used by `spacexr::run.CSIDE()` or
  `spacexr::run.CSIDE.intercept()`. If `NULL`, row names of `design` or
  all spatial spots are used where applicable.

- mode:

  C-SIDE mode. `"auto"` dispatches to `"general"`, `"regions"`,
  `"single"`, or `"intercept"` based on the supplied design inputs.

- assay:

  Assay used in `srt`. If `NULL`, the default assay is used.

- layer:

  Assay layer used as the spatial expression source.

- celltypes:

  Optional cell types passed to C-SIDE as `cell_types`.

- features:

  Optional features retained in the normalized result table. C-SIDE
  itself still applies its own gene filtering through backend parameters
  such as `gene_threshold`.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used to store detailed C-SIDE results in `srt@tools`.

- store_results:

  Whether to store detailed results in `srt@tools`. Compatibility
  summary metadata is updated regardless; existing detailed results are
  retained when `FALSE` and do not represent the current run.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional named parameters passed to the selected C-SIDE backend,
  such as `cell_type_threshold`, `gene_threshold`, `doublet_mode`,
  `cell_type_specific`, or `params_to_test`. When using the stored
  result from
  [`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md)
  with `rctd_mode = "full"`, `doublet_mode` defaults to `FALSE` unless
  explicitly supplied.

## Value

A `Seurat` object with C-SIDE summary metadata and detailed results
stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
thisutils::check_r("dmcable/spacexr", verbose = FALSE)
data(visium_human_pancreas_sub)
data(panc8_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
spatial$region <- cut(
  spatial$x,
  breaks = stats::quantile(spatial$x, probs = seq(0, 1, length.out = 4)),
  include.lowest = TRUE,
  labels = c("left", "middle", "right")
)
reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
  c("ductal", "alpha", "beta")]
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
features_use <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(spatial)
), 300)
spatial <- RunRCTD(
  srt = spatial,
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
  verbose = FALSE
)
spatial <- spatial[, rownames(spatial@tools$RCTD$weights)]
spatial <- RunCSIDE(
  srt = spatial,
  group.by = "region",
  celltypes = c("ductal", "alpha"),
  features = features_use,
  gene_threshold = 0.00005,
  cell_type_threshold = 5,
  verbose = FALSE
)
SpatialSpotPlot(
  spatial,
  group.by = "CSIDE_n_sig",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
