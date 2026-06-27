# Run SPOTlight spatial deconvolution

Estimate spot-level cell type proportions from a spatial `Seurat` object
using a single-cell `Seurat` reference and the optional `SPOTlight`
package.

## Usage

``` r
RunSPOTlight(
  srt,
  reference,
  reference_label = "celltype",
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  mgs = NULL,
  marker_top_n = 100,
  marker_min_logfc = 0,
  gene_id = "gene",
  group_id = "cluster",
  weight_id = "weight",
  min_prop = 0.01,
  scale = TRUE,
  prefix = "SPOTlight",
  tool_name = "SPOTlight",
  store_results = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

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

- mgs:

  Optional marker-gene table passed to `SPOTlight`. It must contain
  columns named by `gene_id`, `group_id`, and `weight_id`. If `NULL`, a
  simple group-vs-rest marker table is generated from the reference
  expression matrix.

- marker_top_n:

  Number of automatically generated marker genes retained per reference
  cell type.

- marker_min_logfc:

  Minimum group-vs-rest log2 fold-change used when generating markers
  automatically.

- gene_id, group_id, weight_id:

  Column names in `mgs` for gene IDs, cell type labels, and marker
  weights.

- min_prop:

  Minimum cell-type proportion passed to `SPOTlight`.

- scale:

  Whether `SPOTlight` scales expression internally.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to
  [`SPOTlight::SPOTlight()`](https://rdrr.io/pkg/SPOTlight/man/SPOTlight.html).

## Value

A `Seurat` object with `SPOTlight` proportion columns in metadata and
dominant cell type summaries. When `store_results = TRUE`, detailed
results are stored in `srt@tools[[tool_name]]`.

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
spotlight_weights <- data.frame(
  SPOTlight_prop_Ductal = seq(0.70, 0.20, length.out = ncol(spatial)),
  SPOTlight_prop_Endocrine = seq(0.20, 0.70, length.out = ncol(spatial)),
  SPOTlight_prop_Immune = 0.10,
  row.names = colnames(spatial)
)
spotlight_weights <- spotlight_weights / rowSums(spotlight_weights)
spatial <- Seurat::AddMetaData(spatial, spotlight_weights)
spatial$SPOTlight_dominant_type <- sub(
  "^SPOTlight_prop_",
  "",
  colnames(spotlight_weights)[max.col(spotlight_weights)]
)
spatial$SPOTlight_max_prop <- apply(spotlight_weights, 1, max)

SpatialSpotPlot(
  spatial,
  group.by = "SPOTlight_dominant_type",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)

if (requireNamespace("scatterpie", quietly = TRUE)) {
  SpatialSpotPlot(
    spatial,
    group.by = "SPOTlight_dominant_type",
    plot_type = "pie",
    overlay_image = FALSE,
    coord.cols = c("x", "y")
  )
}


if (
  requireNamespace("SPOTlight", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
data(pancreas_sub)
features_use <- head(intersect(rownames(spatial), rownames(pancreas_sub)), 300)
spatial <- RunSPOTlight(
  srt = spatial,
  reference = pancreas_sub,
  reference_label = "CellType",
  assay = "Spatial",
  reference_assay = "RNA",
  features = features_use,
  marker_top_n = 20,
  verbose = FALSE
)

spotlight_cols <- grep(
  "^SPOTlight_prop_",
  colnames(spatial@meta.data),
  value = TRUE
)
SpatialSpotPlot(
  spatial,
  group.by = spotlight_cols[1:min(3, length(spotlight_cols))],
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
}
```
