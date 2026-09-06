# Run CytoSPACE spatial assignment

Run CytoSPACE spatial assignment

## Usage

``` r
RunCytoSPACE(
  srt,
  reference,
  reference_label,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  cell_fractions = NULL,
  n_cells_per_spot = NULL,
  mean_cell_numbers = 5,
  scRNA_max_transcripts_per_cell = 1500,
  sampling_method = "duplicates",
  seed = 1,
  prefix = "CytoSPACE",
  store_results = TRUE,
  verbose = TRUE,
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("raw", "legacy_display"),
  backend = c("cpp", "r"),
  max_dense_gib = 8
)
```

## Arguments

- srt:

  A `Seurat` object.

- reference:

  Reference `Seurat` object containing annotated single cells.

- reference_label:

  Metadata column in `reference` with cell type labels.

- assay:

  Assay to use. `NULL` uses the default assay.

- reference_assay:

  Assay used in `reference`.

- layer, reference_layer:

  Assay layers used for spatial and reference expression.

- features:

  Features used for assignment. If `NULL`, shared features are used.

- cell_fractions:

  Optional cell-type fractions. Provide a named numeric vector, one-row
  matrix/data.frame, or a spot-by-cell-type matrix/data.frame.
  Spot-level rows are aggregated to the global composition used by the
  default CytoSPACE assignment workflow.

- n_cells_per_spot:

  Optional number of cells assigned to each spatial spot. If `NULL`,
  counts are estimated from spatial RNA reads with `mean_cell_numbers`.

- mean_cell_numbers:

  Mean number of cells per spot. Default `5`, matching the CytoSPACE
  Visium default.

- scRNA_max_transcripts_per_cell:

  Maximum reference transcripts per cell before assignment. Default
  `1500`, matching CytoSPACE.

- sampling_method:

  Sampling method. Only `"duplicates"` is supported in the package
  runtime.

- seed:

  Random seed used for deterministic reference downsampling and
  duplicate sampling.

- prefix:

  Prefix for metadata columns.

- store_results:

  Whether to store detailed assignment results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- image:

  Optional Seurat image used for spatial coordinates.

- coord.cols:

  Metadata coordinate columns used when no image is selected.

- coordinate_space:

  Coordinate space used for assignment locations. The default is raw
  acquisition coordinates. Use `"legacy_display"` explicitly to
  reproduce the display-scaled locations used before scop 0.9.0.

- backend:

  Numerical backend used to estimate cell-type fractions when
  `cell_fractions` is not supplied. `"cpp"` fuses normalization,
  reference centroid construction, correlation, and weighted
  aggregation; `"r"` keeps the reference implementation. Spot assignment
  uses C++ in both cases.

- max_dense_gib:

  Maximum estimated GiB allowed for dense expression working matrices.

## Value

A `Seurat` object with CytoSPACE metadata columns and detailed results
stored in `srt@tools[["CytoSPACE"]]`.

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
reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
  c("ductal", "alpha", "beta")]
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
features_use <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(spatial)
), 300)
spatial <- RunCytoSPACE(
  srt = spatial,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = features_use,
  mean_cell_numbers = 1,
  coord.cols = c("x", "y"),
  verbose = FALSE
)
SpatialSpotPlot(
  spatial,
  group.by = "CytoSPACE_dominant_type",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
```
