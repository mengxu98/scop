# Run CytoSPACE spatial assignment

Assign reference single cells to spatial transcriptomics spots using a
native R/C++ implementation of the default CytoSPACE spot-level
workflow.

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
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  Reference `Seurat` object containing annotated single cells.

- reference_label:

  Metadata column in `reference` with cell type labels.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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

## Value

A `Seurat` object with CytoSPACE metadata columns and detailed results
stored in `srt@tools[["CytoSPACE"]]`.

## Examples

``` r
data(visium_human_pancreas_sub)
data(pancreas_sub)
pancreas_human <- ConvertHomologs(
  pancreas_sub,
  species_from = "Mus_musculus",
  species_to = "Homo_sapiens",
  verbose = FALSE
)
spatial <- RunCytoSPACE(
  visium_human_pancreas_sub,
  reference = pancreas_human,
  reference_label = "CellType",
  mean_cell_numbers = 1,
  verbose = FALSE
)

SpatialDimPlot(
  visium_human_pancreas_sub,
  group.by = "coda_label",
  theme_use = "theme_scop"
)
#> Error in SpatialDimPlot(visium_human_pancreas_sub, group.by = "coda_label",     theme_use = "theme_scop"): unused argument (theme_use = "theme_scop")

SpatialDimPlot(
  spatial,
  group.by = "CytoSPACE_dominant_type",
  theme_use = "theme_scop"
)
#> Error in SpatialDimPlot(spatial, group.by = "CytoSPACE_dominant_type",     theme_use = "theme_scop"): unused argument (theme_use = "theme_scop")
```
