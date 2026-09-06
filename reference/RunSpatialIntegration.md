# Run multi-sample spatial integration

Integrate multi-slice or multi-sample spatial transcriptomics data with
an optional spatial backend and store standardized embeddings, domains,
and aligned coordinates in a `Seurat` object.

## Usage

``` r
RunSpatialIntegration(
  object,
  method = "PRECAST",
  sample.by = NULL,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  features = NULL,
  image = NULL,
  reduction.name = NULL,
  cluster_colname = NULL,
  tool_name = "SpatialIntegration",
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
)
```

## Arguments

- object:

  A merged spatial `Seurat` object or a list of spatial `Seurat`
  objects.

- method:

  Spatial integration backend.

- sample.by:

  Metadata column identifying samples for a merged `Seurat` object. For
  list input, list names are copied into this column.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer used for expression values.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- features:

  Features to score. If `NULL`, current variable features are used; if
  no variable features are present, all assay features are used.

- image:

  Optional image name or named character vector mapping sample names to
  image names. By default, resolve the unique image covering each sample
  independently. Every input sample must remain represented. For list
  input, image names refer to each original list element, before Seurat
  renames duplicate image keys during merging.

- reduction.name:

  Name of the integrated embedding reduction. If `NULL`, a
  method-specific name is used.

- cluster_colname:

  Metadata column used for spatial domain labels. If `NULL`, a
  method-specific name is used.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store the full result in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- coordinate_space:

  Coordinate system used for integration distances and
  aligned-coordinate input. The default is raw acquisition coordinates;
  `"legacy_display"` remains an explicit compatibility option.

- ...:

  Additional backend-specific arguments.

## Value

A `Seurat` object with spatial integration results stored in metadata,
reductions, and `srt@tools[[tool_name]]`.

## Examples

``` r
data(visium_human_pancreas_pair_sub)
SpatialIntegrationPlot(
  visium_human_pancreas_pair_sub,
  plot_type = "spatial",
  group.by = "domain",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)

SpatialIntegrationPlot(
  visium_human_pancreas_pair_sub,
  plot_type = "embedding"
)


if (FALSE) { # \dontrun{
check_r("PRECAST", verbose = FALSE)
spatial <- RunSpatialIntegration(
  visium_human_pancreas_pair_sub,
  method = "PRECAST", sample.by = "sample", assay = "Spatial",
  coord.cols = c("x", "y"),
  features = rownames(visium_human_pancreas_pair_sub)[1:100],
  verbose = FALSE
)
} # }
```
