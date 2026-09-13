# Run multi-sample spatial integration

Integrate spatial transcriptomics samples and identify shared domains
using PRECAST.

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

  Spatial integration backend. The current stable backend is
  `"PRECAST"`.

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

## Details

Provide spatial samples with shared genes and raw coordinates. Use
`sample.by` to identify samples and `image` to select their spatial
images. PRECAST integrates expression and domain representations; it is
not a physical image-registration method.
