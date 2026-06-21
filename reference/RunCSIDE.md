# Run C-SIDE spatial differential expression

Run `spacexr` C-SIDE after RCTD to test cell type-specific spatial or
condition-aware differential expression.

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

  Named numeric vector used by
  [`spacexr::run.CSIDE.single()`](https://rdrr.io/pkg/spacexr/man/run.CSIDE.single.html).
  Names must match spatial spot names.

- group.by:

  Metadata column used to build C-SIDE regions when `region_list` is not
  supplied.

- condition.by:

  Binary metadata column converted to a 0/1 explanatory variable for
  [`spacexr::run.CSIDE.single()`](https://rdrr.io/pkg/spacexr/man/run.CSIDE.single.html).

- design:

  Numeric design matrix used by
  [`spacexr::run.CSIDE()`](https://rdrr.io/pkg/spacexr/man/run.CSIDE.html).
  Row names must match `barcodes` or spatial spot names.

- region_list:

  Named list of barcode vectors used by
  [`spacexr::run.CSIDE.regions()`](https://rdrr.io/pkg/spacexr/man/run.CSIDE.regions.html).

- barcodes:

  Barcodes used by
  [`spacexr::run.CSIDE()`](https://rdrr.io/pkg/spacexr/man/run.CSIDE.html)
  or `spacexr::run.CSIDE.intercept()`. If `NULL`, row names of `design`
  or all spatial spots are used where applicable.

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

  Whether to store detailed RCTD results in `srt@tools`.

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
library(scop)

data(visium_human_pancreas_sub)
data(panc8_sub)

spatial <- RunRCTD(
  srt = visium_human_pancreas_sub,
  reference = panc8_sub,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  rctd_mode = "full",
  max_cores = 1
)

spatial$region <- ifelse(spatial$col > stats::median(spatial$col), "right", "left")
spatial <- RunCSIDE(
  spatial,
  group.by = "region",
  celltypes = c("alpha", "beta"),
  gene_threshold = 0.00005,
  cell_type_threshold = 125
)

head(spatial@tools$CSIDE$result_table)
SpatialSpotPlot(spatial, group.by = "CSIDE_n_sig")

cside_genes <- head(spatial@tools$CSIDE$result_table$feature, 2)
SpatialSpotPlot(spatial, features = cside_genes, assay = "Spatial")
} # }
```
