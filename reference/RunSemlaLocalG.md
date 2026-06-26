# Run semla local G spatial autocorrelation

Use
[`semla::RunLocalG()`](https://spatial-research.github.io/semla/reference/local-G.html)
on a Staffli-enabled Seurat object. Results are written by semla to
metadata or to an assay, depending on `store_in_metadata`.

## Usage

``` r
RunSemlaLocalG(
  srt,
  features,
  alternative = NULL,
  store_in_metadata = TRUE,
  assay_name = "GiScores",
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object with spatial image data.

- features:

  Features passed to
  [`semla::RunLocalG()`](https://spatial-research.github.io/semla/reference/local-G.html).

- alternative:

  Alternative hypothesis passed to semla. Use `NULL` to keep semla's
  default behavior.

- store_in_metadata:

  Whether semla should store results in metadata.

- assay_name:

  Assay name used by semla when `store_in_metadata = FALSE`.

- image_type:

  Image scale used by
  [`semla::UpdateSeuratForSemla()`](https://spatial-research.github.io/semla/reference/UpdateSeuratForSemla.html)
  when the object does not already contain a Staffli object.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to semla.

## Value

A `Seurat` object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
features <- rownames(spatial)[1:3]

spatial <- RunSemlaLocalG(
  spatial,
  features = features,
  store_in_metadata = TRUE
)

local_g_cols <- grep(features[1], colnames(spatial[[]]), value = TRUE)
head(spatial[[]][local_g_cols])
SpatialSpotPlot(spatial, group.by = local_g_cols[1])
} # }
```
