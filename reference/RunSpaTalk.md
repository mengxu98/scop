# Run official SpaTalk spatial communication analysis

Run the official SpaTalk R backend for single-cell spatial data or spot
data with a single-cell reference. Results are stored as a compact
method bundle and added to SCOP's unified CCC table.

## Usage

``` r
RunSpaTalk(
  srt,
  group.by,
  mode = c("auto", "single_cell", "spot"),
  reference = NULL,
  reference.group.by = NULL,
  species = c("Human", "Mouse"),
  assay = NULL,
  layer = "counts",
  reference.assay = NULL,
  reference.layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  deconvolution = c("NNLM", "RCTD", "none"),
  rctd.tool = "RCTD",
  result.name = "default",
  store.object = c("minimal", "full"),
  overwrite = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A \`Seurat\` spatial object.

- group.by:

  Metadata column containing cell labels for single-cell mode.

- mode:

  Analysis mode. \`"auto"\` selects spot mode when \`reference\` is
  supplied and single-cell mode otherwise.

- reference:

  Optional single-cell \`Seurat\` reference for spot mode.

- reference.group.by:

  Cell-type column in \`reference\`.

- species:

  SpaTalk species database.

- assay, layer:

  Spatial assay and raw-count layer.

- reference.assay, reference.layer:

  Reference assay and raw-count layer.

- image:

  Explicit spatial image. Required for multi-image objects.

- coord.cols:

  Metadata coordinate columns used when no image is present.

- deconvolution:

  Spot deconvolution mode. \`"NNLM"\` runs SpaTalk method 1; \`"RCTD"\`
  reuses \`srt@tools\[\[rctd.tool\]\]\$weights\` from \[RunRCTD()\] (or
  an explicitly supplied \`dec_result\`); \`"none"\` requires a named
  external \`dec_result\` matrix in \`...\`. Imported matrices use
  SpaTalk method 2.

- rctd.tool:

  Tool key containing stored RCTD weights when \`deconvolution =
  "RCTD"\`.

- result.name:

  Stored result name.

- store.object:

  Whether to store the full upstream SpaTalk object.

- overwrite:

  Whether to replace an existing named result.

- backend:

  Backend used only for SCOP CCC result aggregation.

- verbose:

  Whether to print progress messages.

- ...:

  Named arguments forwarded to compatible official SpaTalk steps.

## Value

The input \`Seurat\` object with \`srt@tools\$SpaTalk\` and unified CCC
results updated.

## References

\<https://github.com/ZJUFanLab/SpaTalk\>

## Examples

``` r
if (FALSE) { # \dontrun{
available <- unname(unlist(thisutils::check_r(
  c("linxihui/NNLM", "ZJUFanLab/SpaTalk"),
  verbose = FALSE
)))
if (length(available) == 2L && all(available)) {
  spatial <- RunSpaTalk(
    spatial,
    group.by = "celltype",
    mode = "single_cell",
    image = "slice1"
  )
}
} # }
```
