# Run doublet-calling with Scrublet

Run doublet-calling with Scrublet

## Usage

``` r
db_Scrublet(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- data_type:

  Optional precomputed result from
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  for the input assay. Primarily used internally to avoid repeated scans
  of the same count matrix across nested QC calls.

- ...:

  Additional arguments to be passed to
  [scrublet.Scrublet](https://github.com/swolock/scrublet).

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- db_Scrublet(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.Scrublet_class"
)

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.Scrublet_score"
)
} # }
```
