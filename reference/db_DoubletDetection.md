# Run doublet-calling with DoubletDetection

Run doublet-calling with DoubletDetection

## Usage

``` r
db_DoubletDetection(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  cores = 1,
  data_type = NULL,
  ...
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

- cores:

  The number of CPU cores to use for `doubletdetection`. Default is `1`.

- data_type:

  Optional precomputed result from
  [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  for the input assay. Primarily used internally to avoid repeated scans
  of the same count matrix across nested QC calls.

- ...:

  Additional arguments to be passed to
  [doubletdetection.BoostClassifier](https://github.com/JonathanShor/DoubletDetection).

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- db_DoubletDetection(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.DoubletDetection_class"
)

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.DoubletDetection_score"
)
} # }
```
