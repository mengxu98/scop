# Run doublet-calling for single cell RNA-seq data.

Run doublet-calling for single cell RNA-seq data.

## Usage

``` r
RunDoubletCalling(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  db_method = "scDblFinder",
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

- db_method:

  Method used for doublet-calling. Can be one of `"scDblFinder"`,
  `"Scrublet"`, `"DoubletDetection"`, `"scds_cxds"`, `"scds_bcds"`,
  `"scds_hybrid"`.

- data_type:

  Optional precomputed result from
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  for the input assay. Primarily used internally to avoid repeated scans
  of the same count matrix across nested QC calls.

- ...:

  Additional arguments to be passed to the corresponding doublet-calling
  method.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Returns a Seurat object with the doublet prediction results and
prediction scores stored in the meta.data.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scDblFinder_class"
)

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scDblFinder_score"
)
```
