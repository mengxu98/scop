# Run cell cycle scoring

Estimate cell cycle state with Seurat gene-set scoring,
[`scran::cyclone()`](https://rdrr.io/pkg/scran/man/cyclone.html), or
[tricycle](https://bioconductor.org/packages/tricycle).

## Usage

``` r
RunCellCycle(
  srt,
  method = c("Seurat", "cyclone", "tricycle"),
  assay = NULL,
  layer = "counts",
  species = "Homo_sapiens",
  name = "CellCycle",
  phase_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- method:

  Cell cycle estimation method. One of `"Seurat"`, `"cyclone"`, or
  `"tricycle"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Data layer used by `cyclone` and `tricycle`. Default is `"counts"`.

- species:

  Latin names for animals, i.e., `"Homo_sapiens"`, `"Mus_musculus"`

- name:

  Prefix for metadata columns and tricycle reduction names. Default is
  `"CellCycle"`.

- phase_col:

  Optional metadata column used to store the final phase call, for
  example `"Phase"`. Default is `NULL`, which avoids writing a
  compatibility phase column.

- overwrite:

  Whether to overwrite existing output columns. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the selected method.

## Value

A `Seurat` object with cell cycle metadata and, for `tricycle`, a
tricycle embedding reduction.

## Examples

``` r
data(pancreas_sub)
srt <- pancreas_sub[, 1:80]

if (requireNamespace("scran", quietly = TRUE)) {
  srt <- RunCellCycle(
    srt,
    method = "cyclone",
    species = "Mus_musculus",
    name = "Cyclone"
  )
}
#> ℹ [2026-06-28 14:51:51] Start cell cycle scoring
#> 'select()' returned 1:many mapping between keys and columns
#> ℹ [2026-06-28 15:01:45] Map input feature names to ENSEMBL IDs with org.Mm.eg.db for scran::cyclone
#> ✔ [2026-06-28 15:01:55] Cell cycle scoring completed

if (requireNamespace("tricycle", quietly = TRUE)) {
  srt <- RunCellCycle(
    srt,
    method = "tricycle",
    species = "Mus_musculus",
    name = "Tricycle"
  )
  if ("Cyclone_cyclone_Phase" %in% colnames(srt@meta.data)) {
    CellDimPlot(
      srt,
      reduction = "Tricycle_tricycleEmbedding",
      group.by = "Cyclone_cyclone_Phase"
    )
  }
  FeatureDimPlot(
    srt,
    reduction = "Tricycle_tricycleEmbedding",
    features = "Tricycle_tricyclePosition"
  )
}
#> ℹ [2026-06-28 15:01:55] Start cell cycle scoring
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘scale.data’ is empty
#> Warning: 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning: 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
#> No custom reference projection matrix provided. The ref learned from mouse Neuroshpere data will be used.
#> The number of projection genes found in the new data is 485.
#> ✔ [2026-06-28 15:02:22] Cell cycle scoring completed
```
