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

  A `Seurat` object.

- method:

  Cell cycle estimation method. One of `"Seurat"`, `"cyclone"`, or
  `"tricycle"`.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Data layer used by `cyclone` and `tricycle`.

- species:

  Latin names for animals, i.e., `"Homo_sapiens"`, `"Mus_musculus"`

- name:

  Prefix for metadata columns and tricycle reduction names.

- phase_col:

  Optional metadata column used to store the final phase call, for
  example `"Phase"`. Default is `NULL`, which avoids writing a
  compatibility phase column.

- overwrite:

  Whether to overwrite existing output columns.

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

srt <- RunCellCycle(
  srt,
  method = "cyclone",
  species = "Mus_musculus",
  name = "Cyclone"
)
#> ℹ [2026-09-06 21:53:08] Start cell cycle scoring
#> 'select()' returned 1:many mapping between keys and columns
#> ℹ [2026-09-06 21:53:08] Map input feature names to ENSEMBL IDs with org.Mm.eg.db for scran::cyclone
#> ✔ [2026-09-06 21:53:17] Cell cycle scoring completed

srt <- RunCellCycle(
  srt,
  method = "tricycle",
  species = "Mus_musculus",
  name = "Tricycle"
)
#> ℹ [2026-09-06 21:53:17] Start cell cycle scoring
#> No custom reference projection matrix provided. The ref learned from mouse Neuroshpere data will be used.
#> The number of projection genes found in the new data is 485.
#> ✔ [2026-09-06 21:53:17] Cell cycle scoring completed
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
```
