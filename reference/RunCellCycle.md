# Run cell cycle scoring

Estimate cell cycle state with Seurat gene-set scoring,
[`scran::cyclone()`](https://rdrr.io/pkg/scran/man/cyclone.html), or
[tricycle](https://bioconductor.org/packages/tricycle).

## Usage

``` r
RunCellCycle(
  object,
  method = c("Seurat", "cyclone", "tricycle"),
  assay = NULL,
  layer = "counts",
  species = "Homo_sapiens",
  name = "CellCycle",
  phase_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...,
  srt = NULL
)
```

## Arguments

- object:

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

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with cell cycle metadata and, for `tricycle`, a
tricycle embedding reduction.

## Examples

``` r
data(pancreas_sub)
srt <- pancreas_sub[, 1:80]

if (FALSE) { # \dontrun{
srt <- RunCellCycle(
  srt,
  method = "cyclone",
  species = "Mus_musculus",
  name = "Cyclone"
)

srt <- RunCellCycle(
  srt,
  method = "tricycle",
  species = "Mus_musculus",
  name = "Tricycle"
)
CellDimPlot(
  srt,
  reduction = "Tricycle_tricycleEmbedding",
  group.by = "Cyclone_cyclone_Phase"
)
FeatureDimPlot(
  srt,
  reduction = "Tricycle_tricycleEmbedding",
  features = "Tricycle_tricyclePosition"
)
} # }
```
