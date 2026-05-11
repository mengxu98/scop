# Run DoRothEA transcription factor activity inference

Run DoRothEA transcription factor activity inference

## Usage

``` r
RunDorothea(
  srt,
  assay = NULL,
  layer = "data",
  species = c("Homo_sapiens", "Mus_musculus"),
  confidence = c("A", "B", "C"),
  regulons = NULL,
  method = c("ulm", "viper", "wmean"),
  minsize = 5,
  options = list(),
  assay_name = "dorothea",
  new_assay = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the expression matrix.

- species:

  Species used to select bundled DoRothEA regulons.

- confidence:

  DoRothEA confidence levels to keep.

- regulons:

  Optional regulon table with `tf`, `target`, `mor`, and `confidence`
  columns. If `NULL`, bundled `dorothea_hs` or `dorothea_mm` data are
  loaded from the `dorothea` package.

- method:

  Activity inference backend from `decoupleR`.

- minsize:

  Minimum regulon size passed to `decoupleR`.

- options:

  Additional named options passed to the selected `decoupleR` function.

- assay_name:

  Name of the assay used to store TF activity scores.

- new_assay:

  Whether to store TF activity scores as a new assay.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with DoRothEA results stored in
`srt@tools[["Dorothea"]]`.
