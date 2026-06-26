# Run scMalignantFinder cancer cell state scoring

Score cancer cell state gene sets with `scMalignantFinder` AUCell
utilities and append the resulting activity scores to Seurat metadata.

## Usage

``` r
RunscMalignantStates(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "data",
  cells = NULL,
  gene_sets,
  norm_type = FALSE,
  prefix = "scMalignantState_",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- adata:

  Optional Python AnnData object.

- h5ad:

  Optional path to an `.h5ad` file.

- assay:

  Assay used when `srt` is supplied. Default is `"RNA"`.

- layer:

  Layer used when `srt` is supplied. Default is `"data"`.

- cells:

  Optional cells to run. If supplied with `srt`, results are appended to
  these cells and other cells receive `NA`.

- gene_sets:

  Path to a `.gmt` file containing cancer cell state gene sets, such as
  `Malignant_MPs.Gavish_2023.gmt` from the scMalignantFinder resources.

- norm_type:

  Passed to `scMalignantFinder`. Use `TRUE` for raw counts that should
  be library-size normalized; use `FALSE` for already normalized input.
  Default is `FALSE`.

- prefix:

  Optional prefix for output metadata columns. Default preserves the
  original `scMalignantFinder` column names.

- return_seurat:

  Whether to return a Seurat object when `srt` is supplied. If `FALSE`,
  returns a data frame of predictions.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with cancer-state AUCell scores, or a data frame when
`return_seurat = FALSE`.

## Examples

``` r
if (FALSE) {
srt <- RunscMalignantStates(
  srt,
  gene_sets = "path/to/Malignant_MPs.Gavish_2023.gmt"
)
}
```
