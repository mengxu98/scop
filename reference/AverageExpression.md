# Average expression by cell group

Seurat already performs this operation as sparse matrix multiplication.
scop preserves its complete parameter and object contract.

## Usage

``` r
AverageExpression(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  layer = "data",
  slot = lifecycle::deprecated(),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Seurat object

- assays:

  Which assays to use. Default is all assays

- features:

  Features to analyze. Default is all features in the assay

- return.seurat:

  Whether to return the data as a Seurat object. Default is FALSE

- group.by:

  Category (or vector of categories) for grouping (e.g, ident,
  replicate, celltype); 'ident' by default To use multiple categories,
  specify a vector, such as c('ident', 'replicate', 'celltype')

- add.ident:

  (Deprecated). Place an additional label on each cell prior to
  pseudobulking

- layer:

  Layer(s) to use; if multiple layers are given, assumed to follow the
  order of 'assays' (if specified) or object's assays

- slot:

  (Deprecated). Slots(s) to use

- verbose:

  Print messages and show progress bar

- ...:

  Arguments to be passed to methods such as
  [`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)

## Value

A named list of averaged matrices or a Seurat object.
