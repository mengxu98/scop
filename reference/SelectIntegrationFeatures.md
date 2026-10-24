# Select integration features with exact vectorized ranking

Uses a vectorized rank calculation when every object already has unique
variable features in its default assay. Other branches delegate to
Seurat.

## Usage

``` r
SelectIntegrationFeatures(
  object.list,
  nfeatures = 2000,
  assay = NULL,
  verbose = TRUE,
  fvf.nfeatures = 2000,
  ...
)
```

## Arguments

- object.list:

  List of Seurat objects.

- nfeatures:

  Number of integration features to return.

- assay:

  Assay names, one per object.

- verbose:

  Show Seurat progress messages on delegated paths.

- fvf.nfeatures:

  Features requested if Seurat must calculate them.

- ...:

  Passed to Seurat on delegated paths.

## Value

A character vector of integration features.
