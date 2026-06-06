# Rename features for the Seurat object

Rename features for the Seurat object

## Usage

``` r
RenameFeatures(srt, newnames = NULL, assays = NULL, verbose = TRUE)
```

## Arguments

- srt:

  A Seurat object.

- newnames:

  A vector with the same length of features in Seurat object, or
  characters named with old features.

- assays:

  Assays to rename.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
head(rownames(panc8_sub))
# Simply convert genes from human to mouse and preprocess the data
genenames <- make.unique(
  thisutils::capitalize(rownames(panc8_sub),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_rename <- RenameFeatures(
  panc8_sub,
  newnames = genenames
)
head(rownames(panc8_rename))
```
