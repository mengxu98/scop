# Rename features for the Seurat object

Rename features for the Seurat object

## Usage

``` r
RenameFeatures(srt, newnames = NULL, assays = NULL)
```

## Arguments

- srt:

  A Seurat object.

- newnames:

  A vector with the same length of features in Seurat object, or
  characters named with old features.

- assays:

  Assays to rename.

## Examples

``` r
data(panc8_sub)
head(rownames(panc8_sub))
#> [1] "A1CF"   "A4GALT" "AAAS"   "AACS"   "AADAC"  "AADAT" 
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
#> ℹ [2025-11-19 14:33:43] Rename features for the assay: RNA
head(rownames(panc8_rename))
#> [1] "A1cf"   "A4galt" "Aaas"   "Aacs"   "Aadac"  "Aadat" 
```
