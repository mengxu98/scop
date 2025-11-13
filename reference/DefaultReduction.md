# Find the default reduction name in a Seurat object

Find the default reduction name in a Seurat object

## Usage

``` r
DefaultReduction(srt, pattern = NULL, min_dim = 2, max_distance = 0.1)
```

## Arguments

- srt:

  A Seurat object.

- pattern:

  Character string containing a regular expression to search for.

- min_dim:

  Minimum dimension threshold.

- max_distance:

  Maximum distance allowed for a match.

## Value

Default reduction name.

## Examples

``` r
data(pancreas_sub)
names(pancreas_sub@reductions)
#> [1] "X_pca"  "X_umap"
DefaultReduction(pancreas_sub)
#> [1] "X_umap"

# Searches for matches to "pca"
DefaultReduction(pancreas_sub, pattern = "pca")
#> [1] "X_pca"

# Searches for approximate matches to "pc"
DefaultReduction(pancreas_sub, pattern = "pc")
#> [1] "X_pca"
```
