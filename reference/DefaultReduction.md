# Find the default reduction name in a Seurat object

Find the default reduction name in a Seurat object

## Usage

``` r
DefaultReduction(
  srt,
  pattern = NULL,
  min_dim = 2,
  max_distance = 0.1,
  verbose = TRUE
)
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

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Default reduction name.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
names(pancreas_sub@reductions)

DefaultReduction(pancreas_sub)

DefaultReduction(pancreas_sub, pattern = "pca")

DefaultReduction(pancreas_sub, pattern = "umap")
```
