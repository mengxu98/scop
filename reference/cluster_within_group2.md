# Cluster within group

Cluster within group

## Usage

``` r
cluster_within_group2(mat, factor)
```

## Arguments

- mat:

  A matrix of data

- factor:

  A factor

## Value

A dendrogram with ordered leaves

## Examples

``` r
mat <- matrix(rnorm(100), 10, 10)
factor <- factor(rep(1:2, each = 5))
dend <- cluster_within_group2(mat, factor)
dend
#> 'dendrogram' with 2 branches and 10 members total, at height 7.730525 
plot(dend)
```
