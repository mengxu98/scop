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
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + dendextend   1.19.1 
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got dendextend 1.19.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.05 MB)
#> ✔ Installed dendextend 1.19.1  (73ms)
#> ✔ 1 pkg + 19 deps: kept 19, added 1, dld 1 (5.05 MB) [1.9s]
dend
#> 'dendrogram' with 2 branches and 10 members total, at height 7.730525 
plot(dend)
```
