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
#> ◌ [2025-11-13 12:44:29] Installing: dendextend...
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + dendextend   1.19.1 
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Installed dendextend 1.19.1  (1.1s)
#> ✔ 1 pkg + 19 deps: kept 19, added 1 [2.3s]
#> ✔ [2025-11-13 12:44:31] dendextend installed successfully
dend
#> 'dendrogram' with 2 branches and 10 members total, at height 7.730525 
plot(dend)
```
