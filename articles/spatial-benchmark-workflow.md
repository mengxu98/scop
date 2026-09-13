# Spatial domain benchmark workflow

This article compares spatial domain clustering methods with
\[RunSpatialBenchmark()\]. Each method runs in an isolated R process;
results include ARI, NMI, purity, runtime, and peak memory against a
gold-standard label column. Visualize the result with
\[SpatialBenchmarkPlot()\].

Supported producers include `BayesSpace`, `BANKSY`, and `SmoothClust`.
This workflow is separate from scRNA batch integration benchmarking
(\[RunIntegrationBenchmark()\]).

## Load data and define a gold standard

This workflow uses a small simulated spatial object with three regions
and region-specific expression signals. It demonstrates the benchmark
interface; replace the simulated labels with tissue annotations for real
analyses.

``` r

library(scop)
#>           ⬢          .        ⬡             ⬢     .
#>                      _____ _________  ____
#>                     / ___// ___/ __ ./ __ .
#>                    (__  )/ /__/ /_/ / /_/ /
#>                   /____/ .___/.____/ .___/
#>                                   /_/
#>       ⬢               .      ⬡        .          ⬢
#> ------------------------------------------------------------
#> Version: 0.9.2 (2026-09-12 update)
#> Website: https://mengxu98.github.io/scop/
#> 
#> Python environment initialization is disabled
#> To enable it, set: options(scop_env_init = TRUE)
#> 
#> The message can be suppressed by: 
#>   suppressPackageStartupMessages(library(scop))
#>   or options(log_message.verbose = FALSE)
#> ------------------------------------------------------------

set.seed(42)
n_spots <- 36L
n_genes <- 60L
grid <- expand.grid(x = 1:6, y = 1:6)
spatial_domain <- factor(ifelse(grid$x <= 2, "left", ifelse(grid$x <= 4, "middle", "right")))
counts <- matrix(rpois(n_genes * n_spots, 2),
  nrow = n_genes,
  dimnames = list(paste0("Gene", seq_len(n_genes)), paste0("spot", seq_len(n_spots)))
)
signal <- split(seq_len(n_genes), rep(c("left", "middle", "right"), each = 20))
for (label in names(signal)) {
  counts[signal[[label]], spatial_domain == label] <-
    counts[signal[[label]], spatial_domain == label] + rpois(sum(spatial_domain == label) * 20, 5)
}
spatial <- SeuratObject::CreateSeuratObject(counts)
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
spatial$x <- grid$x
spatial$y <- grid$y
spatial$gold_domain <- spatial_domain
table(spatial$gold_domain)
#> 
#>   left middle  right 
#>     12     12     12
```

## Run the benchmark

Methods share the same spot set. Tune per-method arguments through
`method_params` when backends need different PCA counts or layers.

``` r

bench <- RunSpatialBenchmark(
  spatial,
  gold_standard = "gold_domain",
  method_params = list(
    BayesSpace = list(n.PCs = 5, n.HVGs = 20, coord.cols = c("x", "y")),
    BANKSY = list(layer = "counts", coord.cols = c("x", "y")),
    SmoothClust = list(layer = "counts", min_spots = 1, coord.cols = c("x", "y"))
  ),
  verbose = FALSE
)

bench$summary
#>        method      workflow ARI NMI purity runtime_s baseline_memory_mb
#> 1  BayesSpace SpatialDomain   1   1      1    11.338           617.2969
#> 2      BANKSY SpatialDomain  NA  NA     NA     3.566           617.1719
#> 3 SmoothClust SpatialDomain   1   1      1     3.268           617.4141
#>   peak_memory_mb memory_delta_mb n_evaluated n_clusters  status
#> 1      1311.0234        693.7266          36          3 success
#> 2       886.4414        269.2695           0         NA  failed
#> 3       881.0508        263.6367          36          3 success
#>                               error
#> 1                                  
#> 2 Not enough neighbors in data set!
#> 3
```

## Inspect the result object

[`RunSpatialBenchmark()`](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md)
returns a `spatial_benchmark_result` list:

- `$summary`: quality and resource metrics per method;
- `$predictions`: aligned domain labels per spot;
- `$metrics`: long metric table for custom plots.

``` r

str(bench$summary, max.level = 1)
#> 'data.frame':    3 obs. of  13 variables:
#>  $ method            : chr  "BayesSpace" "BANKSY" "SmoothClust"
#>  $ workflow          : chr  "SpatialDomain" "SpatialDomain" "SpatialDomain"
#>  $ ARI               : num  1 NA 1
#>  $ NMI               : num  1 NA 1
#>  $ purity            : num  1 NA 1
#>  $ runtime_s         : num  11.34 3.57 3.27
#>  $ baseline_memory_mb: num  617 617 617
#>  $ peak_memory_mb    : num  1311 886 881
#>  $ memory_delta_mb   : num  694 269 264
#>  $ n_evaluated       : int  36 0 36
#>  $ n_clusters        : int  3 NA 3
#>  $ status            : chr  "success" "failed" "success"
#>  $ error             : chr  "" "Not enough neighbors in data set!" ""
head(bench$predictions)
#>   spot_id gold_standard     method prediction
#> 1   spot1          left BayesSpace          1
#> 2   spot2          left BayesSpace          1
#> 3   spot3        middle BayesSpace          3
#> 4   spot4        middle BayesSpace          3
#> 5   spot5         right BayesSpace          2
#> 6   spot6         right BayesSpace          2
```

## Plot benchmark views

[`SpatialBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md)
defaults to a publication-oriented overview that combines clustering
agreement with runtime and peak-memory efficiency.

``` r

SpatialBenchmarkPlot(data = bench, plot_type = "overview")
```

![](spatial-benchmark-workflow_files/figure-html/plot-overview-1.png)

Other views:

``` r

SpatialBenchmarkPlot(data = bench, plot_type = "quality")
SpatialBenchmarkPlot(data = bench, plot_type = "efficiency")
SpatialBenchmarkPlot(data = bench, plot_type = "heatmap")
```

## What to report

A spatial benchmark report should include:

- the gold-standard column and number of domains;
- methods compared and any `method_params` overrides;
- `$summary` with ARI, NMI, purity, runtime, and peak memory;
- the overview or quality plot;
- failed or unavailable methods (status column) when backends are
  missing.

Agreement metrics are higher when better; runtime and peak memory are
lower when better.
