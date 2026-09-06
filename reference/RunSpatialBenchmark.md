# Benchmark spatial domain clustering methods

Run supported spatial domain clustering methods from the same immutable
input, compare their aligned labels with a gold standard, and record
classification quality, elapsed time, and sampled peak process-tree
memory. Each method runs in an isolated R process so one failed backend
does not corrupt the input or prevent the remaining methods from being
assessed.

## Usage

``` r
RunSpatialBenchmark(
  srt,
  gold_standard,
  methods = NULL,
  method_params = list(),
  n_clusters = NULL,
  metrics = c("ARI", "NMI"),
  keep_objects = FALSE,
  install_missing = FALSE,
  seed = 1,
  timeout = Inf,
  poll_interval = 0.1,
  verbose = TRUE
)
```

## Arguments

- srt:

  A spatial `Seurat` object.

- gold_standard:

  Either one metadata column in `srt` or a named vector whose names
  match the spot names in `srt`.

- methods:

  Spatial domain methods to benchmark. `NULL` uses every benchmarked
  producer (`BayesSpace`, `BANKSY`, and `SmoothClust`). Method names may
  be written with or without the `Run` prefix.

- method_params:

  Named list of per-method argument lists. Arguments are passed to the
  corresponding SCOP producer, never directly to its backend.

- n_clusters:

  Optional common number of domains. When `NULL`, methods that require a
  cluster count use the number of non-missing gold-standard classes. A
  method-specific `q` or `n_clusters` in `method_params` wins.

- metrics:

  Quality metrics selected by default when plotting. Supported values
  are `"ARI"`, `"NMI"`, and `"purity"`. All three are retained in the
  result summary.

- keep_objects:

  Whether to keep each successful method's full producer result. The
  default keeps only aligned labels and compact run metadata.

- install_missing:

  Whether missing optional backends may enter their producer's normal
  `check_r()` installation path. The default records them as unavailable
  without changing the R library.

- seed:

  Seed set inside each isolated method process and forwarded to
  producers with a public `seed` argument unless overridden in
  `method_params`.

- timeout:

  Maximum wall time in seconds for each isolated run, including process
  start and result serialization. `Inf` disables the timeout.

- poll_interval:

  Seconds between process-tree memory samples.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `spatial_benchmark_result` object (also inherits `benchmark_result`).
Use `result$summary` for the quality table,
[`SpatialBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md)
for visualization, and `$predictions` for the aligned labels.

## See also

[SpatialBenchmarkPlot](https://mengxu98.github.io/scop/reference/SpatialBenchmarkPlot.md),
[RunIntegrationBenchmark](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
visium_human_pancreas_sub$gold_domain <- factor(
  paste0("domain_", (seq_len(ncol(visium_human_pancreas_sub)) - 1) %% 3 + 1)
)
bench <- RunSpatialBenchmark(
  visium_human_pancreas_sub,
  gold_standard = "gold_domain",
  method_params = list(
    BayesSpace = list(n.PCs = 5, n.HVGs = 200),
    BANKSY = list(layer = "counts"),
    SmoothClust = list(layer = "counts", min_spots = 1)
  )
)
bench
SpatialBenchmarkPlot(data = bench)
} # }
```
