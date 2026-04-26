# Milo differential abundance wrapper

Method-specific Milo wrapper used by
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
Stores both group-level and neighborhood-level result lists. When miloR
is available, graph-ready neighborhood node/edge data are also stored
under `details$milo_graph_data` for visualization.

## Usage

``` r
RunProportionTestMilo(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  milo_k = 20L,
  milo_d = 30L,
  n_bootstrap = 500,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Cell-type metadata column.

- split.by:

  Condition metadata column.

- sample.by:

  Sample metadata column.

- comparison:

  Optional comparisons to perform.

- milo_k:

  k used for Milo graph building.

- milo_d:

  Number of dimensions used by Milo.

- n_bootstrap:

  Bootstrap iterations for summary/fallback.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
