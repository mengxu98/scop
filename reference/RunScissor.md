# Run Scissor phenotype-associated cell selection

Run Scissor phenotype-associated cell selection

## Usage

``` r
RunScissor(
  srt,
  bulk_dataset,
  phenotype = NULL,
  condition.by = NULL,
  positive = NULL,
  assay = NULL,
  bulk_assay = "counts",
  layer = "counts",
  features = NULL,
  family = c("gaussian", "binomial", "cox"),
  backend = c("cpp", "r"),
  alpha = NULL,
  cutoff = 0.2,
  tag = NULL,
  graph = NULL,
  dims = 1:10,
  nfeatures = 2000,
  seed = 123,
  prefix = "Scissor",
  tool_name = "Scissor",
  store_inputs = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object containing single-cell expression data.

- bulk_dataset:

  A bulk expression matrix-like object or a `SummarizedExperiment`. Rows
  are genes and columns are bulk samples.

- phenotype:

  Phenotype annotation for bulk samples. For `family = "binomial"`,
  character or factor values are converted to 0/1 according to
  `positive`. If `NULL` and `bulk_dataset` is a `SummarizedExperiment`,
  `condition.by` is used.

- condition.by:

  Column in `colData(bulk_dataset)` used as phenotype when
  `phenotype = NULL`.

- positive:

  Positive phenotype level for binomial Scissor. If `NULL`, the second
  sorted level is used.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- bulk_assay:

  Assay name used when `bulk_dataset` is a `SummarizedExperiment`.

- layer:

  Assay layer used from `srt`.

- features:

  Optional genes used before intersecting bulk and single-cell features.

- family:

  Regression family passed to Scissor.

- backend:

  Scissor backend. `"r"` calls the upstream package path and `"cpp"`
  uses the optimized `scop` path. Legacy aliases `"original"` and
  `"scop"` are accepted as `"r"` and `"cpp"`, respectively.

- alpha:

  Scissor alpha search values. If `NULL`, Scissor's default grid is
  used.

- cutoff:

  Maximum selected-cell fraction used by Scissor's alpha search.

- tag:

  Optional phenotype labels for printed/stored summaries.

- graph:

  Existing Seurat graph used by the `"cpp"` backend. If `NULL`, an assay
  SNN graph or
  [`standard_scop()`](https://mengxu98.github.io/scop/reference/standard_scop.md)
  graph is reused when present, otherwise a temporary Scissor-like graph
  is built.

- dims:

  PCA dimensions used when a temporary graph is built.

- nfeatures:

  Number of variable features used when a temporary graph is built.

- seed:

  Random seed used by Scissor's alpha loop.

- prefix:

  Prefix for metadata column names.

- tool_name:

  Name of the `srt@tools` entry.

- store_inputs:

  Whether to store Scissor regression inputs in `srt@tools`. Default is
  `FALSE` to avoid large objects.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with Scissor status and coefficient columns in
metadata and a Scissor result bundle in `srt@tools[[tool_name]]`.

## References

Sun, D. et al. Identifying phenotype-associated subpopulations by
integrating bulk and single-cell sequencing data. *Nature Biotechnology*
(2021).
[doi:10.1038/s41587-021-01091-3](https://doi.org/10.1038/s41587-021-01091-3)

## See also

[ScissorPlot](https://mengxu98.github.io/scop/reference/ScissorPlot.md)

## Examples

``` r
data(panc8_sub)
data(islet_bulk)
panc8_sub <- standard_scop(panc8_sub, verbose = FALSE)
#> ℹ [2026-06-24 19:06:29] Skip `log1p()` because `layer = data` is not "counts"
panc8_sub <- RunScissor(
  panc8_sub,
  bulk_dataset = islet_bulk,
  condition.by = "condition",
  positive = "bfa",
  family = "binomial",
  features = head(
    intersect(
      rownames(panc8_sub),
      rownames(SummarizedExperiment::assay(islet_bulk, "counts"))
    ), 1000
  ),
  alpha = 0.2,
  cutoff = 0.5
)
#> ℹ [2026-06-24 19:06:39] Scissor alpha 0.2 selected 4 positive and 445 negative cells (28.062%)
#> ✔ [2026-06-24 19:06:39] Scissor stored 4 Scissor+ and 445 Scissor- cells

ScissorPlot(
  panc8_sub,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
```
