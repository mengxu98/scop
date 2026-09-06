# Annotate single cells using SingleR

Annotate single cells using SingleR

## Usage

``` r
RunSingleR(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_assay = "RNA",
  ref_assay = "RNA",
  genes = "de",
  de.method = "wilcox",
  sd.thresh = 1,
  de.n = NULL,
  aggr.ref = FALSE,
  aggr.args = list(),
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  prune = TRUE,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- query_group:

  Column name in the `srt_query` metadata that represents the cell
  grouping.

- ref_group:

  Column name in the `srt_ref` metadata that represents the cell
  grouping.

- query_assay:

  Assay to be used for the query data. Default is the default assay of
  the `srt_query` object.

- ref_assay:

  Assay to be used for the reference data. Default is the default assay
  of the `srt_ref` object.

- genes:

  A string containing `"de"`, indicating that markers should be
  calculated from `ref`. For back compatibility, other string values are
  allowed but will be ignored with a deprecation warning.

  Alternatively, if `ref` is *not* a list, `genes` can be either:

  - A list of lists of character vectors containing DE genes between
    pairs of labels.

  - A list of character vectors containing marker genes for each label.

  If `ref` *is* a list, `genes` can be a list of length equal to `ref`.
  Each element of the list should be one of the two above choices
  described for non-list `ref`, containing markers for labels in the
  corresponding entry of `ref`.

- de.method:

  String specifying how DE genes should be detected between pairs of
  labels. Defaults to `"classic"`, which sorts genes by the log-fold
  changes and takes the top `de.n`. Other options are `"wilcox"` and
  `"t"`, see Details. Ignored if `genes` is a list of markers/DE genes.

- sd.thresh:

  Deprecated and ignored.

- de.n:

  An integer scalar specifying the number of DE genes to use when
  `genes="de"`. If `de.method="classic"`, defaults to
  `500 * (2/3) ^ log2(N)` where `N` is the number of unique labels.
  Otherwise, defaults to 10. Ignored if `genes` is a list of markers/DE
  genes.

- aggr.ref, aggr.args:

  Arguments controlling the aggregation of the references prior to
  annotation, see
  [`trainSingleR`](https://rdrr.io/pkg/SingleR/man/trainSingleR.html).

- quantile:

  "quantile" parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- fine.tune:

  `"fine.tune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- tune.thresh:

  `"tune.thresh"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- prune:

  `"prune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

An annotate `Seurat` object. The annotation results are stored in the
`singler_annotation` column of the meta data, and the corresponding
scores are stored in the `singler_score` column.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)

## Examples

``` r
data(panc8_sub)
keep <- panc8_sub$celltype %in% c("ductal", "alpha", "beta")
panc8_sub <- Seurat::NormalizeData(panc8_sub[, keep], verbose = FALSE)
reference <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
query <- RunSingleR(
  srt_query = query, srt_ref = reference,
  ref_group = "celltype", genes = "de", cores = 1, verbose = FALSE
)
#> ℹ [2026-09-06 22:34:09] Data type is log-normalized
#> ℹ [2026-09-06 22:34:09] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-09-06 22:34:10] Data type is log-normalized
#> ℹ [2026-09-06 22:34:10] Detected `srt_ref` data type: "log_normalized_counts"
query <- RunStandardWorkflow(query, verbose = FALSE, linear_reduction_dims = 10)
#> ℹ [2026-09-06 22:34:11] Skip `log1p()` because `layer = data` is not "counts"
CellDimPlot(
  query, group.by = "singler_annotation",
  label = TRUE, legend.position = "bottom"
)

FeatureStatPlot(
  query, stat.by = "singler_score",
  group.by = "singler_annotation"
)
```
