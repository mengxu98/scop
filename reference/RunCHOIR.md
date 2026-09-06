# Run CHOIR clustering

Runs the optional CHOIR backend on a single-modality `Seurat` object.
CHOIR builds and prunes a hierarchical clustering tree using
random-forest classifiers and permutation tests to identify
statistically distinct cell populations.

The upstream CHOIR records remain available in `srt@misc[[key]]`. This
wrapper also writes a stable cluster column to cell metadata and a
lightweight summary to `srt@tools[[tool_name]]`.

## Usage

``` r
RunCHOIR(
  srt,
  assay = NULL,
  layer = NULL,
  key = "CHOIR",
  cluster_colname = "CHOIR_cluster",
  tool_name = "CHOIR",
  alpha = 0.05,
  p_adjust = c("bonferroni", "fdr", "none"),
  feature_set = c("var", "all"),
  exclude_features = NULL,
  n_iterations = 100,
  n_trees = 50,
  min_accuracy = 0.5,
  max_clusters = "auto",
  normalization_method = c("none", "SCTransform"),
  batch.by = NULL,
  batch_correction_method = NULL,
  reduction = NULL,
  var_features = NULL,
  atac = FALSE,
  n_cores = 1,
  seed = 1,
  store_tool = TRUE,
  verbose = TRUE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay used by CHOIR. If `NULL`, the default assay is used.

- layer:

  Assay layer used by CHOIR. If `NULL`, `"data"` is used for `RNA` and
  `sketch` assays, and `"scale.data"` is used for `SCT` and `integrated`
  assays. Other assays require an explicit layer.

- key:

  Name used by CHOIR to store its complete results in `srt@misc`.

- cluster_colname:

  Metadata column used for the final CHOIR clusters.

- tool_name:

  Name used to store the SCOP result summary in `srt@tools`.

- alpha:

  Significance level used for CHOIR permutation tests.

- p_adjust:

  Multiple-testing correction used by CHOIR.

- feature_set:

  Whether CHOIR random forests use variable or all features.

- exclude_features:

  Features excluded from CHOIR random forests.

- n_iterations:

  Number of bootstrap iterations for each permutation test.

- n_trees:

  Number of trees in each random forest.

- min_accuracy:

  Minimum classifier accuracy required to keep clusters separate.

- max_clusters:

  Must be `"auto"`. Numeric limits are rejected because the pinned
  upstream backend can fail to terminate when the number of clusters
  plateaus.

- normalization_method:

  Normalization performed inside CHOIR. Use `"none"` for previously
  normalized data or `"SCTransform"` with a counts layer. The pinned
  backend does not support `"SCTransform"` for Seurat v5 `Assay5`
  objects.

- batch.by:

  Optional metadata column containing batch labels.

- batch_correction_method:

  Batch correction performed by CHOIR. If `NULL`, `"Harmony"` is
  selected when `batch.by` is supplied and `"none"` otherwise.

- reduction:

  Optional existing dimensional reduction supplied to CHOIR. This can be
  a reduction name in `srt` or a cell-by-dimension matrix.

- var_features:

  Features associated with `reduction`. If `NULL` and a reduction is
  supplied, variable features from `assay` are used.

- atac:

  Whether the selected assay contains ATAC-seq data.

- n_cores:

  Number of cores used by CHOIR. The pinned backend supports macOS and
  Linux; Windows execution is rejected before installation.

- seed:

  Random seed passed to CHOIR.

- store_tool:

  Whether to store a lightweight result summary in
  `srt@tools[[tool_name]]`.

- verbose:

  Whether to print progress messages.

- overwrite:

  Whether to replace existing CHOIR metadata, reduction, and `misc`
  entries, plus the `tools` entry when `store_tool = TRUE`.

- ...:

  Additional named arguments passed to the installed CHOIR entry point.
  Unsupported arguments produce an error rather than being silently
  ignored.

## Value

A `Seurat` object containing CHOIR clusters in `cluster_colname`,
complete upstream records in `srt@misc[[key]]`, and, when
`store_tool = TRUE`, a lightweight summary in `srt@tools[[tool_name]]`.

## References

Sant, C. et al. CHOIR improves significance-based detection of cell
types and states from single-cell data. *Nature Genetics* 57, 1309-1319
(2025).
[doi:10.1038/s41588-025-02148-8](https://doi.org/10.1038/s41588-025-02148-8)

## Examples

``` r
if (FALSE) { # \dontrun{
if (check_r("corceslab/CHOIR", verbose = FALSE)) {
  data(pancreas_sub)
  pancreas_sub <- Seurat::NormalizeData(pancreas_sub, verbose = FALSE)
  pancreas_sub <- RunCHOIR(
    pancreas_sub,
    assay = "RNA",
    n_cores = 2,
    verbose = FALSE
  )
  CellDimPlot(pancreas_sub, group.by = "CHOIR_cluster")
}
} # }
```
