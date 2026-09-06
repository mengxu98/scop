# Run FWP feature-weight phenotype scoring

Run FWP feature-weight phenotype scoring

## Usage

``` r
RunFWP(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nfeatures = 2000,
  phenotype.by = NULL,
  positive = NULL,
  weights = NULL,
  score.name = "FWP_Score",
  tool_name = "FWP",
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object or expression matrix with genes in rows and cells in
  columns.

- assay:

  Assay used for Seurat input.

- layer:

  Layer used for Seurat input.

- features:

  Features used for scoring. If `NULL`, the most variable genes are
  selected.

- nfeatures:

  Number of variable genes selected when `features = NULL`.

- phenotype.by:

  Metadata column used to train feature weights.

- positive:

  Positive phenotype label. If `NULL`, the last ordered label is used.

- weights:

  Optional named vector of precomputed feature weights.

- score.name:

  Metadata column for the FWP score.

- tool_name:

  Name used in `srt@tools`.

- verbose:

  Whether to print progress messages.

## Value

A modified `Seurat` object or a result bundle for matrix input.

## References

Wang Q, Song JJ, Zhang F. Feature-weight based measurement of cancerous
transcriptome using cohort-wide and sample-specific information.
Cellular Oncology, 2024. doi:10.1007/s13402-023-00879-6.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:08:21] Start standard processing workflow...
#> ℹ [2026-09-06 22:08:22] Checking a list of <Seurat>...
#> ! [2026-09-06 22:08:22] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:08:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:08:22] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:08:22] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:08:22] Number of available HVF: 2000
#> ℹ [2026-09-06 22:08:22] Finished check
#> ℹ [2026-09-06 22:08:22] Perform `ScaleData()`
#> ℹ [2026-09-06 22:08:22] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:08:23] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:08:23] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:08:23] Reorder clusters...
#> ℹ [2026-09-06 22:08:23] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:08:23] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:08:30] Standard processing workflow completed
pancreas_sub$IsEndocrine <- pancreas_sub$CellType == "Endocrine"
pancreas_sub <- RunFWP(
  pancreas_sub,
  phenotype.by = "IsEndocrine",
  nfeatures = 300
)
#> ℹ [2026-09-06 22:08:30] Run FWP scoring with 300 features
FeatureDimPlot(
  pancreas_sub,
  features = "FWP_Score"
) + CellDimPlot(
  pancreas_sub,
  group.by = "IsEndocrine",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
```
