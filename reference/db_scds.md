# Run doublet-calling with scds

Run doublet-calling with scds

## Usage

``` r
db_scds(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  method = c("hybrid", "cxds", "bcds"),
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- method:

  The method to be used for doublet-calling. Options are `"hybrid"`,
  `"cxds"`, or `"bcds"`.

- ...:

  Additional arguments to be passed to
  [`scds::cxds_bcds_hybrid()`](https://rdrr.io/pkg/scds/man/cxds_bcds_hybrid.html).

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> Registered S3 method overwritten by 'pROC':
#>   method   from            
#>   plot.roc spatstat.explore
#> Error in xgboost(mm, nrounds = nmax, tree_method = "hist", nthread = 2,     early_stopping_rounds = 2, subsample = 0.5, objective = "binary:logistic",     verbose = 0): argument "y" is missing, with no default
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scds_hybrid_class"
)
#> Error in CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scds_hybrid_class"): "db.scds_hybrid_class" is not in the meta.data of srt object

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scds_hybrid_score"
)
#> Error in FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scds_hybrid_score"): There are no valid features present.
```
