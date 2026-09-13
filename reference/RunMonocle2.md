# Run Monocle2 analysis

Run Monocle2 analysis

## Usage

``` r
RunMonocle2(
  object,
  assay = NULL,
  layer = "counts",
  group.by = NULL,
  expressionFamily = "negbinomial.size",
  features = NULL,
  feature_type = "HVF",
  disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit",
  max_components = 2,
  reduction_method = "DDRTree",
  norm_method = "log",
  residualModelFormulaStr = NULL,
  pseudo_expr = 1,
  root_state = NULL,
  backend = NULL,
  n_neighbors = NULL,
  ddrtree_maxIter = NULL,
  ddrtree_ncenter = NULL,
  ddrtree_tol = NULL,
  show_plot = FALSE,
  xlab = NULL,
  ylab = NULL,
  seed = 11,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- group.by:

  Metadata column(s) used to color cells.

- expressionFamily:

  The distribution family to use for modeling gene expression.

- features:

  Features to use. Defaults to NULL, in which case features were
  determined by `feature_type`.

- feature_type:

  The type of features to use in the analysis. Possible values are "HVF"
  for highly variable features or "Disp" for features selected based on
  dispersion.

- disp_filter:

  Filter to use when `feature_type` is "Disp".

- max_components:

  The maximum number of dimensions to use for dimensionality reduction.

- reduction_method:

  The dimensionality reduction method to use. Possible values are
  `"DDRTree"`, `"ICA"`, `"tSNE"`, `"SimplePPT"`, `"L1-graph"`,
  `"SGL-tree"`.

- norm_method:

  The normalization method to use. Possible values are `"log"` and
  `"none"`.

- residualModelFormulaStr:

  A model formula specifying the effects to subtract.

- pseudo_expr:

  Amount to increase expression values before dimensionality reduction.

- root_state:

  The state to use as the root of the trajectory. If NULL, the
  trajectory is rooted at the first state after an initial ordering.

- backend:

  Deprecated and ignored. Monocle now accelerates cell ordering natively
  with C++, so the separate `"cpp"` backend has been removed.

- n_neighbors:

  Deprecated and ignored.

- ddrtree_maxIter:

  Optional maximum iteration count passed to `DDRTree::DDRTree()` when
  `reduction_method = "DDRTree"`. Lower values can speed up exploratory
  runs but may reduce agreement with the default Monocle2 trajectory.

- ddrtree_ncenter:

  Optional number of DDRTree centers. This can change trajectory
  topology and is intended for advanced exploratory use.

- ddrtree_tol:

  Optional convergence tolerance passed to `DDRTree::DDRTree()`.

- show_plot:

  Whether to print diagnostic plots during the run.

- xlab:

  Plot labels.

- ylab:

  Plot labels.

- seed:

  Random seed.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A Seurat object with the Monocle2 analysis results added to the @tools
slot.

## See also

[RunSlingshot](https://mengxu98.github.io/scop/reference/RunSlingshot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)

## Examples

``` r
if (interactive()) {
  data(pancreas_sub)
  pancreas_sub <- RunStandardWorkflow(pancreas_sub)
  pancreas_sub <- RunMonocle2(
    pancreas_sub,
    group.by = "SubCellType"
  )
  names(pancreas_sub@tools$Monocle2)
  trajectory <- pancreas_sub@tools$Monocle2$trajectory

  p1 <- CellDimPlot(
    pancreas_sub,
    group.by = "Monocle2_State",
    reduction = "DDRTree",
    label = TRUE,
    theme_use = "theme_blank"
  )
  p1

  p1 + trajectory

  FeatureDimPlot(
    pancreas_sub,
    features = "Monocle2_Pseudotime",
    reduction = "UMAP",
    theme_use = "theme_blank"
  )

  pancreas_sub <- RunMonocle2(
    pancreas_sub,
    feature_type = "Disp",
    disp_filter = "mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit"
  )
  trajectory <- pancreas_sub@tools$Monocle2$trajectory
  p2 <- CellDimPlot(
    pancreas_sub,
    group.by = "Monocle2_State",
    reduction = "DDRTree",
    label = TRUE,
    theme_use = "theme_blank"
  )
  p2

  p2 + trajectory

  FeatureDimPlot(
    pancreas_sub,
    features = "Monocle2_Pseudotime",
    reduction = "UMAP",
    theme_use = "theme_blank"
  )
}
```
