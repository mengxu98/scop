# Run Monocle2 analysis

Run Monocle2 analysis

## Usage

``` r
RunMonocle2(
  srt,
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
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `"counts"`.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- expressionFamily:

  The distribution family to use for modeling gene expression. Default
  is `"negbinomial.size"`.

- features:

  A character vector of features to use. Defaults to NULL, in which case
  features were determined by `feature_type`.

- feature_type:

  The type of features to use in the analysis. Possible values are "HVF"
  for highly variable features or "Disp" for features selected based on
  dispersion. Default is `"HVF"`.

- disp_filter:

  A string specifying the filter to use when `feature_type` is "Disp".
  Default is
  `"mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit"`.

- max_components:

  The maximum number of dimensions to use for dimensionality reduction.
  Default is `2`.

- reduction_method:

  The dimensionality reduction method to use. Possible values are
  `"DDRTree"`, `"ICA"`, `"tSNE"`, `"SimplePPT"`, `"L1-graph"`,
  `"SGL-tree"`. Default is `"DDRTree"`.

- norm_method:

  The normalization method to use. Possible values are `"log"` and
  `"none"`. Default is `"log"`.

- residualModelFormulaStr:

  A model formula specifying the effects to subtract. Default is NULL.

- pseudo_expr:

  Amount to increase expression values before dimensionality reduction.
  Default is 1.

- root_state:

  The state to use as the root of the trajectory. If NULL, will prompt
  for user input.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

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
  pancreas_sub <- standard_scop(pancreas_sub)
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
