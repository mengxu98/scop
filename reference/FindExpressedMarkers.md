# Find Expressed Markers

Find Expressed Markers

## Usage

``` r
FindExpressedMarkers(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  min.expression = 0,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  fc.name = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  only.pos = FALSE,
  min.cells.group = 3,
  min.cells.feature = 3,
  norm.method = "LogNormalize",
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object

- ident.1:

  Identity class to define markers for; pass an object of class `phylo`
  or 'clustertree' to find markers for a node in a cluster tree; passing
  'clustertree' requires
  [`BuildClusterTree`](https://satijalab.org/seurat/reference/BuildClusterTree.html)
  to have been run

- ident.2:

  A second identity class for comparison; if `NULL`, use all other cells
  for comparison; if an object of class `phylo` or 'clustertree' is
  passed to `ident.1`, must pass a node to find markers for

- cells.1:

  Vector of cell names belonging to group 1

- cells.2:

  Vector of cell names belonging to group 2

- features:

  A character vector or a named list of features to plot. Features can
  be gene names in Assay or names of numeric columns in meta.data.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `data`.

- min.expression:

  The min.expression used.

- test.use:

  Denotes which test to use. Available options are:

  - "wilcox" : Identifies differentially expressed genes between two
    groups of cells using a Wilcoxon Rank Sum test (default); will use a
    fast implementation by Presto if installed

  - "wilcox_limma" : Identifies differentially expressed genes between
    two groups of cells using the limma implementation of the Wilcoxon
    Rank Sum test; set this option to reproduce results from Seurat v4

  - "bimod" : Likelihood-ratio test for single cell gene expression,
    (McDavid et al., Bioinformatics, 2013)

  - "roc" : Identifies 'markers' of gene expression using ROC analysis.
    For each gene, evaluates (using AUC) a classifier built on that gene
    alone, to classify between two groups of cells. An AUC value of 1
    means that expression values for this gene alone can perfectly
    classify the two groupings (i.e. Each of the cells in cells.1
    exhibit a higher level than each of the cells in cells.2). An AUC
    value of 0 also means there is perfect classification, but in the
    other direction. A value of 0.5 implies that the gene has no
    predictive power to classify the two groups. Returns a 'predictive
    power' (abs(AUC-0.5) \* 2) ranked matrix of putative differentially
    expressed genes.

  - "t" : Identify differentially expressed genes between two groups of
    cells using the Student's t-test.

  - "negbinom" : Identifies differentially expressed genes between two
    groups of cells using a negative binomial generalized linear model.
    Use only for UMI-based datasets

  - "poisson" : Identifies differentially expressed genes between two
    groups of cells using a poisson generalized linear model. Use only
    for UMI-based datasets

  - "LR" : Uses a logistic regression framework to determine
    differentially expressed genes. Constructs a logistic regression
    model predicting group membership based on each feature individually
    and compares this to a null model with a likelihood ratio test.

  - "MAST" : Identifies differentially expressed genes between two
    groups of cells using a hurdle model tailored to scRNA-seq data.
    Utilizes the MAST package to run the DE testing.

  - "DESeq2" : Identifies differentially expressed genes between two
    groups of cells based on a model using DESeq2 which uses a negative
    binomial distribution (Love et al, Genome Biology, 2014).This test
    does not support pre-filtering of genes based on average difference
    (or percent detection rate) between cell groups. However, genes may
    be pre-filtered based on their minimum detection rate (min.pct)
    across both cell groups. To use this method, please install DESeq2,
    using the instructions at
    https://bioconductor.org/packages/release/bioc/html/DESeq2.html

- logfc.threshold:

  Limit testing to genes which show, on average, at least X-fold
  difference (log-scale) between the two groups of cells. Default is 0.1
  Increasing logfc.threshold speeds up the function, but can miss weaker
  signals. If the `slot` parameter is "scale.data" no filtering is
  performed.

- base:

  The base with respect to which logarithms are computed.

- pseudocount.use:

  Pseudocount to add to averaged expression values when calculating
  logFC. 1 by default.

- mean.fxn:

  Function to use for fold change or average difference calculation. The
  default depends on the the value of `fc.slot`:

  - "counts" : difference in the log of the mean counts, with
    pseudocount.

  - "data" : difference in the log of the average exponentiated data,
    with pseudocount. This adjusts for differences in sequencing depth
    between cells, and assumes that "data" has been log-normalized.

  - "scale.data" : difference in the means of scale.data.

- fc.name:

  Name of the fold change, average difference, or custom function column
  in the output data.frame. If NULL, the fold change column will be
  named according to the logarithm base (eg, "avg_log2FC"), or if using
  the scale.data slot "avg_diff".

- min.pct:

  only test genes that are detected in a minimum fraction of min.pct
  cells in either of the two populations. Meant to speed up the function
  by not testing genes that are very infrequently expressed. Default is
  0.01

- min.diff.pct:

  only test genes that show a minimum difference in the fraction of
  detection between the two groups. Set to -Inf by default

- max.cells.per.ident:

  Down sample each identity class to a max number. Default is no
  downsampling. Not activated by default (set to Inf)

- latent.vars:

  Variables to test, used only when `test.use` is one of 'LR',
  'negbinom', 'poisson', or 'MAST'

- only.pos:

  Only return positive markers (FALSE by default)

- min.cells.group:

  Minimum number of cells in one of the groups

- min.cells.feature:

  Minimum number of cells expressing the feature in at least one of the
  two groups, currently only used for poisson and negative binomial
  tests

- norm.method:

  Normalization method for fold change calculation when `slot` is
  “`data`”

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Arguments passed to other methods and to specific DE methods

## See also

\[RunDEtest\]

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
markers <- FindExpressedMarkers(
  pancreas_sub,
  cells.1 = SeuratObject::WhichCells(
    pancreas_sub,
    expression = Phase == "G2M"
  )
)
#> ⠙ [2026-01-13 09:00:40] Running for 1 [1/6145] ■                               …
#> ⠹ [2026-01-13 09:00:40] Running for 3140 [3140/6145] ■■■■■■■■■■■■■■■■          …
#> ✔ [2026-01-13 09:00:40] Completed 6145 tasks in 3.3s
#> 
head(markers)
#>               p_val avg_log2FC pct.1 pct.2    p_val_adj
#> Hmgb2  1.964385e-37  2.1638592 1.000 0.526 3.142624e-33
#> Tuba1b 3.652521e-33  1.8756996 0.986 0.801 5.843303e-29
#> Ran    2.573020e-31  1.2068807 1.000 0.946 4.116317e-27
#> H2afx  9.845758e-31  1.4231361 0.971 0.519 1.575124e-26
#> Ptma   3.541885e-30  0.9021115 1.000 0.999 5.666307e-26
#> Tubb5  1.550857e-29  1.4861586 1.000 0.940 2.481062e-25

FeatureStatPlot(
  pancreas_sub,
  rownames(markers)[1],
  group.by = "Phase",
  add_point = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
```
