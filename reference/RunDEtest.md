# Differential gene test

This function utilizes the Seurat package to perform a differential
expression (DE) test on gene expression data. Users have the flexibility
to specify custom cell groups, marker types, and various options for DE
analysis.

## Usage

``` r
RunDEtest(
  srt,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  markers_type = c("all", "paired", "conserved", "disturbed"),
  grouping.var = NULL,
  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  p.adjust.method = "bonferroni",
  layer = "data",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  A grouping variable in the dataset to define the groups or conditions
  for the differential test. If not provided, the function uses the
  "active.ident" variable in the Seurat object.

- group1:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the first group. If both group.by and group1 are provided,
  group1 takes precedence.

- group2:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the second group. This parameter is only used when group.by
  or group1 is provided.

- cells1:

  A vector of cell IDs specifying the cells that belong to group1. If
  provided, group1 is ignored.

- cells2:

  A vector of cell IDs specifying the cells that belong to group2. This
  parameter is only used when cells1 is provided.

- features:

  A vector of gene names specifying the features to consider for the
  differential test. If not provided, all features in the dataset are
  considered.

- markers_type:

  A character value specifying the type of markers to find. Possible
  values are "all", "paired", "conserved", and "disturbed".

- grouping.var:

  A character value specifying the grouping variable for finding
  conserved or disturbed markers. This parameter is only used when
  markers_type is "conserved" or "disturbed".

- meta.method:

  A character value specifying the method to use for combining p-values
  in the conserved markers test. Possible values are "maximump",
  "minimump", "wilkinsonp", "meanp", "sump", and "votep".

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

- only.pos:

  Only return positive markers (FALSE by default)

- fc.threshold:

  A numeric value used to filter genes for testing based on their
  average fold change between/among the two groups. Default is `1.5`.

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

- min.cells.feature:

  Minimum number of cells expressing the feature in at least one of the
  two groups, currently only used for poisson and negative binomial
  tests

- min.cells.group:

  Minimum number of cells in one of the groups

- norm.method:

  Normalization method for fold change calculation when layer is 'data'.
  Default is `"LogNormalize"`.

- p.adjust.method:

  A character value specifying the method to use for adjusting p-values.
  Default is `"bonferroni"`.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Assay to use in differential expression testing

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- ...:

  Additional arguments to pass to the
  [Seurat::FindMarkers](https://satijalab.org/seurat/reference/FindMarkers.html)
  function.

## See also

[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:36:04] Start standard processing workflow...
#> ℹ [2026-04-02 16:36:05] Checking a list of <Seurat>...
#> ! [2026-04-02 16:36:05] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:05] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:36:07] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:36:08] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:36:08] Number of available HVF: 2000
#> ℹ [2026-04-02 16:36:08] Finished check
#> ℹ [2026-04-02 16:36:08] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:36:08] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:36:12] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:36:12] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:36:13] Reorder clusters...
#> ℹ [2026-04-02 16:36:13] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:36:13] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:36:13] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:36:13] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:36:13] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:36:16] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:36:19] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType"
)
#> ℹ [2026-04-02 16:36:24] Data type is log-normalized
#> ℹ [2026-04-02 16:36:24] Start differential expression test
#> ℹ [2026-04-02 16:36:24] Find all markers(wilcox) among [1] 8 groups...
#> ℹ [2026-04-02 16:36:24] Using 1 core
#> ⠙ [2026-04-02 16:36:24] Running for Ductal [1/8] ■■■■■                         …
#> ✔ [2026-04-02 16:36:24] Completed 8 tasks in 1s
#> 
#> ℹ [2026-04-02 16:36:24] Building results
#> ! [2026-04-02 16:36:24] Found 8 failed results
#> ℹ [2026-04-02 16:36:25] ✖ Error details:
#> ℹ                       ✖ "Ductal": The total size of the 3 globals exported for future expression (‘FUN()’) is 556.27 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.94 MiB of class ‘function’), ‘data.use’ (4.32 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Ngn3-high-EP": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.18 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.54 MiB of class ‘function’), ‘data.use’ (2.63 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Beta": The total size of the 3 globals exported for future expression (‘FUN()’) is 553.99 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (550.43 MiB of class ‘function’), ‘data.use’ (3.56 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Ngn3-low-EP": The total size of the 3 globals exported for future expression (‘FUN()’) is 553.31 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.99 MiB of class ‘function’), ‘data.use’ (3.33 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Pre-endocrine": The total size of the 3 globals exported for future expression (‘FUN()’) is 554.90 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.03 MiB of class ‘function’), ‘data.use’ (3.87 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Alpha": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.68 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.89 MiB of class ‘function’), ‘data.use’ (2.79 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Epsilon": The total size of the 3 globals exported for future expression (‘FUN()’) is 554.39 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (550.70 MiB of class ‘function’), ‘data.use’ (3.68 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "Delta": The total size of the 3 globals exported for future expression (‘FUN()’) is 555.60 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (551.53 MiB of class ‘function’), ‘data.use’ (4.08 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> Error in `[.data.frame`(AllMarkers, , "group1"): undefined columns selected
AllMarkers <- dplyr::filter(
  pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht1 <- GroupHeatmap(
  pancreas_sub,
  features = AllMarkers$gene,
  feature_split = AllMarkers$group1,
  group.by = "SubCellType"
)
#> Error: object 'AllMarkers' not found
ht1$plot
#> Error: object 'ht1' not found

TopMarkers <- AllMarkers |>
  dplyr::group_by(gene) |>
  dplyr::top_n(1, avg_log2FC) |>
  dplyr::group_by(group1) |>
  dplyr::top_n(3, avg_log2FC)
#> Error: object 'AllMarkers' not found
ht2 <- GroupHeatmap(
  pancreas_sub,
  features = TopMarkers$gene,
  feature_split = TopMarkers$group1,
  group.by = "SubCellType",
  show_row_names = TRUE
)
#> Error: object 'TopMarkers' not found
ht2$plot
#> Error: object 'ht2' not found

pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  markers_type = "paired",
  cores = 2
)
#> ℹ [2026-04-02 16:36:30] Data type is log-normalized
#> ℹ [2026-04-02 16:36:30] Start differential expression test
#> ℹ [2026-04-02 16:36:30] Find paired markers(wilcox) among [1] 8 groups...
#> ℹ [2026-04-02 16:36:30] Using 2 cores
#> ⠙ [2026-04-02 16:36:30] Running for 1... [28/56] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-04-02 16:36:30] Completed 56 tasks in 4.9s
#> 
#> ℹ [2026-04-02 16:36:30] Building results
#> ! [2026-04-02 16:36:30] Found 56 failed results
#> ℹ [2026-04-02 16:36:35] ✖ Error details:
#> ℹ                       ✖ "1": The total size of the 3 globals exported for future expression (‘FUN()’) is 549.36 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (547.30 MiB of class ‘function’), ‘data.use’ (2.06 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "2": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.98 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.38 MiB of class ‘function’), ‘data.use’ (612.93 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "3": The total size of the 3 globals exported for future expression (‘FUN()’) is 550.56 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.11 MiB of class ‘function’), ‘data.use’ (2.44 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "4": The total size of the 3 globals exported for future expression (‘FUN()’) is 549.02 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (547.10 MiB of class ‘function’), ‘data.use’ (1.92 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "5": The total size of the 3 globals exported for future expression (‘FUN()’) is 549.84 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (547.63 MiB of class ‘function’), ‘data.use’ (2.21 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "6": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.59 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.12 MiB of class ‘function’), ‘data.use’ (1.47 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "7": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.81 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.62 MiB of class ‘function’), ‘data.use’ (1.19 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "8": The total size of the 3 globals exported for future expression (‘FUN()’) is 552.52 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.41 MiB of class ‘function’), ‘data.use’ (3.10 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "9": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.53 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.09 MiB of class ‘function’), ‘data.use’ (1.44 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "10": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.94 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.34 MiB of class ‘function’), ‘data.use’ (615.56 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "11": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.47 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.72 MiB of class ‘function’), ‘data.use’ (766.07 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "12": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.71 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.52 MiB of class ‘function’), ‘data.use’ (1.19 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "13": The total size of the 3 globals exported for future expression (‘FUN()’) is 548.27 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.57 MiB of class ‘function’), ‘data.use’ (1.70 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "14": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.55 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.77 MiB of class ‘function’), ‘data.use’ (800.79 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "15": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.66 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.84 MiB of class ‘function’), ‘data.use’ (2.82 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "16": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.83 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.60 MiB of class ‘function’), ‘data.use’ (1.22 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "17": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.69 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.84 MiB of class ‘function’), ‘data.use’ (867.58 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "18": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.14 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.50 MiB of class ‘function’), ‘data.use’ (655.45 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "19": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.98 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.37 MiB of class ‘function’), ‘data.use’ (1.61 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "20": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.76 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.57 MiB of class ‘function’), ‘data.use’ (1.19 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "21": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.68 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.18 MiB of class ‘function’), ‘data.use’ (514.65 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "22": The total size of the 3 globals exported for future expression (‘FUN()’) is 550.90 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.32 MiB of class ‘function’), ‘data.use’ (2.58 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "23": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.57 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.76 MiB of class ‘function’), ‘data.use’ (837.77 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "24": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.80 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.91 MiB of class ‘function’), ‘data.use’ (914.63 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "25": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.00 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (543.73 MiB of class ‘function’), ‘data.use’ (277.48 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "26": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.25 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.88 MiB of class ‘function’), ‘data.use’ (1.37 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "27": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.67 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.82 MiB of class ‘function’), ‘data.use’ (862.69 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "28": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.10 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.44 MiB of class ‘function’), ‘data.use’ (679.07 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "29": The total size of the 3 globals exported for future expression (‘FUN()’) is 552.52 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.42 MiB of class ‘function’), ‘data.use’ (3.10 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "30": The total size of the 3 globals exported for future expression (‘FUN()’) is 552.25 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.24 MiB of class ‘function’), ‘data.use’ (3.00 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "31": The total size of the 3 globals exported for future expression (‘FUN()’) is 549.37 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (547.33 MiB of class ‘function’), ‘data.use’ (2.04 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "32": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.55 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.77 MiB of class ‘function’), ‘data.use’ (2.78 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "33": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.78 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.60 MiB of class ‘function’), ‘data.use’ (1.19 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "34": The total size of the 3 globals exported for future expression (‘FUN()’) is 548.07 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.44 MiB of class ‘function’), ‘data.use’ (1.63 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "35": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.32 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.30 MiB of class ‘function’), ‘data.use’ (1.02 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "36": The total size of the 3 globals exported for future expression (‘FUN()’) is 548.60 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.79 MiB of class ‘function’), ‘data.use’ (1.81 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "37": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.43 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (546.01 MiB of class ‘function’), ‘data.use’ (1.42 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "38": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.51 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.74 MiB of class ‘function’), ‘data.use’ (789.51 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "39": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.67 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.15 MiB of class ‘function’), ‘data.use’ (532.83 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "40": The total size of the 3 globals exported for future expression (‘FUN()’) is 548.99 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (547.05 MiB of class ‘function’), ‘data.use’ (1.94 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "41": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.16 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.83 MiB of class ‘function’), ‘data.use’ (1.33 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "42": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.09 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.47 MiB of class ‘function’), ‘data.use’ (640.31 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "43": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.98 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.04 MiB of class ‘function’), ‘data.use’ (967.59 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "44": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.62 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.48 MiB of class ‘function’), ‘data.use’ (1.14 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "45": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.16 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.50 MiB of class ‘function’), ‘data.use’ (671.99 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "46": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.93 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (549.01 MiB of class ‘function’), ‘data.use’ (2.91 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "47": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.19 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.50 MiB of class ‘function’), ‘data.use’ (700.98 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "48": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.63 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.47 MiB of class ‘function’), ‘data.use’ (1.16 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "49": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.81 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.28 MiB of class ‘function’), ‘data.use’ (548.50 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "50": The total size of the 3 globals exported for future expression (‘FUN()’) is 547.02 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.72 MiB of class ‘function’), ‘data.use’ (1.29 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "51": The total size of the 3 globals exported for future expression (‘FUN()’) is 545.64 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.80 MiB of class ‘function’), ‘data.use’ (851.09 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "52": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.56 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (544.07 MiB of class ‘function’), ‘data.use’ (494.69 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "53": The total size of the 3 globals exported for future expression (‘FUN()’) is 551.65 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (548.82 MiB of class ‘function’), ‘data.use’ (2.83 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "54": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.30 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.24 MiB of class ‘function’), ‘data.use’ (1.06 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "55": The total size of the 3 globals exported for future expression (‘FUN()’) is 546.42 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (545.32 MiB of class ‘function’), ‘data.use’ (1.10 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "56": The total size of the 3 globals exported for future expression (‘FUN()’) is 544.13 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (543.79 MiB of class ‘function’), ‘data.use’ (345.04 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> Error in `[.data.frame`(PairedMarkers, , "group1"): undefined columns selected
PairedMarkers <- dplyr::filter(
  pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht3 <- GroupHeatmap(
  pancreas_sub,
  features = PairedMarkers$gene,
  feature_split = PairedMarkers$group1,
  group.by = "SubCellType"
)
#> Error: object 'PairedMarkers' not found
ht3$plot
#> Error: object 'ht3' not found

data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-04-02 16:36:36] Run integration workflow...
#> ℹ [2026-04-02 16:36:36] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-02 16:36:37] Checking a list of <Seurat>...
#> ! [2026-04-02 16:36:37] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:38] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-04-02 16:36:39] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:40] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-04-02 16:36:40] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:41] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-04-02 16:36:42] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:42] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:43] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-04-02 16:36:43] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:36:43] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:45] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-02 16:36:45] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:36:45] Number of available HVF: 2000
#> ℹ [2026-04-02 16:36:46] Finished check
#> ℹ [2026-04-02 16:36:49] Perform Uncorrected integration
#> ℹ [2026-04-02 16:36:51] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:36:52] Perform "pca" linear dimension reduction
#> ℹ [2026-04-02 16:36:57] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-02 16:36:57] Reorder clusters...
#> ℹ [2026-04-02 16:36:58] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:36:58] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:36:58] Error when performing `Seurat::FindClusters()`. Skip this step
#> ℹ [2026-04-02 16:36:58] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ℹ [2026-04-02 16:37:01] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:50)
#> ✔ [2026-04-02 16:37:06] Uncorrected integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("celltype", "tech")
)


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-04-02 16:37:14] Data type is log-normalized
#> ℹ [2026-04-02 16:37:14] Start differential expression test
#> ℹ [2026-04-02 16:37:14] Find conserved markers(wilcox) among [1] 13 groups...
#> ℹ [2026-04-02 16:37:14] Using 2 cores
#> ⠙ [2026-04-02 16:37:14] Running for delta... [7/13] ■■■■■■■■■■■■■■■■■          …
#> ✔ [2026-04-02 16:37:14] Completed 13 tasks in 2.1s
#> 
#> ℹ [2026-04-02 16:37:14] Building results
#> ! [2026-04-02 16:37:14] Found 11 failed results
#> ℹ [2026-04-02 16:37:16] ✖ Error details:
#> ℹ                       ✖ "delta": The total size of the 3 globals exported for future expression (‘FUN()’) is 931.24 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (930.07 MiB of class ‘function’), ‘data.use’ (1.17 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "gamma": The total size of the 3 globals exported for future expression (‘FUN()’) is 930.07 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (929.27 MiB of class ‘function’), ‘data.use’ (819.72 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "acinar": The total size of the 3 globals exported for future expression (‘FUN()’) is 930.05 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (929.26 MiB of class ‘function’), ‘data.use’ (808.84 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "beta": The total size of the 3 globals exported for future expression (‘FUN()’) is 933.09 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (931.32 MiB of class ‘function’), ‘data.use’ (1.77 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "alpha": The total size of the 3 globals exported for future expression (‘FUN()’) is 940.96 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (936.61 MiB of class ‘function’), ‘data.use’ (4.35 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "ductal": The total size of the 3 globals exported for future expression (‘FUN()’) is 939.55 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (935.64 MiB of class ‘function’), ‘data.use’ (3.91 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "mast": The total size of the 3 globals exported for future expression (‘FUN()’) is 935.19 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (932.73 MiB of class ‘function’), ‘data.use’ (2.46 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "macrophage": The total size of the 3 globals exported for future expression (‘FUN()’) is 932.87 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (931.15 MiB of class ‘function’), ‘data.use’ (1.72 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "quiescent-stellate": The total size of the 3 globals exported for future expression (‘FUN()’) is 932.87 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (931.15 MiB of class ‘function’), ‘data.use’ (1.71 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "endothelial": The total size of the 3 globals exported for future expression (‘FUN()’) is 939.33 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (935.51 MiB of class ‘function’), ‘data.use’ (3.81 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "schwann": The total size of the 3 globals exported for future expression (‘FUN()’) is 940.11 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (936.03 MiB of class ‘function’), ‘data.use’ (4.08 MiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> Error in x[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj",     "gene", "group1", "group2")]: incorrect number of dimensions
ConservedMarkers1 <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht4 <- GroupHeatmap(
  panc8_sub,
  layer = "data",
  features = ConservedMarkers1$gene,
  feature_split = ConservedMarkers1$group1,
  group.by = "tech",
  split.by = "celltype",
  within_groups = TRUE
)
#> Error: object 'ConservedMarkers1' not found
ht4$plot
#> Error: object 'ht4' not found

panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "tech",
  grouping.var = "celltype",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-04-02 16:37:22] Data type is log-normalized
#> ℹ [2026-04-02 16:37:22] Start differential expression test
#> ℹ [2026-04-02 16:37:22] Find conserved markers(wilcox) among [1] 5 groups...
#> ℹ [2026-04-02 16:37:22] Using 2 cores
#> ⠙ [2026-04-02 16:37:22] Running for celseq... [3/5] ■■■■■■■■■■■■■■■■■■■        …
#> ✔ [2026-04-02 16:37:22] Completed 5 tasks in 1.4s
#> 
#> ℹ [2026-04-02 16:37:22] Building results
#> ! [2026-04-02 16:37:22] Found 5 failed results
#> ℹ [2026-04-02 16:37:24] ✖ Error details:
#> ℹ                       ✖ "celseq": The total size of the 3 globals exported for future expression (‘FUN()’) is 929.77 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (929.11 MiB of class ‘function’), ‘data.use’ (680.59 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "celseq2": The total size of the 3 globals exported for future expression (‘FUN()’) is 930.64 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (929.71 MiB of class ‘function’), ‘data.use’ (956.01 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "smartseq2": The total size of the 3 globals exported for future expression (‘FUN()’) is 927.87 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (927.79 MiB of class ‘function’), ‘data.use’ (79.49 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "fluidigmc1": The total size of the 3 globals exported for future expression (‘FUN()’) is 929.21 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (928.74 MiB of class ‘function’), ‘data.use’ (480.19 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> ℹ                       ✖ "indrop": The total size of the 3 globals exported for future expression (‘FUN()’) is 929.10 MiB. This exceeds the maximum allowed size 500.00 MiB per plan() argument 'maxSizeOfObjects'. This limit is set to protect against transfering too large objects to parallel workers by mistake, which may not be intended and could be costly. See help("future.globals.maxSize", package = "future") for how to adjust or remove the default threshold via an R option There are three globals: ‘FUN’ (928.67 MiB of class ‘function’), ‘data.use’ (440.96 KiB of class ‘S4’) and ‘j’ (133 bytes of class ‘numeric’)
#> Error in x[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj",     "gene", "group1", "group2")]: incorrect number of dimensions
ConservedMarkers2 <- dplyr::filter(
  panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht4 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = ConservedMarkers2$gene,
  feature_split = ConservedMarkers2$group1,
  group.by = "tech",
  split.by = "celltype"
)
#> Error: object 'ConservedMarkers2' not found
ht4$plot
#> Error: object 'ht4' not found

panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "disturbed",
  cores = 2
)
#> ℹ [2026-04-02 16:37:30] Data type is log-normalized
#> ℹ [2026-04-02 16:37:30] Start differential expression test
#> ℹ [2026-04-02 16:37:30] Find disturbed markers(wilcox) among [1] 13 groups...
#> ℹ [2026-04-02 16:37:30] Using 2 cores
#> ⠙ [2026-04-02 16:37:30] Running for delta... [7/13] ■■■■■■■■■■■■■■■■■          …
#> ✔ [2026-04-02 16:37:30] Completed 13 tasks in 46.6s
#> 
#> ℹ [2026-04-02 16:37:30] Building results
#> ! [2026-04-02 16:37:30] Found 9 failed results
#> ℹ [2026-04-02 16:38:17] ✖ Error details:
#> ℹ                       ✖ "delta": undefined columns selected
#> ℹ                       ✖ "gamma": undefined columns selected
#> ℹ                       ✖ "acinar": undefined columns selected
#> ℹ                       ✖ "beta": undefined columns selected
#> ℹ                       ✖ "alpha": undefined columns selected
#> ℹ                       ✖ "ductal": undefined columns selected
#> ℹ                       ✖ "mast": undefined columns selected
#> ℹ                       ✖ "macrophage": undefined columns selected
#> ℹ                       ✖ "quiescent-stellate": undefined columns selected
#> Error in `[.data.frame`(DisturbedMarkers, , "group1"): undefined columns selected
DisturbedMarkers <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht5 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers$gene,
  feature_split = DisturbedMarkers$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> Error: object 'DisturbedMarkers' not found
ht5$plot
#> Error: object 'ht5' not found

gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'DisturbedMarkers' not found
DisturbedMarkers_specific <- DisturbedMarkers[
  DisturbedMarkers$gene %in% gene_specific,
]
#> Error: object 'DisturbedMarkers' not found
ht6 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> Error: object 'DisturbedMarkers_specific' not found
ht6$plot
#> Error: object 'ht6' not found

ht7 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  aggregate_fun = function(x) mean(expm1(x)) + 1,
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  grouping.var = "tech",
  numerator = "smartseq2"
)
#> Error: object 'DisturbedMarkers_specific' not found
ht7$plot
#> Error: object 'ht7' not found
```
