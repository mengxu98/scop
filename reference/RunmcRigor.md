# Run mcRigor metacell partition assessment

Run mcRigor metacell partition assessment

## Usage

``` r
RunmcRigor(
  srt,
  cell_membership = NULL,
  metacell.by = NULL,
  mode = c("detect", "optimize"),
  tgamma = NULL,
  gamma_names = NULL,
  assay_type = c("RNA", "ATAC"),
  Gammas = NULL,
  aggregate_method = c("mean", "sum", "geom"),
  output_file = NULL,
  Nrep = 1,
  gene_filter = 0.1,
  feature_use = 2000,
  cor_method = c("pearson", "spearman"),
  prePro = TRUE,
  test_cutoff = 0.01,
  thre_smooth = TRUE,
  thre_bw = 1/6,
  D_bw = 10,
  optim_method = c("tradeoff", "dub_rate_large", "dub_rate_small"),
  weight = 0.5,
  dub_rate = 0.1,
  draw = FALSE,
  pur_metric = NULL,
  check_purity = TRUE,
  fields = NULL,
  step_save = FALSE,
  prefix = "mcRigor",
  tool_name = "mcRigor",
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object containing the original single-cell data.

- cell_membership:

  A data frame or matrix with cells in rows and one or more metacell
  partitions in columns. Row names should be cell names. If row names
  are missing and the row count equals `ncol(srt)`, cells are matched in
  `colnames(srt)` order.

- metacell.by:

  Metadata column(s) in `srt` used as metacell partitions when
  `cell_membership = NULL`.

- mode:

  mcRigor task. `"detect"` calls `mcRigor_DETECT()` for one partition;
  `"optimize"` calls `mcRigor_OPTIMIZE()` across candidate partitions.

- tgamma:

  Target partition/gamma for `"detect"`. Can be a membership column name
  or the numeric gamma label used by mcRigor. If `NULL`, the first
  membership column is used.

- gamma_names:

  Optional gamma labels for membership columns. mcRigor requires
  numeric-like column labels; non-numeric labels are mapped internally
  to `1:ncol(cell_membership)` and recorded in the stored result.

- assay_type:

  Assay type passed to mcRigor.

- Gammas:

  Candidate gamma labels for `"optimize"`. Can use original membership
  column names or mapped mcRigor gamma labels.

- aggregate_method:

  Metacell aggregation method passed to mcRigor.

- output_file:

  Optional path where mcRigor writes the `TabMC` RDS file. If `NULL`, a
  temporary file is used to avoid creating files in the working
  directory.

- Nrep:

  Number of permutation repetitions used by mcRigor.

- gene_filter, feature_use, cor_method, prePro, test_cutoff,
  thre_smooth, thre_bw:

  Parameters forwarded to mcRigor.

- D_bw, optim_method, weight, dub_rate:

  Optimization parameters forwarded to `mcRigor_OPTIMIZE()`.

- draw, pur_metric, check_purity, fields, step_save:

  Plotting, purity, and intermediate-save parameters forwarded to
  mcRigor.

- prefix:

  Prefix for metadata columns written to `srt`.

- tool_name:

  Name of the `srt@tools` entry used to store results.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with mcRigor metadata and a result list stored in
`srt@tools[[tool_name]]`.

## References

Liu, P. and Li, J.J. (2024). mcRigor: a statistical method to enhance
the rigor of metacell partitioning in single-cell data analysis.
*bioRxiv*.
[doi:10.1101/2024.10.30.621093](https://doi.org/10.1101/2024.10.30.621093)

## Examples

``` r
data(pancreas_sub)
set.seed(11)
pancreas_sub <- standard_scop(
  pancreas_sub,
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:20,
  nonlinear_reduction_dims = 2,
  verbose = FALSE
)
#> ℹ [2026-06-24 19:18:07] Skip `log1p()` because `layer = data` is not "counts"
mc <- RunMetaCell(
  pancreas_sub,
  method = "supercell",
  gamma = 25
)
#> ℹ [2026-06-24 19:18:19] Running SuperCell with gamma = 25, k.knn = 5 on 1000 cells
#> Error in loadNamespace(name): there is no package called ‘SuperCell’

membership <- data.frame(
  Metacell = mc@misc[["cell_membership"]],
  row.names = names(mc@misc[["cell_membership"]])
)
#> Error: object 'mc' not found

pancreas_sub <- RunmcRigor(
  mc@misc[["original_srt"]],
  cell_membership = membership,
  Nrep = 1,
  feature_use = 100,
  draw = FALSE
)
#> Error: object 'mc' not found

table(pancreas_sub$mcRigor_status)
#> Error in x[[i, drop = TRUE]]: ‘mcRigor_status’ not found in this Seurat object
#>  

CellDimPlot(
  pancreas_sub,
  group.by = "mcRigor_metacell"
)
#> Error in CellDimPlot(pancreas_sub, group.by = "mcRigor_metacell"): "mcRigor_metacell" is not in the meta.data of srt object

CellDimPlot(
  pancreas_sub,
  group.by = "mcRigor_status"
)
#> Error in CellDimPlot(pancreas_sub, group.by = "mcRigor_status"): "mcRigor_status" is not in the meta.data of srt object
```
