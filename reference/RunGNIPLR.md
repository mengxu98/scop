# Infer gene regulatory networks with GNIPLR

Infer gene regulatory networks with GNIPLR

## Usage

``` r
RunGNIPLR(object, ...)

# S3 method for class 'Seurat'
RunGNIPLR(
  object,
  assay = NULL,
  layer = "counts",
  targets = NULL,
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "gniplr",
  envname = "scop_env",
  conda = "auto",
  force = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
RunGNIPLR(object, ...)

# Default S3 method
RunGNIPLR(
  object,
  targets = NULL,
  genes_in = c("rows", "columns"),
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "gniplr",
  envname = "scop_env",
  conda = "auto",
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object or expression matrix.

- ...:

  Additional backend-specific arguments.

- assay:

  Assay used when \`object\` is a Seurat object.

- layer:

  Assay layer used when \`object\` is a Seurat object.

- targets:

  Optional target genes. If \`NULL\`, all genes are considered.

- correlation_threshold:

  Relative correlation filter used before lagged regression. A gene pair
  is tested when its absolute correlation is at least this fraction of
  the maximum absolute correlation for the regulator.

- lasso_degree:

  Polynomial degree used by the LASSO projection step.

- lasso_alpha:

  LASSO regularization strength used by the projection step.

- max_lag:

  Maximum lag used by Granger-style lagged regression. Values above
  \`3\` are capped at \`3\` to match the GNIPLR reference
  implementation.

- backend:

  Runtime backend. Supports \`"cpp"\` and \`"python"\`.

- max_edges_per_target:

  Maximum incoming regulator edges retained per target. The default
  \`Inf\` keeps all positive-importance links.

- output_file:

  Optional path where the adjacency table is written.

- work_dir:

  Working directory used by the Python backend.

- prefix:

  Prefix for temporary files.

- envname:

  Python environment used by the Python backend.

- conda:

  Conda-compatible executable used by the Python backend.

- force:

  Whether to rebuild existing \`output_file\`.

- verbose:

  Whether to print progress messages.

- genes_in:

  Matrix orientation for matrix inputs. \`"rows"\` means genes x cells;
  \`"columns"\` means cells x genes.

## Value

A data frame with columns \`TF\`, \`target\`, \`importance\`, and
\`pvalue\`. The original GNIPLR p-value matrix is stored in
\`attr(result, "grn_matrix")\` when the network is newly inferred.

## References

Zhang Y, Chang X, Liu X. Inference of gene regulatory networks using
pseudo-time series data. Bioinformatics. 2021;37(16):2423-2431.
doi:10.1093/bioinformatics/btab099.

## Examples

``` r
data(pancreas_sub)
expr <- GetAssayData5(
  pancreas_sub,
  assay = SeuratObject::DefaultAssay(pancreas_sub),
  layer = "counts"
)
expr <- as.matrix(expr[, seq_len(8)])
expr <- expr[
  names(sort(apply(expr, 1, stats::var), decreasing = TRUE))[seq_len(5)],
]

grn <- RunGNIPLR(
  expr,
  genes_in = "rows",
  correlation_threshold = 0,
  lasso_degree = 1,
  max_lag = 1,
  max_edges_per_target = 2
)
#> ℹ [2026-09-06 22:09:00] Running GNIPLR with `backend = cpp` on 5 genes and 8 cells
head(grn)
#>     TF target importance      pvalue
#> 3 Iapp   Nnat  2.2454349 0.005682836
#> 9 Rbp4   Spp1  1.6205718 0.023956766
#> 4  Pyy   Nnat  1.1742331 0.066952521
#> 7 Iapp   Rbp4  1.0851237 0.082200850
#> 8 Nnat   Rbp4  1.0308208 0.093149207
#> 5 Nnat    Pyy  0.9471533 0.112939714
attr(grn, "grn_matrix")[1:3, 1:3]
#>           Iapp       Pyy        Nnat
#> Iapp 0.0000000 0.9390915 0.005682836
#> Pyy  0.8998447 0.0000000 0.066952521
#> Nnat 0.3772449 0.1129397 0.000000000
```
