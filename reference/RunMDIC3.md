# Run MDIC3 cell-cell communication inference

Run MDIC3 cell-cell communication inference

## Usage

``` r
RunMDIC3(object, ...)

# S3 method for class 'Seurat'
RunMDIC3(
  object,
  group.by,
  grn = NULL,
  grn_method = c("grnboost2", "gniplr", "correlation"),
  grn_backend = c("cpp", "python"),
  assay = NULL,
  layer = "data",
  backend = c("cpp", "python"),
  envname = "scop_env",
  conda = "auto",
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
RunMDIC3(object, ...)

# Default S3 method
RunMDIC3(
  object,
  labels,
  grn = NULL,
  grn_method = c("grnboost2", "gniplr", "correlation"),
  grn_backend = c("cpp", "python"),
  genes_in = c("rows", "columns"),
  backend = c("cpp", "python"),
  envname = "scop_env",
  conda = "auto",
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object or a gene-by-cell expression matrix.

- ...:

  Additional arguments passed to \[RunGRNBoost2()\] or \[RunGNIPLR()\]
  when those methods are selected through \`grn_method\`.

- group.by:

  Metadata column used as cell labels when \`object\` is a Seurat
  object.

- grn:

  Gene-by-gene GRN matrix aligned to expression genes, or a data frame
  with columns \`TF\`, \`target\`, and \`importance\`.

- grn_method:

  GRN inference method used when \`grn = NULL\`. \`"grnboost2"\` calls
  \[RunGRNBoost2()\]; \`"gniplr"\` calls \[RunGNIPLR()\];
  \`"correlation"\` uses absolute expression correlation as a
  lightweight GRN approximation.

- grn_backend:

  Backend passed to \[RunGRNBoost2()\] when \`grn_method =
  "grnboost2"\`.

- assay:

  Assay used when \`object\` is a Seurat object.

- layer:

  Assay layer used when \`object\` is a Seurat object.

- backend:

  Runtime backend for the MDIC3 scoring step. Supports \`"cpp"\` and
  \`"python"\`.

- envname:

  Python environment used when \`backend = "python"\`.

- conda:

  Conda-compatible executable used when \`backend = "python"\`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- labels:

  Cell labels used when \`object\` is a matrix.

- genes_in:

  Matrix orientation for matrix inputs. \`"rows"\` means genes x cells;
  \`"columns"\` means cells x genes.

## Value

A Seurat object with results stored in \`object@tools\[\["MDIC3"\]\]\`,
or a list for matrix input.

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
labels <- as.character(pancreas_sub$SubCellType[seq_len(8)])

mdic3 <- RunMDIC3(
  expr,
  labels = labels,
  grn_method = "gniplr",
  correlation_threshold = 0,
  lasso_degree = 1,
  max_lag = 1,
  max_edges_per_target = 2
)
#> ℹ [2026-09-06 22:28:17] Running GNIPLR with `backend = cpp` on 5 genes and 8 cells
#> ℹ [2026-09-06 22:28:17] Running MDIC3 with `backend = cpp` and `grn_method = gniplr` on 5 genes and 8 cells
mdic3$celltype_communication
#>                Ductal Ngn3-high-EP      Beta Ngn3-low-EP
#> Ductal       0.000000    0.0000000 1.7308490   0.0000000
#> Ngn3-high-EP 0.000000    0.0000000 0.5554415   0.0000000
#> Beta         1.153515    0.9555135 2.2147613   0.6659227
#> Ngn3-low-EP  0.000000    0.0000000 0.9453457   0.0000000
```
