# Run metacell partitioning for single-cell data

Run metacell partitioning for single-cell data

## Usage

``` r
RunMetaCell(
  object,
  method = c("supercell", "seacells", "metacell"),
  assay = NULL,
  layer = "counts",
  gamma = 20,
  group.by = NULL,
  envname = NULL,
  conda = "auto",
  prefix = "Metacell",
  tool_name = "Metacell",
  verbose = TRUE,
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- method:

  Metacell construction method. One of `"supercell"`, `"seacells"`, or
  `"metacell"`.

- assay:

  Assay to use for metacell construction.

- layer:

  Assay layer used to extract the count matrix.

- gamma:

  Metacell granularity parameter. For SuperCell, larger values produce
  fewer metacells (typical range 10–50). For MetaCell, this is the K
  parameter controlling the number of metacells. For SEACells, the
  comparable parameter is passed via `...` (e.g. `n_metacells`).

- group.by:

  Optional metadata column used to build metacells within each group
  independently (e.g. by sample or cell type), preventing metacells from
  crossing group boundaries.

- envname:

  Python environment name (SEACells only). Passed to
  [`reticulate::use_condaenv()`](https://rstudio.github.io/reticulate/reference/use_python.html)
  when `method = "seacells"`.

- conda:

  Conda executable path (SEACells only).

- prefix:

  Prefix for metadata columns written to `srt`.

- tool_name:

  Name of the `srt@tools` entry.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the underlying metacell method.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A metacell-level `Seurat` object. The original single-cell Seurat is
stored in `@misc[["original_srt"]]` and the cell-to-metacell membership
vector in `@misc[["cell_membership"]]`. The returned object can be
passed directly to any scop function
([`RunStandardWorkflow()`](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md),
[`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
etc.).

## References

Baran Y, Bercovich A, Sebe-Pedros A, et al. (2019). MetaCell: analysis
of single-cell RNA-seq data using K-nn graph partitions. Genome Biology,
20(1), 206. doi:10.1186/s13059-019-1812-2

Bilous M, Tran L, Cianciaruso C, et al. (2022). SuperCell: a versatile
tool for single-cell data analysis. Genome Biology, 23(1), 227.
doi:10.1186/s13059-022-02778-5

Persad S, Choo ZN, Dien C, et al. (2023). SEACells infers
transcriptional and epigenomic cellular states from single-cell genomics
data. Nature Biotechnology, 41(11), 1658-1667.
doi:10.1038/s41587-023-01716-9

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(
  pancreas_sub,
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:20,
  nonlinear_reduction_dims = 2,
  verbose = FALSE
)
#> ℹ [2026-09-20 22:39:21] Skip `log1p()` because `layer = data` is not "counts"

mc1 <- RunMetaCell(
  pancreas_sub,
  method = "supercell",
  gamma = 20
)
#> ℹ [2026-09-20 22:39:28] Running SuperCell with gamma = 20, k.knn = 5 on 1000 cells
#> ℹ [2026-09-20 22:39:29] `RunMetaCell()` ("supercell") built 50 metacells from 1000 cells
#> ℹ [2026-09-20 22:39:29] Metacell size summary: min 5, median 16.5, mean 20, max 56 cells
#> ✔ [2026-09-20 22:39:29] `RunMetaCell()` returned metacell Seurat with 50 metacells. Original cells in `@misc[["original_srt"]]`

MetaCellPlot(mc1, group.by = "CellType")


mc2 <- RunMetaCell(
  pancreas_sub,
  method = "metacell",
  gamma = 20
)
#> ℹ [2026-09-20 22:39:31] Running MetaCell-style KNN partitioning with k = 20 on 1000 cells
#> ℹ [2026-09-20 22:39:31] `RunMetaCell()` ("metacell") built 8 metacells from 1000 cells
#> ℹ [2026-09-20 22:39:31] Metacell size summary: min 16, median 150, mean 125, max 215 cells
#> ✔ [2026-09-20 22:39:31] `RunMetaCell()` returned metacell Seurat with 8 metacells. Original cells in `@misc[["original_srt"]]`

MetaCellPlot(mc2, group.by = "CellType")
```
