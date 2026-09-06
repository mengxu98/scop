# Run scPagwas

Run the optional `scPagwas` package from `scop` without bundling LD,
pathway, or block-annotation resources. The wrapper validates required
GWAS columns, normalizes output paths, and records provenance in Seurat
tools or a result attribute. It prefers the upstream `scPagwas_main2`
runner and applies Seurat 5 compatibility to local function copies
without modifying the installed backend namespace.

## Usage

``` r
RunscPagwas(
  srt = NULL,
  single_data = NULL,
  gwas_data,
  group.by = NULL,
  singlecell = TRUE,
  celltype = TRUE,
  assay = NULL,
  block_annotation = c("hg38", "hg37", "custom"),
  output.dirs = tempdir(),
  cleanup_soar = TRUE,
  return_seurat = !is.null(srt) || inherits(single_data, "Seurat"),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  Optional Seurat object used as single-cell input.

- single_data:

  Optional Seurat object or path to a Seurat `.rds` file.

- gwas_data:

  GWAS summary statistics as a data frame or delimited text file.
  Required columns are `chrom`, `pos`, `rsid`, `se`, `beta`, and `maf`.

- group.by:

  Optional Seurat metadata column used to set cell identities.

- singlecell:

  Whether to calculate single-cell results.

- celltype:

  Whether to calculate cell-type results.

- assay:

  Assay used by `scPagwas`. Defaults to the active assay for a Seurat
  object and to `"RNA"` for other inputs.

- block_annotation:

  Genome build for bundled upstream annotations (`"hg38"` or `"hg37"`)
  or a custom annotation path.

- output.dirs:

  Output directory passed to `scPagwas`.

- cleanup_soar:

  Deprecated compatibility argument. SOAR cleanup is managed by the
  upstream `scPagwas` backend and is ignored by `scop`.

- return_seurat:

  Whether to return a Seurat object when one is available.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the upstream `scPagwas` function after
  filtering by its formal arguments.

## Value

A Seurat object or upstream result list.

## See also

[`PlotscPagwas()`](https://mengxu98.github.io/scop/reference/PlotscPagwas.md)
