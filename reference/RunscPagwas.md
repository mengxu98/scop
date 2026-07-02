# Run scPagwas/scPaGWAS

Run the optional `scPagwas` package from `scop` without bundling LD,
pathway, or block-annotation resources. The wrapper validates required
GWAS columns, normalizes output paths, and records provenance in Seurat
tools or a result attribute.

## Usage

``` r
RunscPagwas(
  srt = NULL,
  single_data = NULL,
  gwas_data,
  celltype_meta = NULL,
  block_annotation = c("hg38", "hg37", "custom"),
  output.dirs = tempdir(),
  cleanup_soar = TRUE,
  return_seurat = !is.null(srt) || inherits(single_data, "Seurat"),
  verbose = TRUE,
  ...
)

RunscPaGWAS(
  srt = NULL,
  single_data = NULL,
  gwas_data,
  celltype_meta = NULL,
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

  Optional Seurat object or path to an `.rds` file used by `scPagwas`.

- gwas_data:

  GWAS summary statistics as a data frame.

- celltype_meta:

  Optional Seurat metadata column used to set identities.

- block_annotation:

  Genome build for bundled upstream annotations (`"hg38"` or `"hg37"`)
  or a custom annotation path.

- output.dirs:

  Output directory passed to `scPagwas`.

- cleanup_soar:

  Whether to remove SOAR objects after the run when the optional SOAR
  package is available.

- return_seurat:

  Whether to return a Seurat object when one is available.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the upstream `scPagwas` function after
  filtering by its formal arguments.

## Value

A Seurat object or upstream result list.
