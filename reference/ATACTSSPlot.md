# Plot ATAC TSS enrichment profiles

Plot previously computed per-position TSS enrichment with
Signac::TSSPlot. A scalar TSS enrichment score alone is insufficient.

## Usage

``` r
ATACTSSPlot(object, assay = NULL, group.by = NULL, idents = NULL,
  cells = NULL)
```

## Arguments

- object:

  A Seurat object containing a ChromatinAssay.

- assay:

  Name of the ChromatinAssay; defaults to the active assay.

- group.by:

  One metadata grouping column, or NULL for identities.

- idents:

  Group labels to retain, or NULL for all groups.

- cells:

  Cell names in the selected assay; NULL uses all assay cells.

## Value

A ggplot object. The input object is not modified.

## Details

First calculate positional enrichment with
`Signac::TSSEnrichment(..., fast = FALSE)` or
`RunATACQC(..., tss.positions = tss, fast = FALSE)`. Explicit TSS
positions force recalculation by RunATACQC even when imported scores
exist. This function never silently recomputes enrichment or invents a
curve from scalar scores. Only selected-assay cells are included. Signac
1.17 marks TSSPlot deprecated; the upstream lifecycle warning is
retained.

## See also

[`RunATACQC`](https://mengxu98.github.io/scop/reference/RunATACQC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
object <- RunATACQC(object, tss.positions = tss, fast = FALSE)
ATACTSSPlot(object, assay = "peaks", group.by = "celltype")
} # }
```
