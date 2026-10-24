# Analyse cytokine responses or immune-cell polarization

Experimental reconstruction of published IREA score, hypergeometric and
cosine-projection analyses. Numerical equivalence to the portal has not
been established for all modes.

## Usage

``` r
RunIREA(object, ...)

# Default S3 method
RunIREA(
  object,
  reference,
  contrast = NULL,
  analysis = c("cytokine_response", "cell_polarization"),
  method = c("score", "hypergeometric"),
  gene_diff_cutoff = 0.25,
  fdr_scope = c("all_contrasts", "selected_contrast"),
  ...
)

# S3 method for class 'Seurat'
RunIREA(
  object,
  reference,
  group.by,
  case,
  control,
  assay = "RNA",
  layer = "data",
  ...
)
```

## Arguments

- object:

  Character gene vector, named numeric gene contrast vector, numeric
  one-column matrix (genes in rows), data frame with genes in its first
  column and numeric contrast columns, path to a .txt/.csv/.xls/.xlsx
  table, or a Seurat object. Positive contrasts mean higher expression
  in the case. Factors and text contrast columns are rejected.
  Missing/empty gene names and nonfinite contrasts are omitted;
  duplicated contrast gene names are rejected.

- ...:

  Arguments passed to the dispatched method.

- reference:

  An `irea_reference` from
  [`PrepareDB()`](https://mengxu98.github.io/scop/reference/PrepareDB.md)
  or
  [`PrepareIREA()`](https://mengxu98.github.io/scop/reference/PrepareIREA.md).

- contrast:

  Column name to analyse in a table with multiple contrasts.

- analysis:

  `"cytokine_response"` or `"cell_polarization"`.

- method:

  `"score"` or `"hypergeometric"` for gene lists; numeric contrasts and
  Seurat inputs use cosine projection.

- gene_diff_cutoff:

  Nonnegative minimum mean reference-gene expression for projection.
  Genes must exceed this value to be retained.

- fdr_scope:

  `"all_contrasts"` adjusts BH jointly across reference terms and all
  supplied table columns before returning the selected contrast.
  `"selected_contrast"` adjusts its terms only. Vectors and Seurat
  inputs have one contrast. This option does not change gene-list
  calculations.

- group.by, case, control:

  Metadata column and distinct group names for Seurat.

- assay, layer:

  Seurat assay and exact expression layer name. Split layers must be
  joined first to compare cells across those layers.

## Value

An `irea_result` list, or the Seurat object with that result stored in
`object@tools$IREA`. Other tools are preserved. The result contains:

- `table`: one row per reference term, sorted by term. `effect` is the
  mean target-minus-PBS summed expression (score) or cosine projection
  difference, or the fraction of query genes in a significant signature
  (hypergeometric). `p_value` is a two-sided asymptotic Wilcoxon P value
  with continuity/tie correction, or an upper-tail hypergeometric P
  value; `fdr` is BH-adjusted. Empty target/control samples yield NA
  P/FDR. `status` uses FDR \< 0.05. Score/projection include `n_target`
  and `n_control`; hypergeometric includes `overlap`, `set_size`,
  `query_size` and `universe_size`. Polarization adds `radar_score`
  between 0 and 1: zero if no positive effect has FDR \< 0.05, otherwise
  each positive effect divided by the largest positive effect, with
  negative effects set to zero.

- `matched_genes`, `input_genes`: gene identifiers in input order after
  matching/filtering and before mapping, respectively.

- `gene_mapping`: human/mouse signature pairs for human input, otherwise
  NULL.

- `parameters`: analysis choices, exact selected contrast and FDR
  family; Seurat results additionally record grouping, case/control and
  assay/layer.

- `reference`: source paths, corresponding checksums and provenance.

- `validation`: experimental numerical-validation status.

## Details

Score and projection compare target reference cells against PBS-treated
reference cells, including polarization. Their P values describe
reference-cell distributions, not biological-replicate inference about
the user's samples. Reference RNA data and the selected input layer must
already be normalized. Human symbols are mapped to mouse reference genes
using the signature tables; ambiguous and unmapped symbols are omitted.

## References

Cui, Ang; Huang, Teddy; Li, Shuqiang; Ma, Aileen; Perez, Jorge L.;
Sander, Chris; Keskin, Derin B.; Wu, Catherine J.; Fraenkel, Ernest;
Hacohen, Nir. Dictionary of immune responses to cytokines at single-cell
resolution. Nature 625, 377-384 (2024). doi:10.1038/s41586-023-06816-9.

## See also

[IREAPlot](https://mengxu98.github.io/scop/reference/IREAPlot.md),
[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
reference <- PrepareDB(
  db = "IREA_Macrophage", species = "Homo_sapiens"
)[["Homo_sapiens"]][["IREA_Macrophage"]]
data(panc8_sub)
counts <- SeuratObject::LayerData(panc8_sub, assay = "RNA", layer = "counts")
macrophage <- which(panc8_sub$celltype == "macrophage")
genes <- intersect(
  c("ISG15", "IFIT3", "BST2"),
  rownames(counts)[Matrix::rowSums(counts[, macrophage, drop = FALSE]) > 0]
)
result <- RunIREA(genes, reference = reference)
IREAPlot(result)
} # }
```
