# Run copy-number alteration inference

Run expression-based single-cell or spatial copy-number alteration
backends and store the results in a unified SCOP schema.

## Usage

``` r
RunCNV(
  srt,
  method = c("copykat", "fastCNV", "scevan", "infercnv"),
  assay = NULL,
  layer = "counts",
  group.by = NULL,
  reference.by = NULL,
  reference = NULL,
  genome = c("hg38", "hg19", "mm10"),
  gene_order = NULL,
  sample.by = NULL,
  output_dir = NULL,
  prefix = "CNV",
  tool_name = "CNV",
  store_matrix = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- method:

  CNA/CNV backend. First-stage expression backends are `"copykat"`,
  `"fastCNV"`, `"scevan"`, and `"infercnv"`.

- assay:

  Assay to use. Defaults to `DefaultAssay(srt)`.

- layer:

  Assay layer used as the expression matrix.

- group.by:

  Optional metadata column forwarded to supported backends and stored as
  cell annotation.

- reference.by:

  Metadata column identifying reference/normal cells. Required for
  `"infercnv"` and `"fastCNV"`.

- reference:

  Reference labels in `reference.by`.

- genome:

  Reference genome label.

- gene_order:

  Gene coordinate table or a path to one. The table should contain gene,
  chromosome, start, and end columns.

- sample.by:

  Optional sample metadata column.

- output_dir:

  Optional backend output directory.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used for `srt@tools`.

- store_matrix:

  Whether to store the normalized CNV matrix in
  `srt@tools[[tool_name]]`.

- verbose:

  Whether to print progress messages.

- ...:

  Additional parameters forwarded to the selected backend.

## Value

A `Seurat` object with CNV metadata columns and a result bundle in
`srt@tools[[tool_name]]`.

## Details

`Numbat` and `CaSpER` are not exposed in this first interface because
they require allele-aware inputs.

## See also

[`CNVPlot`](https://mengxu98.github.io/scop/reference/CNVPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# copykat uses raw counts and can infer diploid/aneuploid cells directly.
srt <- RunCNV(
  srt,
  method = "copykat",
  genome = "hg38"
)

# fastCNV and inferCNV require a normal/reference cell annotation.
srt <- RunCNV(
  srt,
  method = "fastCNV",
  reference.by = "celltype",
  reference = "Normal",
  genome = "hg38"
)

gene_order <- data.frame(
  gene = rownames(srt),
  chr = "chr1",
  start = seq_len(nrow(srt)) * 1000,
  end = seq_len(nrow(srt)) * 1000 + 999
)
srt <- RunCNV(
  srt,
  method = "infercnv",
  reference.by = "celltype",
  reference = "Normal",
  gene_order = gene_order
)

CNVPlot(srt, plot_type = "heatmap", group.by = "CNV_prediction")
CNVPlot(srt, plot_type = "dim", value = "CNV_prediction")
} # }
```
