# Run copy-number alteration inference

Run expression-based single-cell or spatial copy-number alteration
backends and store the results in a unified SCOP schema.

## Usage

``` r
RunCNV(
  srt,
  method = c("copykat", "fastCNV", "scevan", "infercnv", "numbat", "casper"),
  assay = NULL,
  layer = "counts",
  group.by = NULL,
  reference.by = NULL,
  reference = NULL,
  genome = c("hg38", "hg19", "mm10"),
  gene_order = NULL,
  sample.by = NULL,
  allele_counts = NULL,
  reference_counts = NULL,
  loh = NULL,
  loh_name_mapping = NULL,
  cytoband = NULL,
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

  CNA/CNV backend. Supported backends are `"copykat"`, `"fastCNV"`,
  `"scevan"`, `"infercnv"`, `"numbat"`, and `"casper"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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
  chromosome, start, and end columns. If `NULL`, SCOP tries to resolve
  these columns from assay feature metadata.

- sample.by:

  Optional sample metadata column.

- allele_counts:

  Allele count table for `"numbat"`. This is forwarded to
  `numbat::run_numbat()` as `df_allele`.

- reference_counts:

  Reference expression profile for `"numbat"`. This is forwarded to
  `numbat::run_numbat()` as `lambdas_ref`.

- loh:

  B-allele frequency/LOH signal for `"casper"`.

- loh_name_mapping:

  Optional CaSpER LOH-to-cell mapping table.

- cytoband:

  Cytoband table for `"casper"`.

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

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters forwarded to the selected backend.

## Value

A `Seurat` object with CNV metadata columns and a result bundle in
`srt@tools[[tool_name]]`.

## See also

[CNVPlot](https://mengxu98.github.io/scop/reference/CNVPlot.md)

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

# Numbat and CaSpER require allele-aware inputs from matched DNA/allele
# preprocessing workflows.
srt <- RunCNV(
  srt,
  method = "numbat",
  allele_counts = df_allele,
  reference_counts = lambdas_ref,
  genome = "hg38"
)
srt <- RunCNV(
  srt,
  method = "casper",
  reference.by = "celltype",
  reference = "Normal",
  gene_order = gene_order,
  loh = baf_signal,
  cytoband = cytoband_hg38
)

CNVPlot(srt, plot_type = "heatmap", group.by = "CNV_prediction")
CNVPlot(srt, plot_type = "dim", value = "CNV_prediction")
} # }
```
