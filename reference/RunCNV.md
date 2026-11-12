# Run copy-number alteration inference

Run copy-number alteration inference

## Usage

``` r
RunCNV(
  object,
  method = c("copykat", "fastCNV", "scevan", "infercnv", "numbat"),
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
  output_dir = NULL,
  prefix = "CNV",
  tool_name = "CNV",
  store_matrix = TRUE,
  verbose = TRUE,
  srt = NULL,
  ...
)
```

## Source

The example uses expression and cell annotations from
[GSE70630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70630).

## Arguments

- object:

  A `Seurat` object.

- method:

  CNA/CNV method. Supported methods are `"copykat"`, `"fastCNV"`,
  `"scevan"`, `"infercnv"`, and `"numbat"`.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer used as the expression matrix.

- group.by:

  Optional metadata column forwarded to supported methods and stored as
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
  [`numbat::run_numbat()`](https://kharchenkolab.github.io/numbat/reference/run_numbat.html)
  as `df_allele`.

- reference_counts:

  Reference expression profile for `"numbat"`. This is forwarded to
  [`numbat::run_numbat()`](https://kharchenkolab.github.io/numbat/reference/run_numbat.html)
  as `lambdas_ref`.

- output_dir:

  Optional output directory.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used for `srt@tools`.

- store_matrix:

  Whether to store the normalized CNV matrix in
  `srt@tools[[tool_name]]`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.

## Value

A `Seurat` object with CNV metadata columns and a result bundle in
`srt@tools[[tool_name]]`.

## See also

[CNVPlot](https://mengxu98.github.io/scop/reference/CNVPlot.md)

## Examples

``` r
base_url <- paste0(
  "https://raw.githubusercontent.com/broadinstitute/infercnv/",
  "1def6c930335d00ddf251474abb0a31a1c74b0da/inst/extdata/"
)
files <- c("oligodendroglioma_expression_downsampled.counts.matrix.gz",
           "oligodendroglioma_annotations_downsampled.txt")
counts_file <- tempfile(fileext = ".gz")
labels_file <- tempfile()
utils::download.file(paste0(base_url, files[1]), counts_file, mode = "wb")
utils::download.file(paste0(base_url, files[2]), labels_file, mode = "wb")
counts <- as.matrix(utils::read.delim(counts_file, row.names = 1, check.names = FALSE))
labels <- utils::read.delim(labels_file, header = FALSE, row.names = 1)
cohort <- labels[colnames(counts), 1]
colnames(counts) <- paste0("cell_", colnames(counts))
srt <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
srt$group <- factor(
  ifelse(startsWith(cohort, "malignant_"), sub("^malignant_", "Tumor ", cohort), "Reference"),
  levels = c("Reference", "Tumor 93", "Tumor 97", "Tumor MGH36", "Tumor MGH53")
)
set.seed(42)
srt <- RunCNV(srt, method = "copykat", reference.by = "group", reference = "Reference")
#> ℹ [2026-11-12 03:42:33] Run copykat CNV backend with 10338 features and 184 cells
#> ℹ [2026-11-12 03:42:33] Attached copykat to the search path so its lazy data objects are resolvable by the backend
#> [1] "running copykat v1.2.5"
#> [1] "step1: read and filter data ..."
#> [1] "10338 genes, 184 cells in raw data"
#> [1] "9450 genes past LOW.DR filtering"
#> [1] "step 2: annotations gene coordinates ..."
#> [1] "start annotation ..."
#> [1] "step 3: smoothing data with dlm ..."
#> [1] "step 4: measuring baselines ..."
#> [1] "42 known normal cells found in dataset"
#> [1] "run with known normal..."
#> [1] "baseline is from known input"
#> [1] "step 5: segmentation..."
#> [1] "step 6: convert to genomic bins..."
#> [1] "step 7: re-adjust baseline ..."
#> [1] "step 8: final prediction ..."
#> [1] "step 9: saving results..."
#> [1] "step 10: ploting heatmap ..."
#> Time difference of 1.140061 mins
#> ✔ [2026-11-12 03:43:41] copykat CNV results stored in `srt@tools[[CNV]]`
CNVPlot(srt, group.by = "group", heatmap_limits = "auto")
```
