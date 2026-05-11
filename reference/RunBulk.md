# Unified bulk analysis entrypoint

`RunBulk()` provides one bulk-strategy entrypoint for Seurat-centered
workflows. Methods are selected with one character vector, so users do
not need to choose a separate module argument.

## Usage

``` r
RunBulk(
  srt = NULL,
  bulk_se = NULL,
  method = "de_edgeR_qlf",
  sample.by = NULL,
  condition.by = NULL,
  group.by = NULL,
  ref_srt = NULL,
  celltype.by = NULL,
  assay = NULL,
  layer = "counts",
  bulk_assay = "counts",
  ref_assay = NULL,
  ref_layer = "counts",
  condition1 = NULL,
  condition2 = NULL,
  de_markers_type = "single",
  markers_type = NULL,
  run_enrichment = FALSE,
  run_gsea = FALSE,
  enrichment_args = list(),
  gsea_args = list(),
  de_args = list(),
  deconv_args = list(),
  csde_args = list(),
  method_args = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object used for pseudobulk mode.

- bulk_se:

  A `SummarizedExperiment` object used for true bulk mode.

- method:

  Character vector specifying methods to run. Supported values are
  `"de_limma_voom"`, `"de_edgeR_qlf"`, `"de_DESeq2"`, `"de_dream"`,
  `"deconv_MuSiC"`, `"deconv_BisqueRNA"`, `"deconv_BayesPrism"`, and
  `"csde_TOAST"`.

- sample.by:

  Metadata column containing sample IDs in pseudobulk mode.

- condition.by:

  Metadata column containing condition labels. Required when
  `.arg method` includes `"de"` or `"csde"`. Optional for
  deconvolution-only runs.

- group.by:

  Optional metadata column used to stratify DE by subgroup.

- ref_srt:

  Optional `Seurat` reference for deconvolution/CSDE. If omitted in
  pseudobulk mode, `srt` is used as reference.

- celltype.by:

  Metadata column in `ref_srt` defining reference cell types. Required
  when deconvolution or CSDE is requested.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Layer name in `srt` used for pseudobulk counts.

- bulk_assay:

  Assay name in `bulk_se` used as counts matrix.

- ref_assay:

  Assay name in `ref_srt` for reference profiles.

- ref_layer:

  Layer name in `ref_srt` for reference counts.

- condition1:

  First condition label for pairwise comparison.

- condition2:

  Second condition label for pairwise comparison.

- de_markers_type:

  DE comparison mode. One of `"single"` (one pairwise comparison),
  `"all"` (each condition vs the rest), or `"paired"` (all pairwise
  condition comparisons).

- markers_type:

  Alias of `.arg de_markers_type` for consistency with
  [RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md).
  If provided, it overrides `.arg de_markers_type`.

- run_enrichment:

  Whether to run
  [RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md)
  from bulk DE results.

- run_gsea:

  Whether to run
  [RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md) from
  bulk DE results.

- enrichment_args:

  Named list of additional parameters for
  [RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md).

- gsea_args:

  Named list of additional parameters for
  [RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md).

- de_args:

  Named list of additional parameters for DE method functions.

- deconv_args:

  Named list of additional parameters for deconvolution method
  functions. Current deconvolution wrappers use `backend = "internal"`
  by default and record the computational engine in the result bundle.

- csde_args:

  Named list of additional parameters for CSDE method functions. Current
  `csde_TOAST` uses `backend = "limma_interaction"` by default and
  records the computational engine in the result bundle.

- method_args:

  Named list of method-specific parameters. Names should be RunBulk
  method names, such as `"de_edgeR_qlf"` or `"csde_TOAST"`.
  Method-specific values override `.arg de_args`, `.arg deconv_args`,
  `.arg csde_args`, and `.arg ...`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters forwarded to method functions.

## Value

Returns a `Seurat` object in pseudobulk mode with results in
`srt@tools[["Bulk"]]`. If `.arg run_enrichment` or `.arg run_gsea` is
`TRUE`, successful pathway results are also stored in the standard
`srt@tools[["Enrichment_Bulk_wilcox"]]` and
`srt@tools[["GSEA_Bulk_wilcox"]]` slots used by
[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md)
and [GSEAPlot](https://mengxu98.github.io/scop/reference/GSEAPlot.md).
True bulk mode returns a `SummarizedExperiment` object with results in
`metadata(bulk_se)[["Bulk"]]`.

## Examples

``` r
# ----------------------------
# Example 1: Seurat -> pseudobulk (DE) using scop built-in dataset
# ----------------------------
if (requireNamespace("edgeR", quietly = TRUE)) {
  data("pancreas_sub", package = "scop")
  srt <- pancreas_sub

  # Build sample/condition labels for pseudobulk aggregation.
  srt$sample_id <- paste0(
    "sample_",
    (seq_len(ncol(srt)) - 1L) %% 8L + 1L
  )
  srt$condition <- ifelse(
    srt$sample_id %in% paste0("sample_", 1:4),
    "A",
    "B"
  )

  srt <- RunBulk(
    srt = srt,
    method = "de_edgeR_qlf",
    sample.by = "sample_id",
    condition.by = "condition",
    verbose = FALSE
  )
  srt@tools$Bulk$status$de$status
}
#> calcNormFactors has been renamed to normLibSizes
#> [1] "success"

# ----------------------------
# Example 2: pure bulk (SummarizedExperiment -> DE)
# ----------------------------
if (requireNamespace("edgeR", quietly = TRUE) &&
  requireNamespace("SummarizedExperiment", quietly = TRUE) &&
  requireNamespace("S4Vectors", quietly = TRUE)) {
  data("pancreas_sub", package = "scop")
  srt <- pancreas_sub
  srt$sample_id <- paste0(
    "sample_",
    (seq_len(ncol(srt)) - 1L) %% 8L + 1L
  )
  srt$condition <- ifelse(
    srt$sample_id %in% paste0("sample_", 1:4),
    "A",
    "B"
  )

  counts <- GetAssayData5(srt, layer = "counts")
  sid <- as.character(srt$sample_id)
  sid_levels <- unique(sid)
  bulk_counts <- do.call(cbind, lapply(sid_levels, function(x) {
    Matrix::rowSums(counts[, sid == x, drop = FALSE])
  }))
  colnames(bulk_counts) <- sid_levels
  rownames(bulk_counts) <- rownames(counts)

  sample_condition <- tapply(
    as.character(srt$condition),
    sid,
    function(x) unique(x)[1]
  )
  bulk_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(bulk_counts)),
    colData = S4Vectors::DataFrame(
      condition = as.character(sample_condition[colnames(bulk_counts)]),
      row.names = colnames(bulk_counts)
    )
  )

  bulk_se <- RunBulk(
    bulk_se = bulk_se,
    method = "de_edgeR_qlf",
    condition.by = "condition",
    verbose = FALSE
  )
  S4Vectors::metadata(bulk_se)$Bulk$status$de$status
}
#> calcNormFactors has been renamed to normLibSizes
#> [1] "success"

# ----------------------------
# Example 3: deconvolution + CSDE in one call
# ----------------------------
if (requireNamespace("limma", quietly = TRUE)) {
  data("pancreas_sub", package = "scop")
  srt <- pancreas_sub
  srt$sample_id <- paste0(
    "sample_",
    (seq_len(ncol(srt)) - 1L) %% 8L + 1L
  )
  srt$condition <- ifelse(
    srt$sample_id %in% paste0("sample_", 1:4),
    "A",
    "B"
  )
  if ("CellType" %in% colnames(srt@meta.data)) {
    srt$celltype <- as.character(srt$CellType)
  } else {
    srt$celltype <- as.character(SeuratObject::Idents(srt))
  }

  srt <- RunBulk(
    srt = srt,
    method = c("deconv_MuSiC", "csde_TOAST"),
    sample.by = "sample_id",
    condition.by = "condition",
    condition1 = "A",
    condition2 = "B",
    celltype.by = "celltype",
    verbose = FALSE
  )
  srt@tools$Bulk$status$deconv$status
  srt@tools$Bulk$status$csde$status
}
#> [1] "success"
```
