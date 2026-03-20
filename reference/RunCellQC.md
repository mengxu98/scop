# Run cell-level quality control for single cell RNA-seq data.

Run cell-level quality control for single cell RNA-seq data.

## Usage

``` r
RunCellQC(
  srt,
  assay = "RNA",
  split.by = NULL,
  group.by = NULL,
  return_filtered = FALSE,
  qc_metrics = c("doublets", "decontX", "outlier", "umi", "gene", "mito", "ribo",
    "ribo_mito_ratio", "species"),
  db_method = "scDblFinder",
  db_rate = NULL,
  db_coefficient = 0.01,
  decontX_threshold = NULL,
  decontX_batch = NULL,
  decontX_background = NULL,
  decontX_background_assay = NULL,
  decontX_bg_batch = NULL,
  decontX_assay_name = "decontXcounts",
  decontX_store_assay = FALSE,
  decontX_round_counts = TRUE,
  decontX_args = list(),
  outlier_threshold = c("log10_nCount:lower:2.5", "log10_nCount:higher:5",
    "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5"),
  outlier_n = 1,
  UMI_threshold = 3000,
  gene_threshold = 1000,
  mito_threshold = 20,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  ribo_threshold = 50,
  ribo_pattern = c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$",
    "rp[sl]\\d+\\w{0,1}\\d*$"),
  ribo_gene = NULL,
  ribo_mito_ratio_range = c(1, Inf),
  species = NULL,
  species_gene_prefix = NULL,
  species_percent = 95,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- group.by:

  Group labels passed to
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  when `"decontX"` is included in `qc_metrics`. Can be `NULL`, a
  meta.data column name, or a vector aligned to cells. Default is
  `NULL`.

- return_filtered:

  Logical indicating whether to return a cell-filtered Seurat object.
  Default is `FALSE`.

- qc_metrics:

  A character vector specifying the quality control metrics to be
  applied. Available metrics are `"doublets"`, `"decontX"`, `"outlier"`,
  `"umi"`, `"gene"`, `"mito"`, `"ribo"`, `"ribo_mito_ratio"`, and
  `"species"`. Default is
  `c("doublets", "decontX", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio", "species")`.

- db_method:

  Method used for doublet-calling. Can be one of `"scDblFinder"`,
  `"Scrublet"`, `"DoubletDetection"`, `"scds_cxds"`, `"scds_bcds"`,
  `"scds_hybrid"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- db_coefficient:

  The coefficient used to calculate the doublet rate. Default is `0.01`.
  Doublet rate is calculated as `ncol(srt) / 1000 * db_coefficient`.

- decontX_threshold:

  Optional contamination threshold used to filter cells after running
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md).
  Cells with `decontX_contamination` greater than this value are marked
  as failed in `decontX_qc`. Default is `NULL`, which computes decontX
  results without filtering cells by contamination.

- decontX_batch:

  Batch labels passed to
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  when `"decontX"` is included in `qc_metrics`. Default is `NULL`.

- decontX_background:

  Optional background / empty-droplet input passed to
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  when `"decontX"` is included in `qc_metrics`. Default is `NULL`.

- decontX_background_assay:

  Assay name used when `decontX_background` is a `Seurat` object or
  `SingleCellExperiment`. Default is `NULL`.

- decontX_bg_batch:

  Batch labels for `decontX_background` passed to
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md).
  Default is `NULL`.

- decontX_assay_name:

  Name of the assay used to store decontaminated counts from
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md).
  Default is `"decontXcounts"`.

- decontX_store_assay:

  Whether to store decontaminated counts as a new assay when running
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md).
  Default is `FALSE`.

- decontX_round_counts:

  Whether to round decontaminated counts before creating the assay in
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md).
  Default is `TRUE`.

- decontX_args:

  A named list of additional advanced arguments passed to
  [`RunDecontX()`](https://mengxu98.github.io/scop/reference/RunDecontX.md)
  when `"decontX"` is included in `qc_metrics`. Explicit `decontX_*`
  parameters are preferred for common options and take precedence when
  both are supplied. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- outlier_threshold:

  A character vector specifying the outlier threshold. Default is
  `c("log10_nCount:lower:2.5", "log10_nCount:higher:5", "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5")`.

- outlier_n:

  Minimum number of outlier metrics that meet the conditions for
  determining outlier cells. Default is `1`.

- UMI_threshold:

  UMI number threshold. Cells that exceed this threshold will be
  considered as kept. Default is `3000`.

- gene_threshold:

  Gene number threshold. Cells that exceed this threshold will be
  considered as kept. Default is `1000`.

- mito_threshold:

  Percentage of UMI counts of mitochondrial genes. Cells that exceed
  this threshold will be considered as discarded. Default is `20`.

- mito_pattern:

  Regex patterns to match the mitochondrial genes. Default is
  `c("MT-", "Mt-", "mt-")`.

- mito_gene:

  A defined mitochondrial genes. If features provided, will ignore the
  `mito_pattern` matching. Default is `NULL`.

- ribo_threshold:

  Percentage of UMI counts of ribosomal genes. Cells that exceed this
  threshold will be considered as discarded. Default is `50`.

- ribo_pattern:

  Regex patterns to match the ribosomal genes. Default is
  `c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$", "rp[sl]\\d+\\w{0,1}\\d*$")`.

- ribo_gene:

  A defined ribosomal genes. If features provided, will ignore the
  `ribo_pattern` matching. Default is `NULL`.

- ribo_mito_ratio_range:

  A numeric vector specifying the range of ribosomal/mitochondrial gene
  expression ratios for ribo_mito_ratio outlier cells. Default is
  `c(1, Inf)`.

- species:

  Species used as the suffix of the QC metrics. The first is the species
  of interest. Default is `NULL`.

- species_gene_prefix:

  Species gene prefix used to calculate QC metrics for each species.
  Default is `NULL`.

- species_percent:

  Percentage of UMI counts of the first species. Cells that exceed this
  threshold will be considered as kept. Default is `95`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

Returns Seurat object with the QC results stored in the meta.data layer.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-20 09:12:37] Start standard scop workflow...
#> ℹ [2026-03-20 09:12:37] Checking a list of <Seurat>...
#> ! [2026-03-20 09:12:37] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:12:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:12:39] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:12:40] Use the separate HVF from `srt_list`
#> ℹ [2026-03-20 09:12:40] Number of available HVF: 2000
#> ℹ [2026-03-20 09:12:40] Finished check
#> ℹ [2026-03-20 09:12:41] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-20 09:12:41] Perform pca linear dimension reduction
#> ℹ [2026-03-20 09:12:42] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-20 09:12:42] Reorder clusters...
#> ℹ [2026-03-20 09:12:42] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-20 09:12:42] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-20 09:12:46] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-20 09:12:50] Run scop standard workflow completed
pancreas_sub <- RunCellQC(pancreas_sub)
#> ℹ [2026-03-20 09:12:50] Data type is raw counts
#> ℹ [2026-03-20 09:12:51] Data type is raw counts
#> ℹ [2026-03-20 09:12:51] Data type is raw counts
#> ℹ [2026-03-20 09:16:55] Data type is raw counts
#> ℹ [2026-03-20 09:19:36] Running decontX
#> ℹ [2026-03-20 09:19:48] decontX contamination (median/mean/max): 0.0136 / 0.1628 / 0.7465
#> ℹ [2026-03-20 09:19:48] decontX assay stored as decontXcounts
#> ✔ [2026-03-20 09:19:48] decontX decontamination completed
#> ℹ [2026-03-20 09:19:48] decontX contamination estimates stored; no cells filtered because `decontX_threshold` is "NULL".
#> ✔ [2026-03-20 09:19:49] ● Total cells: 1000
#> ✔                       ◉ 957 cells remained
#> ✔                       ◯ 43 cells filtered out:
#> ✔                       ◯   20 potential doublets
#> ✔                       ◯   0 high-contamination cells
#> ✔                       ◯   23 outlier cells
#> ✔                       ◯   0 low-UMI cells
#> ✔                       ◯   0 low-gene cells
#> ✔                       ◯   0 high-mito cells
#> ✔                       ◯   0 high-ribo cells
#> ✔                       ◯   0 ribo_mito_ratio outlier cells
#> ✔                       ◯   0 species-contaminated cells

CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "db_qc", "decontX_qc", "outlier_qc",
    "umi_qc", "gene_qc",
    "mito_qc", "ribo_qc",
    "ribo_mito_ratio_qc",
    "species_qc"
  ),
  plot_type = "upset",
  stat_level = "Fail"
)
#> ! [2026-03-20 09:19:49] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
```
