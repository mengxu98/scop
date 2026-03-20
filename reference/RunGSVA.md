# Perform Gene Set Variation Analysis (GSVA)

Perform Gene Set Variation Analysis (GSVA)

## Usage

``` r
RunGSVA(
  srt = NULL,
  assay = NULL,
  group.by = NULL,
  layer = "data",
  assay_name = "GSVA",
  new_assay = TRUE,
  db = "GO_BP",
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  method = c("gsva", "ssgsea", "zscore", "plage"),
  kcdf = c("Gaussian", "Poisson"),
  abs.ranking = FALSE,
  min.sz = 10,
  max.sz = Inf,
  mx.diff = TRUE,
  tau = 1,
  ssgsea.norm = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing the results of differential expression
  analysis (RunDEtest). If specified, the genes and groups will be
  extracted from the Seurat object automatically. If not specified, the
  `geneID` and `geneID_groups` arguments must be provided.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- group.by:

  Name of metadata column to group cells by for averaging expression. If
  provided, expression will be averaged within each group before GSVA
  analysis (cell-type level). If `NULL`, GSVA is performed on each cell
  individually (single-cell level).

- layer:

  Data layer to use when `group.by = NULL`. Usually `"data"` for
  normalized or `"counts"` for count matrix. Default is `"data"`.

- assay_name:

  Name of the assay to store GSVA scores when `group.by = NULL` and
  `new_assay = TRUE`. Default is `"GSVA"`.

- new_assay:

  Whether to create a new assay for GSVA scores when `group.by = NULL`.
  Default is `TRUE`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF"`.

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the gene annotation
  databases to be retrieved. Default is `"latest"`.

- db_combine:

  Whether to combine multiple databases into one. If `TRUE`, all
  database specified by `db` will be combined as one named "Combined".

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- TERM2GENE:

  A data frame specifying the gene-term mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the gene IDs.

- TERM2NAME:

  A data frame specifying the term-name mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the corresponding term names.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- unlimited_db:

  A character vector specifying the names of databases that do not have
  size restrictions.

- method:

  The method to use for GSVA. Options are `"gsva"`, `"ssgsea"`,
  `"zscore"`, or `"plage"`. Default is `"gsva"`.

- kcdf:

  The kernel cumulative distribution function used for GSVA. Options are
  `"Gaussian"` (for continuous data) or `"Poisson"` (for count data).
  Default is `"Gaussian"`.

- abs.ranking:

  Logical indicating whether to use absolute ranking for GSVA. Default
  is `FALSE`.

- min.sz:

  Minimum size of gene sets to be included in the analysis. Default is
  `10`.

- max.sz:

  Maximum size of gene sets to be included in the analysis. Default is
  `Inf`.

- mx.diff:

  Logical indicating whether to use the maximum difference method.
  Default is `TRUE`.

- tau:

  Exponent for the GSVA method. Default is `1`.

- ssgsea.norm:

  Logical indicating whether to normalize SSGSEA scores. Default is
  `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Returns the modified `Seurat` object. When `group.by` is provided, GSVA
scores are stored in the `tools` slot. When `group.by = NULL`, scores
are stored in a new assay (if `new_assay = TRUE`) and in the `tools`
slot.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-20 09:33:11] Start standard scop workflow...
#> ℹ [2026-03-20 09:33:11] Checking a list of <Seurat>...
#> ! [2026-03-20 09:33:11] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:33:11] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:33:13] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:33:14] Use the separate HVF from `srt_list`
#> ℹ [2026-03-20 09:33:14] Number of available HVF: 2000
#> ℹ [2026-03-20 09:33:14] Finished check
#> ℹ [2026-03-20 09:33:15] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-20 09:33:15] Perform pca linear dimension reduction
#> ℹ [2026-03-20 09:33:16] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-20 09:33:16] Reorder clusters...
#> ℹ [2026-03-20 09:33:16] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-20 09:33:16] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-20 09:33:21] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-20 09:33:25] Run scop standard workflow completed

pancreas_sub <- RunGSVA(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)
#> ℹ [2026-03-20 09:33:25] Start GSVA analysis
#> ℹ [2026-03-20 09:33:25] Averaging expression by "CellType" ...
#> ℹ [2026-03-20 09:33:25] Aggregated expression matrix: 15998 genes x 5 groups
#> ℹ [2026-03-20 09:33:25] Species: "Mus_musculus"
#> ℹ [2026-03-20 09:33:25] Loading cached: GO_BP version: 3.22.0 nterm:15169 created: 2026-03-20 08:29:24
#> ℹ [2026-03-20 09:33:27] Processing database: "GO_BP" ...
#> ℹ [2026-03-20 09:33:28] Initial overlap: 11182 genes out of 15998 expression genes and 16088 genes in gene sets
#> ℹ [2026-03-20 09:33:31] Running GSVA for 5668 gene sets ...
#> ℹ GSVA version 2.4.8
#> ℹ Searching for rows with constant values
#> ! 2 rows with constant values throughout the columns
#> ! Rows with constant values are discarded
#> ℹ Calculating GSVA ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Row-wise ECDF estimation with Gaussian kernels
#> ℹ Calculating row ECDFs
#> ℹ Calculating column ranks
#> ℹ GSVA dense (classical) algorithm
#> ℹ Calculating GSVA scores
#> ✔ Calculations finished
#> ℹ [2026-03-20 09:34:58] GSVA results stored in `tools` slot: "GSVA_CellType_gsva"
#> ✔ [2026-03-20 09:34:58] GSVA analysis done
ht <- GSVAPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "heatmap",
  topTerm = 10,
  width = 1,
  height = 2
)
```
