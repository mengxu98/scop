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
  store_metadata = NULL,
  db = "GO_BP",
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  method = c("gsva", "ssgsea", "zscore", "plage"),
  backend = c("cpp", "r"),
  cpp_chunk_size = NULL,
  kcdf = c("Gaussian", "Poisson", "none"),
  abs.ranking = FALSE,
  min.sz = 10,
  max.sz = Inf,
  mx.diff = TRUE,
  tau = 1,
  ssgsea.norm = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object or `SummarizedExperiment` object containing the
  results of differential expression analysis
  ([`RunDEtest()`](https://mengxu98.github.io/scop/reference/RunDEtest.md)).
  If specified, the genes and groups will be extracted from the object
  automatically. If not specified, the `geneID` and `geneID_groups`
  arguments must be provided.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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

- store_metadata:

  Whether to also store single-cell GSVA scores in `meta.data`. When
  `NULL`, custom `features` or `TERM2GENE` input is stored in
  `meta.data` by default, while database-derived results stay assay-only
  when `new_assay = TRUE`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"`.
  MSigDB subcollections can be requested as `"MSigDB_<collection>"`,
  such as `"MSigDB_H"` for human Hallmark and `"MSigDB_MH"` for mouse
  Hallmark. Note: `"CytoTRACE2"` is species-independent and downloads
  pre-trained model data required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

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

- features:

  A named list of feature lists for custom enrichment gene sets. If
  provided, it takes precedence over `TERM2GENE` and `db`.

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
  `"zscore"`, or `"plage"`. Multiple methods can be supplied at once; in
  single-cell mode they will be stored in method-suffixed assays such as
  `"GSVA_gsva"` and `"GSVA_ssgsea"`. Default is `"gsva"`.

- backend:

  Scoring backend. `"cpp"` is the default and supports all current
  `method` values. `"r"` uses the original
  [`GSVA::gsva()`](https://rdrr.io/pkg/GSVA/man/gsva.html)
  implementation. `"cpp"` supports `method = "ssgsea"`,
  `method = "zscore"`, `method = "plage"`, and `method = "gsva"` with
  `kcdf = "Gaussian"`, `kcdf = "Poisson"`, or `kcdf = "none"`. PLAGE
  scores are oriented to have non-negative dot product with the gene set
  mean z-score so SVD signs are deterministic.

- cpp_chunk_size:

  Optional cell chunk size for C++ GSVA kernels. `NULL` or `"auto"`
  automatically chunks large matrices to reduce peak dense intermediate
  memory; positive values set the chunk size manually.

- kcdf:

  The kernel cumulative distribution function used for GSVA. Options are
  `"Gaussian"` (for continuous data), `"Poisson"` (for count data), or
  `"none"` (skip kernel estimation and use ranks directly). When
  omitted, `backend = "cpp"` with `method = "gsva"` uses `"none"` for
  faster single-cell scoring; explicit `"Gaussian"` or `"Poisson"`
  values are still honored. Other backends and methods default to
  `"Gaussian"`.

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

- ...:

  Passed to other functions.

## Value

Returns the modified `Seurat` object. When `group.by` is provided, GSVA
scores are stored in the `tools` slot. When `group.by = NULL`, scores
are stored in the `tools` slot, optionally in a new assay, and
optionally in `meta.data` for direct use with
[`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
and
[`FeatureStatPlot()`](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md).

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-28 10:29:18] Start standard processing workflow...
#> ℹ [2026-06-28 10:29:19] Checking a list of <Seurat>...
#> ! [2026-06-28 10:29:19] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-28 10:29:19] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 10:29:19] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 10:29:19] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 10:29:20] Number of available HVF: 2000
#> ℹ [2026-06-28 10:29:20] Finished check
#> ℹ [2026-06-28 10:29:20] Perform `ScaleData()`
#> ℹ [2026-06-28 10:29:20] Perform pca linear dimension reduction
#> ℹ [2026-06-28 10:29:20] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-28 10:29:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-28 10:29:21] Reorder clusters...
#> ℹ [2026-06-28 10:29:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-28 10:29:21] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-28 10:29:28] Standard processing workflow completed

pancreas_sub <- RunGSVA(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)
#> ℹ [2026-06-28 10:29:28] Start GSVA analysis
#> ℹ [2026-06-28 10:29:28] Start GSVA analysis
#> ℹ [2026-06-28 10:29:28] Species: "Mus_musculus"
#> ℹ [2026-06-28 10:29:28] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-28 09:42:53
#> ℹ [2026-06-28 10:29:29] Averaging expression by "CellType" ...
#> ℹ [2026-06-28 10:29:29] Aggregated expression matrix: 15998 genes x 5 groups
#> ℹ [2026-06-28 10:29:29] Processing database: "GO_BP" ...
#> ℹ [2026-06-28 10:29:30] Initial overlap: 11277 genes out of 15998 expression genes and 16594 genes in gene sets
#> ℹ [2026-06-28 10:29:30] Running GSVA for 5633 gene sets ...
#> ℹ 47266 nonzeros (less than 2^31) and 16.17% sparsity
#> ℹ [2026-06-28 10:29:34] GSVA results stored in `tools` slot: "GSVA_CellType_gsva"
#> ✔ [2026-06-28 10:29:34] GSVA analysis done
#> ℹ [2026-06-28 10:29:34] Start GSVA analysis
#> ℹ [2026-06-28 10:29:34] Species: "Mus_musculus"
#> ℹ [2026-06-28 10:29:34] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-28 09:42:53
#> ℹ [2026-06-28 10:29:36] Averaging expression by "CellType" ...
#> ℹ [2026-06-28 10:29:36] Aggregated expression matrix: 15998 genes x 5 groups
#> ℹ [2026-06-28 10:29:36] Processing database: "GO_BP" ...
#> ℹ [2026-06-28 10:29:37] Initial overlap: 11277 genes out of 15998 expression genes and 16594 genes in gene sets
#> ℹ [2026-06-28 10:29:37] Running GSVA for 5633 gene sets ...
#> ℹ [2026-06-28 10:29:39] GSVA results stored in `tools` slot: "GSVA_CellType_ssgsea"
#> ✔ [2026-06-28 10:29:39] GSVA analysis done
#> ℹ [2026-06-28 10:29:39] Start GSVA analysis
#> ℹ [2026-06-28 10:29:39] Species: "Mus_musculus"
#> ℹ [2026-06-28 10:29:39] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-28 09:42:53
#> ℹ [2026-06-28 10:29:40] Averaging expression by "CellType" ...
#> ℹ [2026-06-28 10:29:40] Aggregated expression matrix: 15998 genes x 5 groups
#> ℹ [2026-06-28 10:29:40] Processing database: "GO_BP" ...
#> ℹ [2026-06-28 10:29:41] Initial overlap: 11277 genes out of 15998 expression genes and 16594 genes in gene sets
#> ℹ [2026-06-28 10:29:42] Running GSVA for 5633 gene sets ...
#> ℹ [2026-06-28 10:29:43] GSVA results stored in `tools` slot: "GSVA_CellType_zscore"
#> ✔ [2026-06-28 10:29:43] GSVA analysis done
#> ℹ [2026-06-28 10:29:43] Start GSVA analysis
#> ℹ [2026-06-28 10:29:43] Species: "Mus_musculus"
#> ℹ [2026-06-28 10:29:43] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-06-28 09:42:53
#> ℹ [2026-06-28 10:29:44] Averaging expression by "CellType" ...
#> ℹ [2026-06-28 10:29:44] Aggregated expression matrix: 15998 genes x 5 groups
#> ℹ [2026-06-28 10:29:44] Processing database: "GO_BP" ...
#> ℹ [2026-06-28 10:29:46] Initial overlap: 11277 genes out of 15998 expression genes and 16594 genes in gene sets
#> ℹ [2026-06-28 10:29:46] Running GSVA for 5633 gene sets ...
#> ℹ [2026-06-28 10:29:52] GSVA results stored in `tools` slot: "GSVA_CellType_plage"
#> ✔ [2026-06-28 10:29:52] GSVA analysis done
ht <- GSVAPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "heatmap",
  topTerm = 10,
  width = 1,
  height = 2
)
#> ! [2026-06-28 10:29:52] Multiple GSVA results found for "CellType". Using "GSVA_CellType_gsva"
#> Warning: Data is of class matrix. Coercing to dgCMatrix.

features_all <- rownames(pancreas_sub)
pancreas_sub <- RunGSVA(
  pancreas_sub,
  features = list(
    A = features_all[1:20],
    B = features_all[21:40]
  ),
  method = c("gsva", "ssgsea")
)
#> ℹ [2026-06-28 10:29:52] Start GSVA analysis
#> ℹ [2026-06-28 10:29:52] Start GSVA analysis
#> ℹ [2026-06-28 10:29:52] Single-cell GSVA mode: using expression matrix directly ...
#> ℹ [2026-06-28 10:29:52] Expression matrix: 15998 genes x 1000 cells
#> ℹ [2026-06-28 10:29:52] Processing database: "custom" ...
#> ℹ [2026-06-28 10:29:52] Initial overlap: 40 genes out of 15998 expression genes and 40 genes in gene sets
#> ℹ [2026-06-28 10:29:52] Running GSVA for 2 gene sets ...
#> ℹ 6830 nonzeros (less than 2^31) and 82.92% sparsity
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> ℹ [2026-06-28 10:29:53] GSVA results stored in assay "GSVA_gsva", meta.data, and tools slot "GSVA_cell_gsva"
#> ✔ [2026-06-28 10:29:53] GSVA analysis done
#> ℹ [2026-06-28 10:29:53] Start GSVA analysis
#> ℹ [2026-06-28 10:29:53] Single-cell GSVA mode: using expression matrix directly ...
#> ℹ [2026-06-28 10:29:53] Expression matrix: 15998 genes x 1000 cells
#> ℹ [2026-06-28 10:29:53] Processing database: "custom" ...
#> ℹ [2026-06-28 10:29:53] Initial overlap: 40 genes out of 15998 expression genes and 40 genes in gene sets
#> ℹ [2026-06-28 10:29:53] Running GSVA for 2 gene sets ...
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> ℹ [2026-06-28 10:29:53] GSVA results stored in assay "GSVA_ssgsea", meta.data, and tools slot "GSVA_cell_ssgsea"
#> ✔ [2026-06-28 10:29:53] GSVA analysis done
#> Warning: Key ‘gsvagsva_’ taken, using ‘gsva_’ instead
FeatureDimPlot(
  pancreas_sub,
  features = "GSVA_gsva_A",
  add_density = TRUE
)

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("GSVA_gsva_A", "GSVA_ssgsea_A"),
  group.by = "CellType",
  plot.by = "feature",
  plot_type = "violin",
  stack = TRUE,
  flip = TRUE
)
#> ℹ [2026-06-28 10:29:53] Setting `group.by` to "Features" as `plot.by` is set to "feature"
```
