# Run metabolism pathway scoring

Run metabolism pathway scoring

## Usage

``` r
RunMetabolism(
  srt,
  assay = NULL,
  group.by = NULL,
  layer = "counts",
  db = c("KEGG", "REACTOME"),
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  method = c("AUCell", "GSVA", "ssGSEA", "VISION"),
  minGSSize = 10,
  maxGSSize = 500,
  assay_name = "METABOLISM",
  new_assay = TRUE,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Assay to use as expression matrix. Default is `DefaultAssay(srt)`.

- group.by:

  Name of metadata column to group cells by. If `NULL`, single-cell
  scoring. If provided, expression is averaged by group before scoring
  (cell-type level).

- layer:

  Data layer to use, usually `"counts"` for count matrix.

- db:

  Databases to use for metabolism pathways. One or both of `"KEGG"`,
  `"REACTOME"`.

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

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- method:

  Scoring method, one of `"AUCell"`, `"GSVA"`, `"ssGSEA"`, `"VISION"`.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- assay_name:

  Name of the assay to store metabolism scores when `new_assay = TRUE`.
  Default is `"METABOLISM"`.

- new_assay:

  Whether to create a new assay for metabolism scores when
  `group.by = NULL`. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Returns a `Seurat` object. When `group.by = NULL`, stores scores in
assay `assay_name` and tools. When `group.by` is provided, stores in
tools slot `Metabolism_<group.by>_<method>` for
[MetabolismPlot](https://mengxu98.github.io/scop/reference/MetabolismPlot.md).

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-20 09:37:35] Start standard scop workflow...
#> ℹ [2026-03-20 09:37:36] Checking a list of <Seurat>...
#> ! [2026-03-20 09:37:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:37:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:37:38] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:37:39] Use the separate HVF from `srt_list`
#> ℹ [2026-03-20 09:37:39] Number of available HVF: 2000
#> ℹ [2026-03-20 09:37:39] Finished check
#> ℹ [2026-03-20 09:37:39] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-20 09:37:40] Perform pca linear dimension reduction
#> ℹ [2026-03-20 09:37:41] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-20 09:37:41] Reorder clusters...
#> ℹ [2026-03-20 09:37:41] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-20 09:37:41] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-20 09:37:45] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-20 09:37:50] Run scop standard workflow completed
pancreas_sub <- RunMetabolism(
  pancreas_sub,
  assay = "RNA",
  layer = "counts",
  db = c("KEGG", "REACTOME"),
  group.by = "CellType",
  species = "Mus_musculus",
  method = "AUCell"
)
#> ℹ [2026-03-20 09:37:50] Start metabolism pathway scoring
#> ℹ [2026-03-20 09:37:51] Data type is raw counts
#> ℹ [2026-03-20 09:37:51] Averaging expression by "CellType" ...
#> ℹ [2026-03-20 09:37:51] Aggregated expression: 15998 genes x 5 groups
#> ℹ [2026-03-20 09:37:51] Species: "Mus_musculus"
#> ℹ [2026-03-20 09:37:51] Loading cached: KEGG version: Release 117.0+/03-20, Mar 26 nterm:366 created: 2026-03-20 09:07:07
#> ℹ [2026-03-20 09:37:51] Loading cached: Reactome version: 1.95.0 nterm:1815 created: 2026-03-20 09:07:57
#> ℹ [2026-03-20 09:37:52] Total metabolism gene sets to score: 109
#> ✔ [2026-03-20 09:37:52] Metabolism scores stored in tools slot "Metabolism_CellType_AUCell"
ht <- MetabolismPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "heatmap",
  topTerm = 10,
  width = 1,
  height = 2
)
```
