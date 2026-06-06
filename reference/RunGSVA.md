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
  kcdf = c("Gaussian", "Poisson"),
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
  `kcdf = "Gaussian"` or `kcdf = "Poisson"`. PLAGE scores are oriented
  to have non-negative dot product with the gene set mean z-score so SVD
  signs are deterministic.

- cpp_chunk_size:

  Optional cell chunk size for C++ GSVA kernels. `NULL` or `"auto"`
  automatically chunks large matrices to reduce peak dense intermediate
  memory; positive values set the chunk size manually.

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

pancreas_sub <- RunGSVA(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)
ht <- GSVAPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "heatmap",
  topTerm = 10,
  width = 1,
  height = 2
)

features_all <- rownames(pancreas_sub)
pancreas_sub <- RunGSVA(
  pancreas_sub,
  features = list(
    A = features_all[1:20],
    B = features_all[21:40]
  ),
  method = c("gsva", "ssgsea")
)
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
```
