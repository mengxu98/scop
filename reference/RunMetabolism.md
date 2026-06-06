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
  biomart = NULL,
  max_tries = 5,
  use_preparedb = TRUE,
  method = c("AUCell", "GSVA", "ssGSEA", "VISION"),
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  cpp_chunk_size = NULL,
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
  `"REACTOME"`. `"Reactome"` is also accepted and treated identically to
  `"REACTOME"`. When `use_preparedb = TRUE`, gene sets are built via
  [PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md).

- species:

  Species of the input data. The scMetabolism gene sets contain human
  gene symbols. When `species` is not `"Homo_sapiens"` and
  `convert_species` is `TRUE`,
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  is used to map human genes to the target species via biomaRt homolog
  tables. Default is `"Homo_sapiens"`.

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

  Whether to convert human gene symbols from the scMetabolism gene sets
  to the target species using
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md).
  When `TRUE` (default), genes are mapped via cross-species orthologs
  from Ensembl BioMart. When `FALSE`, only case-insensitive direct
  symbol matching is used.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- biomart:

  BioMart database name passed to
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md).
  Default `NULL` uses `"ensembl"`. Other options: `"protists_mart"`,
  `"fungi_mart"`, `"plants_mart"`.

- max_tries:

  Maximum retry attempts for biomaRt connections in
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md).
  Default is `5`.

- use_preparedb:

  When `TRUE`, gene sets are built via
  [PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)
  which provides species-aware gene mapping via BioMart and
  KEGG/Reactome databases. This automatically handles gene symbol
  conversion for non-human species (e.g., `species = "Mus_musculus"` →
  mouse gene symbols in metabolism pathways). When `FALSE`, raw
  scMetabolism GMT files are downloaded and genes are matched
  case-insensitively with optional
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  supplementation when `convert_species = TRUE`. genes and approximates
  zero ties, `"topk"` ranks only genes that can contribute to AUCell
  AUC, and `"full"` ranks all genes.

- method:

  Scoring method, one of `"AUCell"`, `"GSVA"`, `"ssGSEA"`, `"VISION"`.

- backend:

  Scoring backend. `"cpp"` is the default for supported methods. `"r"`
  uses the original R package implementation. `"cpp"` currently supports
  `method = "AUCell"`, `method = "GSVA"`, and `method = "ssGSEA"`.
  `method = "VISION"` falls back to `"r"` when `backend` is not
  explicitly set. AUCell C++ scores may differ from the R backend when
  tied expression values are randomly ranked.

- cpp_strategy:

  C++ AUCell ranking strategy. `"sparse"` ranks non-zero

- cpp_chunk_size:

  Optional cell chunk size for C++ GSVA kernels. `NULL` or `"auto"`
  automatically chunks large matrices to reduce peak dense intermediate
  memory; positive values set the chunk size manually.

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
pancreas_sub <- RunMetabolism(
  pancreas_sub,
  assay = "RNA",
  layer = "counts",
  db = c("KEGG", "REACTOME"),
  group.by = "CellType",
  species = "Mus_musculus",
  method = "AUCell"
)
ht <- MetabolismPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "heatmap",
  topTerm = 10,
  width = 1,
  height = 2
)
```
