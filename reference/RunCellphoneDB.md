# Run CellphoneDB analysis

Run CellphoneDB analysis

## Usage

``` r
RunCellphoneDB(
  srt,
  group.by,
  species = c("Homo_sapiens", "Mus_musculus"),
  assay = NULL,
  layer = "data",
  counts_data = c("hgnc_symbol", "ensembl", "gene_name"),
  method = c("statistical_analysis", "analysis", "degs_analysis"),
  cpdb_file_path = NULL,
  convert_to_human = TRUE,
  gene_id_from_IDtype = "symbol",
  microenvs = NULL,
  active_tfs = NULL,
  degs = NULL,
  score_interactions = TRUE,
  threshold = 0.1,
  pvalue = 0.05,
  iterations = 1000,
  cores = 1,
  separator = "|",
  debug = FALSE,
  debug_seed = 42,
  result_precision = 3,
  subsampling = FALSE,
  subsampling_log = FALSE,
  subsampling_num_pc = 100,
  subsampling_num_cells = 1000,
  output_path = NULL,
  output_suffix = NULL,
  keep_output = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- species:

  Species of the input data. CellphoneDB is human-centric; when a
  non-human species is supplied and `convert_to_human = TRUE`, gene
  symbols will be converted to human ortholog symbols before running the
  analysis.

- assay:

  Assay to use.

- layer:

  Layer to use as CellphoneDB input. Default is `"data"` because
  CellphoneDB generally expects normalized expression.

- counts_data:

  Gene identifier type used by CellphoneDB. One of `"hgnc_symbol"`,
  `"ensembl"`, or `"gene_name"`.

- method:

  CellphoneDB analysis mode. One of `"statistical_analysis"`,
  `"analysis"`, or `"degs_analysis"`.

- cpdb_file_path:

  Path to `cellphonedb.zip`. If `NULL`, the wrapper will try common
  locations and download it automatically when not found.

- convert_to_human:

  Whether to convert non-human gene symbols to human ortholog symbols
  before running CellphoneDB.

- gene_id_from_IDtype:

  Input gene identifier type passed to
  [GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)
  when `convert_to_human = TRUE`.

- microenvs:

  Optional CellphoneDB microenvironment specification. Can be a
  `data.frame` or a file path.

- active_tfs:

  Optional CellphoneDB active TF specification. Can be a `data.frame` or
  a file path.

- degs:

  Optional DEG table for `method = "degs_analysis"`. Can be a
  `data.frame` or a file path.

- score_interactions:

  Whether to compute CellphoneDB interaction scores.

- threshold:

  Minimum fraction of cells expressing a gene for the gene to be
  considered.

- pvalue:

  P-value threshold used by CellphoneDB.

- iterations:

  Number of permutations for statistical analysis.

- cores:

  Number of cores used by CellphoneDB.

- separator:

  Separator used by CellphoneDB for sender/receiver pair columns.

- debug:

  Whether to ask CellphoneDB to save intermediate debug tables.

- debug_seed:

  Random seed used by CellphoneDB.

- result_precision:

  Result precision used by CellphoneDB.

- subsampling:

  Whether to enable CellphoneDB subsampling.

- subsampling_log:

  Whether to log-transform during CellphoneDB subsampling.

- subsampling_num_pc:

  Number of PCs used during CellphoneDB subsampling.

- subsampling_num_cells:

  Number of cells retained during CellphoneDB subsampling.

- output_path:

  Optional directory to keep CellphoneDB output files.

- output_suffix:

  Optional output suffix passed to CellphoneDB.

- keep_output:

  Whether to keep temporary output files when `output_path` is not
  supplied.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with results stored in `srt@tools[["CellphoneDB"]]`.
