# Run SCENIC gene regulatory network analysis

Run SCENIC gene regulatory network analysis

## Usage

``` r
RunSCENIC(
  srt,
  assay = NULL,
  layer = "counts",
  ranking_dbs = NULL,
  motif_annotations = NULL,
  regulators = NULL,
  targets = NULL,
  work_dir = "scenic_output/",
  species = c("Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"),
  genome = NULL,
  data_dir = NULL,
  prefix = "scenic",
  min_expr_cells = 3,
  min_regulon_size = 10,
  backend = c("cpp", "python"),
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  cores = 1,
  seed = 1234,
  force = FALSE,
  assay_name = "scenic",
  tool_name = "SCENIC",
  return_seurat = TRUE,
  envname = NULL,
  conda = "auto",
  prepare_env = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the count matrix.

- ranking_dbs:

  Character vector of cisTarget ranking feather files. If `NULL`, the
  gene-based v10 cisTarget ranking databases are prepared from
  `species`.

- motif_annotations:

  Motif annotation table used by `scenic ctx`. If `NULL`, the v10
  motif2tf table is prepared from `species`.

- regulators:

  Transcription factors used as candidate regulators in GRNBoost2. This
  can be a character vector of gene names or one text file. If `NULL`,
  the cisTarget TF list is prepared from `species`.

- targets:

  Optional target genes used to restrict the GRN output. This can be a
  character vector of gene names or one text file. Regulator expression
  is still kept as predictor input for GRNBoost2, and the adjacency
  table passed to `scenic ctx` is filtered to these targets.

- work_dir:

  Directory used for SCENIC input and output files.

- species:

  Species used to select cisTarget reference files when `ranking_dbs`,
  `motif_annotations`, or `regulators` is `NULL`. Supported values are
  `"Homo_sapiens"`, `"Mus_musculus"`, and `"Drosophila_melanogaster"`.

- genome:

  Genome build used to select cisTarget reference files when automatic
  references are prepared. Human supports `"hg38"` (default) and
  `"hg19"`. Mouse and fly currently use `"mm10"` and `"dm6"`,
  respectively.

- data_dir:

  Directory used to cache automatically prepared SCENIC reference files.
  If `NULL`, files are stored under
  `tools::R_user_dir("scop", "data")/SCENIC/<species>`.

- prefix:

  Prefix for SCENIC output files.

- min_expr_cells:

  Minimum number of cells where a gene must be detected before
  GRNBoost2. To run SCENIC on metacells, first create a metacell-level
  object with
  [`RunMetaCell()`](https://mengxu98.github.io/scop/reference/RunMetaCell.md)
  and pass that object to `RunSCENIC()`.

- min_regulon_size:

  Minimum regulon size kept after `scenic ctx`.

- backend:

  SCENIC backend. `"cpp"` uses the R/C++ path and `"python"` uses the
  Python pySCENIC path. The selected backend controls GRN, cisTarget
  pruning, and AUCell scoring together.

- n_rounds:

  Number of boosting rounds used by GRNBoost2.

- learning_rate:

  Learning rate used by GRNBoost2.

- max_depth:

  Maximum tree depth used by GRNBoost2.

- max_features:

  Fraction of features sampled by GRNBoost2.

- subsample:

  Row subsampling fraction used by GRNBoost2.

- early_stop_window_length:

  Early-stopping window used by GRNBoost2.

- cores:

  Number of workers used by GRNBoost2, `scenic ctx`, and AUCell scoring.
  If multicore execution is not supported, this is automatically reduced
  to one core.

- seed:

  Random seed used by GRNBoost2 and Seurat overclustering.

- force:

  Whether to rebuild existing SCENIC outputs.

- assay_name:

  Name of the assay used to store regulon activity scores.

- tool_name:

  Name of the `srt@tools` entry.

- return_seurat:

  Whether to return the modified Seurat object. If `FALSE`, a result
  list is returned.

- envname:

  Python environment used for SCENIC. If `NULL`, the isolated
  `"scenic_env"` environment is used.

- conda:

  The path or command name of a conda-compatible executable.

- prepare_env:

  Whether to prepare and configure the SCENIC Python environment before
  running.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with SCENIC results, or a result list when
`return_seurat = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSCENIC(
  pancreas_sub,
  species = "Mus_musculus"
)
} # }
```
