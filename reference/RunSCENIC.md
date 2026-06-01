# Run SCENIC gene regulatory network analysis

Run SCENIC from a Seurat object. The GRN and cisTarget steps use a
metacell count matrix by default, while AUCell scores are calculated on
the original single-cell count matrix.

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
  data_dir = NULL,
  prefix = "scenic",
  metacell = TRUE,
  metacell.by = NULL,
  metacell_resolution = NULL,
  metacell_target = NULL,
  metacell_resolution_candidates = c(0.5, 1, 2, 5, 10, 20, 30, 40, 50, 75, 100),
  metacell_reduction = "pca",
  metacell_dims = 1:30,
  min_expr_cells = 3,
  min_regulon_size = 10,
  cores = 1,
  aucell_batch_size = 500,
  aucell_backend = c("r", "cpp"),
  aucell_cpp_strategy = c("full", "sparse", "topk"),
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
  `motif_annotations`, or `regulators` is `NULL`. Supported values
  include `"Homo_sapiens"`, `"Mus_musculus"`,
  `"Drosophila_melanogaster"` and aliases such as `"human"`, `"mouse"`,
  and `"fly"`.

- data_dir:

  Directory used to cache automatically prepared SCENIC reference files.
  If `NULL`, files are stored under
  `tools::R_user_dir("scop", "data")/SCENIC/<species>`.

- prefix:

  Prefix for SCENIC output files.

- metacell:

  Whether to build a metacell count matrix for GRNBoost2.

- metacell.by:

  Optional metadata column(s) used to keep metacells within groups such
  as samples or cell types.

- metacell_resolution:

  Resolution passed to
  [`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)
  for overclustering. If `NULL`, candidate resolutions are scanned and
  the one closest to `metacell_target` is used.

- metacell_target:

  Target number of metacells used when `metacell_resolution = NULL`. If
  `NULL`, a default target is chosen from the number of cells.

- metacell_resolution_candidates:

  Candidate resolutions scanned when `metacell_resolution = NULL`.

- metacell_reduction:

  Reduction used to build the metacell neighbor graph. The default
  `"pca"` keeps the original behavior and recomputes PCA from the
  selected assay. To use an already batch-corrected embedding such as
  Harmony, run it before `RunSCENIC()` and set
  `metacell_reduction = "Harmony"`. Only the metacell grouping uses this
  reduction; GRNBoost2 still uses raw count sums per metacell.

- metacell_dims:

  Dimensions used for metacell overclustering.

- min_expr_cells:

  Minimum number of cells or metacells where a gene must be detected
  before GRNBoost2.

- min_regulon_size:

  Minimum regulon size kept after `scenic ctx`.

- cores:

  Number of workers used by GRNBoost2, `scenic ctx`, and AUCell batch
  scoring. If multicore execution is not supported, this is
  automatically reduced to one core.

- aucell_batch_size:

  Number of cells scored in each AUCell batch.

- aucell_backend:

  Backend used for AUCell regulon activity scoring. `"r"` uses the
  AUCell package. `"cpp"` uses the package C++ gene-set scoring
  implementation.

- aucell_cpp_strategy:

  C++ AUCell ranking strategy passed to the package gene-set scoring
  backend.

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
