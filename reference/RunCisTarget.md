# Run cisTarget motif enrichment on a GRN adjacency table

cisTarget performs motif enrichment analysis on gene regulatory network
adjacency tables. For each transcription factor, it identifies target
genes whose regulatory regions are enriched for the TF's binding motifs,
producing regulons (TF + enriched target gene sets).

## Usage

``` r
RunCisTarget(
  adj,
  species = c("Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"),
  backend = c("python", "r"),
  ranking_dbs = NULL,
  motif_annotations = NULL,
  expression_mtx = NULL,
  ctx_output = NULL,
  min_regulon_size = 10,
  gmt_output = NULL,
  txt_output = NULL,
  work_dir = NULL,
  prefix = "cisTarget",
  data_dir = NULL,
  envname = NULL,
  conda = "auto",
  prepare_env = TRUE,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- adj:

  A data frame with columns `TF`, `target`, and optionally `importance`,
  or a path to a TSV adjacency file. This is typically the output of
  [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md),
  [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md),
  or
  [`RunSCENIC()`](https://mengxu98.github.io/scop/reference/RunSCENIC.md)
  step 1.

- species:

  Species used to select cisTarget reference files when `ranking_dbs`,
  `motif_annotations`, or `regulators` is `NULL`. Supported values are
  `"Homo_sapiens"`, `"Mus_musculus"`, and `"Drosophila_melanogaster"`.

- backend:

  cisTarget runtime backend.

  `"python"`

  :   Uses the pySCENIC/ctxcore Python pipeline. This is the most tested
      backend and produces results identical to official SCENIC.

  `"r"`

  :   Uses the `RcisTarget` Bioconductor package. Requires the cisTarget
      ranking databases and motif annotations to be available in the
      same format as the Python backend (feather files and motif2tf
      table).

- ranking_dbs:

  Character vector of cisTarget ranking feather files.

- motif_annotations:

  Motif annotation table path (motif2tf).

- expression_mtx:

  Optional expression matrix path (CSV). When `NULL` and
  `backend = "python"`, the expression matrix is reconstructed from the
  unique genes in `adj` as a minimal stub.

- ctx_output:

  Optional output file path for the cisTarget result.

- min_regulon_size:

  Minimum number of target genes per regulon.

- gmt_output:

  Optional output path for the regulon GMT file.

- txt_output:

  Optional output path for the regulon TXT file.

- work_dir:

  Working directory used by Python backend.

- prefix:

  Prefix for output files.

- data_dir:

  Directory used to cache automatically prepared SCENIC reference files.
  If `NULL`, files are stored under
  `tools::R_user_dir("scop", "data")/SCENIC/<species>`.

- envname:

  Python environment name (Python backend only).

- conda:

  Conda-compatible executable (Python backend only).

- prepare_env:

  Whether to prepare the Python environment.

- cores:

  Number of workers.

- force:

  Whether to rebuild existing outputs.

- verbose:

  Whether to print progress messages.

- ...:

  Additional arguments passed to the backend.

## Value

A named list with components `regulons` (list of gene vectors),
`ctx_file` (path to raw cisTarget output), `gmt_file`, and `txt_file`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)

# First run GRN inference
grn <- RunGRNBoost2(
  pancreas_sub,
  regulators = c("Neurod1", "Arx", "Pax6"),
  backend = "cpp"
)

# Then run cisTarget (Python backend)
regulons <- RunCisTarget(
  grn,
  species = "Mus_musculus",
  backend = "python"
)
} # }
```
