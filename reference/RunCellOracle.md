# Run official CellOracle perturbation analysis

Fit cluster-specific regulatory models and simulate gene perturbations
using the unmodified CellOracle 0.20.0 Python implementation.

## Usage

``` r
RunCellOracle(
  object,
  perturbation,
  group.by,
  base_grn = NULL,
  assay = "RNA",
  layer = "counts",
  features = NULL,
  reduction = NULL,
  stage_args = list(),
  model_dir = NULL,
  output_dir = NULL,
  store_expression = FALSE,
  tool_name = "CellOracle",
  envname = "celloracle_env",
  conda = "auto",
  docker_container = NULL,
  cores = 1L,
  seed = 1234L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with non-log RNA counts, clusters and a
  two-dimensional embedding. Perform cell/gene selection before calling
  this function.

- perturbation:

  Named numeric vector of non-negative expression values. Zero simulates
  knockout. Multiple genes specify one joint perturbation.

- group.by:

  Metadata column defining the clusters used for GRN fitting.

- base_grn:

  Official TF-info data.frame with `peak_id`, `gene_short_name` and TF
  columns, a path to its parquet file, or a named list mapping target
  genes to character vectors of candidate TFs (the official `TFdict`
  format). Required for a new fit; do not pass SCOP
  `TF/target/importance` edge tables. No expression-derived GRN is
  generated automatically.

- assay:

  Assay containing raw counts. Default is `"RNA"`.

- layer:

  Non-log count layer. Default is `"counts"`. Fractional counts are
  accepted, matching the official raw-count import (including the
  official demo input). Log-transformed, scaled, integrated and SCT
  residual data should not be supplied; their provenance cannot be
  inferred from numeric values alone.

- features:

  Genes to retain, in the requested order. `NULL` uses all genes in the
  assay. CellOracle recommends selecting approximately 1,000–3,000 genes
  before import; this wrapper does not silently select or normalize
  genes.

- reduction:

  Reduction containing the two-dimensional embedding. `NULL` uses
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md).
  The first two columns are used without recomputing.

- stage_args:

  Named list of named argument lists for official methods: `pca`
  (`perform_PCA`), `imputation` (`knn_imputation`), `links`
  (`get_links`), `filter` (`filter_links`), `fit`
  (`fit_GRN_for_simulation`), `simulation` (`simulate_shift`),
  `transition` (`estimate_transition_prob`), `embedding`
  (`calculate_embedding_shift`), `grid` (`calculate_p_mass`), and
  `mass_filter` (`calculate_mass_filter`). Unspecified arguments retain
  upstream defaults. SCOP controls the perturbation, cluster-specific
  fitting, randomization reference, worker arguments and seed. In
  CellOracle 0.20.0, `calculate_p_mass` internally fixes `smooth=0.8`
  and `n_jobs=-1`; the wrapper preserves that upstream behavior. Do not
  supply these two ineffective grid arguments. For small inputs
  explicitly choose suitable `k`, `n_pca_dims`, transition `n_neighbors`
  and grid `n_neighbors`; no silent clipping occurs.

- model_dir:

  Directory containing `model.celloracle.oracle` and
  `network.celloracle.links` saved by a previous run. Reuse requires
  identical cells, genes, raw counts, cluster labels and embedding.
  Supply `base_grn=NULL` and no fitting-stage arguments when reusing a
  model.

- output_dir:

  Optional new directory for official Oracle and Links HDF5 files.
  Existing directories are never overwritten. These files retain the
  full model and allow subsequent perturbations without refitting.

- store_expression:

  Store cell-by-gene `imputed_count`, `simulated_count` and `delta_X`
  matrices in the returned object. Default `FALSE` limits size.

- tool_name:

  Unused name for the result in `object@tools`. Existing entries are
  rejected; use a different name for another perturbation.

- envname:

  Isolated environment prepared with
  `PrepareEnv(modules = "celloracle")`. Default `"celloracle_env"`.

- conda:

  Conda-compatible executable, passed to
  [`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md).

- docker_container:

  Optional name or ID of an already running Linux Docker container whose
  `python` has CellOracle 0.20.0 installed. This allows Windows R to use
  the same official runner without a local Python installation. SCOP
  transfers inputs/results through a temporary directory; `envname` and
  `conda` are unused. Model directories always refer to the R host. The
  function does not pull images, install packages or start containers.

- cores:

  Positive worker count for official methods that honor it. The official
  grid routine may use all cores; see `stage_args`.

- seed:

  Non-negative integer seed. NumPy and standard Python RNG states are
  restored on exit; the upstream transition routine also seeds Numba.
  Docker launches its child with `PYTHONHASHSEED=seed` to make upstream
  set-derived predictor ordering reproducible across processes. Local
  Python inherits its interpreter's hash seed; separate reference runs
  must use the same hash seed to compare bagging coefficients exactly.

- verbose:

  Whether to display SCOP progress messages. Upstream diagnostic output
  and warnings remain visible. Docker output is also retained in the
  result's `diagnostics` character vector.

## Value

A Seurat object with `object@tools[[tool_name]]` containing `cell_shift`
(one row per input cell, with original coordinates, displacement and
randomized displacement), `grid_flow` (official grid vectors, mass and
exclusion mask), `gene_effect` (mean signed and mean absolute `delta_X`
across cells), `links` and `filtered_links` (upstream edge tables per
cluster), `parameters`, `cells`, `genes`, and optionally `expression`
and `model_dir`. Expression matrices have input cell and selected gene
dimnames. No live Python object is stored.

## Details

Requires official `celloracle==0.20.0`, either in a Linux/macOS Python
environment or in `docker_container`. Native Windows Python is not
supported. Prepare the environment or container explicitly before
analysis. The raw import performs CellOracle's own log transformation;
no additional library-size normalization is introduced. PCA, imputation,
network inference, edge filtering, coefficient fitting, simulation and
projection all run in CellOracle. Network edge p-values are not
perturbation-effect p-values. Results are predictions, not experimental
knockout measurements.

## References

Kamimoto et al. (2023). Dissecting cell identity via network inference
and in silico gene perturbation. Nature.
[doi:10.1038/s41586-022-05688-9](https://doi.org/10.1038/s41586-022-05688-9)
.

## See also

[`CellOraclePlot()`](https://mengxu98.github.io/scop/reference/CellOraclePlot.md),
[`PrepareEnv()`](https://mengxu98.github.io/scop/reference/PrepareEnv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# On Linux/macOS; prepare dependencies explicitly before analysis.
PrepareEnv(modules = "celloracle")
# On Windows use docker_container = "celloracle" with an already prepared,
# running Linux container instead of preparing local Python.
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
# Supply an official TF-info parquet prepared for the matching species.
pancreas_sub <- RunCellOracle(
  pancreas_sub, perturbation = c(Pdx1 = 0), group.by = "CellType",
  base_grn = "mouse_promoter_base_GRN.parquet",
  stage_args = list(
    imputation = list(k = 20L, n_pca_dims = 20L),
    transition = list(n_neighbors = 50L),
    grid = list(n_neighbors = 50L)
  )
)
CellOraclePlot(pancreas_sub, plot_type = "quiver")
} # }
```
