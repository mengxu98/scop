# Run SCENICPlus-style eGRN analysis

Build a SCENICPlus-style eRegulon result from paired RNA and chromatin
assays. The native backend is a `scop` implementation and is not an
exact re-run of the official SCENIC+ Python workflow.

## Usage

``` r
RunSCENICPlus(
  srt,
  rna_assay = "RNA",
  atac_assay = "peaks",
  rna_layer = "counts",
  atac_layer = "counts",
  backend = c("cpp", "python"),
  grn_method = c("grnboost2", "regdiffusion", "genie3"),
  regulators = NULL,
  rna_expr = NULL,
  atac_expr = NULL,
  tf_gene_prior = NULL,
  region_gene_prior = NULL,
  tf_region_prior = NULL,
  triplets_prior = NULL,
  eregulons_prior = NULL,
  auc_expr = NULL,
  auc_rankings = NULL,
  region_gene_window = 250000,
  region_gene_search_space = NULL,
  region_gene_method = c("gbm", "correlation"),
  region_gene_n_rounds = 500,
  max_region_gene = 5,
  max_tf_region = 500,
  min_eregulon_size = 5,
  egrn_rho_threshold = 0.05,
  egrn_quantiles = c(0.85, 0.9, 0.95),
  egrn_top_n_region_gene = c(5L, 10L, 15L),
  egrn_min_target_genes = 10L,
  group.by = NULL,
  grn_top_targets = Inf,
  max_grn_targets = Inf,
  grn_target_scope = c("all", "region_gene"),
  grn_n_rounds = 5000,
  grn_learning_rate = 0.01,
  grn_max_depth = 3,
  grn_max_features = 0.1,
  grn_subsample = 0.9,
  grn_early_stop_window_length = 0,
  seed = 666,
  cores = 1,
  assay_name = "scenicplus",
  tool_name = "SCENICPlus",
  python_result_dir = NULL,
  scplus_object = NULL,
  envname = "scenicplus_env",
  conda = "auto",
  prepare_env = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing RNA and chromatin assays.

- rna_assay:

  RNA assay name.

- atac_assay:

  Chromatin assay name.

- rna_layer:

  RNA count layer used for TF-gene inference and AUC scoring.

- atac_layer:

  ATAC count layer used for peak-gene correlations.

- backend:

  Runtime backend. `"cpp"` uses the native package workflow. `"python"`
  checks the official Python environment but does not silently
  substitute the native result.

- grn_method:

  GRN method used for TF-gene inference. The native backend can call
  [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md)
  or
  [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md);
  RegDiffusion is supported only through the Python backend.

- regulators:

  Candidate transcription factors. If `NULL`, motif TF names are used
  when available; otherwise all RNA genes are considered.

- rna_expr:

  Optional expression matrix used for native TF-gene and region-gene
  inference. Rows should be genes and columns should be cells; a
  cell-by-gene matrix is accepted and transposed when cell names make
  the orientation unambiguous.

- atac_expr:

  Optional accessibility matrix used for native region-gene inference.
  Rows should be regions and columns should be cells; a cell-by-region
  matrix is accepted and transposed when cell names make the orientation
  unambiguous.

- tf_gene_prior:

  Optional precomputed TF-gene table with columns `TF`, `target`, and
  `importance`. When supplied with `backend = "cpp"`, native GRN
  inference is skipped and the table is used as the TF-gene evidence
  layer.

- region_gene_prior:

  Optional precomputed region-gene table with columns `region`, `gene`,
  and `score`. When supplied with `backend = "cpp"`, native peak-to-gene
  correlation is skipped.

- tf_region_prior:

  Optional precomputed TF-region table with columns `TF`, `region`, and
  `score`. When supplied with `backend = "cpp"`, existing motif
  incidence is not required.

- triplets_prior:

  Optional precomputed TF-region-gene table with columns `TF`, `region`,
  `gene`, and `score`. When supplied with `backend = "cpp"`, triplets
  are read directly instead of assembled from the three edge layers.

- eregulons_prior:

  Optional precomputed eRegulon table with columns `regulon` and
  `target`. When supplied with `backend = "cpp"`, eRegulon names and
  gene sets are taken from this table.

- auc_expr:

  Optional expression matrix used only for eRegulon AUCell scoring in
  the native backend. Rows should be genes and columns should be cells;
  a cell-by-gene matrix is accepted and transposed when cell names make
  the orientation unambiguous.

- auc_rankings:

  Optional precomputed 0-based ranking matrix used only for eRegulon
  AUCell scoring in the native backend. Rows should be cells and columns
  should be genes; a gene-by-cell matrix is accepted and transposed when
  cell names make the orientation unambiguous. This is mainly useful for
  official SCENIC+ parity tests because SCENIC+ uses seeded random
  tie-breaking.

- region_gene_window:

  Maximum peak-to-gene distance in base pairs.

- region_gene_search_space:

  Optional search-space table with columns `region`/`Name` and
  `gene`/`Gene`. When supplied, native region-gene scoring is limited to
  these candidate pairs.

- region_gene_method:

  Native region-gene scoring method. `"gbm"` follows the official
  SCENIC+ region-to-gene strategy more closely by combining boosting
  feature importance with Spearman correlation; `"correlation"` keeps
  the older correlation-only approximation.

- region_gene_n_rounds:

  Number of boosting rounds for native region-gene importance scoring.

- max_region_gene:

  Maximum region-gene links retained per gene.

- max_tf_region:

  Maximum motif-supported regions retained per TF.

- min_eregulon_size:

  Minimum genes per eRegulon retained for AUC scoring.

- egrn_rho_threshold:

  Absolute Spearman correlation threshold used to split TF-gene and
  region-gene links into activating/repressing eGRN modules.

- egrn_quantiles:

  Region-gene importance quantiles used to build candidate eModules
  before TF-gene leading-edge filtering.

- egrn_top_n_region_gene:

  Top-N region-gene links per gene used as additional candidate eModules
  before TF-gene leading-edge filtering.

- egrn_min_target_genes:

  Minimum leading-edge target genes required for an eRegulon module. The
  default follows the official SCENIC+ CLI workflow.

- group.by:

  Optional metadata column used to calculate eRegulon specificity
  scores.

- grn_top_targets:

  Maximum TF-gene edges retained per target. The default `Inf` keeps all
  positive links to match arboreto GRNBoost2 output.

- max_grn_targets:

  Maximum genes used as TF-gene GRN targets in the native backend. With
  `grn_target_scope = "region_gene"`, genes are ranked by strongest
  absolute peak-to-gene correlation. Set to `Inf` to keep the full
  selected target universe.

- grn_target_scope:

  TF-gene target universe for the native backend. `"all"` follows
  official SCENIC+ TF-to-gene inference by scoring every RNA gene;
  `"region_gene"` restricts TF-gene inference to genes with peak links.

- grn_n_rounds:

  Number of native tree boosting rounds for TF-gene inference. The
  default follows arboreto `SGBM_KWARGS`.

- grn_learning_rate:

  Native tree boosting learning rate.

- grn_max_depth:

  Maximum depth of each native TF-gene regression tree.

- grn_max_features:

  Fraction of candidate TFs sampled at each native tree split.

- grn_subsample:

  Fraction of cells sampled for each native boosting round.

- grn_early_stop_window_length:

  Out-of-bag improvement window used for native GRNBoost2 early
  stopping. The RunSCENICPlus default is `0`, which disables early
  stopping to match official SCENIC+ `SGBM_KWARGS`.

- seed:

  Random seed used by native stochastic boosting backends. The default
  matches official SCENIC+ command-line wrappers.

- cores:

  Number of workers used by native TF-gene GRN inference.

- assay_name:

  Assay used to store eRegulon AUC scores.

- tool_name:

  Name of the `srt@tools` result entry.

- python_result_dir:

  Directory containing official SCENIC+ outputs already exported as
  `tf_gene`, `region_gene`, `tf_region`, `triplets`, `eregulons`, and
  `auc` tables. Used only for `backend = "python"`.

- scplus_object:

  Optional official SCENIC+ Python object path (`.pkl`/`.pickle`) to
  extract into `python_result_dir` before standardization. Used only for
  `backend = "python"`.

- envname:

  Python environment used when `backend = "python"`.

- conda:

  Conda-compatible executable used when `backend = "python"`.

- prepare_env:

  Whether to prepare the Python environment for `backend = "python"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with SCENICPlus-style results.
