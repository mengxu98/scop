# Run an approximate SCENICPlus-compatible eGRN analysis

Build an approximate eRegulon result from paired RNA and chromatin
assays. The C++ workflow is not an exact re-run of the official SCENIC+
Python workflow and requires explicit opt-in.

## Usage

``` r
RunSCENICPlus(
  srt,
  rna_assay = "RNA",
  atac_assay = "peaks",
  rna_layer = "counts",
  atac_layer = "counts",
  backend = c("python", "cpp"),
  allow_approximate = FALSE,
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
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

  Runtime backend. `"python"` reads or extracts official SCENIC+ result
  tables and remains the default. `"cpp"` runs the approximate package
  workflow.

- allow_approximate:

  Whether to allow the approximate C++ workflow. This must be `TRUE`
  when `backend = "cpp"`.

- grn_method:

  GRN method used for TF-gene inference. The C++ workflow can call
  [`RunGRNBoost2()`](https://mengxu98.github.io/scop/reference/RunGRNBoost2.md),
  [`RunGENIE3()`](https://mengxu98.github.io/scop/reference/RunGENIE3.md),
  or
  [`RunGNIPLR()`](https://mengxu98.github.io/scop/reference/RunGNIPLR.md);
  RegDiffusion and GNIPLR are supported through the Python backend.

- regulators:

  Candidate transcription factors. If `NULL`, motif TF names are used
  when available; otherwise all RNA genes are considered.

- rna_expr:

  Optional expression matrix used for C++ TF-gene and region-gene
  inference. Rows should be genes and columns should be cells; a
  cell-by-gene matrix is accepted and transposed when cell names make
  the orientation unambiguous.

- atac_expr:

  Optional accessibility matrix used for C++ region-gene inference. Rows
  should be regions and columns should be cells; a cell-by-region matrix
  is accepted and transposed when cell names make the orientation
  unambiguous.

- tf_gene_prior:

  Optional precomputed TF-gene table with columns `TF`, `target`, and
  `importance`. When supplied with `backend = "cpp"`, C++ GRN inference
  is skipped and the table is used as the TF-gene evidence layer.

- region_gene_prior:

  Optional precomputed region-gene table with columns `region`, `gene`,
  and `score`. When supplied with `backend = "cpp"`, C++ peak-to-gene
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
  the C++ workflow. Rows should be genes and columns should be cells; a
  cell-by-gene matrix is accepted and transposed when cell names make
  the orientation unambiguous.

- auc_rankings:

  Optional precomputed 0-based ranking matrix used only for eRegulon
  AUCell scoring in the C++ workflow. Rows should be cells and columns
  should be genes; a gene-by-cell matrix is accepted and transposed when
  cell names make the orientation unambiguous. This is mainly useful for
  official SCENIC+ parity tests because SCENIC+ uses seeded random
  tie-breaking.

- region_gene_window:

  Maximum peak-to-gene distance in base pairs.

- region_gene_search_space:

  Optional search-space table with columns `region`/`Name` and
  `gene`/`Gene`. When supplied, C++ region-gene scoring is limited to
  these candidate pairs.

- region_gene_method:

  C++ region-gene scoring method. `"gbm"` follows the official SCENIC+
  region-to-gene strategy more closely by combining boosting feature
  importance with Spearman correlation; `"correlation"` keeps the older
  correlation-only approximation.

- region_gene_n_rounds:

  Number of boosting rounds for C++ region-gene importance scoring.

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

  Maximum genes used as TF-gene GRN targets in the C++ workflow. With
  `grn_target_scope = "region_gene"`, genes are ranked by strongest
  absolute peak-to-gene correlation. Set to `Inf` to keep the full
  selected target universe.

- grn_target_scope:

  TF-gene target universe for the C++ workflow. `"all"` follows official
  SCENIC+ TF-to-gene inference by scoring every RNA gene;
  `"region_gene"` restricts TF-gene inference to genes with peak links.

- grn_n_rounds:

  Number of C++ tree boosting rounds for TF-gene inference. The default
  follows arboreto `SGBM_KWARGS`.

- grn_learning_rate:

  C++ tree boosting learning rate.

- grn_max_depth:

  Maximum depth of each C++ TF-gene regression tree.

- grn_max_features:

  Fraction of candidate TFs sampled at each C++ tree split.

- grn_subsample:

  Fraction of cells sampled for each C++ boosting round.

- grn_early_stop_window_length:

  Out-of-bag improvement window used for C++ GRNBoost2 early stopping.
  The RunSCENICPlus default is `0`, which disables early stopping to
  match official SCENIC+ `SGBM_KWARGS`.

- seed:

  Random seed used by the C++ stochastic boosting code. The default
  matches official SCENIC+ command-line wrappers.

- cores:

  Number of workers used by C++ TF-gene GRN inference.

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

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with SCENICPlus-style results.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunSCENICPlus(
  pancreas_sub,
  backend = "python",
  python_result_dir = "test/scenicplus"
)
scenicplus_dot <- SCENICPlusPlot(
  pancreas_sub,
  group.by = "CellType"
)
example_regulons <- unique(scenicplus_dot$top_table$regulon)[1:2]
FeatureDimPlot(
  pancreas_sub,
  features = example_regulons,
  assay = "scenicplus",
  reduction = "StandardpcaUMAP2D"
)
SCENICPlusPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "eregulon_dim",
  features = example_regulons
)
} # }
```
