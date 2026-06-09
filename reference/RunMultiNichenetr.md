# Run MultiNicheNet analysis

Run MultiNicheNet analysis

## Usage

``` r
RunMultiNichenetr(
  srt,
  group.by,
  sample.by,
  condition.by,
  condition_oi,
  condition_reference,
  receiver_celltypes,
  sender_celltypes = NULL,
  assay = NULL,
  sample_agnostic = FALSE,
  contrast_tbl = NULL,
  batches = NULL,
  covariates = NULL,
  species = c("Homo_sapiens", "Mus_musculus"),
  lr_network = NULL,
  ligand_target_matrix = NULL,
  fraction_cutoff = 0.05,
  min_cells = 10,
  empirical_pval = TRUE,
  top_n_interactions = 250,
  backend = c("cpp", "r"),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Metadata column defining cell types.

- sample.by:

  Metadata column defining biological samples.

- condition.by:

  Metadata column defining conditions.

- condition_oi:

  Condition of interest.

- condition_reference:

  Reference condition.

- receiver_celltypes:

  Receiver cell types of interest.

- sender_celltypes:

  Sender cell types of interest. Default is all available cell types.

- assay:

  Assay to use.

- sample_agnostic:

  Whether to use the sample-agnostic MultiNicheNet wrapper.

- contrast_tbl:

  Optional contrast table passed to MultiNicheNet. If `NULL`, a simple
  contrast table will be created from `condition_oi` and
  `condition_reference`.

- batches:

  Optional metadata column(s) used as batches.

- covariates:

  Optional metadata column(s) used as covariates.

- species:

  Species for default NicheNet prior model loading.

- lr_network:

  Optional ligand-receptor prior model or path to an `.rds` file.

- ligand_target_matrix:

  Optional ligand-target prior model or path to an `.rds` file.

- fraction_cutoff:

  Minimum expression fraction cutoff used by MultiNicheNet.

- min_cells:

  Minimum number of cells per sample-celltype combination.

- empirical_pval:

  Whether to use empirical p-values.

- top_n_interactions:

  Number of top prioritized interactions kept in the standardized long
  table.

- backend:

  Backend used for post-processing aggregation. Upstream MultiNicheNet
  inference is unchanged.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with standardized MultiNicheNet results stored in
`srt@tools[["MultiNichenetr"]]`.
