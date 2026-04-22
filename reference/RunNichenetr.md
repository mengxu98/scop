# Run NicheNet analysis

Run NicheNet analysis

## Usage

``` r
RunNichenetr(
  srt,
  group.by,
  receiver,
  sender = "all",
  condition.by = NULL,
  condition_oi = NULL,
  condition_reference = NULL,
  mode = c("aggregate", "aggregate_cluster_de", "custom"),
  assay = NULL,
  expression_pct = 0.1,
  geneset = NULL,
  background_expressed_genes = NULL,
  top_n_ligands = 30,
  top_n_targets = 200,
  species = c("Homo_sapiens", "Mus_musculus"),
  lr_network = NULL,
  ligand_target_matrix = NULL,
  weighted_networks = NULL,
  cutoff_visualization = 0.33,
  lfc_cutoff = 0.25,
  use_sender_agnostic_background = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Metadata column defining cell types.

- receiver:

  Receiver cell type(s).

- sender:

  Sender cell type(s). Use `"all"` to use all non-receiver cell types.

- condition.by:

  Metadata column defining conditions.

- condition_oi:

  Condition of interest.

- condition_reference:

  Reference condition.

- mode:

  NicheNet wrapper mode. Supported values are `"aggregate"`,
  `"aggregate_cluster_de"`, and `"custom"`.

- assay:

  Assay to use.

- expression_pct:

  Minimum fraction of cells expressing a gene.

- geneset:

  Optional target gene set for `mode = "custom"`.

- background_expressed_genes:

  Optional background genes for `mode = "custom"`.

- top_n_ligands:

  Number of top ligands to keep in standardized summaries.

- top_n_targets:

  Number of top targets per ligand to keep in standardized summaries.

- species:

  Species for default NicheNet prior model loading.

- lr_network:

  Optional ligand-receptor prior model or path to an `.rds` file.

- ligand_target_matrix:

  Optional ligand-target prior model or path to an `.rds` file.

- weighted_networks:

  Optional weighted network list or path to an `.rds` file.

- cutoff_visualization:

  Cutoff passed to `nichenetr::prepare_ligand_target_visualization()`.

- lfc_cutoff:

  Log fold change cutoff for DE analysis in `"aggregate"` and
  `"aggregate_cluster_de"` modes. Default is `0.25`.

- use_sender_agnostic_background:

  Whether to use all sender cell types when `sender = "all"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with standardized NicheNet results stored in
`srt@tools[["Nichenetr"]]`.
